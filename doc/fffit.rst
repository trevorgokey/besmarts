Force Field Fitting
===================

The most interesting application, and largely the application that inspired
this package, is force field fitting with automatic chemical perception. What
this means is that we want the code to fit the force field parameters, but also
find the best parameters that give the best fit as well.

This example tries to split the parameter b4 from OpenFF 2.1.0 using a single
molecule. The force field can be found `here <https://github.com/openforcefield/openff-forcefields/blob/main/openforcefields/offxml/openff-2.1.0.offxml>`_.
The SMILES is "C1=CC=C(C(=O)Cl)O1" although we need a fully hydrogenated,
indexed SMILES in order for the code to work (to match the xyz coordinates
produced by a DFT optimization). The physical objectives to fit to are the
geometry and forces from a DFT calculation using B3LYP-D3BJ/DZVP in Psi4. This
example will first fit the physical parameters before splitting, and then try
to find new bond parameters that are able to further improve the fit (only
targeting b4). The physical parameter fits will only adjust the bonds, and so
angles, torsions, non-bonded, etc. are all unchanged throughout the run.
Furthermore, only the equilibrium lengths of the bonds will be modified.

.. code-block:: python

>>> from typing import Dict, List, Tuple
>>> from besmarts.mechanics import fits
>>> from besmarts.mechanics import smirnoff_xml
>>> from besmarts.mechanics import smirnoff_models
>>> from besmarts.mechanics import molecular_models as mm
>>> from besmarts.mechanics import optimizers_scipy
>>> from besmarts.mechanics import fits
>>> from besmarts.core import graphs
>>> from besmarts.core import topology
>>> from besmarts.core import perception
>>> from besmarts.core import arrays
>>> from besmarts.core import assignments
>>> from besmarts.core import optimization
>>> from besmarts.assign import hierarchy_assign_rdkit
>>> from besmarts.codecs import codec_rdkit
>>> from besmarts.core import compute
>>> from besmarts.core import configs
>>> 
>>> configs.processors = 1
>>> configs.remote_compute_enable = False
>>> configs.workqueue_port = 54321
>>> 
>>> xyz_positions = """11
>>> 
>>>   C -1.448194 -0.849408  0.168489
>>>   C -1.594013  0.503187 -0.016781
>>>   C -0.273976  1.026226 -0.135035
>>>   C  0.580644 -0.047164 -0.013031
>>>   C  2.034612 -0.068609 -0.059252
>>>   O  2.728097  0.901087 -0.219099
>>>  Cl  2.762146 -1.707341  0.146556
>>>   O -0.138973 -1.200446  0.173518
>>>   H -2.152268 -1.658361  0.306090
>>>   H -2.527430  1.048099 -0.061800
>>>   H  0.029352  2.052732 -0.289658
>>> """
>>> xyz_grad = """11
>>> 
>>>   C      0.49755    0.17370   -0.04115
>>>   C     -0.00884    0.07632   -0.01031
>>>   C      0.20074   -0.69547    0.09073
>>>   C     -0.02955    1.24848   -0.17483
>>>   C      0.55229   -1.91119    0.25039
>>>   O     -0.15948    0.65794   -0.08724
>>>  Cl     -0.33030    0.82983   -0.10559
>>>   O     -0.73720   -0.66864    0.11909
>>>   H     -0.11502    0.11021   -0.01168
>>>   H     -0.00691    0.04737   -0.00649
>>>   H      0.02566   -0.05163    0.00657
>>> """
>>> def load_xyz(flist, indices=None) -> assignments.graph_db_row:
>>>     s = 0
>>>     N = None
>>>     gdr = assignments.graph_db_row()
>>>     for f in flist:
>>>         lines = f.split('\n')
>>>         if N is None:
>>>             N = int(lines[0].split()[0])
>>>         if indices is None:
>>>             indices = [*range(1,N+1)]
>>>         assert N == int(lines[0].split()[0])
>>>         for chunk in arrays.batched(lines, N+2):
>>>             if chunk and chunk[0]:
>>>                 sel = {}
>>>                 for i, line in enumerate(chunk, -2):
>>>                     if i >= 0:
>>>                         sym, x, y, z = line.split()[:4]
>>>                         sel[indices[i],] = list(map(float, (x, y, z)))
>>> 
>>>                 gdc = assignments.graph_db_column()
>>>                 gdc.selections.update(sel)
>>>                 gdr.columns[s] = gdc
>>>                 s += 1
>>>     return gdr
>>> 
>>> def make():
>>>     smi = "[C:1]1([H:9])=[C:2]([H:10])[C:3]([H:11])=[C:4]([C:5](=[O:6])[Cl:7])[O:8]1"
>>>     s = xyz_positions
>>>     g = xyz_grad
>>>     d  = {
>>>         smi: [
>>>             {
>>>                 assignments.POSITIONS: s,
>>>                 assignments.GRADIENTS: g,
>>>             },
>>>         ],
>>>     }
>>>     return d
>>> 
>>> def new_gdb(f: Dict[str, List[str]]) -> assignments.graph_db:
>>>     gcd = codec_rdkit.graph_codec_rdkit()
>>>     gdb = assignments.graph_db()
>>> 
>>>     ne = 0
>>>     for smi, fn_dict in f.items():
>>> 
>>>         g = gcd.smiles_decode(smi)
>>>         gid = assignments.graph_db_add_graph(gdb, smi, g)
>>> 
>>>         gdb.graphs[gid] = g
>>>         gdb.smiles[gid] = smi
>>>         gdb.selections[topology.index_of(topology.atom)] = {
>>>             gid: {k: v for k, v in enumerate(graphs.graph_atoms(g))}
>>>         }
>>>         gde = assignments.graph_db_entry()
>>>         gdb.entries[len(gdb.entries)] = gde
>>>         for rid, rdata in enumerate(fn_dict):
>>>             tid = assignments.POSITIONS
>>>             gdt = assignments.graph_db_table(topology.atom)
>>>             gdg = assignments.graph_db_graph()
>>>             gdt.graphs[gid] = gdg
>>>             fn = rdata[tid]
>>>             # indices=dict(sorted([(j, x) for j, x in enumerate(g.nodes, 1)], key=lambda x: x[1]))
>>>             indices = None
>>>             gdr = load_xyz([fn], indices=indices)
>>>             gdg.rows[0] = gdr
>>>             gde.tables[tid] = gdt
>>>             tid = assignments.GRADIENTS
>>>             if tid in rdata:
>>>                 gdt = assignments.graph_db_table(topology.atom)
>>>                 gdg = assignments.graph_db_graph()
>>>                 gdt.graphs[gid] = gdg
>>>                 fn = rdata[tid]
>>>                 # indices=dict(sorted([(j, x) for j, x in enumerate(g.nodes)], key=lambda x: x[1]))
>>>                 gdr = load_xyz([fn], indices=indices)
>>>                 gdg.rows[0] = gdr
>>>                 gde.tables[tid] = gdt
>>>                 gx = [x for y in gdr[0].selections.values() for x in y]
>>>                 gdt.values.extend(gx)
>>>             tid = assignments.ENERGY
>>>             if tid in rdata:
>>>                 gdt = assignments.graph_db_table(topology.null)
>>>                 fn = rdata[tid]
>>>                 ene = [*map(float,
>>>                     [x for x in open(fn).read().split('\n') if x]
>>>                 )]
>>>                 gdt.values.extend(ene)
>>>                 gde.tables[tid] = gdt
>>>     return gdb
>>> 
>>> def run(d, ff_fn):
>>>     # build the dataset and input ff
>>>     gcd = codec_rdkit.graph_codec_rdkit()
>>>     labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
>>>     pcp = perception.perception_model(gcd, labeler)
>>>     csys = smirnoff_models.smirnoff_load(ff_fn, pcp)
>>>     gdb = new_gdb(d)
>>>     psys = fits.gdb_to_physical_systems(gdb, csys)
>>>     models = {0: ["b4"]}
>>>     strat = fits.forcefield_optimization_strategy_default(csys, models=models)
>>>     co = fits.chemical_objective
>>> 
>>>     fit_models = [0]
>>>     final = fits.objective_tier()
>>>     final.objectives = {
>>>         0: fits.objective_config_position(
>>>                 fits.graph_db_address(
>>>                     eid=[0],
>>>                 ),
>>>                 scale=1
>>>         ),
>>>         1: fits.objective_config_gradient(
>>>                 fits.graph_db_address(
>>>                     eid=[0],
>>>                 ),
>>>                 scale=1
>>>         ),
>>>     }
>>>     final.objectives[0].verbose = 2
>>>     final.objectives[1].verbose = 2
>>>     # final.key_filter = lambda x: x[0] in fit_models and x[1] == 'l'
>>>     final.fit_models = fit_models
>>>     final.fit_symbols = ["l"]
>>> 
>>>     final.method = "L-BFGS-B" 
>>> 
>>>     onestep = fits.objective_tier()
>>>     onestep.objectives = final.objectives
>>>     onestep.step_limit = 2
>>>     onestep.accept = 3
>>>     # onestep.key_filter = lambda x: x[0] in fit_models and x[1] == 'l'
>>>     onestep.fit_models = fit_models
>>>     onestep.fit_symbols = ["l"]
>>>     onestep.method = "L-BFGS-B" 
>>> 
>>>     tiers = [onestep] # have a default
>>> 
>>>     initial = final
>>> 
>>>     kv0 = mm.chemical_system_iter_keys(csys)
>>>     newcsys, (P0, P), (C0, C) = fits.ff_optimize(
>>>         csys,
>>>         gdb,
>>>         psys,
>>>         strat,
>>>         co,
>>>         initial,
>>>         tiers,
>>>         final
>>>     )
>>> 
>>>     print("Modified parameters:")
>>>     kv = mm.chemical_system_iter_keys(newcsys)
>>>     for k, v in kv.items():
>>>         v0 = kv0.get(k)
>>>         if v0 is not None:
>>>             dv = v-v0
>>>             if abs(dv) > 1e-7:
>>>                 print(f"{str(k):20s} | New: {v:12.6g} Ref {v0:12.6g} Diff {dv:12.6g}")
>>>         else:
>>>             print(f"{str(k):20s} + New: {v:12.6g}")
>>>     print("Initial objectives:")
>>>     X0 = P0 + C0
>>>     X = P + C
>>>     print(f"Total= {X0:15.8g} Physical {P0:15.8g} Chemical {C0:15.8g}")
>>>     print("Final objectives:")
>>>     print(f"Total= {X:15.8g} Physical {P:15.8g} Chemical {C:15.8g}")
>>>     print("Differences:")
>>>     print(f"Total= {100*(X-X0)/X0:14.2f}% Physical {100*(P-P0)/P0:14.2f}% Chemical {100*(C-C0)/C0:14.2f}%")
>>> 
>>> run(make(), "openff-2.1.0.offxml")

A few important parameters need some explanation. The `onestep` objective tier
is a filtering device to prevent wasting time on trying to perform costly fits
on parameters that are not promising. The `onestep.step_limit` inidicates only
two fitting steps will be done, and `onestep.accept` indicates that the top 3
candidates will be passed on to the `final` tier. In this tier, a full fit is
performed, and the best parameter is accepted and incorporated into the
parameter set. Notice that we only fit bonds (model 0 in `fit_models`, and we
only try to split on b4 as defined by the `models` dictionary that is passed to
the `forcefield_optimization_strategy` class. We also indicate that we only want
to fit to equilibrium lengths as given by the parameter term symbol "l" (one
could also include "k" to also fit spring force constants).

Now for the a few snippets from the output:


-- Snip 1 --
.. code-block::

    ### BESMARTS chemical perception on the following assignments ###
    Model:
    Tree:
     0   0 Bonds  
     1   0  b4   [#6X3:1]-[#6X3:2] k: [540.3345953498] l: [1.466199291912]
     1   0  b6   [#6X3:1]=[#6X3:2] k: [898.589948525] l: [1.382361687103]
     1   0  b17  [#6X3:1]-[#8X2:2] k: [598.9859275918] l: [1.357746519746]
     1   0  b21  [#6:1]=[#8X1+0,#8X2+1:2] k: [1527.019744047] l: [1.221668642702]
     1   0  b70  [#6:1]-[#17:2] k: [368.4266150848] l: [1.722215272811]
     1   0  b85  [#6X3:1]-[#1:2] k: [775.3853383846] l: [1.081823673944]

-- Snip 2 --
.. code-block::

    >>> Initial Objective 1.6989e+05
    >>> Final Objective        32494
    >>> Percent change       -80.874%
    2024-12-17 12:35:18.639912 Computing chemical objective
    2024-12-17 12:35:18.782179 C0=1175.14392020606
    2024-12-17 12:35:18.782213 Initial objective: X=      33669.2 P=      32494.1 C=      1175.14
    (0, 'l', 'b4', 0)    | New:      1.43878 Ref       1.4662 Diff   -0.0274194
    (0, 'l', 'b6', 0)    | New:      1.37761 Ref      1.38236 Diff  -0.00475042
    (0, 'l', 'b17', 0)   | New:      1.36398 Ref      1.35775 Diff   0.00623027
    (0, 'l', 'b21', 0)   | New:      1.19576 Ref      1.22167 Diff   -0.0259073
    (0, 'l', 'b70', 0)   | New:      1.81793 Ref      1.72222 Diff    0.0957145
    (0, 'l', 'b85', 0)   | New:      1.08405 Ref      1.08182 Diff   0.00222776



-- Snip 3 --
.. code-block::

    >>> Initial Objective      32494
    >>> Final Objective        31987
    >>> Percent change       -1.5608%
    2024-12-17 12:35:21.879206 Macro objective: P=      31986.9 C=      1176.22 DX=     -506.089

-- Snip 4 --
.. code-block::

     0   0 Bonds  
     1   0  b4   [#6X3:1]-[#6X3:2] k: [540.3345953498] l: [1.4442513728922]
     2   0   B92  [#6X3:1]@;-[#6X3:2] k: [540.3345953498] l: [1.4333615336318999]


-- Snip 5 --
.. code-block::

    Modified parameters:
    (0, 'k', 'B92', 0)   + New:      540.335
    (0, 'l', 'b4', 0)    | New:      1.44425 Ref       1.4662 Diff   -0.0219478
    (0, 'l', 'b6', 0)    | New:      1.37741 Ref      1.38236 Diff  -0.00495234
    (0, 'l', 'b17', 0)   | New:      1.36453 Ref      1.35775 Diff    0.0067858
    (0, 'l', 'b21', 0)   | New:      1.19645 Ref      1.22167 Diff   -0.0252213
    (0, 'l', 'b70', 0)   | New:      1.82021 Ref      1.72222 Diff    0.0979973
    (0, 'l', 'b85', 0)   | New:      1.08326 Ref      1.08182 Diff   0.00143904
    (0, 'l', 'B92', 0)   + New:      1.43339
    Initial objectives:
    Total=       171069.08 Physical       169893.93 Chemical       1175.1439
    Final objectives:
    Total=       33163.051 Physical       31986.828 Chemical       1176.2232
    Differences:
    Total=         -80.61% Physical         -81.17% Chemical           0.09%

Most of the output is intermediate and diagnostics and is not shown here.
Instead, 4 snippets show the important results of the calculation. Snip 3 shows
that B92 was split from b4, which separates rings versus linear bonds where all
considered bonds are between atoms that are bonded to 3 other atoms including
the bond. In particular, it specialized the somewhat generic
`[#6X3:1]-[#6X3:2]` b4 bond with `[#6X3:1]@;-[#6X3:2]`. The new bond parameter
decreased the bond length to 1.433 A, while the original b4 parameter was
originally at 1.466 A and also decreased to 1.444 A in the final result. One
thing to examine is the objective change after the initial fit, but before any
bond parameters were added. From snip 2, we see that the total objective after
the first fit is 32494, or reduced by 80.874% of the objective produced by
OpenFF 2.1.0. The improvement due to adding the B91 was 31987, or 1.5608% from
the fit objective. Most of the improvement was from fitting the force. It is
important to note that the geometry was not ruined either; the initial geometry
RMSE was 0.105 A in 2.1.0. Fitting before splitting lead to a geometry RMSE of
0.107 A, and after adding the new bond parameter the geometry RMSE increased to
0.109 A, of which both are smaller than the original objective. On the other
hand, the initial force RMSE was 124.28 kJ/mol/A, and reduced to 54.35 kJ/mol/A
after the initial fit and 53.92 kJ/mol/A after the bond was split. We can also
observe that b70 changed the most due to the fit, which corresponded to the C-Cl
bond.

Keep in mind that these results are a simple example and not meant to be
accurate. In particular, the fit was focused solely on b4, which may or may not
be a parameter that needs examination in the first place. The fits fixed all
other degrees of freedom (angles, torsions, etc) in addition to freezing the
bond force constants. Likely better performance can be achieved by searching
more parameter space, and allowing more parameters to be fit during the search. 

As a final note, the setup and results shown here is still undergoing
development. In the future, the setup should be easier and shorter, and the
output should be cleaner and more informative.


