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
>>> configs.processors = 4
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
>>>         # 0: fits.objective_config_position(
>>>         #         fits.graph_db_address(
>>>         #             eid=[0],
>>>         #         ),
>>>         #         scale=1
>>>         # ),
>>>         1: fits.objective_config_gradient(
>>>                 fits.graph_db_address(
>>>                     eid=[0],
>>>                 ),
>>>                 scale=1
>>>         ),
>>>     }
>>>     # final.key_filter = lambda x: x[0] in fit_models and x[1] == 'l'
>>>     final.fit_models = fit_models
>>>     final.fit_symbols = ["l"]
>>> 
>>>     onestep = fits.objective_tier()
>>>     onestep.objectives = final.objectives
>>>     onestep.step_limit = 2
>>>     onestep.accept = 3
>>>     # onestep.key_filter = lambda x: x[0] in fit_models and x[1] == 'l'
>>>     onestep.fit_models = fit_models
>>>     onestep.fit_symbols = ["l"]
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

Now for the output:

.. code-block::

    Initial assignments:
    Model:
    Tree:
     0   0 Bonds  
     1   0  b4   [#6X3:1]-[#6X3:2] k: [540.3345953498] l: [1.466199291912]
     1   0  b6   [#6X3:1]=[#6X3:2] k: [898.589948525] l: [1.382361687103]
     1   0  b17  [#6X3:1]-[#8X2:2] k: [598.9859275918] l: [1.357746519746]
     1   0  b21  [#6:1]=[#8X1+0,#8X2+1:2] k: [1527.019744047] l: [1.221668642702]
     1   0  b70  [#6:1]-[#17:2] k: [368.4266150848] l: [1.722215272811]
     1   0  b85  [#6X3:1]-[#1:2] k: [775.3853383846] l: [1.081823673944]
    Tree:
     0   1 Angles  
    Tree:
     0   2 Torsions  
    Tree:
     0   3 OutOfPlanes  
    Tree:
    Tree:
     0   5 vdW   
    Tree:
    2024-06-02 19:25:07.386449 Computing physical objective
    workspace listening on local host. Remote connections prohibited.
    2024-06-02 19:25:07.615330 Calculating initial obj
    2024-06-02 19:25:07.615374 Starting physical parameter optimization
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.466199 d=              0
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.382362 d=              0
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.357747 d=              0
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.221669 d=              0
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.722215 d=              0
    Setting (0, 'l', 'b85', 0)   from        1.081824 to        1.081824 d=              0
    2024-06-02 19:25:07.616558 2024-06-02 19:25:07.616543 Generating 2 objectives
    2024-06-02 19:25:07.617307 2024-06-02 19:25:07.617296 Starting 26 tasks
    2024-06-02 19:25:11.617824 2024-06-02 19:25:11.617807 Calculating 26 tasks
      0000 | X2=    0.1174278 |g|=       1.1377
      0001 | X2=      1010692 |g|= 1.301469e+07
    >>> X2=      1010692 |g|=1.301469e+07
    2024-06-02 19:25:11.618670 2024-06-02 19:25:11.618664 Done. 26 tasks complete
    RUNNING THE L-BFGS-B CODE
    
               * * *
    
    Machine precision = 2.220D-16
     N =            6     M =           10
    
    At X0         0 variables are exactly at the bounds
    
    At iterate    0    f=  1.01069D+06    |proj g|=  8.85641D+06
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.466199 d=   -1.63901e-07
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.382362 d=  -1.545291e-07
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.357746 d=  -1.517775e-07
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.221669 d=  -1.365658e-07
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.863104 d=      0.1408885
    Setting (0, 'l', 'b85', 0)   from        1.081824 to        2.071849 d=      0.9900255
    2024-06-02 19:25:11.619286 2024-06-02 19:25:11.619278 Generating 2 objectives
    2024-06-02 19:25:11.619539 2024-06-02 19:25:11.619533 Starting 26 tasks
    2024-06-02 19:25:16.620010 2024-06-02 19:25:16.619995 Calculating 26 tasks
      0000 | X2=     3.157511 |g|=     7.094438
      0001 | X2= 5.405538e+07 |g|= 1.234751e+08
    >>> X2= 5.405539e+07 |g|=1.234751e+08
    2024-06-02 19:25:16.620822 2024-06-02 19:25:16.620815 Done. 26 tasks complete
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.466199 d=  -1.182603e-08
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.382362 d=  -1.114982e-08
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.357747 d=  -1.095128e-08
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.221669 d=  -9.853705e-09
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.732381 d=      0.0101656
    Setting (0, 'l', 'b85', 0)   from        1.081824 to        1.153257 d=     0.07143382
    2024-06-02 19:25:16.621137 2024-06-02 19:25:16.621130 Generating 2 objectives
    2024-06-02 19:25:16.621416 2024-06-02 19:25:16.621409 Starting 26 tasks
    2024-06-02 19:25:20.621828 2024-06-02 19:25:20.621812 Calculating 26 tasks
      0000 | X2=    0.1848241 |g|=      1.42369
      0001 | X2=     687962.2 |g|=      7067806
    >>> X2=     687962.4 |g|=     7067806
    2024-06-02 19:25:20.622626 2024-06-02 19:25:20.622618 Done. 26 tasks complete
    
    At iterate    1    f=  6.87962D+05    |proj g|=  1.16371D+06
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.468517 d=    0.002317368
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.356166 d=    -0.02619574
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.333313 d=    -0.02443394
    Setting (0, 'l', 'b21', 0)   from        1.221669 to         1.18993 d=     -0.0317383
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.738999 d=     0.01678331
    Setting (0, 'l', 'b85', 0)   from        1.081824 to          1.1424 d=     0.06057656
    2024-06-02 19:25:20.622977 2024-06-02 19:25:20.622969 Generating 2 objectives
    2024-06-02 19:25:20.623221 2024-06-02 19:25:20.623216 Starting 26 tasks
    2024-06-02 19:25:24.623661 2024-06-02 19:25:24.623640 Calculating 26 tasks
      0000 | X2=    0.1465604 |g|=     6.576907
      0001 | X2=       468257 |g|=      3665671
    >>> X2=     468257.1 |g|=     3665671
    2024-06-02 19:25:24.624264 2024-06-02 19:25:24.624256 Done. 26 tasks complete
    
    At iterate    2    f=  4.68257D+05    |proj g|=  7.85671D+05
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.458463 d=    -0.00773658
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.330903 d=    -0.05145825
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.296144 d=    -0.06160265
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.184767 d=    -0.03690201
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.747941 d=     0.02572617
    Setting (0, 'l', 'b85', 0)   from        1.081824 to        1.129634 d=     0.04781079
    2024-06-02 19:25:24.624655 2024-06-02 19:25:24.624647 Generating 2 objectives
    2024-06-02 19:25:24.624887 2024-06-02 19:25:24.624880 Starting 26 tasks
    2024-06-02 19:25:28.625297 2024-06-02 19:25:28.625280 Calculating 26 tasks
      0000 | X2=   0.02065438 |g|=    0.4807981
      0001 | X2=     343519.3 |g|=      2947268
    >>> X2=     343519.3 |g|=     2947268
    2024-06-02 19:25:28.626131 2024-06-02 19:25:28.626125 Done. 26 tasks complete
    
    At iterate    3    f=  3.43519D+05    |proj g|=  1.43971D+06
    Setting (0, 'l', 'b4', 0)    from        1.466199 to          1.4134 d=    -0.05279882
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.310344 d=    -0.07201766
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.230325 d=     -0.1274216
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.183258 d=    -0.03841096
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.770322 d=     0.04810685
    Setting (0, 'l', 'b85', 0)   from        1.081824 to         1.11173 d=     0.02990661
    2024-06-02 19:25:28.626483 2024-06-02 19:25:28.626476 Generating 2 objectives
    2024-06-02 19:25:28.626727 2024-06-02 19:25:28.626722 Starting 26 tasks
    2024-06-02 19:25:31.627089 2024-06-02 19:25:31.627070 Calculating 26 tasks
      0000 | X2=  0.008378642 |g|=    0.2610441
      0001 | X2=     223516.1 |g|=      1616430
    >>> X2=     223516.1 |g|=     1616430
    2024-06-02 19:25:31.627694 2024-06-02 19:25:31.627686 Done. 26 tasks complete
    
    At iterate    4    f=  2.23516D+05    |proj g|=  1.17292D+06
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.379461 d=    -0.08673801
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.280888 d=     -0.1014738
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.229064 d=     -0.1286827
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.195531 d=    -0.02613809
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.784425 d=     0.06220929
    Setting (0, 'l', 'b85', 0)   from        1.081824 to        1.096121 d=     0.01429765
    2024-06-02 19:25:31.628073 2024-06-02 19:25:31.628063 Generating 2 objectives
    2024-06-02 19:25:31.628297 2024-06-02 19:25:31.628291 Starting 26 tasks
    2024-06-02 19:25:35.628730 2024-06-02 19:25:35.628710 Calculating 26 tasks
      0000 | X2=  0.006316887 |g|=    0.1733382
      0001 | X2=     220114.7 |g|=      2089175
    >>> X2=     220114.7 |g|=     2089175
    2024-06-02 19:25:35.629373 2024-06-02 19:25:35.629364 Done. 26 tasks complete
    
    At iterate    5    f=  2.20115D+05    |proj g|=  1.63291D+06
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.384041 d=    -0.08215835
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.297088 d=    -0.08527407
    Setting (0, 'l', 'b17', 0)   from        1.357747 to         1.22681 d=     -0.1309367
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.188653 d=    -0.03301529
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.782929 d=     0.06071404
    Setting (0, 'l', 'b85', 0)   from        1.081824 to        1.102412 d=     0.02058866
    2024-06-02 19:25:35.629750 2024-06-02 19:25:35.629741 Generating 2 objectives
    2024-06-02 19:25:35.629968 2024-06-02 19:25:35.629961 Starting 26 tasks
    2024-06-02 19:25:39.630388 2024-06-02 19:25:39.630373 Calculating 26 tasks
      0000 | X2=  0.006330228 |g|=    0.1908843
      0001 | X2=     202074.9 |g|=       248332
    >>> X2=     202074.9 |g|=    248332.1
    2024-06-02 19:25:39.631211 2024-06-02 19:25:39.631205 Done. 26 tasks complete
    
    At iterate    6    f=  2.02075D+05    |proj g|=  1.48456D+05
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.382827 d=    -0.08337228
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.298374 d=    -0.08398804
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.228541 d=     -0.1292058
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.188312 d=      -0.033357
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.784104 d=     0.06188911
    Setting (0, 'l', 'b85', 0)   from        1.081824 to        1.102776 d=     0.02095277
    2024-06-02 19:25:39.631561 2024-06-02 19:25:39.631553 Generating 2 objectives
    2024-06-02 19:25:39.631802 2024-06-02 19:25:39.631796 Starting 26 tasks
    2024-06-02 19:25:43.632231 2024-06-02 19:25:43.632213 Calculating 26 tasks
      0000 | X2=  0.006248285 |g|=    0.1914787
      0001 | X2=     201629.8 |g|=     110490.9
    >>> X2=     201629.8 |g|=      110491
    2024-06-02 19:25:43.632839 2024-06-02 19:25:43.632830 Done. 26 tasks complete
    
    At iterate    7    f=  2.01630D+05    |proj g|=  9.40870D+04
    Setting (0, 'l', 'b4', 0)    from        1.466199 to        1.381922 d=    -0.08427753
    Setting (0, 'l', 'b6', 0)    from        1.382362 to        1.298647 d=     -0.0837149
    Setting (0, 'l', 'b17', 0)   from        1.357747 to        1.229711 d=     -0.1280353
    Setting (0, 'l', 'b21', 0)   from        1.221669 to        1.188198 d=    -0.03347093
    Setting (0, 'l', 'b70', 0)   from        1.722215 to        1.786276 d=     0.06406052
    Setting (0, 'l', 'b85', 0)   from        1.081824 to        1.102783 d=     0.02095931
    2024-06-02 19:25:43.633226 2024-06-02 19:25:43.633217 Generating 2 objectives
    2024-06-02 19:25:43.633447 2024-06-02 19:25:43.633441 Starting 26 tasks
    2024-06-02 19:25:47.633893 2024-06-02 19:25:47.633875 Calculating 26 tasks
      0000 | X2=  0.006199739 |g|=     0.190105
      0001 | X2=     201437.2 |g|=     84285.65
    >>> X2=     201437.2 |g|=    84285.65
    2024-06-02 19:25:47.634525 2024-06-02 19:25:47.634516 Done. 26 tasks complete
    
    At iterate    8    f=  2.01437D+05    |proj g|=  6.96895D+04
    
               * * *
    
    Tit   = total number of iterations
    Tnf   = total number of function evaluations
    Tnint = total number of segments explored during Cauchy searches
    Skip  = number of BFGS updates skipped
    Nact  = number of active bounds at final generalized Cauchy point
    Projg = norm of the final projected gradient
    F     = final function value
    
               * * *
    
       N    Tit     Tnf  Tnint  Skip  Nact     Projg        F
        6      8     10     12     0     0   6.969D+04   2.014D+05
      F =   201437.23375790875     
    
    CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH             
    >>> Initial Objective 1.0107e+06
    >>> Final Objective   2.0144e+05
    >>> Percent change       -80.069%
    2024-06-02 19:25:47.714874 Computing chemical objective
    2024-06-02 19:25:47.885242 C0=2014.3723375790876
    2024-06-02 19:25:47.885277 Initial objective: X=       203452 P=       201437 C=      2014.37
    (0, 'l', 'b4', 0)    | New:      1.38192 Ref      1.38192 Diff            0
    (0, 'l', 'b6', 0)    | New:      1.29865 Ref      1.29865 Diff            0
    (0, 'l', 'b17', 0)   | New:      1.22971 Ref      1.22971 Diff            0
    (0, 'l', 'b21', 0)   | New:       1.1882 Ref       1.1882 Diff            0
    (0, 'l', 'b70', 0)   | New:      1.78628 Ref      1.78628 Diff            0
    (0, 'l', 'b85', 0)   | New:      1.10278 Ref      1.10278 Diff            0
    Optimization strategy is building steps...
    2024-06-02 19:25:47.885493 The optimization strategy has the following iterations:
    ->   1:00. op= 1 m=[0] a=[0] b=1->1 d=0->0 n=0->0
         2:00. op=-1 m=[0] a=[0] b=0->0 d=0->0 n=0->0
    Targets for this macro step 1:
    1 (0, 0, 0) b4
    2 (0, 0, 0) b6
    3 (0, 0, 0) b17
    4 (0, 0, 0) b21
    5 (0, 0, 0) b70
    6 (0, 0, 0) b85
    N Targets: 6
    Step tracker for current macro step 1
    
    *******************
     iteration=   0 macro=  1/2 micro=1 X=2.0345e+05 P=2.0144e+05 C=   2014.4 models=0:Bonds
    *******************
    
    2024-06-02 19:25:47.892950 Saving checkpoint to chk.cst.p
    2024-06-02 19:25:47.921788 Collecting SMARTS for b4 and setting to depth=0
    
     == iteration=   1 macro=  1/2 micro=  1/1 operation=1 cluster=b4   N= 2 overlap=[0] bits=1->1 depth=0->0 branch=0->0
    
    Attempting to split b4:
    S0: [#6X3:1]-[#6X3:2]
    Matched N=2
    000001 (0, (2, 3))              [#6H1X3x2r5A+0:2]@;-[#6H1X3x2r5A+0:3]
    000002 (0, (4, 5))              [#6H0X3x2r5A+0:4]!@;-[#6H0X3x0!rA+0:5]
    
    2024-06-02 19:25:47.924682 Union merging=2
    2024-06-02 19:25:47.940748 Union merging=1
    2024-06-02 19:25:47.941078 Union is [#6;H0,H1;X3;x2;r5;A;+0:2]-[#6;H0,H1;X3;x0,x2;!r,r5;A;+0:3]
    2024-06-02 19:25:47.941682 Generating splits
    2024-06-02 19:25:47.942463 Generating single splits
    2024-06-02 19:25:47.943520 Generated 20 splits
    BIT [_H0_____:1]_;_[_______:2]
    BIT [_!H0_____:1]_;_[_______:2]
    BIT [_H1_____:1]_;_[_______:2]
    BIT [_!H1_____:1]_;_[_______:2]
    BIT [_______:1]!@;_[_______:2]
    BIT [_______:1]@;_[_______:2]
    BIT [_______:1]@;_[_______:2]
    BIT [_______:1]!@;_[_______:2]
    BIT [_______:1]_;_[___x0___:2]
    BIT [_______:1]_;_[___!x0___:2]
    BIT [_______:1]_;_[___x2___:2]
    BIT [_______:1]_;_[___!x2___:2]
    BIT [_______:1]_;_[____!r__:2]
    BIT [_______:1]_;_[____r__:2]
    BIT [_______:1]_;_[____r5__:2]
    BIT [_______:1]_;_[____!r5__:2]
    2024-06-02 19:25:47.954355 Building tasks
    workspace listening on local host. Remote connections prohibited.
    Progress: 100.00%        32/32
    Finished: 100.00%        32/32
    Removing workspace ('127.0.0.1', 42773)
    Closing workspace
    2024-06-02 19:25:52.248673 Calculating partitions for hits=2
    workspace listening on local host. Remote connections prohibited.
    Submitting 2 packets of work
    Chunk: 100.00%         2/2
    Finished: 100.00%         2/2
    Removing workspace ('127.0.0.1', 34663)
    Closing workspace
    2024-06-02 19:25:55.535221 Unique hits 2/2
    2024-06-02 19:25:55.535244 Searching atoms done; data=2 hits=2
    2024-06-02 19:25:55.535771 Collecting new candidates
    2024-06-02 19:25:55.535829 Scanning done.
    2024-06-02 19:25:55.535838
    
    
    Generating SMARTS on 2
    Model:
    Tree:
     0   0 Bonds  
     1   0  b4   [#6X3:1]-[#6X3:2] k: [540.3345953498] l: [1.3819217653738018]
     1   0  b6   [#6X3:1]=[#6X3:2] k: [898.589948525] l: [1.2986467914734003]
     1   0  b17  [#6X3:1]-[#8X2:2] k: [598.9859275918] l: [1.2297112213860095]
     1   0  b21  [#6:1]=[#8X1+0,#8X2+1:2] k: [1527.019744047] l: [1.1881977141798903]
     1   0  b70  [#6:1]-[#17:2] k: [368.4266150848] l: [1.786275796097973]
     1   0  b85  [#6X3:1]-[#1:2] k: [775.3853383846] l: [1.1027829873202926]
    Tree:
     0   1 Angles  
    Tree:
     0   2 Torsions  
    Tree:
     0   3 OutOfPlanes  
    Tree:
    Tree:
     0   5 vdW   
    Tree:
    Scoring and filtering 2 candidates for operation=1
    Tier 0: Scoring and filtering 2 candidates for operation=1
    Tier 0: Accepting all candidates so we skip
    Scanning 2 candidates for operation=1
    2024-06-02 19:25:55.549736 Generated 2 x 1 = 2 candidate evalulation tasks
    2024-06-02 19:25:55.549760 Dispatching candidate tasks= 2 in serial
    2024-06-02 19:25:55.549772 Running candidate task 1/2
    workspace listening on local host. Remote connections prohibited.
    2024-06-02 19:25:56.208546 Calculating initial obj
    2024-06-02 19:25:56.208595 Starting physical parameter optimization
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.381922 d=              0
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.298647 d=              0
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.229711 d=              0
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.188198 d=              0
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.786276 d=              0
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102783 d=              0
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.381922 d=              0
    2024-06-02 19:25:56.209729 2024-06-02 19:25:56.209716 Generating 2 objectives
    2024-06-02 19:25:56.210160 2024-06-02 19:25:56.210152 Starting 30 tasks
    2024-06-02 19:26:00.212932 2024-06-02 19:26:00.212917 Calculating 30 tasks
      0000 | X2=  0.006199739 |g|=    0.1915674
      0001 | X2=     201437.2 |g|=      2131355
    >>> X2=     201437.2 |g|=     2131355
    2024-06-02 19:26:00.213879 2024-06-02 19:26:00.213872 Done. 30 tasks complete
    RUNNING THE L-BFGS-B CODE
    
               * * *
    
    Machine precision = 2.220D-16
     N =            7     M =           10
    
    At X0         0 variables are exactly at the bounds
    
    At iterate    0    f=  2.01437D+05    |proj g|=  1.52718D+06
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        2.380856 d=      0.9989342
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.298646 d=  -8.494508e-07
    Setting (0, 'l', 'b17', 0)   from        1.229711 to         1.22971 d=  -8.043597e-07
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.195441 d=    0.007243539
    Setting (0, 'l', 'b70', 0)   from        1.786276 to         1.83186 d=      0.0455842
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102782 d=  -7.213354e-07
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.381921 d=  -9.039213e-07
    2024-06-02 19:26:00.214463 2024-06-02 19:26:00.214455 Generating 2 objectives
    2024-06-02 19:26:00.214747 2024-06-02 19:26:00.214741 Starting 30 tasks
    2024-06-02 19:26:05.226661 2024-06-02 19:26:05.226645 Calculating 30 tasks
      0000 | X2=     2.191041 |g|=     4.455478
      0001 | X2=      8633944 |g|= 2.739841e+07
    >>> X2=      8633946 |g|=2.739841e+07
    2024-06-02 19:26:05.227606 2024-06-02 19:26:05.227599 Done. 30 tasks complete
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.458577 d=     0.07665549
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.298647 d=  -6.518453e-08
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.229711 d=  -6.172437e-08
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.188754 d=   0.0005558494
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.789774 d=    0.003498007
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102783 d=  -5.535331e-08
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.381922 d=  -6.936445e-08
    2024-06-02 19:26:05.228012 2024-06-02 19:26:05.228003 Generating 2 objectives
    2024-06-02 19:26:05.228295 2024-06-02 19:26:05.228290 Starting 30 tasks
    2024-06-02 19:26:10.238915 2024-06-02 19:26:10.238897 Calculating 30 tasks
      0000 | X2=   0.03403521 |g|=     0.933406
      0001 | X2=     142778.9 |g|=      2143720
    >>> X2=     142778.9 |g|=     2143720
    2024-06-02 19:26:10.239619 2024-06-02 19:26:10.239610 Done. 30 tasks complete
    
    At iterate    1    f=  1.42779D+05    |proj g|=  1.21476D+06
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.496818 d=      0.1148964
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.318481 d=     0.01983407
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.240725 d=     0.01101365
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.219009 d=      0.0308113
    Setting (0, 'l', 'b70', 0)   from        1.786276 to         1.79802 d=      0.0117447
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102538 d=  -0.0002447489
    Setting (0, 'l', 'B91', 0)   from        1.381922 to         1.34528 d=    -0.03664141
    2024-06-02 19:26:10.240021 2024-06-02 19:26:10.240012 Generating 2 objectives
    2024-06-02 19:26:10.240274 2024-06-02 19:26:10.240267 Starting 30 tasks
    2024-06-02 19:26:15.272536 2024-06-02 19:26:15.272513 Calculating 30 tasks
      0000 | X2=   0.06321656 |g|=     1.105779
      0001 | X2=     129021.1 |g|=      3407824
    >>> X2=     129021.1 |g|=     3407824
    2024-06-02 19:26:15.273528 2024-06-02 19:26:15.273517 Done. 30 tasks complete
    
    At iterate    2    f=  1.29021D+05    |proj g|=  6.00914D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.480683 d=     0.09876135
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.311205 d=     0.01255816
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.242294 d=     0.01258272
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.202336 d=     0.01413863
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.799835 d=     0.01355879
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.104195 d=    0.001412117
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.344783 d=    -0.03713897
    2024-06-02 19:26:15.273925 2024-06-02 19:26:15.273918 Generating 2 objectives
    2024-06-02 19:26:15.274204 2024-06-02 19:26:15.274199 Starting 30 tasks
    2024-06-02 19:26:20.289341 2024-06-02 19:26:20.289322 Calculating 30 tasks
      0000 | X2=   0.04072391 |g|=    0.2894206
      0001 | X2=     91126.56 |g|=      1048334
    >>> X2=      91126.6 |g|=     1048334
    2024-06-02 19:26:20.290047 2024-06-02 19:26:20.290039 Done. 30 tasks complete
    
    At iterate    3    f=  9.11266D+04    |proj g|=  3.68215D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.477414 d=     0.09549244
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.308033 d=    0.009386385
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.246155 d=     0.01644362
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.198106 d=    0.009908732
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.804497 d=     0.01822149
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102794 d=   1.143844e-05
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.334101 d=    -0.04782061
    2024-06-02 19:26:20.290468 2024-06-02 19:26:20.290460 Generating 2 objectives
    2024-06-02 19:26:20.290717 2024-06-02 19:26:20.290711 Starting 30 tasks
    2024-06-02 19:26:25.321009 2024-06-02 19:26:25.320992 Calculating 30 tasks
      0000 | X2=   0.04073447 |g|=    0.3875965
      0001 | X2=     81438.74 |g|=     683347.1
    >>> X2=     81438.78 |g|=    683347.2
    2024-06-02 19:26:25.321979 2024-06-02 19:26:25.321971 Done. 30 tasks complete
    
    At iterate    4    f=  8.14388D+04    |proj g|=  2.71606D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.475884 d=     0.09396177
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.302319 d=    0.003672524
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.250245 d=     0.02053391
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.196968 d=    0.008770196
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.814084 d=     0.02780814
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.097224 d=   -0.005558942
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.313809 d=    -0.06811315
    2024-06-02 19:26:25.322401 2024-06-02 19:26:25.322392 Generating 2 objectives
    2024-06-02 19:26:25.322693 2024-06-02 19:26:25.322687 Starting 30 tasks
    2024-06-02 19:26:28.323080 2024-06-02 19:26:28.323061 Calculating 30 tasks
      0000 | X2=   0.01546186 |g|=    0.1913408
      0001 | X2=     71978.99 |g|=     613166.8
    >>> X2=        71979 |g|=    613166.8
    2024-06-02 19:26:28.323835 2024-06-02 19:26:28.323824 Done. 30 tasks complete
    
    At iterate    5    f=  7.19790D+04    |proj g|=  4.04255D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.472374 d=      0.0904523
    Setting (0, 'l', 'b6', 0)    from        1.298647 to         1.29622 d=   -0.002426962
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.241563 d=     0.01185217
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.198028 d=    0.009829893
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.825393 d=     0.03911703
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.091226 d=    -0.01155689
    Setting (0, 'l', 'B91', 0)   from        1.381922 to         1.29802 d=    -0.08390154
    2024-06-02 19:26:28.324277 2024-06-02 19:26:28.324267 Generating 2 objectives
    2024-06-02 19:26:28.324529 2024-06-02 19:26:28.324520 Starting 30 tasks
    2024-06-02 19:26:32.331836 2024-06-02 19:26:32.331821 Calculating 30 tasks
      0000 | X2=   0.01392054 |g|=    0.3636936
      0001 | X2=     66325.96 |g|=     329329.5
    >>> X2=     66325.98 |g|=    329329.6
    2024-06-02 19:26:32.332765 2024-06-02 19:26:32.332759 Done. 30 tasks complete
    
    At iterate    6    f=  6.63260D+04    |proj g|=  2.88255D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.474668 d=     0.09274669
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.296168 d=   -0.002478367
    Setting (0, 'l', 'b17', 0)   from        1.229711 to         1.24009 d=      0.0103784
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.200627 d=     0.01242952
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.831228 d=     0.04495185
    Setting (0, 'l', 'b85', 0)   from        1.102783 to         1.08889 d=    -0.01389342
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.295411 d=     -0.0865107
    2024-06-02 19:26:32.333140 2024-06-02 19:26:32.333133 Generating 2 objectives
    2024-06-02 19:26:32.333415 2024-06-02 19:26:32.333410 Starting 30 tasks
    2024-06-02 19:26:37.337889 2024-06-02 19:26:37.337875 Calculating 30 tasks
      0000 | X2=   0.01465877 |g|=    0.3703143
      0001 | X2=     65711.61 |g|=     166508.1
    >>> X2=     65711.63 |g|=    166507.9
    2024-06-02 19:26:37.338822 2024-06-02 19:26:37.338816 Done. 30 tasks complete
    
    At iterate    7    f=  6.57116D+04    |proj g|=  1.43515D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.474238 d=     0.09231586
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.295682 d=   -0.002965275
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.239183 d=    0.009471377
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.200692 d=     0.01249471
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.831539 d=     0.04526368
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.089854 d=    -0.01292852
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.295551 d=      -0.086371
    2024-06-02 19:26:37.339220 2024-06-02 19:26:37.339213 Generating 2 objectives
    2024-06-02 19:26:37.339502 2024-06-02 19:26:37.339497 Starting 30 tasks
    2024-06-02 19:26:42.348906 2024-06-02 19:26:42.348889 Calculating 30 tasks
      0000 | X2=   0.01472726 |g|=    0.3829256
      0001 | X2=     65612.06 |g|=     71994.11
    >>> X2=     65612.08 |g|=    71994.11
    2024-06-02 19:26:42.349858 2024-06-02 19:26:42.349848 Done. 30 tasks complete
    
    At iterate    8    f=  6.56121D+04    |proj g|=  3.16585D+04
    Setting (0, 'l', 'b4', 0)    from        1.381922 to         1.47421 d=     0.09228864
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.295903 d=   -0.002744068
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.239008 d=    0.009297198
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.200526 d=     0.01232822
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.831825 d=     0.04554951
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.089911 d=    -0.01287199
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.295963 d=    -0.08595864
    2024-06-02 19:26:42.350279 2024-06-02 19:26:42.350269 Generating 2 objectives
    2024-06-02 19:26:42.350568 2024-06-02 19:26:42.350563 Starting 30 tasks
    2024-06-02 19:26:47.353206 2024-06-02 19:26:47.353192 Calculating 30 tasks
      0000 | X2=   0.01468777 |g|=    0.3795553
      0001 | X2=     65586.96 |g|=     41687.39
    >>> X2=     65586.97 |g|=    41687.31
    2024-06-02 19:26:47.354153 2024-06-02 19:26:47.354144 Done. 30 tasks complete
    
    At iterate    9    f=  6.55870D+04    |proj g|=  2.71753D+04
    
               * * *
    
    Tit   = total number of iterations
    Tnf   = total number of function evaluations
    Tnint = total number of segments explored during Cauchy searches
    Skip  = number of BFGS updates skipped
    Nact  = number of active bounds at final generalized Cauchy point
    Projg = norm of the final projected gradient
    F     = final function value
    
               * * *
    
       N    Tit     Tnf  Tnint  Skip  Nact     Projg        F
        7      9     11     13     0     0   2.718D+04   6.559D+04
      F =   65586.974331744845     
    
    CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH             
    >>> Initial Objective 2.0144e+05
    >>> Final Objective        65587
    >>> Percent change        -67.44%
    2024-06-02 19:26:47.550376 Running candidate task 2/2
    workspace listening on local host. Remote connections prohibited.
    2024-06-02 19:26:48.241999 Calculating initial obj
    2024-06-02 19:26:48.242064 Starting physical parameter optimization
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.381922 d=              0
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.298647 d=              0
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.229711 d=              0
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.188198 d=              0
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.786276 d=              0
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102783 d=              0
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.381922 d=              0
    2024-06-02 19:26:48.243215 2024-06-02 19:26:48.243202 Generating 2 objectives
    2024-06-02 19:26:48.243674 2024-06-02 19:26:48.243658 Starting 30 tasks
    2024-06-02 19:26:52.265172 2024-06-02 19:26:52.265155 Calculating 30 tasks
      0000 | X2=  0.006199739 |g|=    0.1915674
      0001 | X2=     201437.2 |g|=      2131355
    >>> X2=     201437.2 |g|=     2131355
    2024-06-02 19:26:52.265867 2024-06-02 19:26:52.265858 Done. 30 tasks complete
    RUNNING THE L-BFGS-B CODE
    
               * * *
    
    Machine precision = 2.220D-16
     N =            7     M =           10
    
    At X0         0 variables are exactly at the bounds
    
    At iterate    0    f=  2.01437D+05    |proj g|=  1.52718D+06
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.381921 d=  -9.039213e-07
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.298646 d=  -8.494508e-07
    Setting (0, 'l', 'b17', 0)   from        1.229711 to         1.22971 d=  -8.043597e-07
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.195441 d=    0.007243539
    Setting (0, 'l', 'b70', 0)   from        1.786276 to         1.83186 d=      0.0455842
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102782 d=  -7.213354e-07
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        2.380856 d=      0.9989342
    2024-06-02 19:26:52.266491 2024-06-02 19:26:52.266483 Generating 2 objectives
    2024-06-02 19:26:52.266740 2024-06-02 19:26:52.266734 Starting 30 tasks
    2024-06-02 19:26:57.281407 2024-06-02 19:26:57.281388 Calculating 30 tasks
      0000 | X2=     2.191041 |g|=     4.455478
      0001 | X2=      8633944 |g|= 2.739841e+07
    >>> X2=      8633946 |g|=2.739841e+07
    2024-06-02 19:26:57.282138 2024-06-02 19:26:57.282130 Done. 30 tasks complete
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.381922 d=  -6.936445e-08
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.298647 d=  -6.518453e-08
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.229711 d=  -6.172437e-08
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.188754 d=   0.0005558494
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.789774 d=    0.003498007
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102783 d=  -5.535331e-08
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.458577 d=     0.07665549
    2024-06-02 19:26:57.282489 2024-06-02 19:26:57.282481 Generating 2 objectives
    2024-06-02 19:26:57.282733 2024-06-02 19:26:57.282726 Starting 30 tasks
    2024-06-02 19:27:02.310757 2024-06-02 19:27:02.310738 Calculating 30 tasks
      0000 | X2=   0.03403521 |g|=     0.933406
      0001 | X2=     142778.9 |g|=      2143720
    >>> X2=     142778.9 |g|=     2143720
    2024-06-02 19:27:02.311540 2024-06-02 19:27:02.311528 Done. 30 tasks complete
    
    At iterate    1    f=  1.42779D+05    |proj g|=  1.21476D+06
    Setting (0, 'l', 'b4', 0)    from        1.381922 to         1.34528 d=    -0.03664141
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.318481 d=     0.01983407
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.240725 d=     0.01101365
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.219009 d=      0.0308113
    Setting (0, 'l', 'b70', 0)   from        1.786276 to         1.79802 d=      0.0117447
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102538 d=  -0.0002447489
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.496818 d=      0.1148964
    2024-06-02 19:27:02.312000 2024-06-02 19:27:02.311991 Generating 2 objectives
    2024-06-02 19:27:02.312228 2024-06-02 19:27:02.312221 Starting 30 tasks
    2024-06-02 19:27:06.312632 2024-06-02 19:27:06.312617 Calculating 30 tasks
      0000 | X2=   0.06321656 |g|=     1.105779
      0001 | X2=     129021.1 |g|=      3407824
    >>> X2=     129021.1 |g|=     3407824
    2024-06-02 19:27:06.313570 2024-06-02 19:27:06.313565 Done. 30 tasks complete
    
    At iterate    2    f=  1.29021D+05    |proj g|=  6.00914D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.344783 d=    -0.03713897
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.311205 d=     0.01255816
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.242294 d=     0.01258272
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.202336 d=     0.01413863
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.799835 d=     0.01355879
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.104195 d=    0.001412117
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.480683 d=     0.09876135
    2024-06-02 19:27:06.313952 2024-06-02 19:27:06.313944 Generating 2 objectives
    2024-06-02 19:27:06.314227 2024-06-02 19:27:06.314222 Starting 30 tasks
    2024-06-02 19:27:09.320597 2024-06-02 19:27:09.320583 Calculating 30 tasks
      0000 | X2=   0.04072391 |g|=    0.2894206
      0001 | X2=     91126.56 |g|=      1048334
    >>> X2=      91126.6 |g|=     1048334
    2024-06-02 19:27:09.321490 2024-06-02 19:27:09.321484 Done. 30 tasks complete
    
    At iterate    3    f=  9.11266D+04    |proj g|=  3.68215D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.334101 d=    -0.04782061
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.308033 d=    0.009386385
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.246155 d=     0.01644362
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.198106 d=    0.009908732
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.804497 d=     0.01822149
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.102794 d=   1.143844e-05
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.477414 d=     0.09549244
    2024-06-02 19:27:09.321868 2024-06-02 19:27:09.321860 Generating 2 objectives
    2024-06-02 19:27:09.322133 2024-06-02 19:27:09.322128 Starting 30 tasks
    2024-06-02 19:27:13.357968 2024-06-02 19:27:13.357953 Calculating 30 tasks
      0000 | X2=   0.04073447 |g|=    0.3875965
      0001 | X2=     81438.74 |g|=     683347.1
    >>> X2=     81438.78 |g|=    683347.2
    2024-06-02 19:27:13.358900 2024-06-02 19:27:13.358894 Done. 30 tasks complete
    
    At iterate    4    f=  8.14388D+04    |proj g|=  2.71606D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.313809 d=    -0.06811315
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.302319 d=    0.003672524
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.250245 d=     0.02053391
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.196968 d=    0.008770196
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.814084 d=     0.02780814
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.097224 d=   -0.005558942
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.475884 d=     0.09396177
    2024-06-02 19:27:13.359282 2024-06-02 19:27:13.359275 Generating 2 objectives
    2024-06-02 19:27:13.359556 2024-06-02 19:27:13.359551 Starting 30 tasks
    2024-06-02 19:27:16.382482 2024-06-02 19:27:16.382465 Calculating 30 tasks
      0000 | X2=   0.01546186 |g|=    0.1913408
      0001 | X2=     71978.99 |g|=     613166.8
    >>> X2=        71979 |g|=    613166.8
    2024-06-02 19:27:16.383107 2024-06-02 19:27:16.383099 Done. 30 tasks complete
    
    At iterate    5    f=  7.19790D+04    |proj g|=  4.04255D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to         1.29802 d=    -0.08390154
    Setting (0, 'l', 'b6', 0)    from        1.298647 to         1.29622 d=   -0.002426962
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.241563 d=     0.01185217
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.198028 d=    0.009829893
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.825393 d=     0.03911703
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.091226 d=    -0.01155689
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.472374 d=      0.0904523
    2024-06-02 19:27:16.383492 2024-06-02 19:27:16.383484 Generating 2 objectives
    2024-06-02 19:27:16.383736 2024-06-02 19:27:16.383730 Starting 30 tasks
    2024-06-02 19:27:20.411746 2024-06-02 19:27:20.411731 Calculating 30 tasks
      0000 | X2=   0.01392054 |g|=    0.3636936
      0001 | X2=     66325.96 |g|=     329329.5
    >>> X2=     66325.98 |g|=    329329.6
    2024-06-02 19:27:20.412526 2024-06-02 19:27:20.412516 Done. 30 tasks complete
    
    At iterate    6    f=  6.63260D+04    |proj g|=  2.88255D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.295411 d=     -0.0865107
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.296168 d=   -0.002478367
    Setting (0, 'l', 'b17', 0)   from        1.229711 to         1.24009 d=      0.0103784
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.200627 d=     0.01242952
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.831228 d=     0.04495185
    Setting (0, 'l', 'b85', 0)   from        1.102783 to         1.08889 d=    -0.01389342
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.474668 d=     0.09274669
    2024-06-02 19:27:20.412924 2024-06-02 19:27:20.412916 Generating 2 objectives
    2024-06-02 19:27:20.413140 2024-06-02 19:27:20.413135 Starting 30 tasks
    2024-06-02 19:27:25.465338 2024-06-02 19:27:25.465325 Calculating 30 tasks
      0000 | X2=   0.01465877 |g|=    0.3703143
      0001 | X2=     65711.61 |g|=     166508.1
    >>> X2=     65711.63 |g|=    166507.9
    2024-06-02 19:27:25.466226 2024-06-02 19:27:25.466221 Done. 30 tasks complete
    
    At iterate    7    f=  6.57116D+04    |proj g|=  1.43515D+05
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.295551 d=      -0.086371
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.295682 d=   -0.002965275
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.239183 d=    0.009471377
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.200692 d=     0.01249471
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.831539 d=     0.04526368
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.089854 d=    -0.01292852
    Setting (0, 'l', 'B91', 0)   from        1.381922 to        1.474238 d=     0.09231586
    2024-06-02 19:27:25.466605 2024-06-02 19:27:25.466598 Generating 2 objectives
    2024-06-02 19:27:25.466872 2024-06-02 19:27:25.466868 Starting 30 tasks
    2024-06-02 19:27:30.483839 2024-06-02 19:27:30.483824 Calculating 30 tasks
      0000 | X2=   0.01472726 |g|=    0.3829256
      0001 | X2=     65612.06 |g|=     71994.11
    >>> X2=     65612.08 |g|=    71994.11
    2024-06-02 19:27:30.484783 2024-06-02 19:27:30.484773 Done. 30 tasks complete
    
    At iterate    8    f=  6.56121D+04    |proj g|=  3.16585D+04
    Setting (0, 'l', 'b4', 0)    from        1.381922 to        1.295963 d=    -0.08595864
    Setting (0, 'l', 'b6', 0)    from        1.298647 to        1.295903 d=   -0.002744068
    Setting (0, 'l', 'b17', 0)   from        1.229711 to        1.239008 d=    0.009297198
    Setting (0, 'l', 'b21', 0)   from        1.188198 to        1.200526 d=     0.01232822
    Setting (0, 'l', 'b70', 0)   from        1.786276 to        1.831825 d=     0.04554951
    Setting (0, 'l', 'b85', 0)   from        1.102783 to        1.089911 d=    -0.01287199
    Setting (0, 'l', 'B91', 0)   from        1.381922 to         1.47421 d=     0.09228864
    2024-06-02 19:27:30.485177 2024-06-02 19:27:30.485170 Generating 2 objectives
    2024-06-02 19:27:30.485453 2024-06-02 19:27:30.485448 Starting 30 tasks
    2024-06-02 19:27:35.485947 2024-06-02 19:27:35.485931 Calculating 30 tasks
      0000 | X2=   0.01468777 |g|=    0.3795553
      0001 | X2=     65586.96 |g|=     41687.39
    >>> X2=     65586.97 |g|=    41687.31
    2024-06-02 19:27:35.486852 2024-06-02 19:27:35.486845 Done. 30 tasks complete
    
    At iterate    9    f=  6.55870D+04    |proj g|=  2.71753D+04
    
               * * *
    
    Tit   = total number of iterations
    Tnf   = total number of function evaluations
    Tnint = total number of segments explored during Cauchy searches
    Skip  = number of BFGS updates skipped
    Nact  = number of active bounds at final generalized Cauchy point
    Projg = norm of the final projected gradient
    F     = final function value
    
               * * *
    
       N    Tit     Tnf  Tnint  Skip  Nact     Projg        F
        7      9     11     13     0     0   2.718D+04   6.559D+04
      F =   65586.974331269666     
    
    CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH             
    >>> Initial Objective 2.0144e+05
    >>> Final Objective        65587
    >>> Percent change        -67.44%
    The unfiltered results of the candidate scan N=2 total=2 oper=1:
     Initial objectives:  X= 203451.60610 P= 201437.23376 C= 2014.37234
    Cnd.    1/2 N=      1 K= Y dP=  -135850.25943 dC=       40.63667 d(P+C)=  -135809.62275 d%=    -66.753% b4      [#6H1X3:1]-[#6X3:2]
    Cnd.    2/2 N=      1 K= Y dP=  -135850.25943 dC=       40.63667 d(P+C)=  -135809.62275 d%=    -66.753% b4      [#6H0X3:1]-[#6X3:2]
    
                                                                                                                                       
    Nanostep 1: The filtered results of the candidate scan N=2 total=2 oper=1:
    ->    1 Cnd.    2/2 N=      1 K= Y dP=  -135850.25943 dC=       40.63667 d(P+C)=  -135809.62275 d%=    -66.753% b4      [#6H0X3:1]-[#6X3:2]
          2 Cnd.    1/2 N=      1 K= Y dP=  -135850.25943 dC=       40.63667 d(P+C)=  -135809.62275 d%=    -66.753% b4      [#6H1X3:1]-[#6X3:2]
    Performing 1 operations
    There are 1 nodes returned
    Operations per parameter for this micro:
    Counter({'b4': 1})
    Micro total: 1 should be 1
    Operations per parameter for this macro:
    Counter({'b4': 1})
    Macro total: 1 should be 1
    Only one modification, keeping result
    2024-06-02 19:27:36.688497 Macro objective:         65587 C=      2055.01 DX=      -135810
    Model:
    Tree:
     0   0 Bonds  
     1   0  b4   [#6X3:1]-[#6X3:2] k: [540.3345953498] l: [1.2959631260483113]
     2   0   B91  [#6H0X3:1]-[#6X3:2] k: [540.3345953498] l: [1.4742104085148617]
     1   0  b6   [#6X3:1]=[#6X3:2] k: [898.589948525] l: [1.2959027229929716]
     1   0  b17  [#6X3:1]-[#8X2:2] k: [598.9859275918] l: [1.2390084195590223]
     1   0  b21  [#6:1]=[#8X1+0,#8X2+1:2] k: [1527.019744047] l: [1.2005259371996082]
     1   0  b70  [#6:1]-[#17:2] k: [368.4266150848] l: [1.8318253015203048]
     1   0  b85  [#6X3:1]-[#1:2] k: [775.3853383846] l: [1.089910995441042]
    Tree:
     0   1 Angles  
    Tree:
     0   2 Torsions  
    Tree:
     0   3 OutOfPlanes  
    Tree:
    Tree:
     0   5 vdW   
    Tree:
    There were 1 successful operations
    2024-06-02 19:27:36.696128 Visited {'b4', 'B91'}
    Restarting optimization search
    Targets for this macro step 1:
    1 (0, 0, 0) b4
    2 (0, 0, 0) B91
    3 (0, 0, 0) b6
    4 (0, 0, 0) b17
    5 (0, 0, 0) b21
    6 (0, 0, 0) b70
    7 (0, 0, 0) b85
    N Targets: 7
    Step tracker for current macro step 1
    ((0, 0, 0), 'B91') 1
    
    *******************
     iteration=   1 macro=  1/2 micro=2 X=    67642 P=    65587 C=     2055 models=0:Bonds
    *******************
    
    2024-06-02 19:27:36.763149 Saving checkpoint to chk.cst.p
    2024-06-02 19:27:36.791341 Collecting SMARTS for b4 and setting to depth=0
    
     == iteration=   2 macro=  1/2 micro=  1/2 operation=1 cluster=b4   N= 1 overlap=[0] bits=1->1 depth=0->0 branch=0->0
    
    Attempting to split b4:
    S0: [#6X3:1]-[#6X3:2]
    Matched N=1
    000001 (0, (2, 3))              [#6H1X3x2r5A+0:2]@;-[#6H1X3x2r5A+0:3]
    
    Skipping b4 since all graphs are the same
    2024-06-02 19:27:36.792939 Collecting SMARTS for B91 and setting to depth=0
    
     == iteration=   3 macro=  1/2 micro=  2/2 operation=1 cluster=B91  N= 1 overlap=[0] bits=1->1 depth=0->0 branch=0->0
    
    Attempting to split B91:
    S0: [#6H0X3:1]-[#6X3:2]
    Matched N=1
    000001 (0, (4, 5))              [#6H0X3x2r5A+0:4]!@;-[#6H0X3x0!rA+0:5]
    
    Skipping B91 since all graphs are the same
    2024-06-02 19:27:36.794294 Scanning done.
    2024-06-02 19:27:36.794304
    
    
    Generating SMARTS on 0
    Model:
    Tree:
     0   0 Bonds  
     1   0  b4   [#6X3:1]-[#6X3:2] k: [540.3345953498] l: [1.2959631260483113]
     2   0   B91  [#6H0X3:1]-[#6X3:2] k: [540.3345953498] l: [1.4742104085148617]
     1   0  b6   [#6X3:1]=[#6X3:2] k: [898.589948525] l: [1.2959027229929716]
     1   0  b17  [#6X3:1]-[#8X2:2] k: [598.9859275918] l: [1.2390084195590223]
     1   0  b21  [#6:1]=[#8X1+0,#8X2+1:2] k: [1527.019744047] l: [1.2005259371996082]
     1   0  b70  [#6:1]-[#17:2] k: [368.4266150848] l: [1.8318253015203048]
     1   0  b85  [#6X3:1]-[#1:2] k: [775.3853383846] l: [1.089910995441042]
    Tree:
     0   1 Angles  
    Tree:
     0   2 Torsions  
    Tree:
     0   3 OutOfPlanes  
    Tree:
    Tree:
     0   5 vdW   
    Tree:
    Scoring and filtering 0 candidates for operation=1
    Tier 0: Scoring and filtering 0 candidates for operation=1
    Tier 0: Accepting all candidates so we skip
    Scanning 0 candidates for operation=1
    2024-06-02 19:27:36.811804 Generated 0 x 1 = 0 candidate evalulation tasks
    2024-06-02 19:27:36.811826 Dispatching candidate tasks= 0 in serial
    The unfiltered results of the candidate scan N=0 total=0 oper=1:
     Initial objectives:  X= 67641.98334 P= 65586.97433 C= 2055.00901
    
                                                                                                                                       
    Nanostep 1: The filtered results of the candidate scan N=0 total=0 oper=1:
    There were 0 successful operations
    2024-06-02 19:27:36.811902 Visited set()
    Targets for this macro step 2:
    1 (0, 0, 0) b4
    2 (0, 0, 0) B91
    3 (0, 0, 0) b6
    4 (0, 0, 0) b17
    5 (0, 0, 0) b21
    6 (0, 0, 0) b70
    7 (0, 0, 0) b85
    N Targets: 7
    Step tracker for current macro step 2
    ((0, 0, 0), 'B91') 2
    ((0, 0, 0), 'b4') 2
    
    *******************
     iteration=   3 macro=  2/2 micro=2 X=    67642 P=    65587 C=     2055 models=0:Bonds
    *******************
    
    2024-06-02 19:27:36.850726 Saving checkpoint to chk.cst.p
    2024-06-02 19:27:36.878860 Collecting SMARTS for b4 and setting to depth=0
    
     == iteration=   4 macro=  2/2 micro=  1/2 operation=-1 cluster=b4   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 19:27:36.878954 Collecting SMARTS for B91 and setting to depth=0
    
     == iteration=   5 macro=  2/2 micro=  2/2 operation=-1 cluster=B91  N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 19:27:36.878984 Scanning done.
    2024-06-02 19:27:36.878991
    
    
    Generating SMARTS on 1
    Model:
    Tree:
     0   0 Bonds  
     1   0  b4   [#6X3:1]-[#6X3:2] k: [540.3345953498] l: [1.2959631260483113]
     2   0   B91  [#6H0X3:1]-[#6X3:2] k: [540.3345953498] l: [1.4742104085148617]
     1   0  b6   [#6X3:1]=[#6X3:2] k: [898.589948525] l: [1.2959027229929716]
     1   0  b17  [#6X3:1]-[#8X2:2] k: [598.9859275918] l: [1.2390084195590223]
     1   0  b21  [#6:1]=[#8X1+0,#8X2+1:2] k: [1527.019744047] l: [1.2005259371996082]
     1   0  b70  [#6:1]-[#17:2] k: [368.4266150848] l: [1.8318253015203048]
     1   0  b85  [#6X3:1]-[#1:2] k: [775.3853383846] l: [1.089910995441042]
    Tree:
     0   1 Angles  
    Tree:
     0   2 Torsions  
    Tree:
     0   3 OutOfPlanes  
    Tree:
    Tree:
     0   5 vdW   
    Tree:
    Scoring and filtering 1 candidates for operation=-1
    Tier 0: Scoring and filtering 1 candidates for operation=-1
    Tier 0: Accepting all candidates so we skip
    Scanning 1 candidates for operation=-1
    2024-06-02 19:27:36.896660 Generated 1 x 1 = 1 candidate evalulation tasks
    2024-06-02 19:27:36.896682 Dispatching candidate tasks= 1 in serial
    2024-06-02 19:27:36.896696 Running candidate task 1/1
    workspace listening on local host. Remote connections prohibited.
    2024-06-02 19:27:37.568336 Calculating initial obj
    2024-06-02 19:27:37.568376 Starting physical parameter optimization
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.295963 d=              0
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.295903 d=              0
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.239008 d=              0
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.200526 d=              0
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.831825 d=              0
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.089911 d=              0
    2024-06-02 19:27:37.569330 2024-06-02 19:27:37.569319 Generating 2 objectives
    2024-06-02 19:27:37.569799 2024-06-02 19:27:37.569791 Starting 26 tasks
    2024-06-02 19:27:42.570353 2024-06-02 19:27:42.570338 Calculating 26 tasks
      0000 | X2=   0.02037902 |g|=     3.017232
      0001 | X2=     390173.9 |g|=      5194615
    >>> X2=       390174 |g|=     5194616
    2024-06-02 19:27:42.571204 2024-06-02 19:27:42.571198 Done. 26 tasks complete
    RUNNING THE L-BFGS-B CODE
    
               * * *
    
    Machine precision = 2.220D-16
     N =            6     M =           10
    
    At X0         0 variables are exactly at the bounds
    
    At iterate    0    f=  3.90174D+05    |proj g|=  3.65705D+06
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        2.295962 d=      0.9999991
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.295902 d=   -3.54357e-07
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.239008 d=  -3.387996e-07
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.200526 d=  -3.282768e-07
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.831825 d=  -5.009019e-07
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.091224 d=    0.001313048
    2024-06-02 19:27:42.571810 2024-06-02 19:27:42.571802 Generating 2 objectives
    2024-06-02 19:27:42.572046 2024-06-02 19:27:42.572041 Starting 26 tasks
    2024-06-02 19:27:47.590312 2024-06-02 19:27:47.590298 Calculating 26 tasks
      0000 | X2=     2.744048 |g|=     6.691107
      0001 | X2= 1.715364e+07 |g|= 4.743588e+07
    >>> X2= 1.715364e+07 |g|=4.743589e+07
    2024-06-02 19:27:47.591101 2024-06-02 19:27:47.591095 Done. 26 tasks complete
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.385507 d=     0.08954358
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.295903 d=  -3.173042e-08
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.239008 d=  -3.033735e-08
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.200526 d=   -2.93951e-08
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.831825 d=  -4.485258e-08
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.090029 d=   0.0001175752
    2024-06-02 19:27:47.591415 2024-06-02 19:27:47.591408 Generating 2 objectives
    2024-06-02 19:27:47.591642 2024-06-02 19:27:47.591637 Starting 26 tasks
    2024-06-02 19:27:51.594475 2024-06-02 19:27:51.594460 Calculating 26 tasks
      0000 | X2=  0.004698566 |g|=    0.1295538
      0001 | X2=     226440.9 |g|=      2256879
    >>> X2=     226440.9 |g|=     2256879
    2024-06-02 19:27:51.595304 2024-06-02 19:27:51.595298 Done. 26 tasks complete
    
    At iterate    1    f=  2.26441D+05    |proj g|=  1.60712D+06
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.385148 d=     0.08918453
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.294575 d=   -0.001327634
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.231367 d=   -0.007641334
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.178207 d=    -0.02231884
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.828391 d=   -0.003434781
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.114245 d=     0.02433418
    2024-06-02 19:27:51.595655 2024-06-02 19:27:51.595648 Generating 2 objectives
    2024-06-02 19:27:51.595901 2024-06-02 19:27:51.595895 Starting 26 tasks
    2024-06-02 19:27:55.619240 2024-06-02 19:27:55.619224 Calculating 26 tasks
      0000 | X2=   0.03152384 |g|=    0.1943937
      0001 | X2=     230653.2 |g|=      2873942
    >>> X2=     230653.2 |g|=     2873942
    2024-06-02 19:27:55.619864 2024-06-02 19:27:55.619855 Done. 26 tasks complete
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.385337 d=      0.0893734
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.295273 d=  -0.0006292975
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.235387 d=   -0.003621907
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.189947 d=    -0.01057885
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.830197 d=   -0.001628065
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.101507 d=     0.01159592
    2024-06-02 19:27:55.620211 2024-06-02 19:27:55.620203 Generating 2 objectives
    2024-06-02 19:27:55.620429 2024-06-02 19:27:55.620424 Starting 26 tasks
    2024-06-02 19:27:58.643184 2024-06-02 19:27:58.643170 Calculating 26 tasks
      0000 | X2=  0.006431124 |g|=    0.1703743
      0001 | X2=     208251.2 |g|=     668055.5
    >>> X2=     208251.2 |g|=    668055.6
    2024-06-02 19:27:58.644005 2024-06-02 19:27:58.643999 Done. 26 tasks complete
    
    At iterate    2    f=  2.08251D+05    |proj g|=  4.51798D+05
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.387266 d=     0.09130304
    Setting (0, 'l', 'b6', 0)    from        1.295903 to         1.29819 d=    0.002286926
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.233259 d=   -0.005749202
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.190585 d=    -0.00994134
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.827963 d=   -0.003862074
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.103092 d=     0.01318143
    2024-06-02 19:27:58.644366 2024-06-02 19:27:58.644358 Generating 2 objectives
    2024-06-02 19:27:58.644626 2024-06-02 19:27:58.644620 Starting 26 tasks
    2024-06-02 19:28:02.673210 2024-06-02 19:28:02.673196 Calculating 26 tasks
      0000 | X2=  0.006250314 |g|=    0.1769401
      0001 | X2=     205948.6 |g|=     400558.2
    >>> X2=     205948.6 |g|=    400558.3
    2024-06-02 19:28:02.674034 2024-06-02 19:28:02.674028 Done. 26 tasks complete
    
    At iterate    3    f=  2.05949D+05    |proj g|=  2.04644D+05
    Setting (0, 'l', 'b4', 0)    from        1.295963 to          1.3896 d=     0.09363656
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.301496 d=    0.005592847
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.230654 d=    -0.00835454
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.191728 d=   -0.008798085
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.821041 d=    -0.01078462
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.105606 d=     0.01569525
    2024-06-02 19:28:02.674387 2024-06-02 19:28:02.674380 Generating 2 objectives
    2024-06-02 19:28:02.674654 2024-06-02 19:28:02.674647 Starting 26 tasks
    2024-06-02 19:28:06.715022 2024-06-02 19:28:06.715006 Calculating 26 tasks
      0000 | X2=  0.006111424 |g|=    0.1857443
      0001 | X2=     203922.7 |g|=     250168.6
    >>> X2=     203922.8 |g|=    250168.6
    2024-06-02 19:28:06.715850 2024-06-02 19:28:06.715844 Done. 26 tasks complete
    
    At iterate    4    f=  2.03923D+05    |proj g|=  5.99700D+04
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.389198 d=     0.09323441
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.302949 d=    0.007045885
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.230086 d=   -0.008922628
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.190695 d=   -0.009830959
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.814125 d=    -0.01769988
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.105293 d=     0.01538227
    2024-06-02 19:28:06.716218 2024-06-02 19:28:06.716212 Generating 2 objectives
    2024-06-02 19:28:06.716466 2024-06-02 19:28:06.716460 Starting 26 tasks
    2024-06-02 19:28:10.750233 2024-06-02 19:28:10.750219 Calculating 26 tasks
      0000 | X2=  0.005951266 |g|=    0.1853704
      0001 | X2=     202898.9 |g|=     296696.9
    >>> X2=     202898.9 |g|=    296696.8
    2024-06-02 19:28:10.751057 2024-06-02 19:28:10.751051 Done. 26 tasks complete
    
    At iterate    5    f=  2.02899D+05    |proj g|=  9.89974D+04
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.387727 d=     0.09176425
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.302019 d=    0.006115811
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.229941 d=   -0.009067249
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.190494 d=     -0.0100324
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.799148 d=    -0.03267737
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.105625 d=     0.01571413
    2024-06-02 19:28:10.751429 2024-06-02 19:28:10.751422 Generating 2 objectives
    2024-06-02 19:28:10.751674 2024-06-02 19:28:10.751668 Starting 26 tasks
    2024-06-02 19:28:14.752504 2024-06-02 19:28:14.752489 Calculating 26 tasks
      0000 | X2=  0.006117953 |g|=    0.1938061
      0001 | X2=     201528.9 |g|=     246285.5
    >>> X2=     201528.9 |g|=    246285.4
    2024-06-02 19:28:14.753330 2024-06-02 19:28:14.753324 Done. 26 tasks complete
    
    At iterate    6    f=  2.01529D+05    |proj g|=  8.66989D+04
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.384522 d=     0.08855868
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.300122 d=    0.004219161
    Setting (0, 'l', 'b17', 0)   from        1.239008 to        1.230359 d=   -0.008649701
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.188459 d=    -0.01206694
    Setting (0, 'l', 'b70', 0)   from        1.831825 to        1.793898 d=    -0.03792739
    Setting (0, 'l', 'b85', 0)   from        1.089911 to         1.10302 d=     0.01310884
    2024-06-02 19:28:14.753702 2024-06-02 19:28:14.753695 Generating 2 objectives
    2024-06-02 19:28:14.753948 2024-06-02 19:28:14.753943 Starting 26 tasks
    2024-06-02 19:28:18.764284 2024-06-02 19:28:18.764266 Calculating 26 tasks
      0000 | X2=  0.006083782 |g|=    0.1876833
      0001 | X2=     201159.2 |g|=     162245.1
    >>> X2=     201159.2 |g|=      162245
    2024-06-02 19:28:18.765183 2024-06-02 19:28:18.765174 Done. 26 tasks complete
    
    At iterate    7    f=  2.01159D+05    |proj g|=  8.83989D+04
    Setting (0, 'l', 'b4', 0)    from        1.295963 to        1.384738 d=     0.08877498
    Setting (0, 'l', 'b6', 0)    from        1.295903 to        1.299488 d=    0.003585422
    Setting (0, 'l', 'b17', 0)   from        1.239008 to         1.23046 d=   -0.008548149
    Setting (0, 'l', 'b21', 0)   from        1.200526 to        1.189077 d=    -0.01144935
    Setting (0, 'l', 'b70', 0)   from        1.831825 to         1.79511 d=    -0.03671554
    Setting (0, 'l', 'b85', 0)   from        1.089911 to        1.103434 d=     0.01352321
    2024-06-02 19:28:18.765551 2024-06-02 19:28:18.765544 Generating 2 objectives
    2024-06-02 19:28:18.765781 2024-06-02 19:28:18.765777 Starting 26 tasks
    2024-06-02 19:28:22.792019 2024-06-02 19:28:22.792004 Calculating 26 tasks
      0000 | X2=  0.006159883 |g|=    0.1889322
      0001 | X2=     201068.1 |g|=     15443.45
    >>> X2=     201068.1 |g|=    15443.36
    2024-06-02 19:28:22.792862 2024-06-02 19:28:22.792854 Done. 26 tasks complete
    
    At iterate    8    f=  2.01068D+05    |proj g|=  3.39776D+03
    
               * * *
    
    Tit   = total number of iterations
    Tnf   = total number of function evaluations
    Tnint = total number of segments explored during Cauchy searches
    Skip  = number of BFGS updates skipped
    Nact  = number of active bounds at final generalized Cauchy point
    Projg = norm of the final projected gradient
    F     = final function value
    
               * * *
    
       N    Tit     Tnf  Tnint  Skip  Nact     Projg        F
        6      8     11     12     0     0   3.398D+03   2.011D+05
      F =   201068.14154365813     
    
    CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH             
    >>> Initial Objective 3.9017e+05
    >>> Final Objective   2.0107e+05
    >>> Percent change       -48.467%
    The unfiltered results of the candidate scan N=1 total=1 oper=-1:
     Initial objectives:  X= 67641.98334 P= 65586.97433 C= 2055.00901
    Cnd.    1/1 N=      2 K= N dP=   135481.16721 dC=      -40.63667 d(P+C)=   135440.53054 d%=    200.231% b4      [#6H0X3:1]-[#6X3:2]
    
                                                                                                                                       
    Nanostep 1: The filtered results of the candidate scan N=0 total=1 oper=-1:
    There were 0 successful operations
    2024-06-02 19:28:23.792058 Visited {'B91'}
    Nothing found. Done.
    Start time: 2024-06-02 19:25:07.319094
    End   time: 2024-06-02 19:28:23.820893
    Model:
    Tree:
     0   0 Bonds  
     1   0  b4   [#6X3:1]-[#6X3:2] k: [540.3345953498] l: [1.2959631260483113]
     2   0   B91  [#6H0X3:1]-[#6X3:2] k: [540.3345953498] l: [1.4742104085148617]
     1   0  b6   [#6X3:1]=[#6X3:2] k: [898.589948525] l: [1.2959027229929716]
     1   0  b17  [#6X3:1]-[#8X2:2] k: [598.9859275918] l: [1.2390084195590223]
     1   0  b21  [#6:1]=[#8X1+0,#8X2+1:2] k: [1527.019744047] l: [1.2005259371996082]
     1   0  b70  [#6:1]-[#17:2] k: [368.4266150848] l: [1.8318253015203048]
     1   0  b85  [#6X3:1]-[#1:2] k: [775.3853383846] l: [1.089910995441042]
    Tree:
     0   1 Angles  
    Tree:
     0   2 Torsions  
    Tree:
     0   3 OutOfPlanes  
    Tree:
    Tree:
     0   5 vdW   
    Tree:
    Modified parameters:
    (0, 'k', 'B91', 0)   + New:      540.335
    (0, 'l', 'b4', 0)    | New:      1.29596 Ref       1.4662 Diff    -0.170236
    (0, 'l', 'b6', 0)    | New:       1.2959 Ref      1.38236 Diff    -0.086459
    (0, 'l', 'b17', 0)   | New:      1.23901 Ref      1.35775 Diff    -0.118738
    (0, 'l', 'b21', 0)   | New:      1.20053 Ref      1.22167 Diff   -0.0211427
    (0, 'l', 'b70', 0)   | New:      1.83183 Ref      1.72222 Diff      0.10961
    (0, 'l', 'b85', 0)   | New:      1.08991 Ref      1.08182 Diff   0.00808732
    (0, 'l', 'B91', 0)   + New:      1.47421
    Initial objectives:
    Total=       1012706.6 Physical       1010692.2 Chemical       2014.3723
    Final objectives:
    Total=       67641.983 Physical       65586.974 Chemical        2055.009
    Differences:
    Total=         -93.32% Physical         -93.51% Chemical           2.02%

Most of the output is intermediate and diagnostics. However, at the bottom we
see that the objective improved 93.32%. This shows that a single bond was found
which dropped the objective further. In particular, it specialized the somewhat
generic `[#6X3:1]-[#6X3:2]` b4 bond with `[#6H0X3:1]-[#6X3:2]`. The new bond
increased the bond length to 1.474 A, while the original b4 parameter was
originally at 1.466 A but dropped to 1.296 A in the final result. One thing to
examine is the objective change after the initial fit, but before any bond
parameters were added. From the output, we see that the total objective after
the first fit is 203451.60610, or reduced by 80% of the objective produced by
OpenFF 2.1.0. The improvement due to adding the B91 was 65587, or 67% from the
fit objective (and 93% of OpenFF 2.1.0). Most of the improvement was from
fitting the force. It is important to note that the geometry was not ruined
either; the initial geometry objective was 0.1174278 A^2 in 2.1.0. The
objective is the Cartesian residual sum of squares (RSS) in Angstroms, and so
this objective indicates a relatively good MM-minimized geometry. Fitting lead
to a geometry objective of 0.006199739 A^2, and after adding the new bond
parameter the geometry objective increased to 0.01468777 A^2, of which both are
smaller than the original objective.


As a final note, the setup and results shown here is still undergoing
development. In the future, the setup should be easier and shorter, and the
output should be cleaner and more informative.


