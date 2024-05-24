Force field fitting
===================

The most interesting application, and largely the application that inspired
this package, was force field fitting with automatic chemical perception. What
this means is that we want the code to fit the force field parameters, but also
find the best parameters that give the best fit as well.

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
>>> # configs.processors = 1
>>> xyz_positions = """11
>>> Energy =  -802.9904298568 (+0.195 kcal/mol)
>>>   C -1.44819400 -0.84940800  0.16848900
>>>   C -1.59401300  0.50318700 -0.01678100
>>>   C -0.27397600  1.02622600 -0.13503500
>>>   C  0.58064400 -0.04716400 -0.01303100
>>>   C  2.03461200 -0.06860900 -0.05925200
>>>   O  2.72809700  0.90108700 -0.21909900
>>>  Cl  2.76214600 -1.70734100  0.14655600
>>>   O -0.13897300 -1.20044600  0.17351800
>>>   H -2.15226800 -1.65836100  0.30609000
>>>   H -2.52743000  1.04809900 -0.06180000
>>>   H  0.02935200  2.05273200 -0.28965800
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
>>>             indices = list(range(N))
>>>         assert N == int(lines[0].split()[0])
>>>         for chunk in arrays.batched(lines, N+2):
>>>             if chunk and chunk[0]:
>>>                 sel = {}
>>>                 for i, line in enumerate(chunk, -2):
>>>                     if i >= 0:
>>>                         sym, x, y, z = line.split()[:4]
>>>                         sel[indices[i]] = list(map(float, (x, y, z)))
>>> 
>>>                 gdc = assignments.graph_db_column()
>>>                 gdc.selections.update(sel)
>>>                 gdr.columns[s] = gdc
>>>                 s += 1
>>>     return gdr
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
>>> def new_gdb(f: Dict[str, List[str]]) -> assignments.graph_db:
>>>     gcd = codec_rdkit.graph_codec_rdkit()
>>>     gdb = assignments.graph_db()
>>> 
>>>     ne = 0
>>>     for smi, fn_dict in f.items():
>>> 
>>>         g = gcd.smiles_decode(smi)
>>>         gid = assignments.graph_db_add_graph(gdb, g)
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
>>>             indices=dict(sorted([(j, x) for j, x in enumerate(g.nodes)], key=lambda x: x[1]))
>>>             gdr = load_xyz([fn], indices=list(indices))
>>>             gdg.rows[0] = gdr
>>>             gde.tables[tid] = gdt
>>>             tid = assignments.GRADIENTS
>>>             if tid in rdata:
>>>                 gdt = assignments.graph_db_table(topology.atom)
>>>                 gdg = assignments.graph_db_graph()
>>>                 gdt.graphs[gid] = gdg
>>>                 fn = rdata[tid]
>>>                 indices=dict(sorted([(j, x) for j, x in enumerate(g.nodes)], key=lambda x: x[1]))
>>>                 gdr = load_xyz([fn], indices=list(indices))
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
>>> def run(d, ff_fn):
>>>     # build the dataset and input ff
>>>     gcd = codec_rdkit.graph_codec_rdkit()
>>>     labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
>>>     pcp = perception.perception_model(gcd, labeler)
>>>     csys = smirnoff_models.smirnoff_load(ff_fn, pcp)
>>>     gdb = new_gdb(d)
>>>     psys = fits.gdb_to_physical_systems(gdb, csys)
>>>     models = [0,1]
>>>     fit_models = [0,1]
>>>     strat = fits.forcefield_optimization_strategy_default(csys, models=models)
>>>     co = fits.chemical_objective
>>>     final = None # default would be a full fit
>>>     initial = final # have a default, full fit
>>>     onestep = None # do a forcefield fit with one step 
>>>     final = fits.objective_tier()
>>>     final.objectives = {
>>>         0: fits.objective_config_position(
>>>                 fits.graph_db_address(
>>>                     eid=[0],
>>>                 ),
>>>                 scale=10
>>>         ),
>>>         1: fits.objective_config_gradient(
>>>                 fits.graph_db_address(
>>>                     eid=[0],
>>>                 ),
>>>                 scale=1e-5
>>>         ),
>>>     }
>>>     final.key_filter = lambda x: x[0] in fit_models and x[1] == 'l'
>>>     onestep = fits.objective_tier()
>>>     onestep.objectives = final.objectives
>>>     onestep.step_limit = 2
>>>     onestep.key_filter = lambda x: x[0] in fit_models and x[1] == 'l'
>>>     initial = final
>>>     tiers = [onestep] # have a default
>>>     kv0 = mm.chemical_system_iter_keys(csys)
>>>     newcsys, P, C = fits.ff_optimize(
>>>         csys,
>>>         gdb,
>>>         psys,
>>>         strat,
>>>         co,
>>>         initial,
>>>         tiers,
>>>         final
>>>     )
>>>     fits.print_chemical_system(newcsys)
>>>     print(f"X {P+C:15.8g} P {P:15.8g} C {C:15.8g}")
>>> run(make(), "openff-2.1.0.offxml")



We loaded the OpenFF 2.1.0 FF and fit bonds and angles to the geometry and
gradient of a molecule that was refined using B3LYP-D3BJ/dzvp in Psi4. The
output is much to long to show here. The important part is the initial objective:

.. code-block::

    RUNNING ANTECHAMBER FOR CHARGES
    2024-05-23 18:21:20.972785 Computing chemical objective
    2024-05-23 18:21:21.134459 Computing physical objective
    Started local workspace on ('0.0.0.0', 43477)
    2024-05-23 18:21:21.421784 Building physical systems
    2024-05-23 18:21:21.424306 Calculating initial obj
    2024-05-23 18:21:24.142884 Calculating fit
      0000 | X2=     1.1743 |g|=    5.6886
      0001 | X2=     10.107 |g|=    65.197
    >>> X2=     11.281 |g|=     64.39
      0000 | X2=     9.7226 |g|=    13.239
      0001 | X2=     599.91 |g|=    700.06
    >>> X2=     609.63 |g|=    710.77
      0000 | X2=     1.2118 |g|=    7.1619
      0001 | X2=     7.2125 |g|=    46.401
    >>> X2=     8.4243 |g|=    45.603
      0000 | X2=      0.199 |g|=    2.0676
      0001 | X2=     3.4447 |g|=    15.047
    >>> X2=     3.6437 |g|=    15.366
      0000 | X2=    0.18452 |g|=    2.1156
      0001 | X2=     2.7826 |g|=    11.486
    >>> X2=     2.9671 |g|=    12.639
      0000 | X2=     0.2805 |g|=   0.95383
      0001 | X2=     1.8705 |g|=     7.083
    >>> X2=      2.151 |g|=    7.3985
      0000 | X2=    0.31856 |g|=    1.7494
      0001 | X2=     1.5835 |g|=    9.4642
    >>> X2=     1.9021 |g|=    10.372
      0000 | X2=   0.048051 |g|=   0.80182
      0001 | X2=      1.373 |g|=    3.2839
    >>> X2=      1.421 |g|=    3.3126
      0000 | X2=   0.044377 |g|=   0.75911
      0001 | X2=     1.3004 |g|=    4.5262
    >>> X2=     1.3448 |g|=    4.2591
      0000 | X2=   0.042355 |g|=   0.73586
      0001 | X2=      1.206 |g|=    3.0855
    >>> X2=     1.2483 |g|=    2.6508
      0000 | X2=    0.12102 |g|=    1.6246
      0001 | X2=     1.1662 |g|=   0.70169
    >>> X2=     1.2873 |g|=    1.3192
      0000 | X2=   0.042573 |g|=   0.73299
      0001 | X2=     1.1988 |g|=    2.8028
    >>> X2=     1.2414 |g|=    2.3761
      0000 | X2=   0.042602 |g|=   0.73252
      0001 | X2=     1.1978 |g|=    2.7641
    >>> X2=     1.2404 |g|=    2.3387
    >>> Initial Objective     11.281
    >>> Final Objective       1.2404
    >>> Percent change       -89.004%
    2024-05-23 18:22:07.060736 Initial objective:       1.24045 C=      11.5003
    (0, 'l', 'b4', 0)    | New:      1.35782 Ref       1.4662 Diff    -0.108375
    (0, 'l', 'b6', 0)    | New:      1.30577 Ref      1.38236 Diff   -0.0765937
    (0, 'l', 'b17', 0)   | New:      1.19428 Ref      1.35775 Diff    -0.163462
    (0, 'l', 'b21', 0)   | New:      1.18682 Ref      1.22167 Diff   -0.0348451
    (0, 'l', 'b70', 0)   | New:      1.80143 Ref      1.72222 Diff    0.0792125
    (0, 'l', 'b85', 0)   | New:      1.08916 Ref      1.08182 Diff   0.00733215
    (1, 'l', 'a10', 0)   | New:      2.14613 Ref      2.09145 Diff    0.0546734
    (1, 'l', 'a14', 0)   | New:      2.14641 Ref      2.17754 Diff   -0.0311309
    (1, 'l', 'a29', 0)   | New:      1.76572 Ref      1.88807 Diff    -0.122343

Showing that a single fit is able to fit the gradient much better. At the bottom of the output we see the final
objective:

.. code-block::

      2024-05-23 18:28:01.567502 Calculating fit
      0000 | X2=   0.042602 |g|=   0.73241
      0001 | X2=     1.1978 |g|=    9.1799
    >>> X2=     1.2404 |g|=    9.1236
      0000 | X2=     11.698 |g|=    16.705
      0001 | X2=      137.9 |g|=    240.13
    >>> X2=      149.6 |g|=     244.8
      0000 | X2=    0.08572 |g|=    1.7278
      0001 | X2=    0.95806 |g|=    7.1875
    >>> X2=     1.0438 |g|=    7.5149
      0000 | X2=    0.22711 |g|=    11.757
      0001 | X2=    0.41881 |g|=    2.9737
    >>> X2=    0.64592 |g|=    9.4268
      0000 | X2= 0.00074601 |g|=  0.065677
      0001 | X2=      2.927 |g|= 1.207e+06
    >>> X2=     2.9278 |g|= 1.207e+06
      0000 | X2=          0 |g|=         0
      0001 | X2=    0.84076 |g|=    14.431
    >>> X2=    0.84076 |g|=    14.431
      0000 | X2= 0.00033565 |g|=   0.14047
      0001 | X2=    0.44734 |g|=    5.1086
    >>> X2=    0.44767 |g|=    5.1702
      0000 | X2= 4.1158e-05 |g|=   0.02142
      0001 | X2=    0.43314 |g|=    4.2271
    >>> X2=    0.43318 |g|=    4.2355
      0000 | X2=    0.11432 |g|=     3.489
      0001 | X2=    0.39621 |g|=    1.3572
    >>> X2=    0.51053 |g|=    3.3933
      0000 | X2= 4.1528e-05 |g|=  0.021769
      0001 | X2=     0.4291 |g|=    3.9855
    >>> X2=    0.42914 |g|=    3.9938
      0000 | X2= 4.1589e-05 |g|=    3.1302
      0001 | X2=    0.42855 |g|=    3.9516
    >>> X2=    0.42859 |g|=    6.3183
      0000 | X2= 0.00052464 |g|=   0.22546
      0001 | X2=    0.40583 |g|=    2.2435
    >>> X2=    0.40635 |g|=    2.3086
      0000 | X2=          0 |g|=         0
      0001 | X2=    0.38849 |g|=    1.5621
    >>> X2=    0.38849 |g|=    1.5621
      0000 | X2= 0.00062675 |g|=    0.2775
      0001 | X2=    0.38099 |g|=    1.7393
    >>> X2=    0.38162 |g|=    1.6208
    >>> Initial Objective     1.2404
    >>> Final Objective      0.38162
    >>> Percent change       -69.236%
    2024-05-23 18:28:50.825973 Accepting objective:      0.381616 C=      11.5215 DX=    -0.837628
    (0, 'l', 'b4', 0)    | New:      1.29592 Ref      1.35782 Diff   -0.0619061
    (0, 'l', 'b6', 0)    | New:      1.30216 Ref      1.30577 Diff  -0.00360538
    (0, 'l', 'b17', 0)   | New:      1.19934 Ref      1.19428 Diff   0.00505461
    (0, 'l', 'b21', 0)   | New:      1.19253 Ref      1.18682 Diff   0.00570217
    (0, 'l', 'b70', 0)   | New:      1.80154 Ref      1.80143 Diff  0.000109376
    (0, 'l', 'b85', 0)   | New:      1.08376 Ref      1.08916 Diff  -0.00539522
    (0, 'l', 'B91', 0)   | New:      1.42599 Ref      1.35782 Diff    0.0681625
    (1, 'l', 'a10', 0)   | New:      2.14545 Ref      2.14613 Diff -0.000678516
    (1, 'l', 'a14', 0)   | New:      2.14672 Ref      2.14641 Diff   0.00031012
    (1, 'l', 'a29', 0)   | New:      1.76722 Ref      1.76572 Diff   0.00150007

This shows that a bond was found which dropped the objective further. In particular,
it specialized the somewhat generic `[#6X3:1]-[#6X3:2]` b4 bond with
`[#6H0X3:1]-[#6X3:2]`. The new bond increased the bond length to 1.426 A, while
the original b4 parameter was originally at 1.466 A but dropped to 1.296 A in
the final result.

As a final note, the setup and results shown here is still undergoing
development. In the future, the setup should be easier and shorter, and the
output should be cleaner and more informative.


