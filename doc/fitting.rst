
Fitting
-------

The first thing that needs to be created is the data. The data is a
`graph_topology_db`, and will hold the reference data such as the positions,
gradients, and hessians. For vector quantities such as positions and
gradients, the xyz and sdf formats are supported via their loaders:

```
tbl = file_xyz.file_xyz_to_graph_db_table(fn) ->         graph_db_table
tbl = file_xyz.file_xyz_to_smiles_assignment(fn, smi) -> smiles_assignment
tbl = file_xyz.file_xyz_to_graph_assignment(fn, g) ->    smiles_assignment
tbl = file_sdf.file_sdf_to_graph_db_table(fn) ->         graph_db_table
tbl = file_sdf.file_sdf_to_smiles_assignment(fn, smi) -> smiles_assignment
tbl = file_sdf.file_sdf_to_graph_assignment(fn, g) ->    smiles_assignment
```

These can then be added to a graph_db object:

```
graph_db_add_table(gdb, tbl, assignments.POSITIONS)
```

Once a gdb object is made, we need to define the calculations and the
objectives. Lets assume we want to fit the minima and so we need to 
set up a minimization calculation for the gdb. Now we setup an objective for
positions:

```
addr = {xid: [0], aid: [assignments.POSITIONS]}
obj_cfg_a = objective_cfg(addr, min=True)
obj_cfg_b = objective_cfg(addr, min=False)
```

and would generate a compute for each.

```
calc_a = compute_cfg(obj_cfg_a)
calc_b = compute_cfg(obj_cfg_b)
calcs = merge_calcs({0: calcs_a, 1: calcs_b})
result = refold(submit(set(calcs.values())))
obj = compute_obj(obj_cfg[0], gdb.get(obj_cfg[0].addr), result[0],  
obj.total()
```

This will pull all positions from a gdb list and create a new gdb object, and in this case
with new positions, gradients, and hessians.

This will create two objective objects. The goal here is to prevent this from
calculating the same thing multiple times. Each objective must define a 
calculation and must have a uniqueness property. We should therefore be
able to supply an objective with any calculation as long as it matches the
proper signature. I will also need to merge calculations:

```
{ objective_id: calc_obj}
{
    0: calc_ax,
    1: calc_ay,
    2: calc_bx,
    3: calc_by,
}
```
becomes
```
{ objective_id: calc_obj}
{
    0: calc_axy,
    1: calc_axy,
    2: calc_bxy,
    3: calc_bxy,
}
```
and so each calc will have a return that can be used by the objective. We create
a separate result
```
{ objective_id: result_obj}
{
    0: result_axy,
    1: result_axy,
    2: result_bxy,
    3: result_bxy,
}
```
and now we calculate as
compute(obj[cid], gdb.get(obj[cid].addr), res[cid])


Next, an initial force field must be set up:

```
csys = smirnoff.smirnoff_load(fn)
```

With these in 
