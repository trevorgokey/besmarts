
"""
besmarts.mechanics.optimizers_scipy
"""

import copy
import scipy.optimize

from besmarts.core import geometry
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import objectives


def optimize_positions_scipy(csys, psys: mm.physical_system):

    pos = copy.deepcopy(psys.models[0].positions[0])
    
    args, keys = objectives.array_flatten_assignment(pos.selections)
    for pair in list(zip(keys, args)):
        print(pair)

    jac = objectives.array_geom_gradient
    # jac = None

    result = scipy.optimize.minimize(
        objectives.array_geom_energy,
        args,
        jac=jac,
        args=(keys, csys, psys),
        options={'disp': False, 'gtol': .1}
    )

    n_confs = len(list(pos.selections.values())[0])
    for (c, n, i), v in zip(keys, result.x):
        pos.selections[n][c][i] = v

    return pos

def ff_obj(x, keys, csys, refpsys, ref, refcsys):
    # print("input array:")
    # print(x)
    pos = copy.deepcopy(ref)
    recalc = set([x[0] for x in keys])

    reuse = set(range(len(csys.models))).difference(recalc)
    psys = mm.chemical_system_to_physical_system(csys, pos, ref=refpsys, reuse=reuse)
    optpos = optimize_positions_scipy(csys, psys)

    for k, v in zip(keys, x):
        mm.chemical_system_set_value(csys, k, v)

    obj = 0
    for ic, rxyz in ref[0].selections.items():
        oxyz = optpos.selections[ic]
        for a, b in zip(rxyz, oxyz):
            dx = geometry.array_distance(a, b)
            dx2 = dx*dx
            obj += dx2
            # print(obj, dx2)
    for k in keys:
        newv = mm.chemical_system_get_value(csys, k)
        refv = mm.chemical_system_get_value(refcsys, k)
        print(k, f"| New: {newv} Ref {refv} Diff {newv-refv}")
    print("Total obj:", obj)
    return obj

def optimize_forcefield_scipy(csys, pos):
    """
    Take a csys and a psys, then fit the params according to the flattened list
    The flattened list is the psys values.
    I need to append all psys values into the same list.
    Then, I need to take all modified values and reassign to psys
    """
    kv = mm.chemical_system_iter_keys(csys)
    kv = {k:v for k, v in kv.items() if k[0] in [0,1,2] and k[1] in "kl"}
    keys = list(kv.keys())
    x0 = list(kv.values())
    psys = mm.chemical_system_to_physical_system(csys, pos)
    # x0 = [x*.5 for x in x0]

    result = scipy.optimize.minimize(
        ff_obj,
        x0,
        args=(keys, csys, psys, pos, copy.deepcopy(csys)),
        options={'disp': False, 'gtol': .1}
    )
    return {k:v for k,v in zip(keys, result.x)}
