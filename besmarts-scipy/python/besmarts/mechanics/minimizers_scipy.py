
"""
besmarts.mechanics.minimizers_scipy
"""

import copy
import scipy.optimize

from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import objectives

def minimization_scipy(csys, psys: mm.physical_system):
    # create the system, minimize, return new state

    pos = copy.deepcopy(psys.models[0].positions[0])
    # ff = mm.physical_system_iter_keys([psys], csys)
    
    args, keys = objectives.array_flatten_assignment(pos.selections)

    jac = objectives.array_geom_gradient
    # jac = None

    result = scipy.optimize.minimize(
        objectives.array_geom_energy,
        args,
        jac=jac,
        args=(keys, csys, psys),
        options={'disp': True, 'gtol': .01}
    )

    n_confs = len(list(pos.selections.values())[0])
    for (c, n, i), v in zip(keys, result.x):
        pos.selections[n][c][i] = v

    return pos

