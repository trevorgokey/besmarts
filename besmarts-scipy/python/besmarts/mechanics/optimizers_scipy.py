
"""
besmarts.mechanics.optimizers_scipy
"""

import copy
import scipy.optimize

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
