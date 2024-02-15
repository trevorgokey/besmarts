
from besmarts.core import topology
from besmarts.core import assignments

def derivative_angle(pos, indices):
    """
    need to write out dq/dx...
    """
    pass

def derivative_bond(pos, indices):

    xyz = pos.selections
    selections = {}

    for bond in indices:
        c1, c2 = bond
        selections[bond] = []
        for (x0,y0,z0), (x1,y1,z1) in zip(xyz[c1,], xyz[c2,]):
            r = ((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z0 - z1) ** 2) ** 0.5
            zeros = list([0.0] * 3*len(pos.selections))
            zeros[3*(c1-1)+0] = (x1 - x0)/r
            zeros[3*(c1-1)+1] = (y1 - y0)/r
            zeros[3*(c1-1)+2] = (z1 - z0)/r
            zeros[3*(c2-1)+0] = (x0 - x1)/r
            zeros[3*(c2-1)+1] = (y0 - y1)/r
            zeros[3*(c2-1)+2] = (z0 - z1)/r
            # this will be an IC x XYZ matrix with each IC a key in the dict
            selections[bond].append(zeros)

    return bond_assignment_float(selections)

def geometry_derivative(top, pos):

    if top in [topology.bond, topology.pair]:
        return jacobian_bond(pos)
    

