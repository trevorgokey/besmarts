"""
besmarts.core.geometry

Functions to process the geometry of a molecular graph.
"""

import math

def is_outofplane(combo, edges) -> bool:
    return (
        tuple(sorted((combo[0], combo[1]))) in edges
        and tuple(sorted((combo[1], combo[2]))) in edges
        and tuple(sorted((combo[1], combo[3]))) in edges
    )

def is_torsion(combo, edges) -> bool:
    return (
        tuple(sorted((combo[2], combo[3]))) in edges
        and tuple(sorted((combo[0], combo[1]))) in edges
        and tuple(sorted((combo[1], combo[2]))) in edges
    )

def is_dihedral(combo, edges) -> bool:
    return is_torsion(combo, edges) or is_outofplane(combo, edges)

def measure_distance(xyz1, xyz2):

    result = []

    for (x0,y0,z0), (x1,y1,z1) in zip(xyz1, xyz2):
        r = ((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z0 - z1) ** 2) ** 0.5
        result.append(r)

    return result

def measure_angle(xyz1, xyz2, xyz3):
    result = []
    for (x0,y0,z0), (x1,y1,z1), (x2,y2,z2) in zip(xyz1, xyz2, xyz3):
        rr10 = (x0 - x1),  (y0 - y1), (z0 - z1)
        r10 = sum([x**2 for x in rr10])**.5
        rr10 = [x/r10 for x in rr10]

        rr12 = (x2 - x1),  (y2 - y1), (z2 - z1)
        r12 = sum([x**2 for x in rr12])**.5
        rr12 = [x/r12 for x in rr12]

        proj = sum([a*b for a,b in zip(rr10, rr12)])
        theta = math.acos(proj)

        result.append(theta)

    return result

def measure_dihedral(xyz1, xyz2, xyz3, xyz4):
    result = []
    for (x0,y0,z0), (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) in zip(xyz1, xyz2, xyz3, xyz4):

        rr10 = (x0 - x1),  (y0 - y1), (z0 - z1)
        r10  = sum([x**2 for x in rr10])**.5
        rr10 = [x/r10 for x in rr10]

        rr12 = (x2 - x1),  (y2 - y1), (z2 - z1)
        r12 = sum([x**2 for x in rr12])**.5
        rr12 = [x/r12 for x in rr12]
        
        n1 = (
            (rr10[1]*rr12[2] - rr10[2]*rr12[1]),
            (rr10[2]*rr12[0] - rr10[0]*rr12[2]),
            (rr10[0]*rr12[1] - rr10[1]*rr12[0])
        )
        e1 = sum([x**2 for x in n1])**.5
        n1 = [x/e1 for x in n1]

        rr23 = (x3 - x2),  (y3 - y2), (z3 - z2)
        r23 = sum([x**2 for x in rr23])**.5
        rr23 = [x/r23 for x in rr23]

        n2 = (
            (rr23[1]*rr12[2] - rr23[2]*rr12[1]), 
            (rr23[2]*rr12[0] - rr23[0]*rr12[2]),
            (rr23[0]*rr12[1] - rr23[1]*rr12[0])
        )
        e2 = sum([x**2 for x in n2])**.5
        n2 = [x/e2 for x in n2]

        proj = sum([a*b for a,b in zip(n1, n2)])
        theta = math.acos(proj)

        result.append(theta)

    return result

def bond(x):
    if x[1] < x[0]:
        x = x[::-1]
    return tuple(x)

def pair(x):
    if x[1] < x[0]:
        x = x[::-1]
    return tuple(x)

def angle(x):
    if x[2] < x[0]:
        x = x[::-1]
    return tuple(x)

def torsion(x):
    if x[3] < x[0]:
        x = x[::-1]
    return tuple(x)

def outofplane(x):
    y = sorted((x[0], x[2], x[3]))
    return tuple((y[0], x[1], *y[1:]))
