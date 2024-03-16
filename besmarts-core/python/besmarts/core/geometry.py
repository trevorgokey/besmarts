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

def array_scale(a, s):
    return type(a)((i*s for i in a))

def array_add(a, b):
    return type(a)((i+j for i,j in zip(a,b)))

def array_difference(a, b):
    return type(a)((i-j for i,j in zip(a,b)))

def array_multiply(a, b):
    return type(a)((i*j for i,j in zip(a,b)))

def array_inner_product(a, b):
    return sum((i*j for i,j in zip(a,b)))

def array_cross(a, b):
    return (
        (a[1]*b[2] - a[2]*b[1]),
        (a[2]*b[0] - a[0]*b[2]),
        (a[0]*b[1] - a[1]*b[0])
    )

def array_unit(a, b):
    """
    unit vector from a to b
    """
    return array_scale(array_difference(b, a), 1/array_distance(b, a))

def array_basis(a, b):
    """
    unit vector from a to b and its projection (magnitude)
    """
    r = array_distance(a, b)
    return array_scale(array_difference(b, a), 1/r), r

def array_magnitude(a) -> float:
    return sum([x*x for x in a])**.5


def array_distance(a,b) -> float:
    return sum([x*x for x in array_difference(b,a)])**.5

def measure_distance(xyz1, xyz2):

    result = [[array_distance(a,b)] for a, b in zip(xyz1, xyz2)]

    return result

def jacobian_distance(xyz1, xyz2):

    result = []

    for a, b in zip(xyz1, xyz2):

        r = array_distance(a, b)
        d = array_difference(b, a)
        
        # this will be an IC x XYZ matrix with each IC a key in the dict
        result.append([array_scale(d, -1/r), array_scale(d, 1/r)])

    return result

def jacobian_angle_geometric(xyz1, xyz2, xyz3):
    result = []
    basis = array_basis
    for a, b, c in zip(xyz1, xyz2, xyz3):
        u, u_norm =  basis(b, a)
        v, v_norm =  basis(b, c)
        w =  array_cross(u, v)
        w_norm = sum([x**2 for x in w])**.5
        w = array_scale(w, 1/w_norm)

        # w = w_prime / np.linalg.norm(w_prime)
        term1 = array_scale(array_cross(u, w), 1/ u_norm)
        term2 = array_scale(array_cross(w, v), 1/ v_norm)
        result.append([term1, array_scale(array_add(term1, term2), -1), term2])
    return result

def jacobian_angle(xyz1, xyz2, xyz3):

    return jacobian_angle_geometric(xyz1, xyz2, xyz3)
    result = []

    thetas = measure_angle(xyz1, xyz2, xyz3)
    for a, b, c, theta in zip(xyz1, xyz2, xyz3, thetas):
        theta = theta[0]
        eba, rba = array_basis(b, a) 
        ebc, rbc = array_basis(b, c) 

        cos = math.cos(theta)
        sin = math.sin(theta)
        s1 = array_scale(array_difference(array_scale(eba, cos), ebc), 1/(rba*sin))
        s2 = array_scale(array_difference(array_scale(ebc, cos), eba), 1/(rbc*sin))
        # s1 = array_scale(array_difference(eba, ebc), cos/(rba*sin))
        # s2 = array_scale(array_difference(ebc, eba), cos/(rbc*sin))
        s3a = array_scale(eba, rba - rbc*cos) 
        s3b = array_scale(ebc, rbc - rba*cos) 
        s3 = array_scale(array_add(s3a, s3b), rba*rbc*sin)
        result.append([s1, s3, s2])

    return result

def jacobian_outofplane_geometric(xyz1, xyz2, xyz3, xyz4):
    return jacobian_torsion_geometric(xyz1, xyz2, xyz3, xyz4)

def jacobian_outofplane(xyz1, xyz2, xyz3, xyz4):

    return jacobian_outofplane_geometric(xyz1, xyz2, xyz3, xyz4)
    result = []

    thetas = measure_dihedral(xyz1, xyz2, xyz3, xyz4)
    phi1 = measure_angle(xyz1, xyz2, xyz3)
    phi2 = measure_angle(xyz3, xyz2, xyz4)
    phi3 = measure_angle(xyz1, xyz2, xyz4)

    for a, b, c, d, t, p1, p2, p3 in zip(xyz1, xyz2, xyz3, xyz4, thetas, phi1, phi2, phi3):

        t = t[0]
        p1 = p1[0]
        p2 = p2[0]
        p3 = p3[0]
        eba, rba = array_basis(b, a)
        ebc, rbc = array_basis(b, c)
        ebd, rbd = array_basis(b, d)
        cos = math.cos(t)
        cos1 = math.cos(p1)
        sin1 = math.sin(p1)
        tan = math.tan(t)
        s1a = array_scale(array_cross(eba, ebc), 1/(cos*sin1)) 
        s1b = array_scale(ebd, tan)
        s1 = array_scale(array_difference(s1a, s1b), 1/rbd)

        s2a = array_scale(array_cross(ebc, eba), 1/(cos*sin1)) 
        s2b = array_scale(array_difference(ebc, array_scale(ebd, cos1)), tan/sin1**2)
        s2 = array_scale(array_difference(s2a, s2b), 1/rbc)

        s3a = array_scale(array_cross(eba, ebc), 1/(cos*sin1)) 
        s3b = array_scale(array_difference(ebd, array_scale(ebc, cos1)), tan/sin1**2)
        s3 = array_scale(array_difference(s1a, s1b), 1/rbd)
        
        s4 = array_difference(array_difference(array_scale(s1, -1.0), s2), s3)

        result.append([s1,s2,s3,s4])

    return result


def jacobian_torsion_geometric(xyz1, xyz2, xyz3, xyz4):

    """
    Adapted from Lee-Ping Wang's version in geomeTRIC

    Seems similar to Bright Decius and Cross except sin2 = 1 - dot and other
    tricks to avoid trig
    """
    dot = array_inner_product
    cross = array_cross
    scale = array_scale
    basis = array_basis

    result = []
    for a, b, c, d in zip(xyz1, xyz2, xyz3, xyz4):
        u, u_norm =  basis(b, a)
        w, w_norm =  basis(b, c)
        v, v_norm =  basis(c, d)
        if (1 - dot(u, w)**2) < 1e-6:
            term1 = scale(cross(u, w) , 0)
            term3 = scale(cross(u, w) , 0)
        else:
            term1 = scale(cross(u, w), 1 / (u_norm * (1 - dot(u, w)**2)))
            term3 = scale(cross(u, w), dot(u, w) / (w_norm * (1 - dot(u, w)**2)))
        if (1 - dot(v, w)**2) < 1e-6:
            term2 = scale(cross(v, w) , 0)
            term4 = scale(cross(v, w) , 0)
        else:
            term2 = scale(cross(v, w), 1 / (v_norm * (1 - dot(v, w)**2)))
            term4 = scale(cross(v, w), dot(v, w) / (w_norm * (1 - dot(v, w)**2)))

        s1 = term1
        s2 = array_difference(array_add(scale(term1, -1), term3), term4)
        s3 = array_difference(array_add(term2, term4), term3)
        s4 = scale(term2, -1)

        result.append([s1, s2, s3, s4])
    return result

def jacobian_torsion(xyz1, xyz2, xyz3, xyz4):

    # return jacobian_torsion_geometric(xyz1, xyz2, xyz3, xyz4)

    # works
    result = []

    thetas = measure_dihedral(xyz1, xyz2, xyz3, xyz4)
    phi2 = measure_angle(xyz1, xyz2, xyz3)
    phi3 = measure_angle(xyz2, xyz3, xyz4)

    for a, b, c, d, t, p2, p3 in zip(xyz1, xyz2, xyz3, xyz4, thetas, phi2, phi3):

        t = t[0]
        p2 = p2[0]
        p3 = p3[0]
        eab, rab = array_basis(a, b)
        ebc, rbc = array_basis(b, c)
        sin2 = math.sin(p2)
        sin3 = math.sin(p3)
        s1 = array_scale(array_cross(eab, ebc), -1/(rab*sin2**2))

        cos2 = math.cos(p2)
        cos3 = math.cos(p3)
        s2a = (rbc - rab * cos2)/(rbc*rab*sin2**2)
        s2b = cos3/(rbc*sin3**2)

        edc, rdc = array_basis(d, c)
        ecb, rcb = array_basis(c, b)
        s2 = array_difference(
            array_scale(array_cross(eab, ebc), s2a), 
            array_scale(array_cross(edc, ecb), s2b)
        )

        s1p = array_scale(array_cross(edc, ecb), cos3/(rcb*sin3**2))
        s2 = array_difference(
            array_scale(s1, -(rbc - rab * cos2)/rbc),
            s1p
        )


        s3a = (rcb - rdc * cos3)/(rcb*rdc*sin3**2)
        s3b = cos2/(rcb*sin2**2)
        s3 = array_difference(
            array_scale(array_cross(edc, ecb), s3a), 
            array_scale(array_cross(eab, ebc), s3b)
        )

        # 43 is dc
        s4 = array_scale(array_cross(edc, ecb), -1/(rdc*sin3**2))

        result.append([s1,s2,s3,s4])

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

        result.append([theta])

    return result

def measure_dihedral(xyz1, xyz2, xyz3, xyz4):
    result = []
    for a, b, c, d in zip(xyz1, xyz2, xyz3, xyz4):

        v1 = array_difference(b, a)
        v2 = array_difference(c, b)
        v3 = array_difference(d, c)
        c1 = array_cross(v2, v3)
        c2 = array_cross(v1, v2)
        y = sum(array_multiply(v1, c1)) * array_magnitude(v2)
        x = sum(array_multiply(c1, c2))
        theta = math.atan2(y, x)
        result.append([theta])
        # rr10 = (x0 - x1),  (y0 - y1), (z0 - z1)
        # r10  = sum([x**2 for x in rr10])**.5
        # rr10 = [x/r10 for x in rr10]

        # rr12 = (x2 - x1),  (y2 - y1), (z2 - z1)
        # r12 = sum([x**2 for x in rr12])**.5
        # rr12 = [x/r12 for x in rr12]
        
        # n1 = array_cross(rr10, rr12)
        # # n1 = (
        # #     (rr10[1]*rr12[2] - rr10[2]*rr12[1]),
        # #     (rr10[2]*rr12[0] - rr10[0]*rr12[2]),
        # #     (rr10[0]*rr12[1] - rr10[1]*rr12[0])
        # # )
        # e1 = sum([x**2 for x in n1])**.5
        # n1 = [x/e1 for x in n1]

        # rr23 = (x3 - x2),  (y3 - y2), (z3 - z2)
        # r23 = sum([x**2 for x in rr23])**.5
        # rr23 = [x/r23 for x in rr23]

        # n2 = array_cross(rr23, rr12)
        # # n2 = (
        # #     (rr23[1]*rr12[2] - rr23[2]*rr12[1]), 
        # #     (rr23[2]*rr12[0] - rr23[0]*rr12[2]),
        # #     (rr23[0]*rr12[1] - rr23[1]*rr12[0])
        # # )
        # e2 = sum([x**2 for x in n2])**.5
        # n2 = [x/e2 for x in n2]

        # proj = sum([a*b for a,b in zip(n1, n2)])
        # theta = math.acos(proj)

        # result.append([theta])

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
