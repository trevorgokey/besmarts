"""
besmarts.core.geometry

Functions to process the geometry of a molecular graph.
"""

import math
from besmarts.core import topology

INF = math.inf

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

def array_outer_product(a, b):
    return [x for ai in a for x in array_scale(b,  ai)]

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
    rinv = INF
    if r > 0.0:
        rinv = 1/r
    u = array_scale(array_difference(b, a), rinv)
    return array_scale(array_difference(b, a), rinv), r

def array_magnitude(a) -> float:
    return sum([x*x for x in a])**.5


def array_distance(a,b) -> float:
    return math.sqrt(sum([x*x for x in array_difference(b,a)]))

def array_round(a, b):
    return [round(x, b) for x in a]

def measure_distance(xyz1, xyz2):

    result = [[array_distance(a,b)] for a, b in zip(xyz1, xyz2)]

    return result

def jacobian_distance(xyz1, xyz2):

    result = []

    for a, b in zip(xyz1, xyz2):

        r = array_distance(a, b)
        d = array_difference(b, a)
        
        # this will be an IC x XYZ matrix with each IC a key in the dict
        rinv = INF
        if r > 0.0:
            rinv = 1/r
        result.append([array_scale(d, -rinv), array_scale(d, rinv)])

    return result

def jacobian2_distance_terms(u, lu):

    x = [(u[i]*u[j] - int(i == j))/lu for i in range(3) for j in range(3)]

    return x

def jacobian2_distance(xyz1, xyz2):

    """
    from Bakken & Helgaker 2012 10.1063/1.1515483
    gives 2x2 1x9 arrays per position, (4 3x3)
    """

    result = []

    for a, b in zip(xyz1, xyz2):

        mat = []

        u, lu = array_basis(a, b)
        t = jacobian2_distance_terms(u, lu)

        row = []
        row.append(array_scale(t, -1))
        row.append(array_scale(t,  1))
        mat.append(row)
        
        row = []
        row.append(array_scale(t,  1))
        row.append(array_scale(t, -1))
        mat.append(row)

        result.append(mat)

    return result

def jacobian_angle(xyz1, xyz2, xyz3):
    result = []
    basis = array_basis
    for a, b, c in zip(xyz1, xyz2, xyz3):
        u, u_norm =  basis(b, a)
        v, v_norm =  basis(b, c)
        w =  array_cross(u, v)
        w_norm = sum([x**2 for x in w])**.5
        if w_norm < 1e-9:
            w_norm = 1.0
        w = array_scale(w, 1/w_norm)

        # w = w_prime / np.linalg.norm(w_prime)
        term1 = array_scale(array_cross(u, w), 1/ u_norm)
        term2 = array_scale(array_cross(w, v), 1/ v_norm)
        result.append([term1, array_scale(array_add(term1, term2), -1), term2])
    return result

def jacobian2_angle_terms(uu, uv, vv, lu, lv, cosq, sinq):

    scale = array_scale
    lus = 1/(lu*lu*sinq)
    lvs = 1/(lv*lv*sinq)

    luvs = 1/(lu*lv*sinq)


    t1 = [
        (uv[3*i+j] + uv[3*j+i] - 3*uu[3*i+j]*cosq + int(i==j)*cosq) * lus
            for i in range(3) for j in range(3)
    ]

    t2 = [
        (uv[3*j+i] + uv[3*i+j] - 3*vv[3*i+j]*cosq + int(i==j)*cosq) * lvs
            for i in range(3) for j in range(3)
    ]

    t3 = [
        (uu[3*i+j] + vv[3*j+i] - uv[3*i+j]*cosq - int(i==j)) * luvs
            for i in range(3) for j in range(3)
    ]

    t4 = [
        (uu[3*j+i] + vv[3*i+j] - uv[3*j+i]*cosq - int(i==j)) * luvs
            for i in range(3) for j in range(3)
    ]

    return [t1,t2,t3,t4]

    

def jacobian2_angle(xyz1, xyz2, xyz3):

    """
    gives 3x3 1x9 arrays
    """

    result = []
    basis = array_basis
    dq = jacobian_angle(xyz1, xyz2, xyz3)

    for m, o, n, (dqm, dqo, dqn) in zip(xyz1, xyz2, xyz3, dq):
        u, u_norm =  basis(o, m)
        v, v_norm =  basis(o, n)

        uu = array_outer_product(u, u)
        uv = array_outer_product(u, v)
        vv = array_outer_product(v, v)
        cosq = array_inner_product(u, v)
        sinq = math.sqrt(1 - min(1.0, cosq**2))

        if abs(sinq) < 1e-7:
            cs = 0.0
            t = []
        else:
            cs = -cosq/sinq
            t = jacobian2_angle_terms(uu, uv, vv, u_norm, v_norm, cosq, sinq)

        mat = []
        row = []

        # a=m, b=m [++, 00, +0, 0+]
        row.append(jacobian2_angle_term_reduce([ 1, 0, 0, 0], t, cs, dqm, dqm))
        # a=m, b=o [+-, 0-, +-, 0-]
        row.append(jacobian2_angle_term_reduce([-1, 0,-1, 0], t, cs, dqm, dqo))
        # a=m, b=n [+0, 0+, ++, 00]
        row.append(jacobian2_angle_term_reduce([ 0, 0, 1, 0], t, cs, dqm, dqn))

        mat.append(row)
        row = []
        # a=o, b=m [-+, -0, -0, -+]
        row.append(jacobian2_angle_term_reduce([-1, 0, 0,-1], t, cs, dqo, dqm))
        # a=o, b=o [--, --, --, --]
        row.append(jacobian2_angle_term_reduce([ 1, 1, 1, 1], t, cs, dqo, dqo))
        # a=o, b=n [-0, -+, -+, -0]
        row.append(jacobian2_angle_term_reduce([ 0,-1,-1, 0], t, cs, dqo, dqn))

        mat.append(row)
        row = []
        # a=n, b=m [0+, +0, 00, ++]
        row.append(jacobian2_angle_term_reduce([ 0, 0, 0, 1], t, cs, dqn, dqm))
        # a=n, b=o [0-, +-, 0-, +-]
        row.append(jacobian2_angle_term_reduce([ 0,-1, 0,-1], t, cs, dqn, dqo))
        # a=n, b=n [00, ++, 0+, +0]
        row.append(jacobian2_angle_term_reduce([ 0, 1, 0, 0], t, cs, dqn, dqn))

        mat.append(row)
        result.append(mat)

        # for visualization; order is aa bb ab ba
        # [++, 00, +0, 0+] [+-, 0-, +-, 0-] [+0, 0+, ++, 00]
        # [-+, -0, -0, -+] [--, --, --, --] [-0, -+, -+, -0]
        # [0+, +0, 00, ++] [0-, +-, 0-, +-] [00, ++, 0+, +0]

    return result

def jacobian2_angle_term_reduce(coef, terms, cs, dq1, dq2):

    x = list([0]*9)
    for a, term in zip(coef, terms):
        x = array_add(array_scale(term, a), x)
    t0 = [cs*dq1[i]*dq2[j] for i in range(3) for j in range(3)]
    x = array_add(t0, x)
    return x


def jacobian_outofplane(xyz1, xyz2, xyz3, xyz4):
    return jacobian_torsion(xyz1, xyz2, xyz3, xyz4)

def jacobian2_outofplane(xyz1, xyz2, xyz3, xyz4):
    return jacobian2_torsion(xyz1, xyz2, xyz3, xyz4)

def jacobian_outofplane_v2(xyz1, xyz2, xyz3, xyz4):

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

        rbdinv = float("inf")
        if rbdinv != 0.0:
            rbdinv = 1/rbd
        s1 = array_scale(array_difference(s1a, s1b), rbdinv)

        cossininv = float("inf")
        if (cos*sin1) != 0.0:
            cossininv = 1/(cos*sin1)
        s2a = array_scale(array_cross(ebc, eba), cossininv) 

        tansin = 0.0
        if tan == 0.0:
            tansin = 0.0
        elif sin1 == 0.0:
            tansin = INF
        else:
            tansin = tan/sin1**2
        s2b = array_scale(array_difference(ebc, array_scale(ebd, cos1)), tansin)

        rbcinv = INF
        if rbcinv != 0.0:
            rbcinv = 1/rbc
        s2 = array_scale(array_difference(s2a, s2b), rbcinv)

        s3a = array_scale(array_cross(eba, ebc), cossininv) 
        s3b = array_scale(array_difference(ebd, array_scale(ebc, cos1)), tansin)

        s3 = array_scale(array_difference(s1a, s1b), rbdinv)
        s4 = array_difference(array_difference(array_scale(s1, -1.0), s2), s3)

        result.append([s1,s2,s3,s4])

    return result


def jacobian_torsion(xyz1, xyz2, xyz3, xyz4):

    """
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

def jacobian_torsion_v2(xyz1, xyz2, xyz3, xyz4):

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

def jacobian2_torsion_terms(u, v, w, lu, lv, lw, cosu, cosv):

    """
    """

    """
    Incorporates corrections from psi4numpy and geomTRIC of 10.1063/1.1515483
    
    All of the signs agree with psi4numpy, however follow geomeTRICs zeta
    coefficients. psi4numpy from static analysis of the code shows that
    term7 and term8 will only hit 3 times, when it seems like they should hit 6
    times each. The atoms that hit term8 according to geomeTRIC are
         o   p   n
      o  0  -1   1
      p  1   0  -1
      n -1   1   0

    with the remaining pairs all 0, whereas psi4numpy is
         o   p   n
      o  0   0   0
      p  1   0   0
      n -1   1   0
     
    which doesn't feel right when you consider that the zeta sign functions
    allow swapping of indices with a sign change, which is what happens in the
    geomeTRIC case where eg po = -op.

    That said, they seem to get the same result for an arbitrary torsion, so 
    it's hard to tell what is what. 
    """

    sinu = math.sqrt(1 - min(1.0, cosu**2))
    sinv = math.sqrt(1 - min(1.0, cosv**2))

    if abs(sinu) < 1e-12 or abs(sinv) < 1e-12:
        return [*[[*[0]*9]]*8]

    uxw = array_cross(u, w)
    vxw = array_cross(v, w)

    s2u = sinu * sinu
    s4u = s2u * s2u

    s2v = sinv * sinv
    s4v = s2v * s2v

    c2u = cosu * cosu
    c3u = c2u * cosu

    c2v = cosv * cosv
    c3v = c2v * cosv

    t1 = [
        (uxw[i] * (w[j]*cosu - u[j]) 
       + uxw[j] * (w[i]*cosu - u[i]))/(lu*lu*s4u)
            for i in range(3) for j in range(3)
    ]

    t2 = [
        (vxw[i] * (w[j]*cosv + v[j])
       + vxw[j] * (w[i]*cosv + v[i]))/(lv*lv*s4v)
            for i in range(3) for j in range(3)
    ]

    t3 = [
        (uxw[i] * (w[j] - 2.*u[j]*cosu + w[j]*c2u) 
       + uxw[j] * (w[i] - 2.*u[i]*cosu + w[i]*c2u)) / (2.*lu*lw*s4u)
            for i in range(3) for j in range(3)
    ]

    # psi4numpy says this should be v[j] instead of u[j]
    t4 = [
        (vxw[i] * (w[j] + 2.*v[j]*cosv + w[j]*c2v) 
       + vxw[j] * (w[i] + 2.*v[i]*cosv + w[i]*c2v)) / (2.*lv*lw*s4v)
            for i in range(3) for j in range(3)
    ]

    t5 = [
        (uxw[i] * (u[j] + u[j]*c2u - 3.*w[j]*cosu + w[j]*c3u) 
       + uxw[j] * (u[i] + u[i]*c2u - 3.*w[i]*cosu + w[i]*c3u)) / (2.*lw*lw*s4u) 
            for i in range(3) for j in range(3)
    ]

    t6 = [
        (vxw[i] * (-v[j] - v[j]*c2v - 3.*w[j]*cosv + w[j]*c3v) 
       + vxw[j] * (-v[i] - v[i]*c2v - 3.*w[i]*cosv + w[i]*c3v)) / (2.*lw*lw*s4v)
            for i in range(3) for j in range(3)
    ]

    # and psi4numpy says it is sin2 instead of sin
    t7 = [0 if i == j else
        (j-i)*(-.5)**(abs(j-i)) * (-w[3-i-j]*cosu + u[3-i-j])/(lu*lw*s2u)
            for i in range(3) for j in range(3)
    ]

    t8 = [0 if i == j else 
        (j-i)*(-.5)**(abs(j-i)) * (-w[3-i-j]*cosv - v[3-i-j])/(lv*lw*s2v)
            for i in range(3) for j in range(3) 
    ]

    return [t1, t2, t3, t4, t5, t6, t7, t8]
    
def jacobian2_reshape_to_matrix(b2):
    """
    Give the canonical 3Nx3N form
    """
    c2 = []
    for a in b2:
        for i in range(3):
            c2.append([x for xyz in a for x in xyz[3*i:3*i+3]])
    return c2

def jacobian2_torsion(xyz1, xyz2, xyz3, xyz4):

    """
    Returns 4x4 1x9s
    """
    dot = array_inner_product
    cross = array_cross
    scale = array_scale
    basis = array_basis

    result = []

    for m, o, p, n in zip(xyz1, xyz2, xyz3, xyz4):
        mat = []
        u, ul = basis(o, m)
        w, wl = basis(o, p)
        v, vl = basis(p, n)
        cosu = dot(u, w)
        cosv = -dot(w, v)

        t = jacobian2_torsion_terms(u, v, w, ul, vl, wl, cosu, cosv)

        mat = []
        row = []

        # from the paper, change t6 from aopbop          to apobop
        # from the paper, change t7 from amobop + apobom to amobpo + apobom
        # from the paper, change t8 from anobop + apobom to anobpo + apobon

        #a=m b=m [++, 00, +0p0-, 00p00, 00, 00, 0(+0p0-), 0(00p00)]
        row.append(jacobian2_torsion_term_reduce([ 1, 0, 0, 0, 0, 0, 0, 0], t))
        #a=m b=o [+-, 00, ++p0+, 0-p00, 0-, 0+, +(+-p0+), +(0-p00)]
        row.append(jacobian2_torsion_term_reduce([-1, 0, 1, 0, 0, 0,-1, 0], t))
        #a=m b=p [+0, 0-, +-p00, 0+p0-, 0+, 0-, +(++p00), +(0+p00)]
        row.append(jacobian2_torsion_term_reduce([ 0, 0,-1, 0, 0, 0, 1, 0], t))
        #a=m b=n [+0, 0+, +0p00, 00p0+, 00, 00, +(+0p00), +(00p0+)]
        row.append(jacobian2_torsion_term_reduce([ 0, 0, 0, 0, 0, 0, 0, 0], t))
        mat.append(row)
        
        row = []
        #a=o b=m [-+, 00, -0p--, 00p-0, +0, -0, +(-0p--), +(-0p-0)]
        row.append(jacobian2_torsion_term_reduce([-1, 0, 1, 0, 0, 0, 1, 0], t))
        #a=o b=o [--, 00, -+p-+, 0-p-0, +-, -+, 0(--p-+), 0(--p-+)]
        row.append(jacobian2_torsion_term_reduce([ 1, 0,-2, 0,-1,-1, 0, 0], t))
        #a=o b=p [-0, 0-, --p-0, 0+p--, ++, --, +(-+p-0), +(-+p-0)]
        row.append(jacobian2_torsion_term_reduce([ 0, 0, 1, 1, 1, 1,-1,-1], t))
        #a=o b=n [-0, 0+, -0p-0, 00p-+, +0, -0, +(-0p-0), +(-0p--)]
        row.append(jacobian2_torsion_term_reduce([ 0, 0, 0,-1, 0, 0, 0, 1], t))
        mat.append(row)

        row = []
        #a=p b=m [0+, -0, 00p+-, -0p+0, -0, +0, +(00p+-), +(00p+0)]
        row.append(jacobian2_torsion_term_reduce([ 0, 0,-1, 0, 0, 0,-1, 0], t))
        #a=p b=o [0-, -0, 0+p++, --p+0, --, ++, +(0-p++), +(0-p++)]
        row.append(jacobian2_torsion_term_reduce([ 0, 0, 1, 1, 1, 1, 1, 1], t))
        #a=p b=p [00, --, 0-p+0, -+p+-, -+, +-, 0(0+p+0), 0(0+p+0)]
        row.append(jacobian2_torsion_term_reduce([ 0, 1, 0,-2,-1,-1, 0, 0], t))
        #a=p b=n [00, -+, 00p+0, -0p++, -0, +0, +(00p+0), +(00p+-)]
        row.append(jacobian2_torsion_term_reduce([ 0,-1, 0, 1, 0, 0, 0,-1], t))
        mat.append(row)

        row = []
        #a=n b=m [0+, +0, 00p0-, +0p00, 00, 00, +(00p0-), +(+0p00)]
        row.append(jacobian2_torsion_term_reduce([ 0, 0, 0, 0, 0, 0, 0, 0], t))
        #a=n b=o [0-, +0, 0+p0+, +-p00, 0-, 0+, +(0-p0+), +(+-p0+)]
        row.append(jacobian2_torsion_term_reduce([ 0, 0, 0,-1, 0, 0, 0,-1], t))
        #a=n b=p [00, +-, 0-p00, ++p0-, 0+, 0-, +(0+p00), +(++p00)]
        row.append(jacobian2_torsion_term_reduce([ 0,-1, 0, 1, 0, 0, 0, 1], t))
        #a=n b=n [00, ++, 00p00, +0p0+, 00, 00, 0(00p00), 0(+0p0-)]
        row.append(jacobian2_torsion_term_reduce([ 0, 1, 0, 0, 0, 0, 0, 0], t))
        mat.append(row)

        result.append(mat)

    return result

def jacobian2_torsion_term_reduce(coef, terms):

    x = list([0]*12)
    for a, term in zip(coef, terms):
        x = array_add(array_scale(term, a), x)
    return x


def measure_angle(xyz1, xyz2, xyz3):
    result = []
    for (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) in zip(xyz1, xyz2, xyz3):
        rr10 = (x0 - x1), (y0 - y1), (z0 - z1)
        r10 = sum([x**2 for x in rr10])**.5

        rr12 = (x2 - x1),  (y2 - y1), (z2 - z1)
        r12 = sum([x**2 for x in rr12])**.5

        if r10 == 0.0 or r12 == 0.0:
            result.append([0.0])
            continue

        rr10 = [x/r10 for x in rr10]
        rr12 = [x/r12 for x in rr12]

        proj = sum([a*b for a, b in zip(rr10, rr12)])

        if proj >= 1.0:
            theta = 0.0
        elif proj <= -1.0:
            theta = math.pi
        else:
            theta = math.acos(proj)

        result.append([theta])

    return result


def measure_dihedral(xyz1, xyz2, xyz3, xyz4):
    result = []
    for a, b, c, d in zip(xyz1, xyz2, xyz3, xyz4):

        v1 = array_difference(b, a)
        v2 = array_difference(c, b)
        c2 = array_cross(v1, v2)

        v3 = array_difference(d, c)
        c1 = array_cross(v2, v3)

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



