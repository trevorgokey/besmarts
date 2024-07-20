"""
besmarts.mechanics.hessians
"""

import sys
import copy
import numpy as np
import itertools
from besmarts.core import assignments
from besmarts.core import arrays
from besmarts.core import graphs
from besmarts.core import geometry
from besmarts.mechanics import objectives
from besmarts.mechanics import optimizers_scipy
# from besmarts.mechanics import optimizers_openmm
from besmarts.codecs import codec_rdkit
from besmarts.mechanics import molecular_models as mm
import numpy as np
import sys
import pprint
from pprint import pprint

from besmarts.mechanics import vibration
au2kcal = 627.51/(0.529177**2)

def get_xyz(fn):
    frames = []
    syms = []
    frame = []
    sym = []
    with open(fn) as f:
        for i, line in enumerate(f):
            # print("loading line", i)
            line = line.replace('\n','')
            if i == 0:
                natom = int(line.strip())
                continue
            if i % (natom+2) == 0 and frame:
                frames.append(frame)
                # print("Added frame", len(frames))
                frame = []
                syms.append(sym)
                sym = []
                continue
            elif i % (natom+2) > 1:
                # print("loading xyz", line)
                s,x,y,z = line.split()

                xyz = [float(x), float(y), float(z)]
                sym.append(s)
                frame.append(xyz)
        if frame:
            frames.append(frame)
            syms.append(sym)
    return syms[0], np.array(frames[0])

def load():
    syms, frame = get_xyz(xyzfile)
    g = gcd.smiles_decode(smiles)
    N = len(g.nodes)

def b2matrix(pos, bonds=True, angles=True, torsions=True, outofplanes=True, pairs=True, remove1_3=False):

    N = len(pos.selections)
    i = 0
    weights = []
    if bonds:
        bonds = assignments.graph_assignment_jacobian2_bonds(pos)

    if angles:
        angles = assignments.graph_assignment_jacobian2_angles(pos)

    if torsions:
        torsions = assignments.graph_assignment_jacobian2_torsions(pos)

    if outofplanes:
        outofplanes = assignments.graph_assignment_jacobian2_outofplanes(pos)

    if pairs:
        pairs = assignments.graph_assignment_jacobian2_pairs(pos)
        if angles and remove1_3:
            for a in angles.selections:
                if (a[0], a[2]) in pairs.selections:
                    del pairs.selections[(a[0], a[2])]

    B2 = []
    ics = []

    for valence in filter(None, (bonds, angles, torsions, outofplanes, pairs)):
        ics.extend(valence.selections)
        for ic, confs in valence.selections.items():

            b2 = []
            for conf in confs:
                b2.append(geometry.jacobian2_reshape_to_matrix(conf))
            B2.append(b2)

    return ics, B2
    
def bmatrix(pos, bonds=True, angles=True, torsions=True, outofplanes=True, pairs=True, remove1_3=False):


    def print_ic(ic, i=0):
        for j, (k, v) in enumerate(ic.selections.items(), i+1):
            print(j, k, v)
        return j

    N = len(pos.selections)
    i = 0
    if bonds:
        bonds = assignments.graph_assignment_jacobian_bonds(pos)

    if angles:
        angles = assignments.graph_assignment_jacobian_angles(pos)

    if torsions:
        torsions = assignments.graph_assignment_jacobian_torsions(pos)

    if outofplanes:
        outofplanes = assignments.graph_assignment_jacobian_outofplanes(pos)

    if pairs:
        pairs = assignments.graph_assignment_jacobian_pairs(pos)
        if angles and remove1_3:
            for a in angles.selections:
                if (a[0], a[2]) in pairs.selections:
                    del pairs.selections[(a[0], a[2])]


    conf = 0
    b = {}
    B = []
    ics = []
    for valence in filter(None, (bonds, angles, torsions, outofplanes, pairs)):
        ics.extend(valence.selections)
        for ic, confs in valence.selections.items():
            brow = np.zeros(3*N)
            for i, nid in enumerate(ic):
                d = (nid-1)*3
                brow[d:d+3] += confs[conf][i]
            B.append(brow)

    B = np.array(B)
    return ics, B

def transform(B):

    """
    full inv G
    """
    G = np.dot(B, B.T)
    G = G.round(12)
    u, v = np.linalg.eigh(G)
    u = u.round(12)
    v = v.round(12)
    # print("G eigenvals", u)

    # w = []

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-11:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.outer(vi/ui, vi)
            vv = vv.round(12)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(12)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B) # MxN
    GB = GB.round(12)
    
    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB

def transform(B):

    """
    full inv G
    """
    G = np.dot(B, B.T)
    G = G.round(12)
    u, v = np.linalg.eigh(G)
    u = u.round(12)
    v = v.round(12)
    # print("G eigenvals", u)

    # w = []

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-11:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.outer(vi/ui, vi)
            vv = vv.round(12)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(12)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B) # MxN
    GB = GB.round(12)
    
    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB

def transform_v2(B):

    """
    just use diag of inv G
    """
    G = np.dot(B, B.T)
    G = G.round(12)
    u, v = np.linalg.eigh(G)
    u = u.round(12)
    v = v.round(12)
    # print("G eigenvals", u)

    # w = []

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-10:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.diag(vi*vi / ui)
            vv = vv.round(12)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(12)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B) # MxN
    GB = GB.round(12)
    
    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB


def transform_v3(B):

    """
    just use diag of G
    """
    G = np.dot(B, B.T)
    G = np.diag(np.diag(G))
    G = G.round(12)
    u, v = np.linalg.eigh(G)
    u = u.round(12)
    v = v.round(12)
    # print("G eigenvals", u)

    # w = []

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-10:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.diag(vi*vi / ui)
            vv = vv.round(12)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(12)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B) # MxN
    GB = GB.round(12)
    
    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB

def project_ics(B, H):

    H = np.array(H, dtype=float)
    H = H.round(12)
    
    GB = transform(B)
    Q = np.dot(np.dot(GB, H), GB.T)
    Q = Q.round(12)
    return Q

def project_gradient(B, gx):

    gx = np.array(gx, dtype=float)
    gx = gx.round(12)
    
    GB = transform(B)
    gq = np.dot(GB, gx)
    gq = gq.round(12)
    return gq


def subtract_projection(Bnb, H_qm_x, H_mm_x):

    if not len(Bnb):
        return H_qm_x
    H_mm_nb = project_ics(Bnb, H_mm_x)

    H_qm_x = np.array(H_qm_x, dtype=float)
    H_qm_x = H_qm_x.round(12)

    H_mm_x = np.dot(Bnb.T, np.dot(H_mm_nb, Bnb))
    H_mm_x = np.round(H_mm_x, 12)

    H_qm_nb = project_ics(Bnb, H_qm_x)
    H_qm_nb = H_qm_nb.round(12)

    H_qm_nbx = np.dot(Bnb.T, np.dot(H_qm_nb, Bnb))
    H_qm_nbx = np.round(H_qm_nbx, 12)

    H = H_qm_x 
    # H = H_qm_x - H_qm_nbx - H_mm_x
    # H = H_qm_x - H_mm_x
    return H


def bmatrix_blocked(ba, t, bat):
    """
    """
    gmats = []
    M = 2
    sz = len(ba) + len(t)
    N = len(ba[0])
    
    bat = [x for x in ba]
    bat.extend([x for x in t])
    bat = np.array(bat)

    ba_t = []
    for v in ba:
        ba_t.append(v)
    for v in t:
        ba_t.append([0 for _ in range(N)])
    ba_t = np.array(ba_t)

    t_ba = []
    for v in ba:
        t_ba.append([0 for _ in range(N)])
    for v in t:
        t_ba.append(v)
    t_ba = np.array(t_ba)

    def invert(m):
        n = None
        U, V = np.linalg.eigh(m)
        for u, v in zip(U, V.T):
            if u > 1e-10:
                if n is None:
                    n = np.outer(v, v/u)
                else:
                    n += np.outer(v, v/u)
        return n

    g_ba = invert(np.dot(ba_t, ba_t.T))
    g_t = invert(np.dot(t_ba, t_ba.T))
    g_bat = invert(np.dot(bat, bat.T))

    G = np.vstack((g_bat[0], g_ba[0], g_t[0], [[1.0 for _ in range(len(g_bat[0]))] for _ in range(len(g_bat) + len(g_ba) + len(g_t) - 3)]))
    # for n in len(g_bat) + len(g_ba) + len(g_t) - 3:
    # G = np.vstack((g_bat[0], g_ba[0], g_t[0]))

    # G = g_bat - g_ba - g_t
    np.save("debug.npy", (g_bat, g_ba, g_t))
    g_bat_p = vibration.gram_schmidt(G).T
    # g_bat_p = vibration.gram_schmidt(g_bat).T

    return g_bat_p


# def collect(psys, m, ic_fcs):

#     fc_clusters = {}

#     m_lbls = psys.models[m].labels[0]
#     new_k = {}

#     for ic, assn in m_lbls.items():
#         fc_k = assn['k']
#         if fc_k not in fc_clusters:
#             fc_clusters[fc_k] = []
#         fc_clusters[fc_k].append(ic_fcs[ic])
    
#     for lbl, values in fc_clusters.items():
#         k_key = (m, 'k', lbl, 0)
#         new_k[k_key] = sum(values)/len(values)
    
#     return new_k

# def objective(alpha, ic_qm_fcs, ics, B, csys, psys):

#     # step 3
#     # cluster the BA fcs according to FF labels
#     # psys already labeled
#     qube = {
#         (0, 'k', 'b1', 0): 229.998,
#         (0, 'k', 'b2', 0): 353.322,
#         (1, 'k', 'a1', 0): 46,
#         (1, 'k', 'a2', 0): 38,
#     }
    

#     new_k = {}
#     bond_sel = collect(psys, 0, ic_qm_fcs)
#     new_k.update(bond_sel)
#     angle_sel = collect(psys, 1, ic_qm_fcs)
#     new_k.update(angle_sel)

#     if alpha == 0.0:
#         alpha=1.0
#         new_k = qube
#         print("Using qube values")

#     print("Setting FF k values")
#     for k, v in new_k.items():
#         mm.chemical_system_set_value(csys, k, alpha*v)
#         print(k, alpha*v)

#     psys = mm.chemical_system_to_physical_system(csys, psys.models[0].positions, ref=psys, reuse=[2,3,4,5])

#     # step 4
#     # calculate the MM hessian and get the IC
#     hess_mm = objectives.physical_system_hessian(psys, csys)
#     hess_mm_ic = project_ics(B, hess_mm)
#     hess_mm_ic /= 4.184
#     ic_mm_fcs = dict(zip(ics, np.diag(hess_mm_ic)))
#     print("MM FCs (kcal/A/A)")
#     for k,v in ic_mm_fcs.items():
#         print(f"{str(k):20s} {v:15.8g} diff= {v-ic_qm_fcs[k]:15.8g}")
    
#     new_k_mm = {}
#     new_k_mm.update(collect(psys, 0, ic_mm_fcs))
#     new_k_mm.update(collect(psys, 1, ic_mm_fcs))
    
#     dx = arrays.array_difference([*new_k.values()], [*new_k_mm.values()])
#     x2 = arrays.array_inner_product(dx, dx)

#     return x2*1e-10/len(dx)**.5


def remove_tr(H):
    """
    be lazy and just remove the first 6
    """
    u, v = np.linalg.eigh(H)

    u = np.atleast_2d(u[6:])
    v = v[:,6:]
    return np.dot(v*u, v.T)


def project_onto_ics_from_data(psys, fc_map):

    ic_msm_fcs = {}

    d = fc_map


    for ic, k in d["bonds"].items():
        bond = geometry.bond(tuple(int(x) for x in ic.split()))
        if k < 0:
            k = 0.0
        elif k > 4000:
            k = 4000
        ic_msm_fcs[bond] = k

    for ic, k in d["angles"].items():
        if k < 0:
            k = 0.0
        elif k > 4000:
            k = 4000
        angle = geometry.angle(tuple(int(x) for x in ic.split()))
        ic_msm_fcs[angle] = k
    
    pos = psys.models[0].positions[0]
    return assignments.graph_assignment(pos.smiles, ic_msm_fcs, pos.graph)

        

def hessian_project_onto_ics(csys, psys: mm.physical_model, hess_qm, verbose=False, B=None, shm=None) -> dict:

    # hess_qm in kcal/A/A

    # step 1
    # get the qm hessian
    # given as arg

    # step 2
    # project hessian onto ICs and return the diag
    pos = copy.deepcopy(psys.models[0].positions[0])
    g = pos.graph

    xyz = np.vstack([x[0] for x in pos.selections.values()], dtype=float)
    xyz = xyz.round(12)
    sym = graphs.graph_symbols(pos.graph)
    mass = np.array([[vibration.mass_table[sym[n]]]*3 for n in sym])

    sym = list(sym.values())
    remove1_3 = True
    torsions = False
    outofplanes = False
    pairs = False

    hess_qm_freq, hess_qm_modes = vibration.hessian_modes(hess_qm, sym, xyz, mass, 0, remove=0, stdtr=True, verbose=False)
    omega = np.round(hess_qm_freq, 12)
    # omega_qm = omega
    if verbose:
        print("Ref Hessian Frequencies (cm-1):")
        print(omega)

    if B is None:
        ics, B = bmatrix(pos, torsions=torsions, outofplanes=outofplanes, pairs=pairs,remove1_3=remove1_3) # pairs=False, torsions=False, outofplanes=False)
    else:
        ics, B = B
        B = np.array(B)
    hess_qm_ic = project_ics(B, hess_qm)

    hess_qm_ic = np.array(np.diag(hess_qm_ic))

    if (hess_qm_ic < 0).any():
        print("Warning, negative force constants found")
    if (hess_qm_ic > 4000).any():
        print("Warning, large force constants found")

    hess_qm_ic[hess_qm_ic < 0] = 0.0
    hess_qm_ic[hess_qm_ic > 4000] = 4000

    ic_qm_fcs = dict(zip(ics, hess_qm_ic))

    print("Projected MM Fcs")
    pprint(ic_qm_fcs, sort_dicts=False)

    # u, v = np.linalg.eigh(hess_qm_ic)
    # with np.printoptions(precision=5, linewidth=160, formatter={'float_kind': "{:12.5f}".format}):
    #     print("QM IC eigvals:\n                      ", np.array(u) )
    #     # print("QM IC eigvals:\n", u )
    #     # v[1,:] = 0
    #     # v[6,:] = 0
    #     # v[8,:] = 0
    #     print(v)
    #     print(ics)
    #     pprint(dict(zip([f"{str(x):12s}" for x in ics], v)))
    #     print("Projected")
    #     print(dict(zip(ics, (v[:,3:]*np.reshape(u[3:], (1, -1))).sum(axis=1))))
    #     print("Project transformed")
    #     vv = np.dot(v*u, v.T)
    #     pprint(dict(zip([f"{str(x):12s}" for x in ics], np.diag(vv))))

    # hess_mm = objectives.physical_system_hessian(psys, csys, h=1e-7)
    # hess_mm = np.array(hess_mm)/4.184 * au2kcal

    # if verbose:
    #     print("Computed Hessian Frequencies (cm-1) (pre-proj):")
    #     hess_mm_freq, hess_mm_modes = vibration.hessian_modes(hess_mm, sym, xyz, mass, 1, remove=0, stdtr=True)
    #     omega = np.round(hess_mm_freq, 12)
    #     print(omega)

    # print("Computed Hessian Frequencies (cm-1) (post-proj):")
    # hess_mm_freq, hess_mm_modes, hess_mm_clean = vibration.hessian_modes(hess_mm, sym, xyz, mass, 2, remove=6, stdtr=True)
    # omega = np.round(hess_mm_freq, 12)
    # print(omega)

    # posmm = pos 
    # Bmm = B
    # posmm = optimizers_scipy.optimize_positions_scipy(csys, psys, tol=1e-10) 
    # posmm = optimizers_openmm.optimize_positions_openmm(csys, psys) 
    # ics, Bmm = bmatrix(posmm, torsions=torsions, pairs=pairs,remove1_3=remove1_3) # pairs=False, torsions=False, outofplanes=False)

    # xyz = np.vstack([x[0] for x in posmm.selections.values()], dtype=float)
    # xyz = xyz.round(12)
    # args, keys = objectives.array_flatten_assignment(posmm.selections)
    # hess_mm = objectives.array_geom_hessian(args, keys, csys, psys, h=1e-7)
    # hess_mm = objectives.physical_system_hessian(psys, csys)
    # hess_mm = np.array(hess_mm)/4.184

    # gradx_mm= objectives.array_geom_gradient(args, keys, csys, psys)
    # gradic_mm = project_gradient(Bmm, gradx_mm)
    # gradic_mm /= 4.148
    # ic_grad = dict(zip(ics, gradic_mm))
    # if verbose:
    #     print("IC gradients:")
    #     print(ic_grad)
    #     print("IC gradient sum:", sum(gradic_mm))

    # b2ics, B2 = b2matrix(posmm, remove1_3=True, pairs=pairs)

    # plug in atom index to get matrix index
    # id_map = {v: k for k, v in enumerate(pos.graph.nodes)}

    # hgx = np.zeros_like(hess_mm)

    # for ic, b2 in zip(b2ics, B2):
    #     for ai, a in enumerate(ic):
    #         a = id_map[a] 
    #         for bi, b in enumerate(ic):
    #             b = id_map[b]
    #             for i in range(3):
    #                 for j in range(3):
    #                     b2ab = b2[0][3*ai + i][3*bi + j]
    #                     hgx[3*a + i][3*b + j] += ic_grad[ic] * b2ab
    
    # hess_mm_ic_corrected = project_ics(B, hess_mm + hgx)
    # hgx *= 10.0
    # if verbose:
    #     hess_mm_ic_uncorrected = project_ics(B, hess_mm)
    #     print("Uncorrected MM ICs:")
    #     print(dict(zip(ics, np.diag(hess_mm_ic_uncorrected))))
    # if verbose:
    #     print("corrected MM ICs:")
    #     print(dict(zip(ics, np.diag(hess_mm_ic_corrected))))


    # if verbose:
    #     A = transform(B)
    #     P = np.dot(B, A.T)
    #     hess_mm_ic_corrected = np.dot(P, np.dot(hess_mm_ic_corrected, P))
    #     for i, ic in enumerate(ics):
    #         k0 = hess_mm_ic_uncorrected[i][i]
    #         k1 = hess_mm_ic_corrected[i][i]
    #         print(f"{str(ic):15s} New: {k1:15.7f} Ref: {k0:15.7f} Diff{k1-k0:15.7f}")

     
    # ic_mm_fcs = np.diag(hess_mm_ic_corrected)

    # print("Computed Hessian Frequencies (cm-1) (pre-proj):")
    # hess_mm_freq, hess_mm_modes, hess_mm_clean = vibration.hessian_modes(hess_mm, sym, xyz, mass, 3, remove=0, stdtr=True)
    # omega = np.round(hess_mm_freq, 12)
    # print(omega)
    # hess_mm_freq, hess_mm_modes = vibration.hessian_modes(hess_mm +hgx ,  sym, xyz, mass, 4, remove=6, stdtr=True)
    # if verbose:
    #     print("Computed Hessian Frequencies (cm-1) with correction (min pre-proj):")
    #     print("hess_mm_freq")
    #     print(hess_mm_freq)
    # omega = np.round(hess_mm_freq, 1)
    # f = list(omega[6:])
    # # f = f[:-2] + [f[-1], f[-2]]
    # f = np.array(f)
    # if verbose:
    #     for w in f:
    #         print(f"{w:6.1f}")
    #     print(f"{np.sum(np.abs(omega_qm[6:] - f)):6.1f}\n{np.sum(np.abs(omega_qm[6:] - f))/(len(f)):6.1f}")

    # if verbose:
    #     print("Difference in Computed Hessian Frequencies (cm-1) with correction (min pre-proj):")
    #     omega = np.round(hess_mm_freq, 1)
    #     f = list(omega[6:])
    #     # f = f[:-2] + [f[-1], f[-2]]
    #     f = np.array(f)
    #     for w, w0 in zip(f, omega_qm[6:]):
    #         print(f"{w:6.1f} - {w0:6.1f} = {w-w0:6.1f}")
    #     # print(f"{np.sum(np.abs(omega_qm[6:] - f)):6.1f}\n{np.sum(np.abs(omega_qm[6:] - f))/(len(f)):6.1f}")

    # A = transform(B)
    # P = np.dot(B.T, A)
    # hess_mm2 = np.dot(P, np.dot(hess_mm + hgx, P))
    # # print("Computed Hessian Frequencies (cm-1) (pre-proj):")
    # # hess_mm_freq, hess_mm_modes, hess_mm_clean = vibration.hessian_modes(hess_mm, sym, xyz, mass, 3, remove=0, stdtr=True)
    # # omega = np.round(hess_mm_freq, 12)
    # # print(omega)
    # hess_mm_freq, hess_mm_modes = vibration.hessian_modes(hess_mm2, sym, xyz, mass, 4, remove=6, stdtr=True)
    # if verbose:
    #     print("Computed Hessian Frequencies (cm-1) with sanitzed correction (min pre-proj):")
    #     print("hess_mm_freq")
    #     print(hess_mm_freq)
    #     omega = np.round(hess_mm_freq, 1)
    #     f = list(omega[6:])
    #     # f = f[:-2] + [f[-1], f[-2]]
    #     f = np.array(f)
    #     for w in f:
    #         print(f"{w:6.1f}")
    #     print(f"{np.sum(np.abs(omega_qm[6:] - f)):6.1f}\n{np.sum(np.abs(omega_qm[6:] - f))/(len(f)):6.1f}")

    # if verbose:
    #     print("Difference in Computed Hessian Frequencies (cm-1) with correction (min pre-proj):")
    #     omega = np.round(hess_mm_freq, 1)
    #     f = list(omega[6:])
    #     # f = f[:-2] + [f[-1], f[-2]]
    #     f = np.array(f)
    #     for w, w0 in zip(f, omega_qm[6:]):
    #         print(f"{w:6.1f} - {w0:6.1f} = {w-w0:6.1f}")
    #     # print(f"{np.sum(np.abs(omega_qm[6:] - f)):6.1f}\n{np.sum(np.abs(omega_qm[6:] - f))/(len(f)):6.1f}")

    # if verbose:
    #     hess_mm_freq, hess_mm_modes = vibration.hessian_modes(hess_mm , sym, xyz, mass, 4, remove=6, stdtr=True)
    #     print("Computed Hessian Frequencies (cm-1) without (min pre-proj):")
    #     print("hess_mm_freq")
    #     print(hess_mm_freq)
    #     omega = np.round(hess_mm_freq, 1)
    #     f = list(omega[6:])
    #     # f = f[:-2] + [f[-1], f[-2]]
    #     f = np.array(f)
    #     for w in f:
    #         print(f"{w:6.1f}")
    #     print(f"{np.sum(np.abs(omega_qm[6:] - f)):6.1f}\n{np.sum(np.abs(omega_qm[6:] - f))/(len(f)):6.1f}")

    # if verbose:
    #     print("Difference in Computed Hessian Frequencies (cm-1) without (min pre-proj):")
    #     omega = np.round(hess_mm_freq, 1)
    #     f = list(omega[6:])
    #     # f = f[:-2] + [f[-1], f[-2]]
    #     f = np.array(f)
    #     for w, w0 in zip(f, omega_qm[6:]):
    #         print(f"{w:6.1f} - {w0:6.1f} = {w-w0:6.1f}")
    #     # print(f"{np.sum(np.abs(omega_qm[6:] - f)):6.1f}\n{np.sum(np.abs(omega_qm[6:] - f))/(len(f)):6.1f}")

    # if verbose:
    #     print("Modes:")
    #     with np.printoptions(precision=5, linewidth=160, formatter={'float_kind': "{:12.5f}".format}):
    #         print(omega)
    #         for mode in hess_mm_modes:
    #             print(f"{mode}")
    # # hess_mm_freq, hess_mm_modes, hess_mm_clean = vibration.hessian_modes(hess_mm, sym, xyz, mass, 5, remove=6, stdtr=True)
    # # print("Computed Hessian Frequencies (cm-1) (min post-proj):")
    # # omega = np.round(hess_mm_freq, 12)
    # # print(omega)

    # if False:
    #     Bnb = [B[i] for i, ic in enumerate(ics) if ic in graphs.graph_pairs(g) and i < len(B)] 
    #     # Bnb = [B[i] for i, ic in enumerate(ics) if ic in graphs.graph_torsions(g) and ic in graphs.angles(g) and  i < len(B)] 

    #     if len(Bnb):
    #         Bnb = np.vstack(Bnb)

    #     ics_nonb, B_nonb = bmatrix(pos, torsions=True, bonds=True, pairs=False,remove1_3=True, outofplanes=False) # pairs=False, torsions=False, outofplanes=False)
    #     # ics_nonb, B_nonb = bmatrix(pos, torsions=True, bonds=False, angles=True, pairs=True, remove1_3=False) # pairs=False, torsions=False, outofplanes=False)
    #     hess_qm_nommnb = subtract_projection(Bnb, hess_qm, hess_mm)

    #     hess_qm_nonb = project_ics(B_nonb, hess_qm_nommnb)
    #     print("QM IC FCs after removing MM NB")
    #     ic_qm_fcs = dict(zip(ics_nonb, np.diag(hess_qm_nonb)))
    #     pprint(ic_qm_fcs, sort_dicts=False)
    
    # # hess_mm_ic = project_ics(Bmm, hess_mm)
    # hess_mm_ic = hess_mm_ic_corrected

    # u, v = np.linalg.eigh(hess_mm_ic)
    # if verbose:
    #     with np.printoptions(precision=5, linewidth=160, formatter={'float_kind': "{:12.5f}".format}):
    #         print("MM IC eigvals:")
    #         print(" "*15, u)
    #         for ic, mode in zip(ics, v):
    #             ic = str(ic)
    #             print(f"{ic:15s} {mode}")

    #         # pprint(dict(zip([f"{str(x):12s}" for x in ics], v)))
    #         # print(v[0,2:],u[2:])
    #         # print((v[0,2:]*u[2:]).sum())
    #         # v[1,:] = 0
    #         # v[6,:] = 0
    #         # v[8,:] = 0
    #         # print(v)
    #         # print(ics)
    #         pprint(dict(zip([f"{str(x):12s}" for x in ics], v)), sort_dicts=False)
    #         print("Projected")
    #         print(dict(zip(ics, (v[:,:]*np.reshape(u[:], (1, -1))).sum(axis=1))))
    #         print("Project transformed")
    #         vv = np.dot(v*u, v.T)
    #         # ic_qm_fcs = dict(zip(ics_nonb, np.diag(vv)))
    #         pprint(dict(zip([f"{str(x):12s}" for x in ics], np.diag(vv))), sort_dicts=False)

    # # hess_mm_ic /= 4.184
    # ic_mm_fcs = dict(zip(ics, np.diag(hess_mm_ic)))

    # hmm = [*ic_mm_fcs.values()]
    # hqm = [*ic_qm_fcs.values()]
    # hqhm = arrays.array_inner_product(hmm, hqm)
    # hmhm = arrays.array_inner_product(hmm, hmm)
    # if verbose:
    #     print("Hessian FCs (kcal/A/A)")
    # dx = []
    # for k,v in ic_qm_fcs.items():
    #     d = ic_mm_fcs[k] - v
    #     dx.append(d)
    # x0 = arrays.array_inner_product(dx, dx)
    # if we set the objective to the rss and set deriv to 0 wrt alpha
    # we get this analytic form
    # alpha = hqhm/hmhm
    # alpha = 1.0 #hqhm/hmhm

    # hess_mm_freq, _, _ = vibration.hessian_modes(hess_mm*alpha, sym, xyz, mass, 3, remove=6, stdtr=True)
    # print("Computed Hessian Frequencies (cm-1) (post-scale):")
    # omega = np.round(hess_mm_freq, 12)
    # print(omega)

    # hess_mm_ic = project_ics(B, hess_mm_clean)
    # hess_mm_ic /= 4.184
    # ic_mm_fcs = dict(zip(ics, np.diag(hess_mm_ic)))

    # dx = []
    # for k in ic_mm_fcs:
    #     ic_mm_fcs[k] *= alpha
    #     d = ic_mm_fcs[k] - ic_qm_fcs[k]
    #     dx.append(d)
    #     if verbose:
    #         print(f"{str(k):25s} MM0={ic_mm_fcs[k]/alpha:12.6g} MM={ic_mm_fcs[k]:12.6g} QM={ic_qm_fcs[k]:12.6g} Diff={d:12.6g}")
    # x1 = arrays.array_inner_product(dx, dx)
    # if verbose:
    #     print(f"Alpha {alpha:12.6g} X2_0: {x0:12.6g} X2_1: {x1:12.6g} Diff: {x1-x0:12.6g}")

    # NOTE: returning the QM fcs! useful for penalties
    return assignments.graph_assignment(pos.smiles, ic_qm_fcs, pos.graph)


def hessian_frequencies(g, hess_mm, grad_mm, DL, ics, B, B2):

    sym = graphs.graph_symbols(g)
    mass = np.array([[vibration.mass_table[sym[n]]]*3 for n in sym])

    # ic_grad = dict(zip(ics, project_gradient(B, grad_mm)))

    # plug in atom index to get matrix index
    # id_map = {v: k for k, v in enumerate(g.nodes)}

    hgx = np.zeros_like(hess_mm)

    # for ic, b2 in zip(ics, B2):
    #     for ai, a in enumerate(ic):
    #         a = id_map[a]
    #         for bi, b in enumerate(ic):
    #             b = id_map[b]
    #             for i in range(3):
    #                 for j in range(3):
    #                     b2ab = b2[0][3*ai + i][3*bi + j]
    #                     hgx[3*a + i][3*b + j] += ic_grad[ic] * b2ab

    hess_mm_au = vibration.hessian_transform_mass_weighted(hess_mm + hgx, mass)
    hess_mm_freq = np.diag(np.dot(np.dot(DL.T, hess_mm_au), DL))
    hess_mm_freq = vibration.converteig(hess_mm_freq)
    hess_mm_freq = np.round(hess_mm_freq, 12)

    return hess_mm_freq
