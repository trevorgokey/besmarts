"""
besmarts.mechanics.hessians
"""

import copy
import numpy as np
import pprint

from besmarts.core import assignments
from besmarts.core import graphs
from besmarts.core import geometry
from besmarts.core import configs
from besmarts.mechanics import molecular_models as mm

from besmarts.mechanics import vibration

au2kcal = 627.51/(0.529177**2)

PRECISION = configs.precision


def transform(B):
    """
    full inv G
    """
    B = np.atleast_2d(np.asarray(B))
    G = np.dot(B, B.T)
    G = G.round(PRECISION)
    u, v = np.linalg.eigh(G)
    u = u.round(PRECISION)
    v = v.round(PRECISION)

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-11:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.outer(vi/ui, vi)
            vv = vv.round(PRECISION)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(PRECISION)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B)  # MxN
    GB = GB.round(PRECISION)
    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB


def transform(B):
    """
    full inv G
    """
    B = np.atleast_2d(np.asarray(B))
    G = np.atleast_2d(np.dot(B, B.T))
    G = G.round(PRECISION)
    u, v = np.linalg.eigh(G)
    u = u.round(PRECISION)
    v = v.round(PRECISION)
    # print("G eigenvals", u)

    # w = []

    c = 0.99
    s = 0
    N = sum(u)
    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        s += ui
        if s/N < c:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.outer(vi/ui, vi)
            vv = vv.round(PRECISION)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(PRECISION)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B)  # MxN
    GB = GB.round(PRECISION)

    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB


def transform_v2(B):
    """
    just use diag of inv G
    """
    G = np.dot(B, B.T)
    G = G.round(PRECISION)
    u, v = np.linalg.eigh(G)
    u = u.round(PRECISION)
    v = v.round(PRECISION)
    # print("G eigenvals", u)

    # w = []

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-10:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.diag(vi*vi / ui)
            vv = vv.round(PRECISION)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(PRECISION)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B)  # MxN
    GB = GB.round(PRECISION)
    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB


def transform_v3(B):
    """
    just use diag of G
    """
    G = np.dot(B, B.T)
    G = np.diag(np.diag(G))
    G = G.round(PRECISION)
    u, v = np.linalg.eigh(G)
    u = u.round(PRECISION)
    v = v.round(PRECISION)
    # print("G eigenvals", u)

    # w = []

    ginv = np.zeros_like(G)
    for ui, vi in zip(u[::-1], v.T[::-1]):
        if ui > 1e-10:
            # w.append(vi)
            # ginv += np.outer(vi/ui,vi)
            vv = np.diag(vi*vi / ui)
            vv = vv.round(PRECISION)
            # print(f"VV is u={ui}", vv.min(), vv.mean(), vv.max())
            ginv += vv

    ginv = ginv.round(PRECISION)
    # print("GINV is", ginv.min(), ginv.mean(), ginv.max())
    # w = np.array(w).T

    # print("G inverse")
    # print(ginv)

    GB = np.dot(ginv, B)  # MxN
    GB = GB.round(PRECISION)

    # print("GB is", GB.min(), GB.mean(), GB.max())
    return GB

def project_ics(B, H):

    H = np.array(H, dtype=float)
    H = H.round(PRECISION)

    GB = transform(B)
    Q = np.dot(np.dot(GB, H), GB.T)
    Q = Q.round(PRECISION)
    return Q


def project_gradient(B, gx, shm=None):

    gx = np.array(gx, dtype=float)
    gx = gx.round(PRECISION)

    GB = transform(B)
    gq = np.dot(GB, gx)
    gq = gq.round(PRECISION)
    return gq


def subtract_projection(Bnb, H_qm_x, H_mm_x):

    if not len(Bnb):
        return H_qm_x
    H_mm_nb = project_ics(Bnb, H_mm_x)

    H_qm_x = np.array(H_qm_x, dtype=float)
    H_qm_x = H_qm_x.round(PRECISION)

    H_mm_x = np.dot(Bnb.T, np.dot(H_mm_nb, Bnb))
    H_mm_x = np.round(H_mm_x, PRECISION)

    H_qm_nb = project_ics(Bnb, H_qm_x)
    H_qm_nb = H_qm_nb.round(PRECISION)

    H_qm_nbx = np.dot(Bnb.T, np.dot(H_qm_nb, Bnb))
    H_qm_nbx = np.round(H_qm_nbx, PRECISION)

    H = H_qm_x

    return H


def remove_tr(H):
    """
    be lazy and just remove the first 6
    """
    u, v = np.linalg.eigh(H)

    u = np.atleast_2d(u[6:])
    v = v[:, 6:]
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


def hessian_project_onto_ics(
    csys,
    psys: mm.physical_model,
    hess_qm,
    verbose=False,
    B=None,
    shm=None
) -> dict:

    pos = psys.models[0].positions

    xyz = np.vstack([x[0] for posi in pos for x in posi.selections.values()], dtype=float)
    xyz = xyz.round(PRECISION)
    sym = [s for posi in pos for s in graphs.graph_symbols(posi.graph).values()]
    mass = np.array([[vibration.mass_table[s]]*3 for s in sym])

    # sym = list(sym.values())
    remove1_3 = True
    torsions = True
    outofplanes = True
    pairs = False

    hess_qm_freq, hess_qm_modes = vibration.hessian_modes(
        hess_qm,
        sym,
        xyz,
        mass,
        0,
        remove=0,
        stdtr=True,
        verbose=False
    )
    omega = np.round(hess_qm_freq, PRECISION)

    if verbose:
        print("Ref Hessian Frequencies (cm-1):")
        print(omega)

    if B is None:
        ics, B = assignments.bmatrix(
            pos,
            torsions=torsions,
            outofplanes=outofplanes,
            pairs=pairs,
            remove1_3=remove1_3
        )
    else:
        ics, B = B
        B = np.array(B)

    ic_qm_fcs = {}
    if len(B.shape) > 1 and B.shape[0] > 0 and B.shape[1] > 0:
        hess_qm_ic = project_ics(B, hess_qm)

        hess_qm_ic = np.array(np.diag(hess_qm_ic))

        if (hess_qm_ic < 0).any() and verbose:
            print("Warning, negative force constants found")
        if (hess_qm_ic > 4000).any() and verbose:
            print("Warning, large force constants found")

        hess_qm_ic[hess_qm_ic < 0] = 0.0
        hess_qm_ic[hess_qm_ic > 4000] = 4000

        ic_qm_fcs = dict(zip(ics, hess_qm_ic))

        if verbose:
            print("Projected MM Fcs")
            pprint.pprint(ic_qm_fcs, sort_dicts=False)

    return ic_qm_fcs
    # return assignments.graph_assignment(pos.smiles, ic_qm_fcs, pos.graph)


def hessian_transform(mass, hess_mm, grad_mm, DL, ics, B, B2):

    # sym = graphs.graph_symbols(g)
    # mass = np.array([[vibration.mass_table[sym[n]]]*3 for n in sym])
    mass = np.array(mass)

    hgx = np.zeros_like(hess_mm)

    # needed if we do analytic Hessian (TODO)
    # ic_grad = dict(zip(ics, project_gradient(B, grad_mm)))

    # plug in atom index to get matrix index
    # id_map = {v: k for k, v in enumerate(g.nodes)}

    # for ic, b2 in zip(ics, B2):
    #     for ai, a in enumerate(ic):
    #         a = id_map[a]
    #         for bi, b in enumerate(ic):
    #             b = id_map[b]
    #             for i in range(3):
    #                 for j in range(3):
    #                     b2ab = b2[0][3*ai + i][3*bi + j]
    #                     hgx[3*a + i][3*b + j] += ic_grad[ic] * b2ab

    hess_mm_au = vibration.hessian_transform_mass_weighted(hess_mm - hgx, mass)
    hess_mm_freq = np.dot(np.dot(DL.T, hess_mm_au), DL)
    hess_mm_freq = vibration.converteig(hess_mm_freq)
    hess_mm_freq = np.round(hess_mm_freq, PRECISION)

    return hess_mm_freq


def hessian_frequencies(mass, hess_mm, grad_mm, DL, ics, B, B2):

    return np.diag(hessian_transform(mass, hess_mm, grad_mm, DL, ics, B, B2))


