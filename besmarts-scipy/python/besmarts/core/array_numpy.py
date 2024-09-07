
"""
besmarts.core.geometry_numpy.py
"""

from typing import Dict, List, Tuple
import numpy as np

from besmarts.core import topology
from besmarts.core import assignments

def b_matrix(pos, internals: Dict[topology.structure_topology, List]):
    pass

def g_matrix(pos, internals: Dict[topology.structure_topology, List]):
    pass

def dlc_matrix(pos, internals: Dict[topology.structure_topology, List], cutoff=1e-5):
    pass

def gradient_cartesian_to_internal_array(pos, gx, internals: Dict[topology.structure_topology, List]):
    pass

def gradient_cartesian_to_internal_matrix(pos, gx, internals: Dict[topology.structure_topology, List]):
    pass

def ginverse(G):

    ginv = np.zeros_like(G)
    u, v = np.linalg.eigh(G)
    for ui, vi in zip(u, v):
        if ui > 1e-1:
            ginv += np.outer(vi/ui,vi)
    return ginv

def dlcmatrix(gx, B):
    B = np.array(B)
    G = np.dot(B, B.T)
    ginv = ginverse(G)
    return ginv

def bmatrix_dense(valence, N):

    # N = len(valence.selections)
    conf = 0
    B = []

    for ic, confs in valence.selections.items():
        brow = np.zeros(3*N)
        for i, nid in enumerate(ic):
            d = (nid-1)*3
            brow[d:d+3] += confs[conf][i]
        B.append(brow)
    return dict(zip(valence.selections.keys(), B))

def bmatrix_ic(pos, jac):
    ics = jac(pos)
    return bmatrix_dense(ics, len(pos.selections))

def bmatrix_bonds(pos):
    return bmatrix_ic(pos, assignments.graph_assignment_jacobian_bonds)

def bmatrix_angles(pos):
    return bmatrix_ic(pos, assignments.graph_assignment_jacobian_angles)

def bmatrix_torsions(pos):
    return bmatrix_ic(pos, assignments.graph_assignment_jacobian_torsions)

def bmatrix_outofplanes(pos):
    return bmatrix_ic(pos, assignments.graph_assignment_jacobian_outofplanes)

def bmatrix_pairs(pos):
    return bmatrix_ic(pos, assignments.graph_assignment_jacobian_pairs)

def bmatrix(pos, bonds=True, angles=True, torsions=True, outofplanes=True, pairs=True) -> dict:

    B = {}

    if bonds:
        B.update(bmatrix_bonds(pos))

    if angles:
        B.update(bmatrix_angles(pos))

    if torsions:
        B.update(bmatrix_torsions(pos))

    if outofplanes:
        B.update(bmatrix_outofplanes(pos))

    if pairs:
        B.update(bmatrix_pairs(pos))

    return B

def transform(B):

    B = np.array(B)

    G = np.dot(B, B.T)
    u, v = np.linalg.eigh(G)

    ginv = np.zeros_like(G)
    for ui, vi in zip(u, v):
        if ui > 1e-1:
            ginv += np.outer(vi/ui,vi)

    GB = np.dot(ginv, B)
    
    return GB.tolist()

def transform_dlcs(B):

    """
    just use diag of inv G
    """
    G = np.dot(B, B.T)
    G = G.round(12)
    u, v = np.linalg.eigh(G)
    u = u.round(12)
    v = v.round(12)

    return u, v

def transform_v2(B):

    """
    just use diag of inv G
    """
    B = np.array(B)
    G = np.dot(B, B.T)
    u, v = np.linalg.eigh(G)

    ginv = np.zeros_like(G)
    for ui, vi in zip(u, v):
        if ui > 1e-1:
            ginv += np.diag(vi*vi / ui)

    GB = np.dot(ginv, B) # MxN
    
    return GB.tolist()

def transform_v3(B):

    """
    just use diag of G
    """
    B = np.array(B)
    G = np.dot(B, B.T)
    G = np.diag(np.diag(G))
    u, v = np.linalg.eigh(G)

    ginv = np.zeros_like(G)
    for ui, vi in zip(u, v):
        if ui > 1e-1:
            ginv += np.outer(vi/ui,vi)

    GB = np.dot(ginv, B) # MxN
    
    return GB.tolist()

def project_hessian(B, H):

    GB = transform_v2(B)
    return np.dot(np.dot(GB, H), GB.T)

def project_gradient(B, H):

    GB = transform_v2(B)
    return np.dot(np.dot(GB, H), GB.T)

def dlcmatrix_project_gradients(pos, gx, eps=1e-2):

    B = bmatrix(pos)
    B = np.vstack(list(B.values()))
    B = B.round(12)
    gx = np.atleast_2d(gx).T
    gx = gx.round(12)

    du, dv = transform_dlcs(B)

    gq = []
    for ui, vi in zip(du, dv.T):
        if ui > eps:
            gb = np.dot(np.outer(vi/ui, vi), B)
            gb = gb.round(12)
            gb = np.dot(gb, gx).reshape(-1)
            gb = gb.round(12)
            gq.extend(gb.tolist())
    return gq

