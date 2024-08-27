"""
besmarts.sample

Sample molecules based on internal coordinate SMARTS
"""

from besmarts.core import mapper
from besmarts.core import graphs
import random


def sample_ic(cfg, icgl, n=1):
    """
    """
    d = {}
    for gi, g in enumerate(icgl):
        print(f"{gi:5d}/{len(icgl)}", end="\r")
        mapper.mapper_smarts_extend(cfg, g)
        for ic in g:
            d.setdefault(ic, []).append(gi)
    print()

    ret = set()
    for ic, gids in d.items():
        N = len(gids)
        if N == n:
            ret.update(gids)
        elif len(gids) <= n:
            # warning
            print(f"Warning, sample pool size is {N} but sample size is {n}")
            # print(graphs.structure_print(ic))
            ret.update(gids)
        else:
            ret.update(random.sample(gids, k=n))

    return ret


def sample_bonds(cfg, gl, n=1):
    icgl = [graphs.graph_to_structure_bonds(g) for g in gl]
    ret = set()
    ret.update(sample_ic(cfg, icgl, n=n))
    return ret


def sample_angles(cfg, gl, n=1):
    icgl = [graphs.graph_to_structure_angles(g) for g in gl]
    ret = set()
    ret.update(sample_ic(cfg, icgl, n=n))
    return ret


def sample_torsions(cfg, gl, n=1):
    icgl = [graphs.graph_to_structure_torsions(g) for g in gl]
    ret = set()
    ret.update(sample_ic(cfg, icgl, n=n))
    return ret


def sample_outofplanes(cfg, gl, n=1):
    icgl = [graphs.graph_to_structure_outofplanes(g) for g in gl]
    ret = set()
    ret.update(sample_ic(cfg, icgl, n=n))
    return ret


def sample_valence(cfg, gl, n=1):
    ret = set()

    fn_list = sample_bonds, sample_angles, sample_torsions, sample_outofplanes
    for fn in fn_list:
        x = fn(cfg, gl, n=n)
        ret.update(x)
        print(fn, len(x), len(gl), "Returning", len(ret))

    return ret
