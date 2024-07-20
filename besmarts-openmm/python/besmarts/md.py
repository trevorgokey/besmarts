"""
besmarts.md
"""

import sys
from besmarts.mechanics import smirnoff_models
from besmarts.mechanics import fits
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import optimizers_openmm
from besmarts.codecs import codec_rdkit
from besmarts.assign import hierarchy_assign_rdkit
from besmarts.core import perception


def run(csys, xyz, smi, ts=.002, steps=1000):

    if type(csys) is str:
        gcd = codec_rdkit.graph_codec_rdkit()
        labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
        pcp = perception.perception_model(gcd, labeler)
        csys = smirnoff_models.smirnoff_load(csys, pcp)
    else:
        gcd = csys.perception.gcd

    if type(xyz) is str:
        with open(xyz) as f:
            xyzdata = f.read()

        sym, xyz = optimizers_openmm.parse_xyz(xyzdata)
        pos = fits.xyz_to_graph_assignment(gcd, smi, xyz)
    else:
        pos = xyz

    psys = mm.chemical_system_to_physical_system(csys, [pos])
    optimizers_openmm.molecular_dynamics(psys, ts, steps)


if __name__ == "__main__":
    run(*sys.argv[1:])
