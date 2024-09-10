"""
besmarts.md
"""

import sys
from besmarts.mechanics import smirnoff_models
from besmarts.mechanics import fits
from besmarts.core import assignments
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import optimizers_openmm
from besmarts.codecs import codec_rdkit
from besmarts.assign import hierarchy_assign_rdkit
from besmarts.core import perception


def run(csys, xyz, smi, ts=.002, steps=1000, temperature=278.15):

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
        pos = assignments.xyz_to_graph_assignment(gcd, smi, xyzdata)
    else:
        pos = xyz

    psys = mm.chemical_system_to_physical_system(csys, [pos])
    optimizers_openmm.molecular_dynamics(
        psys,
        float(ts),
        int(steps),
        temperature=float(temperature)
    )


if __name__ == "__main__":
    run(*sys.argv[1:])
