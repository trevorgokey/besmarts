
"""

"""

import sys
import os
import pickle
from typing import Dict, List, Tuple
from besmarts.mechanics import smirnoff_models
from besmarts.mechanics import optimizers_scipy
from besmarts.mechanics import fits
from besmarts.core import graphs
from besmarts.core import topology
from besmarts.core import perception
from besmarts.core import arrays
from besmarts.core import assignments
from besmarts.core import primitives
from besmarts.assign import hierarchy_assign_rdkit
from besmarts.codecs import codec_rdkit

def graph_db_insert_sdf_files(gdb, gcd, fname_pos, fname_grad=None, fname_hessian=None, energy=None):

    pos = None
    gx = None
    hx = None

    pos, extras = gcd.sdf_decode(fname_pos)
    smi = pos.smiles

    ene = extras.pop('qm_energy (kcal/mol)')
    if ene is not None and energy is None:
        energy = ene * 4.184

    if fname_grad:
        with open(fname_grad) as f:
            xyzdata = f.read()
        sym, xyz = parse_xyz(xyzdata)
        gx = xyz_to_graph_assignment(gcd, smi, xyz)

    if fname_hessian:
        with open(fname_grad) as f:
            xyzdata = f.read()
        sym, xyz = parse_xyz(xyzdata)
        gx = xyz_to_graph_assignment(gcd, smi, xyz)

    eid, gid = assignments.graph_db_add_single_molecule_state(gdb, pos, gradients=gx, hessian=hx, energy=energy)
    return eid, gid

def print_xyz(pos, comment="") -> List[str]:
    lines = []
    lines.append(str(len(pos.selections)))
    lines.append(comment)
    for ic, xyz in pos.selections.items():
        n = pos.graph.nodes[ic[0]]
        sym = primitives.element_tr[str(n.primitives['element'].on()[0])]
        x, y, z = xyz[0][:3]
        lines.append(f"{sym:8s} {x:.6f} {y:.6f} {z:.6f}")
    return lines

def run(ff_fname, sdf_fname):
    # build the dataset and input ff
    gdb = assignments.graph_db()
    gcd = codec_rdkit.graph_codec_rdkit()

    labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
    pcp = perception.perception_model(gcd, labeler)
    csys = smirnoff_models.smirnoff_load(ff_fname, pcp)

    
    eid, gid = graph_db_insert_sdf_files(gdb, gcd, sdf_fname)

    if os.path.exists("psys.p"):
        print("Loading cached psys.p")
        with open("psys.p", 'rb') as f:
            psys = pickle.load(f)

    psys = fits.gdb_to_physical_systems(gdb, csys)[eid]
    # psys = psys[6878]
    print("\n".join(print_xyz(psys.models[0].positions[0])))

    pos = optimizers_scipy.optimize_positions_scipy(csys, psys)

    print("\n".join(print_xyz(pos)))

if __name__ == "__main__":
    ff = sys.argv[1]
    sdf = sys.argv[2]
    run(ff, sdf)
