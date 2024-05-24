
"""
"""

from typing import Dict, List, Tuple
from besmarts.mechanics import fits
from besmarts.mechanics import smirnoff_xml
from besmarts.mechanics import smirnoff_models
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import optimizers_scipy
from besmarts.mechanics import fits
from besmarts.core import graphs
from besmarts.core import topology
from besmarts.core import perception
from besmarts.core import arrays
from besmarts.core import assignments
from besmarts.core import optimization
from besmarts.assign import hierarchy_assign_rdkit
from besmarts.codecs import codec_rdkit
from besmarts.core import compute
from besmarts.core import configs
# configs.processors = 1
xyz_positions = """11
Energy =  -802.9904298568 (+0.195 kcal/mol)
  C -1.44819400 -0.84940800  0.16848900
  C -1.59401300  0.50318700 -0.01678100
  C -0.27397600  1.02622600 -0.13503500
  C  0.58064400 -0.04716400 -0.01303100
  C  2.03461200 -0.06860900 -0.05925200
  O  2.72809700  0.90108700 -0.21909900
 Cl  2.76214600 -1.70734100  0.14655600
  O -0.13897300 -1.20044600  0.17351800
  H -2.15226800 -1.65836100  0.30609000
  H -2.52743000  1.04809900 -0.06180000
  H  0.02935200  2.05273200 -0.28965800
"""
xyz_grad = """11

  C      0.49755    0.17370   -0.04115
  C     -0.00884    0.07632   -0.01031
  C      0.20074   -0.69547    0.09073
  C     -0.02955    1.24848   -0.17483
  C      0.55229   -1.91119    0.25039
  O     -0.15948    0.65794   -0.08724
 Cl     -0.33030    0.82983   -0.10559
  O     -0.73720   -0.66864    0.11909
  H     -0.11502    0.11021   -0.01168
  H     -0.00691    0.04737   -0.00649
  H      0.02566   -0.05163    0.00657
"""
def load_xyz(flist, indices=None) -> assignments.graph_db_row:
    s = 0
    N = None
    gdr = assignments.graph_db_row()
    for f in flist:
        lines = f.split('\n')
        if N is None:
            N = int(lines[0].split()[0])
        if indices is None:
            indices = list(range(N))
        assert N == int(lines[0].split()[0])
        for chunk in arrays.batched(lines, N+2):
            if chunk and chunk[0]:
                sel = {}
                for i, line in enumerate(chunk, -2):
                    if i >= 0:
                        sym, x, y, z = line.split()[:4]
                        sel[indices[i]] = list(map(float, (x, y, z)))

                gdc = assignments.graph_db_column()
                gdc.selections.update(sel)
                gdr.columns[s] = gdc
                s += 1
    return gdr
def make():
    smi = "[C:1]1([H:9])=[C:2]([H:10])[C:3]([H:11])=[C:4]([C:5](=[O:6])[Cl:7])[O:8]1"
    s = xyz_positions 
    g = xyz_grad
    d  = {
        smi: [
            {
                assignments.POSITIONS: s,
                assignments.GRADIENTS: g,
            },
        ],
    }
    return d
def new_gdb(f: Dict[str, List[str]]) -> assignments.graph_db:
    gcd = codec_rdkit.graph_codec_rdkit()
    gdb = assignments.graph_db()

    ne = 0
    for smi, fn_dict in f.items():

        g = gcd.smiles_decode(smi)
        gid = assignments.graph_db_add_graph(gdb, g)

        gdb.graphs[gid] = g
        gdb.smiles[gid] = smi
        gdb.selections[topology.index_of(topology.atom)] = {
            gid: {k: v for k, v in enumerate(graphs.graph_atoms(g))}
        }
        gde = assignments.graph_db_entry()
        gdb.entries[len(gdb.entries)] = gde
        for rid, rdata in enumerate(fn_dict):
            tid = assignments.POSITIONS
            gdt = assignments.graph_db_table(topology.atom)
            gdg = assignments.graph_db_graph()
            gdt.graphs[gid] = gdg
            fn = rdata[tid]
            indices=dict(sorted([(j, x) for j, x in enumerate(g.nodes)], key=lambda x: x[1]))
            gdr = load_xyz([fn], indices=list(indices))
            gdg.rows[0] = gdr
            gde.tables[tid] = gdt
            tid = assignments.GRADIENTS
            if tid in rdata:
                gdt = assignments.graph_db_table(topology.atom)
                gdg = assignments.graph_db_graph()
                gdt.graphs[gid] = gdg
                fn = rdata[tid]
                indices=dict(sorted([(j, x) for j, x in enumerate(g.nodes)], key=lambda x: x[1]))
                gdr = load_xyz([fn], indices=list(indices))
                gdg.rows[0] = gdr
                gde.tables[tid] = gdt
                gx = [x for y in gdr[0].selections.values() for x in y]
                gdt.values.extend(gx)
            tid = assignments.ENERGY
            if tid in rdata:
                gdt = assignments.graph_db_table(topology.null)
                fn = rdata[tid]
                ene = [*map(float, 
                    [x for x in open(fn).read().split('\n') if x]
                )]
                gdt.values.extend(ene)
                gde.tables[tid] = gdt
    return gdb
def run(d, ff_fn):
    # build the dataset and input ff
    gcd = codec_rdkit.graph_codec_rdkit()
    labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
    pcp = perception.perception_model(gcd, labeler)
    csys = smirnoff_models.smirnoff_load(ff_fn, pcp)
    gdb = new_gdb(d)
    psys = fits.gdb_to_physical_systems(gdb, csys)
    models = [0,1]
    fit_models = [0,1]
    strat = fits.forcefield_optimization_strategy_default(csys, models=models)
    co = fits.chemical_objective
    final = None # default would be a full fit
    initial = final # have a default, full fit
    onestep = None # do a forcefield fit with one step 
    final = fits.objective_tier()
    final.objectives = {
        0: fits.objective_config_position(
                fits.graph_db_address(
                    eid=[0],
                ),
                scale=10
        ),
        1: fits.objective_config_gradient(
                fits.graph_db_address(
                    eid=[0],
                ),
                scale=1e-5
        ),
    }
    final.key_filter = lambda x: x[0] in fit_models and x[1] == 'l'
    onestep = fits.objective_tier()
    onestep.objectives = final.objectives
    onestep.step_limit = 2
    onestep.key_filter = lambda x: x[0] in fit_models and x[1] == 'l'
    initial = final
    tiers = [onestep] # have a default
    kv0 = mm.chemical_system_iter_keys(csys)
    newcsys, P, C = fits.ff_optimize(
        csys,
        gdb,
        psys,
        strat,
        co,
        initial,
        tiers,
        final
    )
    fits.print_chemical_system(newcsys)
    print(f"X {P+C:15.8g} P {P:15.8g} C {C:15.8g}")
if __name__ == "__main__":
    run(make(), "openff-2.1.0.offxml")
