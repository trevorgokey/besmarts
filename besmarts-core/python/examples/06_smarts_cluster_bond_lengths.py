from besmarts.cluster.cluster_assignment import smiles_assignment_float
from besmarts.core.assignments import smiles_assignment_group_bonds
from besmarts.cluster.cluster_objective import clustering_objective_mean_separation
from besmarts.cluster.cluster_optimization import cluster_means
from besmarts.codecs.codec_rdkit import graph_codec_rdkit
from besmarts.assign import hierarchy_assign_rdkit
from besmarts.core import configs

configs.processors = 1
configs.remote_compute_enable = False
configs.workqueue_port = 59321

gcd = graph_codec_rdkit()
labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
objective = clustering_objective_mean_separation(split_separation=0.1)

smi = "[C:1]([H:3])#[C:2][H:4]"
assns = {
    (1, 2): [1.1],
    (1, 3): [1.3],
    (2, 4): [1.3]
}

sa = smiles_assignment_float(smi, assns)

sag = smiles_assignment_group_bonds([sa])
cst = cluster_means(gcd, labeler, sag, objective=objective)
