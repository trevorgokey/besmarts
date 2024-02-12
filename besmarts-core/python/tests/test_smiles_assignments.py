"""
besmarts.tests.test_smiles_assignments


"""
from besmarts.cluster.cluster_assignment import smiles_assignment_float
from besmarts.core.assignments import smiles_assignment_group_bonds
from besmarts.cluster.cluster_objective import clustering_objective_mean_separation
from besmarts.cluster.cluster_optimization import cluster_means
from besmarts.codecs.codec_rdkit import graph_codec_rdkit

gcd = graph_codec_rdkit()

smi = "[C:1]([H:3])#[C:2][H:4]"
assns = {(1,2): [1.1], (1,3): [1.3], (2,4): [1.3]}
sa  = smiles_assignment_float(smi, assns)

objective = clustering_objective_mean_separation(split_separation=0.1)


sag = smiles_assignment_group_bonds([sa])
# cst = cluster_means(gcd, sag, objective=objective)
