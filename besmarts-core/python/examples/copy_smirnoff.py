
# from besmarts.hierarchy.hierarchy_index import hindex, hentry, hindex_node_depth
# from besmarts.hierarchy.assign_rdkit import assign_smiles_rdkit

# from besmarts.hierarchy.merge import hindex_merge

from besmarts.hierarchy.smirnoff import hindex_smirnoff_load_xml, hindex_smirnoff_save_xml
from besmarts.codecs.codec_rdkit import graph_codec_rdkit
# from besmarts.hierarchy.hindex_iters import hindex_iter_dive

gcd = graph_codec_rdkit()

hA = hindex_smirnoff_load_xml("in.offxml", gcd)
hindex_smirnoff_save_xml(hA, "in.copy.offxml")
# hA = hindex_smirnoff_load_xml("in.copy.offxml", gcd)
# hindex_smirnoff_save_xml(hA, "in.copy.2.offxml")


