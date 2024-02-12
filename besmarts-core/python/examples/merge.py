
# from besmarts.hierarchy.hierarchy_index import hindex, hentry, hindex_node_depth
# from besmarts.hierarchy.hindex_iters import hindex_iter_breadth_first, hindex_add_hindex
# from besmarts.hierarchy.assign_rdkit import assign_smiles_rdkit



# from besmarts.hierarchy.merge import hindex_merge

# from besmarts.hierarchy.smirnoff import hindex_smirnoff_load_xml, hindex_smirnoff_save_xml
# from besmarts.hierarchy import smirnoff
# from besmarts.codecs.codec_rdkit import graph_codec_rdkit
# from besmarts.hierarchy.hindex_iters import hindex_iter_dive

# gcd = graph_codec_rdkit()

# hA = hindex_smirnoff_load_xml("in.offxml", gcd)
# # for e in hindex_iter_dive(hA, hA.root):
# #     s = " " * hA.node_depth(e)
# #     print(s, e.index, e.key)
# hB = hindex_smirnoff_load_xml("A.offxml", gcd)
# # for e in hindex_iter_dive(hB, hB.root):
# #     s = " " * hB.node_depth(e)
# #     print(s, e.index, e.key)

# mols = ["CCC", "CC", "CCCC", "CC=O", "CCO", "CCN", "C"]
# mols = None

# ff = hindex_merge(hA, hB, gcd, mols)

# for e in hindex_iter_dive(ff, ff.root):
#     s = " " * hindex_node_depth(ff, e)
#     if e.key in ff.db:
#         print(s, e.index, e.key)
#         if "A" not in ff.db[e.key]:
#             continue
#         print(s, "   ", ff.db[e.key].get("smarts"))
#         param = ff.db[e.key].get("parameter")
#         if param:
#             print(s, "     :", param.smarts)
#             for term, value in param.terms.items():
#                 print(s, "     :", term, value)
#         param = ff.db[e.key]["A"].get("parameter")
#         if param:
#             print(s, "    A:", param.smarts)
#             for term, value in param.terms.items():
#                 print(s, "    A:", term, value)
#         param = ff.db[e.key]["B"].get("parameter")
#         if param:
#             print(s, "    B:", param.smarts)
#             for term, value in param.terms.items():
#                 print(s, "    B:", term, value)
#     else:
#         print(s, e.index, e.key)

# smirnoff.hindex_smirnoff_save_xml(ff, "merged.offxml")

# ff = smirnoff.hindex_smirnoff_rekey(ff, None)

# smirnoff.hindex_smirnoff_save_xml(ff, "merged.rekey.offxml")

