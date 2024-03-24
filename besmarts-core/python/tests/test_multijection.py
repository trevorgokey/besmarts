"""
tests/test_multijection
"""

from pprint import pprint

from besmarts import graphs, topologies, bijections, multijections

from besmarts.codecs import codec_native as native
from besmarts.codecs.codec_rdkit import graph_codec_rdkit

gcd = graph_codec_rdkit()

topo = topology.atom

ref = graphs.subgraph_to_structure(gcd.smarts_decode("[*:1]~[*](~[*])"), topo)
smarts_graphs = list(map(lambda x: graphs.subgraph_to_structure(x, topo), map(gcd.smarts_decode, [
    "[#6:1]-[#6]-[#1]",
    "[#6:1]-[#6]-[#6]",
    "[#6:1]-[#6]-[#8]",
    "[#6:1]-[#7]-[#1]",
    "[#6:1]-[#6](-[#6])-[#6]",
    "[#6:1]-[#6](-[#6])-[#1]"
])))

FF = multijections.multijection_fit_reference(smarts_graphs, ref, fill=True)
print(gcd.smarts_encode(multijections.multijection_union(FF)))
