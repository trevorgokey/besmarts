from besmarts.codecs import codec_native
from pprint import pprint

graphs = []
for g in codec_native.graph_codec_native_load("./out.bes"):
    pprint(g)
    pprint(g.nodes)
    pprint(g.edges)
    # pprint(g.select)
    graphs.append(g)

codec_native.graph_codec_native_save("./out.bes", graphs)
