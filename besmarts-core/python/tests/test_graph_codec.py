"""
besmarts.tests.test_graph_codec

Run the graph codec and check that can be instantiated correctly
"""

from besmarts.codecs.codec_native import graph_codec_native, primitive_codecs_get

codecs = primitive_codecs_get()
gcd = graph_codec_native(codecs, ['element', 'formal_charge', 'connectivity_ring'], ['bond_order'])

sma = "[!#8!#15!#16!#17!#35!+!X4:1]~[*:2]"
print(sma)
# TODO write the parser for the native codec

# g = gcd.smarts_decode(sma)
#assert gcd.smarts_encode(g) == sma
