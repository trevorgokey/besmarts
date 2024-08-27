"""
tests/test_intvec.py
"""

import io
import unittest
from besmarts.core.codecs import intvec_codec
from besmarts.codecs import codec_native
from timeit import timeit

def build_input():
    inp = """#GRAPH
#ATOM element 
#BOND bond_order
-21 -21     2   
-10 -10    64   
 -9  -9    64   
 -8  -8    64   
 -7  -7    64   
-12 -12    64   
-11 -11    64   
-22 -22     2   
-23 -23     2   
 -6  -6   128   
-18 -18     2   
 -4  -4 65536 
 -5  -5   256   
 -3  -3   256   
 -2  -2    64   
-16 -16     2   
-17 -17     2   
 -1  -1    64   
-13 -13     2   
-14 -14     2   
-15 -15     2   
-19 -19     2   
-20 -20     2   
 10  21     2
  9  10    32
 10  11    32
  8   9    32
  9  20     2
  7   8    32
  8  19     2
  7  12    32
  6   7     2
 11  12    32
 12  23     2
 11  22     2
  6  18     2
  4   6     2
  4   5     4
  3   4     2
  2   3     2
  2  16     2
  2  17     2
  1   2     2
  1  13     2
  1  14     2
  1  15     2
"""
    return io.StringIO(inp)

class test_intvec(unittest.TestCase):
    def test_intvec(self):

        n = 1 

        codecs = codec_native.primitive_codecs_get()
        atom_prims = ["element"]
        bond_prims = ["bond_order"]
        gcd = codec_native.graph_codec_native(codecs, atom_prims, bond_prims)

        inp = build_input()
        g = codec_native.graph_codec_native_read(inp)[0]
        # t = timeit('codec_native.graph_codec_native_read(build_input())[0]', globals=globals(), number=n)/n

        atom_prims = []
        for k, v in g.nodes.items():
            atom_prims.extend(list(v.primitives))
            break

        bond_prims = []
        for k, v in g.edges.items():
            bond_prims.extend(list(v.primitives))
            break


        # gcd.atom_primitives = atom_prims
        # gcd.bond_primitives = bond_prims


        icd = intvec_codec(gcd.primitive_codecs, atom_prims, bond_prims)
        # print(f"SMILES loaded, time={t} sec/call")


        smi_a = gcd.smiles_encode(g)
        # t = timeit('gcd.smiles_encode(g)', globals=globals(), number=n)/n
        # print(smi_a)
        # print(f"SMILES encoded, time={t} sec/call")


        intvec = icd.graph_encode(g)
        # t = timeit('icd.graph_encode(g)', globals=globals(), number=n)/n
        # print(f"intvec Encoded, time={t} sec/call")

        g = icd.graph_decode(intvec)
        # t = timeit('icd.graph_decode(intvec)', globals=globals(), number=n)/n
        # print(g)
        # print(f"intvec Decoded, time={t} sec/call")

        smi_b = gcd.smiles_encode(g)
        # t = timeit('gcd.smiles_encode(g)', globals=globals(), number=n)/n
        # print(smi_b)
        # print(f"SMILES encoded, time={t} sec/call")

        self.assertTrue(smi_a == smi_b)
        # print("SMILES are the same; test passed!")
