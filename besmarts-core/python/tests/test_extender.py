import io
import unittest

from besmarts.codecs import codec_rdkit, codec_native
from besmarts.core import graphs
from besmarts.core import configs
from besmarts.core import mapper
from besmarts.core import topology


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


class test_extender(unittest.TestCase):
    def test_extender(self):
        # check depth of 3 on 21-10-11-22 should have node 6. Interesting test
        # because we need to ensure that it doesn't go around the ring the long
        # way

        buf_in = build_input()
        g = codec_native.graph_codec_native_read(buf_in)[0]

        torsions = [
            graphs.graph_to_structure(g, (21, 10, 11, 22), topology.torsion)
        ]

        cfg = configs.smarts_extender_config(3, 3, True)
        mapper.mapper_smarts_extend(cfg, torsions)
        self.assertTrue(6 in torsions[0].select)
