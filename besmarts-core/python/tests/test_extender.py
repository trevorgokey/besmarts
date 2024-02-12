

import io
from besmarts.codecs import codec_rdkit, codec_native
from besmarts.core import graphs
from besmarts.core import configs
from besmarts.core import mapper

def build_input():
    inp = """#GRAPH
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic formal_charge
#BOND bond_ring bond_order
-21 -21   2   1   2   1   1   1   1
-10 -10  64   2   8   4  16   2   1
 -9  -9  64   2   8   4  16   2   1
 -8  -8  64   2   8   4  16   2   1
 -7  -7  64   1   8   4  16   2   1
-12 -12  64   2   8   4  16   2   1
-11 -11  64   2   8   4  16   2   1
-22 -22   2   1   2   1   1   1   1
-23 -23   2   1   2   1   1   1   1
 -6  -6 128   2   8   1   1   1   1
-18 -18   2   1   2   1   1   1   1
 -4  -4 65536   1   8   1   1   1   1
 -5  -5 256   1   2   1   1   1   1
 -3  -3 256   1   4   1   1   1   1
 -2  -2  64   4  16   1   1   1   1
-16 -16   2   1   2   1   1   1   1
-17 -17   2   1   2   1   1   1   1
 -1  -1  64   8  16   1   1   1   1
-13 -13   2   1   2   1   1   1   1
-14 -14   2   1   2   1   1   1   1
-15 -15   2   1   2   1   1   1   1
-19 -19   2   1   2   1   1   1   1
-20 -20   2   1   2   1   1   1   1
 10  21   1   2
  9  10   2  32
 10  11   2  32
  8   9   2  32
  9  20   1   2
  7   8   2  32
  8  19   1   2
  7  12   2  32
  6   7   1   2
 11  12   2  32
 12  23   1   2
 11  22   1   2
  6  18   1   2
  4   6   1   2
  4   5   1   4
  3   4   1   2
  2   3   1   2
  2  16   1   2
  2  17   1   2
  1   2   1   2
  1  13   1   2
  1  14   1   2
  1  15   1   2
"""
    return io.StringIO(inp) 

# check depth of 3 on 21-10-11-22 should get me 6. Interesting test because
# we need to ensure that it doesn't go around the ring the long way

buf_in = build_input()
g = codec_native.graph_codec_native_read(buf_in)[0]

# sigh... doesn't work like this, but works when i load from tbg
#smiles = '[H:21][c:10]1[c:9]([c:8]([c:7]([c:12]([c:11]1[H:22])[H:23])[N:6]([H:18])[S:4](=[O:5])[O:3][C:2]([H:16])([H:17])[C:1]([H:13])([H:14])[H:15])[H:19])[H:20]'
#gcd = codec_rdkit.graph_codec_rdkit()
#g = gcd.smiles_decode(smiles)
# # codec_native.graph_codec_native_write("g.bg", [g] )

torsions = graphs.graph_to_structure_dihedrals(g)

torsions = [x for x in torsions if x.select == (21, 10, 11, 22)]
assert len(torsions) == 1

cfg = configs.smarts_extender_config(3, 3, True)
mapper.mapper_smarts_extend(cfg, torsions)
assert 6 in torsions[0].select







