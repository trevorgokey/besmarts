
"""
Tests for BESMARTS file format round-tripping
"""

import io
import sys
import unittest

from besmarts.codecs import codec_native
from pprint import pprint

def build_input():
    inp = """#GRAPH
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic chirality formal_charge
#BOND bond_ring bond_order
  1   1  64   8  16   1   1   1   1   1
  2   2  64   1  16   4   8   1   1   1
  3   3  64   8  16   1   1   1   1   1
  4   4 256   1   4   4   8   1   1   1
  5   5  64   4  16   4   8   1   1   1
  6   6  64   2  16   4   8   1   2   1
  7   7 256   2   4   1   1   1   1   1
  8   8  64   1  16  16   8   1   2   1
  9   9 256   1   4   4  16   1   1   1
 10  10  64   2  16   4  16   1   2   1
 11  11  64   8  16   1   1   1   1   1
 12  12  64   2  16   8   8   1   4   1
 13  13  64   4  16   4   8   1   1   1
 14  14  64   2   8   4   8   1   1   1
 15  15  64   1   8   8   8   1   1   1
 16  16  64   1  16  16   8   1   4   1
 17  17  64   4  16   4  16   1   1   1
 18  18  64   1   8   4  16   1   1   1
 19  19 256   1   2   1   1   1   1   1
 20  20  64   4  16   4  16   1   1   1
 21  21  64   2  16   8  16   1   2   1
 22  22  64   2  16   8  16   1   2   1
 23  23  64   4  16   4  16   1   1   1
 24  24  64   4  16   4  16   1   1   1
 25  25  64   2  16   8  16   1   2   1
 26  26  64   4  16   4  16   1   1   1
 27  27  64   1  16   8  16   1   4   1
 28  28  64   8  16   1   1   1   1   1
 29  29  64   4  16   4  16   1   1   1
 30  30  64   1   8   8  16   2   1   1
 31  31 128   1   4   4  16   2   1   1
 32  32  64   1   8   8  16   2   1   1
 33  33 128   1   4   4  16   2   1   1
 34  34  64   1   8   8  16   2   1   1
 35  35  64   4  16   4  16   1   1   1
 36  36  64   1  16   8  16   1   2   1
 37  37  64   8  16   1   1   1   1   1
 38  38  64   1   8   8  16   2   1   1
 39  39  64   4  16   4  16   1   1   1
 40  40  64   2  16   8  16   1   4   1
 41  41  64   4  16   4  16   1   1   1
 42  42  64   4  16   4  16   1   1   1
 43  43  64   2  16   8  16   1   4   1
 44  44  64   2  16   8  16   1   4   1
 45  45  64   4  16   4  16   1   1   1
 46  46  64   2  16   4  16   1   2   1
 47  47 256   2   4   1   1   1   1   1
 48  48  64   1  16   8   8   1   2   1
 49  49  64   8  16   1   1   1   1   1
 50  50  64   1   8   8   8   1   1   1
 51  51  64   2   8   4   8   1   1   1
 52  52  64   2  16   8   8   1   4   1
 53  53 256   1   4   4   8   1   1   1
 54  54  64   1  16   8   8   1   4   1
 55  55 256   2   4   1   1   1   1   1
 56  56  64   2  16   4   8   1   4   1
 57  57  64   8  16   1   1   1   1   1
 58  58  64   1  16  16   8   1   4   1
 59  59 256   1   4   4   8   1   1   1
 60  60  64   2  16   4   8   1   4   1
 61  61 256   2   4   1   1   1   1   1
 62  62  64   4  16   4   8   1   1   1
 63  63  64   1  16   4   8   1   2   1
 64  64  64   8  16   1   1   1   1   1
 65  65  64   4  16   1   1   1   1   1
 66  66 256   2   4   1   1   1   1   1
  1   2   1   2
  2   3   1   2
  2   4   2   2
  2   5   2   2
  5   6   2   2
  6   7   1   2
  6   8   2   2
  8   9   2   2
  8  10   2   2
 10  11   1   2
 10  12   2   2
 12  13   2   2
 13  14   2   2
 14  15   2   4
 15  16   2   2
 16  17   2   2
 16  18   2   2
 18  19   1   4
 18  20   2   2
 20  21   2   2
 21  22   2   2
 22  23   2   2
 23  24   2   2
 24  25   2   2
 25  26   2   2
 25  27   2   2
 27  28   1   2
 27  29   2   2
 29  30   2   2
 30  31   2  32
 30  32   2  32
 32  33   2  32
 33  34   2  32
 34  35   2   2
 35  36   2   2
 36  37   1   2
 34  38   2  32
 38  39   2   2
 39  40   2   2
 40  41   2   2
 41  42   2   2
 42  43   2   2
 43  44   2   2
 44  45   2   2
 45  46   2   2
 46  47   1   2
 46  48   2   2
 48  49   1   2
 48  50   2   2
 50  51   2   4
 51  52   2   2
 52  53   2   2
 52  54   2   2
 54  55   1   2
 54  56   2   2
 56  57   1   2
 56  58   2   2
 58  59   2   2
 58  60   2   2
 60  61   1   2
 60  62   2   2
 62  63   2   2
 63  64   1   2
 63  65   1   2
 65  66   1   2
  4   8   2   2
  9  17   2   2
 12  16   2   2
 15  22   2   2
 21  27   2   2
 26  32   2   2
 31  38   2  32
 36  40   2   2
 36  44   2   2
 43  50   2   2
 48  54   2   2
 53  58   2   2
 59  63   2   2
"""
    return io.StringIO(inp) 

class test_io(unittest.TestCase):
    def test_io(self):
        buf_in = build_input()
        A = codec_native.graph_codec_native_read(buf_in)
        self.assertTrue(len(A) == 1)

        buf_out = io.StringIO()
        codec_native.graph_codec_native_write(buf_out, A)
        B = codec_native.graph_codec_native_read(buf_out)
        self.assertTrue(len(B) == 1)

        # graphs might be iterated differently so just check string length and
        # if indices are present
        self.assertTrue(len(buf_in.getvalue()) == len(buf_out.getvalue()))
        self.assertFalse(set(A[0].nodes).symmetric_difference(B[0].nodes))
        self.assertFalse(set(A[0].edges).symmetric_difference(B[0].edges))

