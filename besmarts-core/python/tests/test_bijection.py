"""
besmarts.tests.test_bijection

Test the mapper modes that add/remove nodes during map
"""
import unittest
from pprint import pprint

from besmarts.codecs import codec_native
from besmarts.core import topology
from besmarts.core import bijections
from besmarts.core import graphs

from besmarts.codecs.codec_rdkit import graph_codec_rdkit

gcd = graph_codec_rdkit()
atom = topology.atom_topology()
bond = topology.bond_topology()

bijection_fn = {
    0: bijections.structure_bijection_mcs,
    1: bijections.structure_bijection_constant_right,
    2: bijections.structure_bijection_constant_left,
    3: bijections.structure_bijection_branch
}

class test_bijections(unittest.TestCase):

    def test_mapping(self):

        checks = [
            ["[*:1]",     "[*:1]-[*]", 0, "[*:1]",     "[*:1]"],
            ["[*:1]",     "[*:1]-[*]", 1, "[*:1]~[*]", "[*:1]-[*]"],
            ["[*:1]",     "[*:1]-[*]", 2, "[*:1]~[*]", "[*:1]-[*]"],
            ["[*:1]",     "[*:1]-[*]", 3, "[*:1]",     "[*:1]"    ],
            ["[*:1]-[*]", "[*:1]",     0, "[*:1]",     "[*:1]"    ],
            ["[*:1]-[*]", "[*:1]",     1, "[*:1]-[*]", "[*:1]~[*]"],
            ["[*:1]-[*]", "[*:1]",     2, "[*:1]",     "[*:1]"    ],
            ["[*:1]-[*]", "[*:1]",     3, "[*:1]-[*]", "[*:1]~[*]"],
        ]
        self.assertTrue(check(checks, atom))

        x = "[*:1]~[*]"
        y = "[*:1](~[*])~[*]~[*]"
        checks = [
            [x,     y, 0, x,  x],
            [x,     y, 1, y,  y],
            [x,     y, 2, y,  y],
            [x,     y, 3, x,  x],
        ]
        self.assertTrue(check(checks, atom))

        a = "[*:1]~[#6]"
        b = "[*:1](~[#6])~[*]~[*]"
        checks = [
            [a,     b, 0, a,  a],
            [a,     b, 1, b,  b],
            [a,     b, 2, b,  b],
            [a,     b, 3, a,  a],
        ]
        self.assertTrue(check(checks, atom))

        a = "[*:1]~[#8]"
        b = "[*:1](~[#6])~[*]~[*]"
        checks = [
            [a,     b, 0, a,  "[*:1]~[*]"],
            [a,     b, 1, "[*:1](~[#8])~[*]~[*]",  b],
            [a,     b, 2, "[*:1](~[#8])~[*]~[*]",  b],
            [a,     b, 3, a,  "[*:1]~[*]"],
        ]
        self.assertTrue(check(checks, atom))

        a = "[*:1]~[#8:2]"
        b = "[*:1](~[#6])~[*:2]~[*]"
        checks = [
            [a,     b, 0, a,  "[*:1]~[*:2]"],
            [a,     b, 1, "[*:1](~[*])~[#8:2]~[*]",  b],
            [a,     b, 2, "[*:1](~[*])~[#8:2]~[*]",  b],
            [a,     b, 3, a,  "[*:1]~[*:2]"],
        ]
        self.assertTrue(check(checks, bond))

def check(checks, topo):

    for gs, hs, add_nodes, tg, th in checks:
        g = graphs.subgraph_to_structure(gcd.smarts_decode(gs), topo)
        h = graphs.subgraph_to_structure(gcd.smarts_decode(hs), topo)

        fn = bijection_fn[add_nodes]
        F: bijections.bijection = fn(g, h)
        print("Mapping", gs, "to", hs, "with", ".".join((fn.__module__,fn.__name__)))
        print("A: ", gcd.smarts_encode(g), "is original")
        print("TA:", gcd.smarts_encode(F.g), "should be", tg)
        print("B: ", gcd.smarts_encode(h), "is original")
        print("TB:", gcd.smarts_encode(F.h), "should be", th)
        print("M: ", F.t)
        if gcd.smarts_encode(F.g) != tg:
            T = mapper.map_to(g, h, add_nodes=add_nodes, fill=True)
            return False
        elif gcd.smarts_encode(F.h) != th:
            T = mapper.map_to(g, h, add_nodes=add_nodes, fill=True)
            return False
        assert gcd.smarts_encode(F.g) == tg
        assert gcd.smarts_encode(F.h) == th
    return True

if __name__ == "__main__":
    unittest.main()

