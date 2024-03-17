
"""
besmarts.tests.test_match

Sanity checks for determining matching (i.e. subset relation)

"""

import unittest
from pprint import pprint

from besmarts.codecs import codec_native
from besmarts.core import topology
from besmarts.core import mapper
from besmarts.core import graphs

from besmarts.codecs.codec_rdkit import graph_codec_rdkit

gcd = graph_codec_rdkit()
topo = topology.atom_topology()

class test_match(unittest.TestCase):

    def test_match(self):
        checks = [
            ["[*:1]", "[*:1]-[*]", False],
            ["[*:1]=[*]", "[*:1]-[*]", False],
            ["[*:1]-[#1]", "[*:1]-[*]", True],
            ["[*:1]", "[*:1]-[#7]", False],
            ["[*:1]", "[X1:1][*]", False],
            ["[!r:1](-[#1])-[#6H3]", "[*:1]-[*]", True],
            ["[!r:1](-[#1])-[#6H3]", "[#1]-[!r:1]-[#6]", True],
            ["[!r:1](-[#1])-[#6H3]", "[#1]-[!r:1]-[#6H3]", True],
            ["[#6H3:1](-[#6H2])-[#6H2]", "[*]-[*:1]-[!#8]", True],
            ["[#6H3:1](-[#6H2])-[#8H1]", "[*]-[*:1]-[!#8]", True],
            ["[#6H3:1](-[#6H2])-[#8H1]", "[!#8]-[*:1]-[!#8]", False],
            ["[#6H3:1](-[#6H2])-[#8H1]-[#1]", "[*]-[*:1]-[!#8]-[#1]", False],
            ["[#6H3:1](-[#6H2])-[#8H1]-[#1]", "[*]-[*:1]-[*]-[#1]", True],
            ["[#1;H0;X1;x0;!r;A;+0:37]!@;-[#6;H2;X4;x0;!r;A;+0](!@;-[#6;H0;X3;x0;!r;A;+0])!@;-[#6;H1;X4;x0;!r;A;+0]", "[*:1]~[*]~[!#6]", False]
        ]
        self.assertTrue(check(checks, topo))

        bonds = [
            ["[#6H3:1](-[#1])(-[#1])(-[#1])-[#6H0:2](-[#6H3])(-[#8H0])-[#6H2]", "[*:1]~[*:2]~[H3]", True]
        ]


        self.assertTrue(check(bonds, topology.bond_topology()))

def check(checks, topo):
    for gs, hs, check in checks:
        g = graphs.subgraph_to_structure(gcd.smarts_decode(gs), topo)
        h = graphs.subgraph_to_structure(gcd.smarts_decode(hs), topo)

        print("Does")
        print(gcd.smarts_encode(g))
        print("belong to cluster")
        print(gcd.smarts_encode(h))

        match = mapper.mapper_match(g, h)
        print("Response:", match, "Answer:", check)
        print()
        if match != check:
            breakpoint()
            match = mapper.mapper_match(g, h)
        assert match == check
    return True
