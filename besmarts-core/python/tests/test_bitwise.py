"""
besmarts.tests.test_bitwise
"""

import unittest
from pprint import pprint

from besmarts.core import topology
from besmarts.core import bijections
from besmarts.core import graphs

from besmarts.codecs.codec_rdkit import graph_codec_rdkit

class test_bitwise(unittest.TestCase):

    def setUp(self):
        self.gcd = graph_codec_rdkit()
        self.gcd.atom_primitives = ("element",)
        self.gcd.bond_primitives = ("bond_order",)

    def test_union(self):

        gcd = self.gcd
        topo = topology.angle
        G = graphs.subgraph_to_structure(gcd.smarts_decode("[#6:1]-[#7:2]-[#8:3]"), topo)
        H = graphs.subgraph_to_structure(gcd.smarts_decode("[#1:1]-[#6:2]-[#8:3]"), topo)
        F = bijections.structure_bijection_mcs(G, H)
        Q = bijections.bijection_union(F)
        self.assertEqual(gcd.smarts_encode(Q), "[#1,#6:1]-[#6,#7:2]-[#8:3]")

    def test_intersection(self):

        gcd = self.gcd
        topo = topology.angle
        G = graphs.subgraph_to_structure(gcd.smarts_decode("[#1,#6:1]-[#7:2]-[#8:3]"), topo)
        H = graphs.subgraph_to_structure(gcd.smarts_decode("[#1:1]-[#7:2]-[#8:3]"), topo)
        F = bijections.structure_bijection_mcs(G, H)
        Q = bijections.bijection_intersection(F)
        self.assertEqual(gcd.smarts_encode(Q), "[#1:1]-[#7:2]-[#8:3]")

    def test_xor(self):

        gcd = self.gcd
        topo = topology.angle
        G = graphs.subgraph_to_structure(gcd.smarts_decode("[#1:1]-,=[#7,#8:2]-,=[#6:3]"), topo)
        H = graphs.subgraph_to_structure(gcd.smarts_decode("[#1:1]-[#7:2]-[#8:3]"), topo)
        F = bijections.structure_bijection_mcs(G, H)
        Q = bijections.bijection_xor(F)
        self.assertEqual(gcd.smarts_encode(Q), "[_:1]=[#8:2]=[#6,#8:3]")

    def test_subtract(self):

        gcd = self.gcd
        topo = topology.angle
        G = graphs.subgraph_to_structure(gcd.smarts_decode("[#1,#6:1]-,=[#8:2]-,=[#6,#8:3]"), topo)
        H = graphs.subgraph_to_structure(gcd.smarts_decode("[#6:1]-[#6:2]-[#6:3]"), topo)
        F = bijections.structure_bijection_mcs(G, H)
        Q = bijections.bijection_subtract(F)
        self.assertEqual(gcd.smarts_encode(Q), "[#1:1]=[#8:2]=[#8:3]")

if __name__ == "__main__":
    unittest.main()

