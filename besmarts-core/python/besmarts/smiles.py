#!/usr/bin/env python3

from besmarts.codecs import codec_rdkit
from besmarts.core import graphs
import sys

smi = sys.argv[1]
gcd = codec_rdkit.graph_codec_rdkit()
g = gcd.smiles_decode(smi)
print(gcd.smiles_encode(graphs.graph_to_subgraph(g, tuple(g.nodes))))

