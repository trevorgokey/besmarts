
"""
Example of combining the bond graphs of propane into a single pattern.

pypy can be used to approximately half the runtime needed for this example.
However, pypy can only be used when rdkit is not used. To do this, make sure all
graphs are decoded and saved first, then load them using the native load function
in this example.
"""

from besmarts.codecs import codec_native
from besmarts.core import graphs
from besmarts.core import mapper
from besmarts.core import configs

# load the default native graph codec. It can encode SMARTS/SMILES
codecs = codec_native.primitive_codecs_get()
atom_primitives = list(codec_native.primitive_codecs_get_atom())
bond_primitives = list(codec_native.primitive_codecs_get_bond())
gcd = codec_native.graph_codec_native(codecs, atom_primitives, bond_primitives)


# load in a pre-decoded propane graph
G = codec_native.graph_codec_native_load("propane.bg")[0]
print("SMARTS of propane:")
print(gcd.smarts_encode(G))

gcd.atom_primitives = ("element", "hydrogen")
gcd.bond_primitives = ("bond_order",)
print("SMARTS of propane (element, hydrogen, and bond order):")
print(gcd.smarts_encode(G))

GnoH = graphs.graph_remove_hydrogen(G)
print("SMARTS of propane no hydrogen:")
print(gcd.smarts_encode(GnoH))

bonds = graphs.graph_to_structure_bonds(G)

print("SMARTS of propane bonds")
for i, bond in enumerate(bonds, 1):
    print(f"{i:2d}", gcd.smarts_encode(bond))

print("SMARTS of unique propane bonds")
for i, bond in enumerate(set(bonds), 1):
    print(f"{i:2d}", gcd.smarts_encode(bond))

# default depth is 0, so it will only union the primary atoms (i.e. the bond)
U = mapper.union_list(bonds)

print("SMARTS union (depth=0):")
print(gcd.smarts_encode(U))

print("Breakdown:")
graphs.structure_print(U)


# extend the bonds to include the immediate neighbors
cfg = configs.smarts_extender_config(1, 1, True)
graphs.structure_extend(cfg, bonds)
U = mapper.union_list(bonds)

print("SMARTS union (depth=1):")
print(gcd.smarts_encode(U))

print("Breakdown:")
graphs.structure_print(U)

