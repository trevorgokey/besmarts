Python Examples
===============

SMILES decoding to structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

>>> from besmarts.codecs import codec_native
>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.core import graphs 
>>> gcd = graph_codec_rdkit()
>>> smi_list = ["O", "C"]
>>> structs = []
>>> for smi in smi_list:
>>>     g = gcd.smiles_decode(smi)
>>>     g_structs = graphs.graph_to_structure_bonds(g)
>>>     structs.extend(g_structs)
>>> # saves to a file
>>> codec_native.graph_codec_native_save("bonds.bes", structs)
>>> # saves to memory, print below
>>> encoded = codec_native.graph_codec_native_encode(structs)
>>> print("\n".join(encoded))
>>> # load the data back from file
>>> structs = codec_native.graph_codec_native_load("bonds.bes")
#GRAPH BOND
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic formal_charge
#BOND bond_ring bond_order
 -1  -1 256   4   4   1   1   1   1
 -2  -2   2   1   2   1   1   1   1
  3   3   2   1   2   1   1   1   1
  1   2   1   2
#GRAPH BOND
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic formal_charge
#BOND bond_ring bond_order
 -1  -1 256   4   4   1   1   1   1
  2   2   2   1   2   1   1   1   1
 -3  -3   2   1   2   1   1   1   1
  1   3   1   2
#GRAPH BOND
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic formal_charge
#BOND bond_ring bond_order
 -1  -1  64  16  16   1   1   1   1
 -2  -2   2   1   2   1   1   1   1
  3   3   2   1   2   1   1   1   1
  4   4   2   1   2   1   1   1   1
  5   5   2   1   2   1   1   1   1
  1   2   1   2
#GRAPH BOND
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic formal_charge
#BOND bond_ring bond_order
 -1  -1  64  16  16   1   1   1   1
  2   2   2   1   2   1   1   1   1
 -3  -3   2   1   2   1   1   1   1
  4   4   2   1   2   1   1   1   1
  5   5   2   1   2   1   1   1   1
  1   3   1   2
#GRAPH BOND
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic formal_charge
#BOND bond_ring bond_order
 -1  -1  64  16  16   1   1   1   1
  2   2   2   1   2   1   1   1   1
  3   3   2   1   2   1   1   1   1
 -4  -4   2   1   2   1   1   1   1
  5   5   2   1   2   1   1   1   1
  1   4   1   2
#GRAPH BOND
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic formal_charge
#BOND bond_ring bond_order
 -1  -1  64  16  16   1   1   1   1
  2   2   2   1   2   1   1   1   1
  3   3   2   1   2   1   1   1   1
  4   4   2   1   2   1   1   1   1
 -5  -5   2   1   2   1   1   1   1
  1   5   1   2

Replace graphs.graph_to_structure_bonds with

- graphs.graph_to_structure_atoms
- graphs.graph_to_structure_angles
- graphs.graph_to_structure_torsions
- graphs.graph_to_structure_outofplans

as needed.


