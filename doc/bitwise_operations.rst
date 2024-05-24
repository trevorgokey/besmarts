
Bitwise Operations
==================

SMILES/SMARTS binary operations
-------------------------------

Binary operations only work on structures, which means they are subgraphs with
topology. To combine two atom SMARTS into a single pattern, use the following:

.. code-block:: python

>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.core import graphs 
>>> from besmarts.core import mapper
>>>
>>> gcd = graph_codec_rdkit()
>>> g = 'CCC=O'
>>> g = gcd.smiles_decode(g)
>>>
>>> A, B = graphs.graph_to_structure_atoms(g)[:2]
>>> C = mapper.union(A, B)
>>> print(gcd.smarts_encode(graphs.structure_remove_unselected(A)))
>>> print(gcd.smarts_encode(graphs.structure_remove_unselected(B)))
>>> print(gcd.smarts_encode(C))
[#6H3X4x0!rA+0:1]
[#6H2X4x0!rA+0:2]
[#6;H2,H3;X4;x0;!r;A;+0:1]
    
The atom structures always start with the smallest subset possible, i.e. the 
single atom. Binary operations only use the subgraph, and so the result will
only have one atom in the SMARTS pattern. We can change this by extending the
subgraphs manually:

.. code-block:: python

>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.core import graphs 
>>> from besmarts.core import mapper
>>> from besmarts.core import configs
>>>
>>> g = 'CCC=O'
>>> g = gcd.smiles_decode(g)
>>>
>>> A, B = graphs.graph_to_structure_atoms(g)[:2]:
>>>
>>> min_depth = 1
>>> max_depth = 1
>>> hydrogen = True
>>> config = configs.smarts_extender_config(min_depth, max_depth, hydrogen)
>>> modified = mapper.mapper_smarts_extend(config, [A, B])
>>>
>>> C = mapper.union(A, B)
>>> print(gcd.smarts_encode(graphs.structure_remove_unselected(A)))
>>> print(gcd.smarts_encode(graphs.structure_remove_unselected(B)))
>>> print(gcd.smarts_encode(C))
[#6H3X4x0!rA+0:1](!@;-[#6H2X4x0!rA+0])(!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
[#6H2X4x0!rA+0:2](!@;-[#6H3X4x0!rA+0])(!@;-[#6H1X3x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
[#6;H2,H3;X4;x0;!r;A;+0:1](!@;-[#6;H2,H3;X4;x0;!r;A;+0])(!@;-[#1H0X1x0!rA+0])(!@;-[#1,#6;H0,H1;X1,X3;x0;!r;A;+0])!@;-[#1H0X1x0!rA+0]

And now A and B were extended to exactly a depth of 1, including hydrogen (if 
present), and so the result will include the additional environment of the 
two atoms.

The other operations of interest are

- mapper.intersection
- mapper.xor
- mapper.subtract

SMARTS iteration
----------------

After combining a list of structures, a common task is iterating the SMARTS
information in the resulting structure. Here is how to iterate the bits of the
previous result stored in `C`:

.. code-block:: python

>>> from besmarts.core import graph_visitors
>>> for bit in graph_visitors.structure_iter_bits(C, skip_ones=True, iter_inverse=True):
>>>     print(gcd.smarts_encode(bit))
[_H2_____:1](_;_[_______])(_;_[_______])(_;_[_______])_;_[_______]
[_!H2_____:1](_;_[_______])(_;_[_______])(_;_[_______])_;_[_______]
[_H3_____:1](_;_[_______])(_;_[_______])(_;_[_______])_;_[_______]
[_!H3_____:1](_;_[_______])(_;_[_______])(_;_[_______])_;_[_______]
[_______:1](_;_[_H2_____])(_;_[_______])(_;_[_______])_;_[_______]
[_______:1](_;_[_!H2_____])(_;_[_______])(_;_[_______])_;_[_______]
[_______:1](_;_[_H3_____])(_;_[_______])(_;_[_______])_;_[_______]
[_______:1](_;_[_!H3_____])(_;_[_______])(_;_[_______])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[#1______])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[!#1______])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[#6______])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[!#6______])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[_H0_____])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[_!H0_____])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[_H1_____])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[_!H1_____])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[__X1____])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[__!X1____])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[__X3____])_;_[_______]
[_______:1](_;_[_______])(_;_[_______])(_;_[__!X3____])_;_[_______]

