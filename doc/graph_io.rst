
Graphs
======

Graphs represent a fully described SMILES or SMARTS pattern as represented by
the SMARTS primitives that the graph codec has produced. From a graph, we examine
its internal coordinates and show how to examine each internal coordinate as a
function of depth, or distance from the internal coordinate. From a graph
perspective, internal coordinates are graph structures, and each internal
coordinate has several properties that are important for mapping.

SMARTS
------

First, we load a SMARTS patterns that is intended to describe a bond internal
coordinate:

.. code-block:: python

>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> gcd = graph_codec_rdkit()
>>> smarts = '[#6:1]~[*:2]=[#8]'
>>> g = gcd.smarts_decode(smarts)
>>> print(gcd.smarts_encode(smarts))
[#6:1]~[*:2]=[#8:3]

We have just created a graph, but specifically a subgraph because the SMARTS
pattern has indexed atoms. All atoms are indexed if at least one atom is
indexed. In this example, the oxygen is index 3 even though no index was
provided in the input. The index is preserved in the graph, as shown below:

.. code-block:: python

>>> g.nodes[1]
bechem(element:  0b1000000 hydrogen: ~0b0 connectivity_total: ~0b0 connectivity_ring: ~0b0 ring_smallest: ~0b0 a
romatic: ~0b0 formal_charge: ~0b0)

which corresponds to the carbon atom in the pattern. *NOTE*: because of this feature, 
the graphs are 1-indexed data structures. 

We can also save the decoded results to a file, and read them later to avoid
decoding twice:

.. code-block:: python

>>> from besmarts.codecs import codec_native
>>> codec_native.graph_codec_native_save("g.bes", [g])
>>> g = native.graph_codec_native_load("g.bes")[0]
>>> print(gcd.smarts_encode(g))
[#6:1]~[*:2]=[#8:3]

Letting one start quickly and from a known point and decouples the perception
from the rest of the code. The load and save functions process list of graphs,
and so a bes file can contain multiple graphs.

SMILES
------

Now we move onto SMILES:

.. code-block:: python

>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> gcd = graph_codec_rdkit()
>>> g = 'CCC=O'
>>> g = gcd.smiles_decode(g)
>>> print(gcd.smiles_encode(g))
C([H])([H])([H])C([H])([H])C([H])=O

Here the graphs generated from the SMARTS code path and the SMILES code path are
the same. The only difference is that the SMARTS was tagged, and was therefore
loaded as a subgraph. It is possible to encode SMILES without hydrogen using
either `gcd.smiles_config.strip_hydrogen` or `graphs.graph_remove_hydrogen`.

Primitives
----------

As shown above, `g.nodes[1]` returned a `bechem` object, which is essentially
a list of primitives that the `gcd` graph codec generated. The rdkit graph codec
has a default set of primitives that it will decode from a SMILES/SMARTS. We can
modify this such that only a specific set of primitives is generated, or we can
provide a specific order of primitives for subsequent SMARTS printing. Alternatively,
we may let the graph codec generate a list of primitives, but then subsequently
operate only on a selection of the primitives. Below, we define the set of primitives
and then only print SMARTS patterns using the element and bond order primitives:

.. code-block:: python

>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.core import graphs
>>> from besmarts.core.primitives import primitive_key
>>> # the default primitives
>>> atom_primitives = tuple(
>>>     (
>>>         "element",
>>>         "hydrogen",
>>>         "connectivity_total",
>>>         "connectivity_ring",
>>>         "ring_smallest",
>>>         "aromatic",
>>>         "formal_charge",
>>>     )
>>> )
>>> bond_primitives = tuple(
>>>     (
>>>         "bond_ring",
>>>         "bond_order",
>>>     )
>>> )
>>> gcd = graph_codec_rdkit(
>>>     atom_primitives=atom_primitives, bond_primitives=bond_primitives
>>> )
>>> g = 'CC'
>>> g = gcd.smiles_decode(g)
>>> graphs.graph_set_primitives_atom(g, [primitive_key.ELEMENT])
>>> graphs.graph_set_primitives_bond(g, [primitive_key.BOND_ORDER])
>>> print(gcd.smarts_encode(g))
[#6](-[#1])(-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#1]
>>> graphs.graph_set_primitives_atom(g, atom_primitives[::-1])
[+0A!rx0X4H3#6](-[+0A!rx0X1H0#1])(-[+0A!rx0X1H0#1])(-[+0A!rx0X1H0#1])-[+0A!rx0X4H3#6](-[+0A!rx0X1H0#1])(-[+0A!rx0X1H0#1])-[+0A!rx0X1H0#1]

This will not delete the other primitives as show above. Subsequent operations,
such as mapping, printing, and hashing, will only use the set primitives.

Graph structures: atoms, bonds, angles, torsions, and out-of-planes
-------------------------------------------------------------------

Loading indexed SMARTS or SMILES patterns will cause the graphs to use the
subgraph data type. To describe a internal coordinate such as a bond, one needs
to represent a subgraph as a structure. Structures as subgraphs associated with
a topology which describes the be. All structures assume that the first tagged
atoms define the topology. For example, if have the SMARTS
"[#6:1]([#6])[#6:2]", we can transform it to a bond structure as follows:

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs
    from besmarts.core import topology

    gcd = graph_codec_rdkit()

    g = "[#6:1]([#6])[#6:2]"
    g = gcd.smarts_decode(g)
    topo = topology.bond_topology()
    bond = graphs.subgraph_to_structure(g, topo)

Note that loading a SMARTS always produces a subgraph data type, and so the
untagged atom in the pattern will also be part of the subgraph. This would
result in a subgraph of 3 atoms.

For SMILES, a subgraph is only returned if there are mapped atoms. Otherwise,
a simple graph object is returned. It is possible to generate the necessary
geometries from the graph, and make structures from them. First, we show how
to identify the indices for each of the structure types:

.. code-block:: python

>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.core import graphs 
>>> gcd = graph_codec_rdkit()
>>> g = 'CC'
>>> g = gcd.smiles_decode(g)
>>> print("Atoms:")
>>> for indices in graphs.graph_atoms(g):
>>>     print(indices)
>>> print("Bonds:")
>>> for indices in graphs.graph_bonds(g):
>>>     print(indices)
>>> print("Angles:")
>>> for indices in graphs.graph_angles(g):
>>>     print(indices)
>>> print("Torsions:")
>>> for indices in graphs.graph_torsions(g):
>>>     print(indices)
>>> print("Out-of-Planes:")
>>> for indices in graphs.graph_outofplanes(g):
>>>     print(indices)
Atoms:
(1,)
(2,)
(3,)
(4,)
(5,)
(6,)
(7,)
(8,)
Bonds:
(1, 2)
(1, 3)
(1, 4)
(1, 5)
(2, 6)
(2, 7)
(2, 8)
Angles:
(2, 1, 3)
(2, 1, 4)
(2, 1, 5)
(3, 1, 4)
(3, 1, 5)
(4, 1, 5)
(1, 2, 6)
(1, 2, 7)
(1, 2, 8)
(6, 2, 7)
(6, 2, 8)
(7, 2, 8)
Torsions:
(3, 1, 2, 6)
(3, 1, 2, 7)
(3, 1, 2, 8)
(4, 1, 2, 6)
(4, 1, 2, 7)
(4, 1, 2, 8)
(5, 1, 2, 6)
(5, 1, 2, 7)
(5, 1, 2, 8)
Out-of-Planes:
(2, 1, 3, 4)
(2, 1, 3, 5)
(2, 1, 4, 5)
(3, 1, 4, 5)
(1, 2, 6, 7)
(1, 2, 6, 8)
(1, 2, 7, 8)
(6, 2, 7, 8)

These indices are sorted in a special order. The first indices describe the atoms
of the structure (bond, angle, etc). Because these topologies are invariant to
certain permtuations (e.g. bond is 1-2 or 2-1), the following functions can be
used to get the canonical ordering:

.. code-block:: python

>>> for indices in graphs.graph_bonds(g):
>>>     geometry.bond(indices[::-1])
>>> for indices in graphs.graph_angles(g):
>>>     geometry.angle(indices[::-1])
>>> for indices in graphs.graph_torsions(g):
>>>     geometry.torsion(indices[::-1])
>>> for indices in graphs.graph_outofplanes(g):
>>>     geometry.outofplane(idx)

*NOTE* The central atom in out-of-planes is always index 1 (0-based).

Next, we get a structure for each type of topology. This will create a new graph
for each internal coordinate. Some atoms are indistinguishable, and so we would
create duplicate graphs. Here, we print on the unique graphs using the `set` 
data structure.

.. code-block:: python

>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.core import graphs 
>>> gcd = graph_codec_rdkit()
>>> g = 'CC'
>>> g = gcd.smiles_decode(g)
>>> print("Atoms:")
>>> for struct in set(graphs.graph_to_structure_atoms(g)):
>>>     print(gcd.smarts_encode(struct))
>>> print("Bonds:")
>>> for struct in set(graphs.graph_to_structure_bonds(g)):
>>>     print(gcd.smarts_encode(struct))
>>> print("Angles:")
>>> for struct in set(graphs.graph_to_structure_angles(g)):
>>>     print(gcd.smarts_encode(struct))
>>> print("Torsions:")
>>> for struct in set(graphs.graph_to_structure_torsions(g)):
>>>     print(gcd.smarts_encode(struct))
>>> print("Out-of-Planes:")
>>> for struct in set(graphs.graph_to_structure_outofplanes(g)):
>>>     print(gcd.smarts_encode(struct))
Atoms:
[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#6H3X4x0!rA+0](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
[#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#6H3X4x0!rA+0](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
Bonds:
[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#6H3X4x0!rA+0:2](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0:3])(!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#6H3X4x0!rA+0](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
Angles:
[#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0:4])(!@;-[#1H0X1x0!rA+0])!@;-[#6H3X4x0!rA+0](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
[#6H3X4x0!rA+0:2](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0:3])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
Torsions:
[#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#6H3X4x0!rA+0:2](!@;-[#1H0X1x0!rA+0:6])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
Out-of-Planes:
[#6H3X4x0!rA+0:2](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0:3])(!@;-[#1H0X1x0!rA+0:4])!@;-[#1H0X1x0!rA+0]
[#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0:4])(!@;-[#1H0X1x0!rA+0:5])!@;-[#6H3X4x0!rA+0](!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]

As expected, there are only 2 unique atoms, bonds, angles, torsions, and out-of-planes for ethane.
A few words about the returned structures. The structures contain all nodes and edges from the original
graph, but only the primary nodes are selected. This results in a SMARTS pattern that extends the entire 
molecule, but only the internal coordinate is tagged. This is the default behavior. Nearly all structure methods
only operate on the selection; the unselected atoms in the graph serve other purposes and are kept for certain functionality such as mapping and searching.
We can prune the graphs to only contain the primary nodes in the internal coordinate: 

.. code-block:: python
   
>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.core import graphs, mapper
>>> gcd = graph_codec_rdkit()
>>> g = 'CC'
>>> g = gcd.smiles_decode(g)
>>> print("Atoms:")
>>> for struct in mapperset(graphs.graph_to_structure_atoms(g)):
>>>     print(gcd.smarts_encode(graphs.structure_up_to_depth(struct, 0)))
>>> print("Bonds:")
>>> for struct in set(graphs.graph_to_structure_bonds(g)):
>>>     print(gcd.smarts_encode(graphs.structure_up_to_depth(struct, 0)))
>>> print("Angles:")
>>> for struct in set(graphs.graph_to_structure_angles(g)):
>>>     print(gcd.smarts_encode(graphs.structure_up_to_depth(struct, 0)))
>>> print("Torsions:")
>>> for struct in set(graphs.graph_to_structure_torsions(g)):
>>>     print(gcd.smarts_encode(graphs.structure_up_to_depth(struct, 0)))
>>> print("Out-of-Planes:")
>>> for struct in set(graphs.graph_to_structure_outofplanes(g)):
>>>     print(gcd.smarts_encode(graphs.structure_up_to_depth(struct, 0)))
Atoms:
[#1H0X1x0!rA+0:3]
[#6H3X4x0!rA+0:1]
Bonds:
[#6H3X4x0!rA+0:1]!@;-[#6H3X4x0!rA+0:2]
[#6H3X4x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
Angles:
[#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0:1]!@;-[#1H0X1x0!rA+0:4]
[#6H3X4x0!rA+0:2]!@;-[#6H3X4x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
Torsions:
[#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0:1]!@;-[#6H3X4x0!rA+0:2]!@;-[#1H0X1x0!rA+0:6]
Out-of-Planes:
[#6H3X4x0!rA+0:2]!@;-[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0:3])!@;-[#1H0X1x0!rA+0:4]
[#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0:1](!@;-[#1H0X1x0!rA+0:4])!@;-[#1H0X1x0!rA+0:5]

Keep in mind that the `graphs.graph_to_structure_*` functions remove nodes and
therefore change the graph. Also, since the `graphs.structure_up_to_depth` function
only operates on the selection, depths greater than 0 will not produce larger graphs. Instead,
one should extend the selection beforehand:

.. code-block:: python

>>> from besmarts.core import mapper, configs
>>> cfg = configs.smarts_extender_config(9,9,True)
>>> for struct in set(graphs.graph_to_structure_atoms(g)):
>>>     # modifies structures in-place (just structure.select)
>>>     mapper.mapper_smarts_extend(cfg, [struct])
>>>     print("Atom structure:")
>>>     for d in range(0, graphs.structure_max_depth(struct)):
>>>         s = graphs.structure_up_to_depth(struct, d)
>>>         print(f"Depth {d}: {gcd.smarts_encode(s)}")
Atom structure:
Depth 0: [#1H0X1x0!rA+0:3]
Depth 1: [#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0]
Depth 2: [#1H0X1x0!rA+0:3]!@;-[#6H3X4x0!rA+0](!@;-[#6H3X4x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]
Atom structure:
Depth 0: [#6H3X4x0!rA+0:1]
Depth 1: [#6H3X4x0!rA+0:1](!@;-[#6H3X4x0!rA+0])(!@;-[#1H0X1x0!rA+0])(!@;-[#1H0X1x0!rA+0])!@;-[#1H0X1x0!rA+0]

At this point we have created a graph from a SMILES and examined the SMARTS
patterns of all generated internal coordinates and showed how to generate
SMARTS as a function of depth.
