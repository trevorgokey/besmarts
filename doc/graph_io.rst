
Graphs
======

Reading and writing SMARTS
--------------------------

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    gcd = graph_codec_rdkit()
    smarts = '[#6:1]~[*:2]=[#8]'
    g = gcd.smarts_decode(smarts)

We have just created a graph, but specifically a subgraph because the SMARTS 
pattern has indexed atoms. The index is preserved in the graph:

.. code-block:: python

    g.nodes[1]

corresponds to the carbon atom in the pattern. *NOTE* because of this feature, 
the graphs are 1-indexed data structures. 

We can also save the decoded results to a file, and read them later
to avoid decoding twice:


.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.codecs import native
    gcd = graph_codec_rdkit()
    g = '[#6:1]~[*:2]=[#8]'
    g = gcd.smarts_decode(g)
    native.graph_codec_native_save("g.bes", [g])

Letting one start quickly and from a known point:

.. code-block:: python

    from besmarts.codecs import native
    g = native.graph_codec_native_load("g.bes")[0]

The load and save functions process list of graphs, and so a bes file can contain
multiple graphs.

Reading and writing SMILES
-----------------------------

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit

    gcd = graph_codec_rdkit()

    g = 'CCC=O'
    g = gcd.smiles_decode(g)
    g = gcd.smiles_encode(g)

Here the graphs generated from the SMARTS code path and the SMILES code path are
the same. The only difference is that the SMARTS was tagged, and was therefore
loaded as a subgraph. The graphs can be saved and loaded the same way shown 
above for the SMARTS case.

Graph structures: atoms, bonds, angles, torsions, and impropers
---------------------------------------------------------------

Loading indexed SMARTS or SMILES patterns will cause the graphs to use the
subgraph data type. To build a graph structure, one needs a subgraph definition
with a topology which describes the primary atoms in the subgraph that represent
the topology. All structures assume that the first tagged atoms define the 
topology. For example, if have the SMARTS "[#6:1]([#6])[#6:2]", we can transform
it to a bond structure as follows:

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs
    from besmarts.core import topology

    gcd = graph_codec_rdkit()

    g = "[#6:1]([#6])[#6:2]"
    g = gcd.smarts_decode(g)
    topology = topology.bond_topology()
    bond = graphs.subgraph_to_structure(g, topology)

Note that loading a SMARTS always produces a subgraph data type, and so the
untagged atom in the pattern will also be part of the subgraph. This would
result in a subgraph of 3 atoms.

For SMILES, a subgraph is only returned if there are mapped atoms. Otherwise,
a simple graph object is returned. It is possible to generate the necessary
geometries from the graph, and make structures from them. First, we show how
to identify the indices for each of the structure types:


.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs 

    g = 'CCC=O'
    g = gcd.smiles_decode(g)

    for indices in graphs.graph_atoms(g):
        print(indices)

    for indices in graphs.graph_bonds(g):
        print(indices)

    for indices in graphs.graph_angles(g):
        print(indices)

    for indices in graphs.graph_dihedrals(g):
        print(indices)

    for indices in graphs.graph_impropers(g):
        print(indices)

These indices are sorted in a special order. The first indices describe the atoms
of the structure (bond, angle, etc). Because these topologies are invarient to
certain permtuations (e.g. bond is 1-2 or 2-1), the following functions can be
used to get the canonical ordering:

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs 
    from besmarts.core import geometry

    g = 'CCC=O'
    g = gcd.smiles_decode(g)


    for indices in graphs.graph_bonds(g):
        for idx in indices:
            print(geometry.bond(idx[::-1]))

    for indices in graphs.graph_angles(g):
        for idx in indices:
            print(geometry.angle(idx[::-1]))

    for indices in graphs.graph_dihedrals(g):
        for idx in indices:
            print(geometry.dihedral(idx[::-1]))

    for indices in graphs.graph_impropers(g):
        for idx in indices:
            print(geometry.improper(idx))

*NOTE* The central atom in impropers is always index 1.

Getting a structure for each type of topology:

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs 

    g = 'CCC=O'
    g = gcd.smiles_decode(g)

    for struct in graphs.graph_to_structure_atoms(g):
        print(gcd.smarts_encode(struct))

    for struct in graphs.graph_to_structure_bonds(g):
        print(gcd.smarts_encode(struct))

    for struct in graphs.graph_to_structure_angles(g):
        print(gcd.smarts_encode(struct))

    for struct in graphs.graph_to_structure_dihedrals(g):
        print(gcd.smarts_encode(struct))

    for struct in graphs.graph_to_structure_impropers(g):
        print(gcd.smarts_encode(struct))


SMILES/SMARTS binary operations
-------------------------------

Binary operations only work on structures, which means they are subgraphs with
topology. To combine two atom SMARTS into a single pattern, use the following:

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs 
    from besmarts.core import mapper

    g = 'CCC=O'
    g = gcd.smiles_decode(g)

    A, B = graphs.graph_to_structure_atoms(g)[:2]:
    C = mapper.union(A, B)
    
The atom structures always start with the smallest subset possible, i.e. the 
single atom. Binary operations only use the subgraph, and so the result will
only have one atom in the SMARTS pattern. We can change this by extending the
subgraphs manually:

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs 
    from besmarts.core import mapper
    from besmarts.core import configs

    g = 'CCC=O'
    g = gcd.smiles_decode(g)

    A, B = graphs.graph_to_structure_atoms(g)[:2]:

    min_depth = 2
    max_depth = 2
    hydrogen = True
    config = configs.smarts_extender_config(min_depth, max_depth, hydrogen)
    modified = mapper.mapper_smarts_extend(config, [A, B])

    C = mapper.union(A, B)

And now A and B were extended to exactly a depth of 2, including hydrogen (if 
present), and so the result will include the additional environment of the 
two atoms.

Now if we want perform a union on many structures, it is best to use a reference
structure. Otherwise, the results can become sensitive to the order in which
the operation was performed. 

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs 
    from besmarts.core import mapper
    from besmarts.core import configs

    g = 'CCC=O'
    g = gcd.smiles_decode(g)

    atoms = graphs.graph_to_structure_atoms(g):

    min_depth = 2
    max_depth = 2
    hydrogen = True
    config = configs.smarts_extender_config(min_depth, max_depth, hydrogen)
    modified = mapper.mapper_smarts_extend(config, bonds)

    domain = atoms[0]
    T = mapper.structure_mapper(domain)
    for atom in atoms:
        T.add(atom)
    C = mapper.mapper_union(T)

Adding another atom to the mapper can potentially change the shape of the graphs
if nodes can be modified. In such case, the structure mapper will keep track of
which nodes were added and additionally add nodes to all mapped structures such
that all mappings are bijective. 

Reading and writing mapped SMILES/SMARTS
-------------------------------------------

Mapping can be an expensive process, and so it is a good idea to save the results
if needed. Here is how to do it:

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import graphs 
    from besmarts.core import mapper
    from besmarts.core import configs

    g = 'CCC=O'
    g = gcd.smiles_decode(g)

    atoms = graphs.graph_to_structure_atoms(g):

    min_depth = 2
    max_depth = 2
    hydrogen = True
    config = configs.smarts_extender_config(min_depth, max_depth, hydrogen)
    modified = mapper.mapper_smarts_extend(config, bonds)

    domain = atoms[0]
    T = mapper.structure_mapper(domain)
    for atom in atoms:
        T.add(atom)

    mapper.mapper_save("g.bem", T)
    T = mapper.mapper_load("g.bem")

See file formats to read more about this files structure.
