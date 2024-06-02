Labeling
========

A fundamental task in the BESMARTS package is labeling the internal coordinates
of a graph with a SMARTS pattern. Such an operation is at the core of the
SMIRNOFF force field specification, which assigns parameters based on SMARTS
patterns. The collection of SMARTS patterns in the OpenFF offxml files are a
flat list of parameters/SMARTS, and priority is given to the bottom of the
list. On this page we show how to form a simple hierarchy and then label a
molecule with it.

Hierarchies
-----------

Here, we form an explicit hierarchy of SMARTS. This can be done because the
BESMARTS package defines the hierarchical relationship using a subset relation
on the SMARTS. Details of this are given in another section. Here, it is easy
to see that the torsions supplied are all a subset of `[*:1]~[*:2]~[*:3]~[*:4]`
since this pattern matches any torsion. The following SMARTS patterns are added underneath the generic pattern,
forming a wide hierarchy.

.. code-block:: python

>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.assign.hierarchy_assign_rdkit import smarts_hierarchy_assignment_rdkit
>>> from besmarts.core import hierarchies, assignments, graphs, trees
>>> import pprint
>>> 
>>> gcd = graph_codec_rdkit()
>>> labeler = smarts_hierarchy_assignment_rdkit()
>>> 
>>> hier = trees.tree_index()
>>> 
>>> # create a root node that has no parent
>>> root = hier.node_add_below(None)
>>> root.name = "p0"
>>> 
>>> hidx = hierarchies.smarts_hierarchy(hier, {})
>>> hidx.smarts[root.index] = "[*:1]~[*:2]~[*:3]~[*:4]"
>>> 
>>> # define and add the SMARTS patterns under the root node
>>> p = [
>>>         "[*:1]~[X4:2]~[X4:3]~[*:4]",
>>>         "[*:1]~[X4:2]~[X3:3]~[*:4]",
>>>         "[*:1]~[X4:2]~[X2:3]~[*:4]",
>>>         "[*:1]~[X3:2]~[X3:3]~[*:4]",
>>>         "[*:1]~[X3:2]~[X2:3]~[*:4]",
>>>         "[*:1]~[X2:2]~[X2:3]~[*:4]",
>>> ]
>>> 
>>> for i, smarts in enumerate(p, 1):
>>>     n = hier.node_add_below(root.index)
>>>     n.name = f"p{i}"
>>>     hidx.smarts[n.index] = smarts
>>> 
>>> # Convert to a torsion hierarchy by decoding the SMARTS. 
>>> hidx = hierarchies.smarts_hierarchy_to_structure_hierarchy_torsions(hidx, gcd)
>>> 
>>> smiles = ["c1ccccc1C([NH]CC=CO)=O"]
>>> g = gcd.smiles_decode(smiles[0])
>>> 
>>> # Label the molecule
>>> lbls: assignments.smiles_assignment_group  = labeler.assign(hidx, gcd, smiles, hidx.topology)
>>> glbl = lbls.assignments[0]
>>>
>>> # make a look up table of the nodes by name
>>> nodes = {n.name: n for n in hidx.index.nodes.values()}
>>>
>>> # set the primitives to print to make it easier to read
>>> graphs.graph_set_primitives_atom(g, ["element", "connectivity_total"])
>>> for ic in graphs.graph_torsions(g):
>>>     lbl = glbl.selections.get(ic)
>>>     n = nodes[lbl]
>>>     smarts = hidx.smarts[n.index]
>>>     subg = graphs.subgraph_remove_unselected(graphs.graph_to_subgraph(g, ic))
>>>     print(f"{str(ic):20s} {lbl} {smarts} <- {gcd.smarts_encode(subg)}")
(6, 1, 2, 15)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:6]@;:[#6X3:1]@;:[#6X3:2]!@;-[#1X1:15]
(14, 1, 2, 15)       p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#1X1:14]!@;-[#6X3:1]@;:[#6X3:2]!@;-[#1X1:15]
(2, 1, 6, 5)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:2]@;:[#6X3:1]@;:[#6X3:6]@;:[#6X3:5]
(2, 1, 6, 7)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:2]@;:[#6X3:1]@;:[#6X3:6]!@;-[#6X3:7]
(3, 2, 1, 6)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:3]@;:[#6X3:2]@;:[#6X3:1]@;:[#6X3:6]
(3, 2, 1, 14)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:3]@;:[#6X3:2]@;:[#6X3:1]!@;-[#1X1:14]
(1, 2, 3, 4)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:1]@;:[#6X3:2]@;:[#6X3:3]@;:[#6X3:4]
(1, 2, 3, 16)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:1]@;:[#6X3:2]@;:[#6X3:3]!@;-[#1X1:16]
(15, 2, 3, 16)       p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#1X1:15]!@;-[#6X3:2]@;:[#6X3:3]!@;-[#1X1:16]
(4, 3, 2, 15)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:4]@;:[#6X3:3]@;:[#6X3:2]!@;-[#1X1:15]
(2, 3, 4, 5)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:2]@;:[#6X3:3]@;:[#6X3:4]@;:[#6X3:5]
(2, 3, 4, 17)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:2]@;:[#6X3:3]@;:[#6X3:4]!@;-[#1X1:17]
(16, 3, 4, 17)       p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#1X1:16]!@;-[#6X3:3]@;:[#6X3:4]!@;-[#1X1:17]
(5, 4, 3, 16)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:5]@;:[#6X3:4]@;:[#6X3:3]!@;-[#1X1:16]
(3, 4, 5, 6)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:3]@;:[#6X3:4]@;:[#6X3:5]@;:[#6X3:6]
(3, 4, 5, 18)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:3]@;:[#6X3:4]@;:[#6X3:5]!@;-[#1X1:18]
(17, 4, 5, 18)       p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#1X1:17]!@;-[#6X3:4]@;:[#6X3:5]!@;-[#1X1:18]
(6, 5, 4, 17)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:6]@;:[#6X3:5]@;:[#6X3:4]!@;-[#1X1:17]
(4, 5, 6, 7)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:4]@;:[#6X3:5]@;:[#6X3:6]!@;-[#6X3:7]
(5, 6, 1, 14)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:5]@;:[#6X3:6]@;:[#6X3:1]!@;-[#1X1:14]
(7, 6, 1, 14)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:7]!@;-[#6X3:6]@;:[#6X3:1]!@;-[#1X1:14]
(1, 6, 5, 4)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:1]@;:[#6X3:6]@;:[#6X3:5]@;:[#6X3:4]
(1, 6, 5, 18)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:1]@;:[#6X3:6]@;:[#6X3:5]!@;-[#1X1:18]
(7, 6, 5, 18)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:7]!@;-[#6X3:6]@;:[#6X3:5]!@;-[#1X1:18]
(1, 6, 7, 8)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:1]@;:[#6X3:6]!@;-[#6X3:7]!@;-[#7X3:8]
(1, 6, 7, 13)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:1]@;:[#6X3:6]!@;-[#6X3:7]!@;=[#8X1:13]
(5, 6, 7, 8)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:5]@;:[#6X3:6]!@;-[#6X3:7]!@;-[#7X3:8]
(5, 6, 7, 13)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:5]@;:[#6X3:6]!@;-[#6X3:7]!@;=[#8X1:13]
(6, 7, 8, 9)         p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:6]!@;-[#6X3:7]!@;-[#7X3:8]!@;-[#6X4:9]
(6, 7, 8, 19)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X3:6]!@;-[#6X3:7]!@;-[#7X3:8]!@;-[#1X1:19]
(13, 7, 8, 19)       p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#8X1:13]!@;=[#6X3:7]!@;-[#7X3:8]!@;-[#1X1:19]
(9, 8, 7, 13)        p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X4:9]!@;-[#7X3:8]!@;-[#6X3:7]!@;=[#8X1:13]
(7, 8, 9, 10)        p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#6X3:7]!@;-[#7X3:8]!@;-[#6X4:9]!@;-[#6X3:10]
(7, 8, 9, 20)        p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#6X3:7]!@;-[#7X3:8]!@;-[#6X4:9]!@;-[#1X1:20]
(7, 8, 9, 21)        p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#6X3:7]!@;-[#7X3:8]!@;-[#6X4:9]!@;-[#1X1:21]
(19, 8, 9, 20)       p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#1X1:19]!@;-[#7X3:8]!@;-[#6X4:9]!@;-[#1X1:20]
(19, 8, 9, 21)       p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#1X1:19]!@;-[#7X3:8]!@;-[#6X4:9]!@;-[#1X1:21]
(10, 9, 8, 19)       p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#6X3:10]!@;-[#6X4:9]!@;-[#7X3:8]!@;-[#1X1:19]
(8, 9, 10, 11)       p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#7X3:8]!@;-[#6X4:9]!@;-[#6X3:10]!@;=[#6X3:11]
(8, 9, 10, 22)       p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#7X3:8]!@;-[#6X4:9]!@;-[#6X3:10]!@;-[#1X1:22]
(20, 9, 10, 22)      p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#1X1:20]!@;-[#6X4:9]!@;-[#6X3:10]!@;-[#1X1:22]
(21, 9, 10, 22)      p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#1X1:21]!@;-[#6X4:9]!@;-[#6X3:10]!@;-[#1X1:22]
(11, 10, 9, 20)      p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#6X3:11]!@;=[#6X3:10]!@;-[#6X4:9]!@;-[#1X1:20]
(11, 10, 9, 21)      p2 [*:1]~[X4:2]~[X3:3]~[*:4] <- [#6X3:11]!@;=[#6X3:10]!@;-[#6X4:9]!@;-[#1X1:21]
(9, 10, 11, 12)      p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X4:9]!@;-[#6X3:10]!@;=[#6X3:11]!@;-[#8X2:12]
(9, 10, 11, 23)      p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#6X4:9]!@;-[#6X3:10]!@;=[#6X3:11]!@;-[#1X1:23]
(22, 10, 11, 23)     p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#1X1:22]!@;-[#6X3:10]!@;=[#6X3:11]!@;-[#1X1:23]
(12, 11, 10, 22)     p4 [*:1]~[X3:2]~[X3:3]~[*:4] <- [#8X2:12]!@;-[#6X3:11]!@;=[#6X3:10]!@;-[#1X1:22]
(10, 11, 12, 24)     p5 [*:1]~[X3:2]~[X2:3]~[*:4] <- [#6X3:10]!@;=[#6X3:11]!@;-[#8X2:12]!@;-[#1X1:24]
(23, 11, 12, 24)     p5 [*:1]~[X3:2]~[X2:3]~[*:4] <- [#1X1:23]!@;-[#6X3:11]!@;-[#8X2:12]!@;-[#1X1:24]

A few notes of importance. The labeler requires each node has its `name` set,
and each name is unique. In this example we set each SMARTS node name to pN.
Although not clear from this example, the default behavior is to process the
input list of SMILES in parallel, which is why the `labeler.assign` function
requires a list of SMILES.

Iterating the hierarchy
-----------------------

Beyond labeling a molecule with a SMARTS hierarchy, the other task of interest
is hierarchy traversal. Using the subset relation of SMARTS patterns to define the hiearchy,
the priority of such a hierarchy is top to bottom, left to right. This is neither depth-first search or breadth-first search,
and so we define it as a "dive":

.. code-block:: python

>>> hierarchies.smarts_hierarchy_print(hidx)
>>> from besmarts.core import tree_iterators
>>> def smarts_hierarchy_print(hidx):
>>>     roots = [
>>>         hidx.index.nodes[x] for x, y in hidx.index.above.items() if y is None
>>>     ]
>>>     for root in roots:
>>>         for e in tree_iterators.tree_iter_dive(hidx.index, root):
>>>             s = " " * trees.tree_index_node_depth(hidx.index, e)
>>>             print("**", e.index, s, e.name, hidx.smarts.get(e.index))
>>> 
>>> smarts_hierarchy_print(hidx)
**  0 p0 [*:1]~[*:2]~[*:3]~[*:4]
**   1 p1 [*:1]~[X4:2]~[X4:3]~[*:4]
**   2 p2 [*:1]~[X4:2]~[X3:3]~[*:4]
**   3 p3 [*:1]~[X4:2]~[X2:3]~[*:4]
**   4 p4 [*:1]~[X3:2]~[X3:3]~[*:4]
**   5 p5 [*:1]~[X3:2]~[X2:3]~[*:4]
**   6 p6 [*:1]~[X2:2]~[X2:3]~[*:4]
** 0  p0 [*:1]~[*:2]~[*:3]~[*:4]
** 1   p1 [*:1]~[X4:2]~[X4:3]~[*:4]
** 2   p2 [*:1]~[X4:2]~[X3:3]~[*:4]
** 3   p3 [*:1]~[X4:2]~[X2:3]~[*:4]
** 4   p4 [*:1]~[X3:2]~[X3:3]~[*:4]
** 5   p5 [*:1]~[X3:2]~[X2:3]~[*:4]
** 6   p6 [*:1]~[X2:2]~[X2:3]~[*:4]

As shown above, a convencience print function is provided and slightly modified
version is given which shows the use of the `tree_iter_dive` function.
