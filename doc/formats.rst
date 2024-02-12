File formats
============

Configuration files
-------------------

Config files are in YAML

Graphs
------

BESMARTS native format (bes)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#GRAPH [structure]

#ATOM <primitive, primitive, ...>

#BOND <primitive, primitive, ...>

<i, j> <value, value, ...>

Where structure is an optional topology definition, and can be:

- ATOM
- BOND
- ANGLE
- DIHEDRAL
- IMPROPER

For #ATOM and #BOND, a list of primitives is specified, which both specified
the order of the primitive values defined below, and which codecs are used.

The remaining lines describe the graph. The pair <i, j> specifies the node or 
edge. The remaining values are the primitive values and are encoded
as two's compliment. If i and j are the same, then the the primitive values
describe a node/atom, otherwise they describe a edge/bond. If i and j are the
same and are negative, then they are part of the subgraph. Note that, for 
structures, the topology assumes that the structure starts with the first N
nodes. For example, an angle structure assumes the first 3 nodes form the angle,
and must be negative as they are part of the subgraph.

A single file can contain multiple graphs, and are delimited by the #GRAPH
directive. The #ATOM and #BOND directives must come before the nodes and edges,
but can be in any order.

BESMARTS native format in JSON (bes.json)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The json version follows the above description, except that keys are

- GRAPH : None or str
- ATOM : List[str]
- BOND : List[str]
- (i,j) : List[int]

and the values correspond to the format described above. Note that this schema 
also applies to the native format.

BESMARTS mapped format (bem)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mapped format is similar to the graph format, except that the node and edge
lines are combined in a single line per graph. Each line, or graph, must have 
the same number of primitive values. Each column assumes that they are mapped.
For example, the lines

-1 -1 2 4 6 8 -2 -2 1 2 4 6 1 2 2 4

-2 -2 1 1 1 1 -1 -1 2 2 2 2 1 2 1 1

indicate two graphs with two nodes, with both nodes as part of the subgraph. 
The nodes have 4 primitves, and the edge has 2 primitives. The mapping of these
graphs is 1->2 and 2->1.

BESMARTS mapped format in JSON (bem.json)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The JSON format for mapped graphs follows a similar definition as the graph
format, but instead of (i, j) tuples, we have

- GRAPH : None or str
- ATOM : List[str]
- BOND : List[str]
- MAP : List[int]

