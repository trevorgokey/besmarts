Data formats
------------

A few definitions:

Graph: the entire "system" which comprises all interactions in the graph

There are two primary data formats. Using the example of positions with the
same data density:

assignments:
    {(select): [[xyz1, xyz2]]}
    # e.g. {(1,): [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]}
graph_db_rows:
    {frame1: {select: [xyz1]}, frame2: {select: [xyz2]}}
    # e.g. 0: {(1,): [1.0, 2.0, 3.0]}, 1: {[4.0, 5.0, 6.0]}}

Some important requirements:
1. select is a tuple for assignments but not rows
2. The data is always a vector of ints|floats|strs

The graph_db_structs is more suitable for fitting force fields because it is
more flexible in how to associate different data. For example, it is more
difficult to associate an arbitrary hessian to an arbitrary conformation. With
the graph_db format, this is accomplished by using the same frame index.


