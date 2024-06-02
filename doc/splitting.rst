
Splitting
=========

The goal of SMARTS clustering is to find a SMARTS pattern that discriminates 
a group of SMARTS from another group. As such, we are interested in finding
a SMARTS pattern which partitions a set of SMARTS patterns. Since partitioning
is an expensive problem all by itself. We offer two ways to look for SMARTS
clusters:

    1. Given a set of SMARTS, generate the partitions which minimize the number of SMARTS edits
    2. Given an explicit partition, find satisfactory SMARTS pattern that minimizes the number of partition edits

In the first case, called numerical search, we are interested in all of the
"adjacent" SMARTS patterns and the partitions that they induce. In the second
case, called analytic search, we are seeking a specific partitioning and want
to know a SMARTS pattern that will descriminate between the two. Note that in
both cases we allow a fuzzy result by trying to minimize the number of edits;
the SMARTS edits in the first case and the partition edits in the second case.

Numerical/Iterative search
--------------------------

The following is the code to perform the numerical search.

.. code-block:: python

>>> from besmarts.core import graphs
>>> from besmarts.core import configs
>>> from besmarts.core import topology
>>> from besmarts.core import splits
>>> from besmarts.core import codecs
>>> from besmarts.core import configs
>>> from besmarts.core import compute
>>> 
>>> from besmarts.core.primitives import primitive_key
>>> 
>>> # use the RDKit plugin
>>> from besmarts.codecs import codec_rdkit
>>> 
>>> configs.processors = 1
>>> 
>>> prims = (primitive_key.ELEMENT, primitive_key.HYDROGEN), (
>>>     primitive_key.BOND_ORDER,
>>> )
>>> 
>>> # use all available primitives
>>> # prims = (None, None)
>>> 
>>> gcd = codec_rdkit.graph_codec_rdkit(*prims)
>>> gcd.smiles_config.strip_hydrogen = False
>>> icd = codecs.intvec_codec(
>>>     gcd.primitive_codecs,
>>>     gcd.atom_primitives,
>>>     gcd.bond_primitives
>>> )
>>> 
>>> smi = "CCO"
>>> 
>>> G = {0: gcd.smiles_decode(smi)}
>>> ic_list = [s for s in graphs.graph_to_structure_bonds(G[0])]
>>> selections = [(i, x) for i in G for x in graphs.graph_bonds(G[i])]
>>> G[0] = icd.graph_encode(G[0])
>>> 
>>> # set these all to 1 to split on neighbors too
>>> branch_min = 0
>>> branch_limit = 0
>>> branch_depth_min = 0
>>> branch_depth = 0
>>> 
>>> bit_depth_min = 1
>>> bit_depth_max = 1
>>> 
>>> splitter = configs.smarts_splitter_config(
>>>     bit_depth_min,
>>>     bit_depth_max,
>>>     branch_min,
>>>     branch_limit,
>>>     branch_depth_min,
>>>     branch_depth,
>>>     unique=False,
>>>     return_matches=True,
>>>     max_splits=0,
>>>     split_general=True,
>>>     split_specific=True,
>>>     unique_complements=False,
>>>     unique_complements_prefer_min=True,
>>> )
>>> 
>>> # for this to work, we need to extend our graphs to at least the depth of S0
>>> extender = configs.smarts_extender_config(branch_depth, branch_depth, True)
>>> graphs.structure_extend(extender, ic_list)
>>> ic_list = [graphs.structure_remove_unselected(g) for g in ic_list]
>>> 
>>> S0 = gcd.smarts_decode("[*:1]~[*:2]")
>>> S0 = graphs.structure(S0.nodes, S0.edges, (1, 2), topology.bond)
>>> 
>>> print("Structures")
>>> for i, f in enumerate(ic_list):
>>>     print(i, gcd.smarts_encode(f))
>>> 
>>> configs.remote_compute_enable = False
>>> 
>>> wq = compute.workqueue_local("127.0.0.1", 63210)
>>> results: splits.split_return_type = splits.split_structures_distributed(splitter, S0, G, selections, wq, icd)
>>> wq.close()
Structures:
0 [#6H3:1]-[#6H2:2]
1 [#6H2:2]-[#8H1:3]
2 [#6H3:1]-[#1H0:4]
3 [#6H3:1]-[#1H0:5]
4 [#6H3:1]-[#1H0:6]
5 [#6H2:2]-[#1H0:7]
6 [#6H2:2]-[#1H0:8]
7 [#8H1:3]-[#1H0:9]
2024-05-23 13:26:37.017948 Generating splits
2024-05-23 13:26:37.018074 Union merging=8
2024-05-23 13:26:37.053296 Union merging=1
2024-05-23 13:26:37.054224 Generating single splits
2024-05-23 13:26:37.054898 Generated 22 splits
BIT [#6_:1]_[__:2]
BIT [!#6_:1]_[__:2]
BIT [#8_:1]_[__:2]
BIT [!#8_:1]_[__:2]
BIT [_H1:1]_[__:2]
BIT [_!H1:1]_[__:2]
BIT [_H2:1]_[__:2]
BIT [_!H2:1]_[__:2]
BIT [_H3:1]_[__:2]
BIT [_!H3:1]_[__:2]
BIT [__:1]_[#1_:2]
BIT [__:1]_[!#1_:2]
BIT [__:1]_[_H0:2]
BIT [__:1]_[_!H0:2]
2024-05-23 13:26:37.063593 Building tasks
Started local workspace on ('127.0.0.1', 42357)
workspace listening on local host. Remote connections prohibited.
2024-05-23 13:26:37.608606 Searching atoms=2 data=8 bit_depth=1/1 b_j=1/28 hits=0            
2024-05-23 13:26:37.614127 Searching atoms=2 data=8 bit_depth=1/1 b_j=3/28 hits=2            
2024-05-23 13:26:37.619337 Searching atoms=2 data=8 bit_depth=1/1 b_j=6/28 hits=2            
2024-05-23 13:26:37.624520 Searching atoms=2 data=8 bit_depth=1/1 b_j=9/28 hits=3            
2024-05-23 13:26:37.629774 Searching atoms=2 data=8 bit_depth=1/1 b_j=12/28 hits=4            
2024-05-23 13:26:37.634958 Searching atoms=2 data=8 bit_depth=1/1 b_j=14/28 hits=4            
2024-05-23 13:26:37.640125 Searching atoms=2 data=8 bit_depth=1/1 b_j=17/28 hits=5            
2024-05-23 13:26:37.645391 Searching atoms=2 data=8 bit_depth=1/1 b_j=20/28 hits=6            
2024-05-23 13:26:37.650726 Searching atoms=2 data=8 bit_depth=1/1 b_j=23/28 hits=7            
2024-05-23 13:26:37.655815 Searching atoms=2 data=8 bit_depth=1/1 b_j=26/28 hits=7            
2024-05-23 13:26:37.660988 Searching atoms=2 data=8 bit_depth=1/1 b_j=28/28 hits=8            
Progress: 100.00%        28/28
Finished: 100.00%        28/28
Removing workspace ('127.0.0.1', 42357)
Closing workspace
2024-05-23 13:26:38.048585 Calculating partitions for hits=8
Started local workspace on ('127.0.0.1', 43381)
workspace listening on local host. Remote connections prohibited.
Submitting 8 packets of work
Chunk: 100.00%         8/8
Finished: 100.00%         8/8
Removing workspace ('127.0.0.1', 43381)
Closing workspace
2024-05-23 13:26:38.631840 Searching atoms done; data=8 hits=8

The primary result is the `result` object returned from
`splits.split_structures_distributed`. Starting from the top, the example first
configures the graph codec and the primitives. For simplicity, we only decode
element, hydrogen, and bond order. This means that splitting will only occur on
this primitives. Using the other primitives will likely produce many redundant
splits because the dataset is quite small (1 molecule). Next, the dataset is
flattened into a dictionary and compressed using a intvec graph codec. This is
potentially required for very large datasets (100K molecules) as the code
maintains the dataset in memory. Next, the splitting configuration is set.
Here, we use a basic, quick search. This will find splits that only differ by 1
bit, and will only examine the two primary atoms of the bond. More options will
be discussed after the results are shown below. Next, before the split function
is called, the SMARTS pattern to be split is defined. Normally, this S0 SMARTS
will be taken from a hierarchy. In this example, we use a catch-all S0
(`[*:1]~[*:2]`) and therefore every bond matches. This is primarily done for
simplicity, otherwise the input graphs (`G`) will need to be pruned such that
only those that match S0 are kept. Lastly, we turn off remote computing; this
has little effect here other than not opening up a listening port. This results
in reduced overhead, of which is wasteful as we will be computing locally. If a
large search is performed, one can turn this to `True` and then run `python -m
besmarts.worker 127.0.0.1 63210`, where the host should be changed if the
worker is running from a different machine.

The output is as follows. The structures are shown at depth 0 which corresponds
to the depth of S0 and the depth defined in the search settings. Next, the
splits are found by first combining all structures and then enumerating all
bits. This resulted in 11 unique bits, and 22 splits since we wanted to find
both general and specific splits. This can be seen by the fact that the bit
`[#6_:1]_[__:2]` was found which would produce a specific split (one atom must
be carbon) versus its general counterpart `[!#6_:1]_[__:2]` (one atom must not
be carbon). If multiple bits were searched, the algorith would combine these
bits to produce new splits. However, this is done by using combinations and
therefore grows exponentially. The output then shows that out of the 28
possible splits, only 8 generated unique partitions. This section tries to fail
as fast as possible, and so does not perform full scans and aims to be a quick
filter. The following output then shows that a full match analysis is done on
the 8 valid splits before the result is returned.

Next, some code here is provided to as an example to examine the results.

.. code-block:: python
>>> # custom processing of results
>>> 
>>> seen = {}
>>> keep = {}
>>> 
>>> print("Results:", len(results.splits))
>>> for j, (Sj, matches, bj) in enumerate(
>>>     zip(results.splits, results.matched_idx, results.shards), 1
>>> ):
>>>     Sj = graphs.structure(Sj.nodes, Sj.edges, Sj.select, results.topology)
>>>     atoms, bits = len(Sj.select), graphs.graph_bits(Sj, maxbits=True)
>>>     matches = tuple(sorted(matches))
>>>     unmatches = tuple(
>>>         sorted([i for i in range(len(ic_list)) if i not in matches])
>>>     )
>>>     entry = tuple(sorted([matches, unmatches]))
>>>     if len(matches) > 0 and len(ic_list) != len(matches):
>>>         if entry in seen:
>>>             if (atoms, bits) < seen[entry]:
>>>                 seen[entry] = (atoms, bits)
>>>                 keep[entry] = j
>>>         else:
>>>             seen[entry] = (atoms, bits)
>>>             keep[entry] = j
>>> 
>>> unique = {}
>>> found = 0
>>> for j, (Sj, matches, bj) in enumerate(
>>>     zip(results.splits, results.matched_idx, results.shards), 1
>>> ):
>>> 
>>>     matches = tuple(matches)
>>>     l = unique.get(matches, list())
>>>     l.append((Sj, bj))
>>>     unique[matches] = l
>>> 
>>> for j, (matches, params) in enumerate(unique.items(), 1):
>>>     matches = tuple(matches)
>>>     found += 1
>>>     if splitter.return_matches:
>>>         print(
>>>             f"{found:4d}",
>>>             f"{j:4d}",
>>>             "match:",
>>>             f"{len(matches):4d}",
>>>             "unmatched:",
>>>             f"{len(ic_list) - len(matches):4d}",
>>>         )
>>>     else:
>>>         print(
>>>             f"{found:4d}",
>>>             f"{j:4d}",
>>>         )
>>>     for k, (Sj, bj) in enumerate(params, 1):
>>> 
>>>         Sj = graphs.structure(Sj.nodes, Sj.edges, Sj.select, results.topology)
>>>         print(f"   {k:2d} Sj:", gcd.smarts_encode(Sj))
>>>         # print(f"   {k:2d} Sj:    ", Sj.nodes)
>>>     if splitter.return_matches:
>>>         print("      ", matches)
>>>         for i, f in enumerate(ic_list):
>>>             if i in matches:
>>>                 print(f"{i:4d}", " -> ", f.select, gcd.smarts_encode(f))
>>>             else:
>>>                 print(f"{i:4d}", f.select, gcd.smarts_encode(f))
>>>     print("####################################")
Results: 8
   1    1 match:    7 unmatched:    1
    1 Sj: [!#6:1]~[*:2]
       (1, 2, 3, 4, 5, 6, 7)
   0 (1, 2) [#6H3:1]-[#6H2:2]
   1  ->  (2, 3) [#6H2:2]-[#8H1:3]
   2  ->  (1, 4) [#6H3:1]-[#1H0:4]
   3  ->  (1, 5) [#6H3:1]-[#1H0:5]
   4  ->  (1, 6) [#6H3:1]-[#1H0:6]
   5  ->  (2, 7) [#6H2:2]-[#1H0:7]
   6  ->  (2, 8) [#6H2:2]-[#1H0:8]
   7  ->  (3, 9) [#8H1:3]-[#1H0:9]
####################################
   2    2 match:    7 unmatched:    1
    1 Sj: [#6:1]~[*:2]
       (0, 1, 2, 3, 4, 5, 6)
   0  ->  (1, 2) [#6H3:1]-[#6H2:2]
   1  ->  (2, 3) [#6H2:2]-[#8H1:3]
   2  ->  (1, 4) [#6H3:1]-[#1H0:4]
   3  ->  (1, 5) [#6H3:1]-[#1H0:5]
   4  ->  (1, 6) [#6H3:1]-[#1H0:6]
   5  ->  (2, 7) [#6H2:2]-[#1H0:7]
   6  ->  (2, 8) [#6H2:2]-[#1H0:8]
   7 (3, 9) [#8H1:3]-[#1H0:9]
####################################
   3    3 match:    2 unmatched:    6
    1 Sj: [#8:1]~[*:2]
    2 Sj: [H1:1]~[*:2]
       (1, 7)
   0 (1, 2) [#6H3:1]-[#6H2:2]
   1  ->  (2, 3) [#6H2:2]-[#8H1:3]
   2 (1, 4) [#6H3:1]-[#1H0:4]
   3 (1, 5) [#6H3:1]-[#1H0:5]
   4 (1, 6) [#6H3:1]-[#1H0:6]
   5 (2, 7) [#6H2:2]-[#1H0:7]
   6 (2, 8) [#6H2:2]-[#1H0:8]
   7  ->  (3, 9) [#8H1:3]-[#1H0:9]
####################################
   4    4 match:    4 unmatched:    4
    1 Sj: [H2:1]~[*:2]
       (0, 1, 5, 6)
   0  ->  (1, 2) [#6H3:1]-[#6H2:2]
   1  ->  (2, 3) [#6H2:2]-[#8H1:3]
   2 (1, 4) [#6H3:1]-[#1H0:4]
   3 (1, 5) [#6H3:1]-[#1H0:5]
   4 (1, 6) [#6H3:1]-[#1H0:6]
   5  ->  (2, 7) [#6H2:2]-[#1H0:7]
   6  ->  (2, 8) [#6H2:2]-[#1H0:8]
   7 (3, 9) [#8H1:3]-[#1H0:9]
####################################
   5    5 match:    4 unmatched:    4
    1 Sj: [H3:1]~[*:2]
       (0, 2, 3, 4)
   0  ->  (1, 2) [#6H3:1]-[#6H2:2]
   1 (2, 3) [#6H2:2]-[#8H1:3]
   2  ->  (1, 4) [#6H3:1]-[#1H0:4]
   3  ->  (1, 5) [#6H3:1]-[#1H0:5]
   4  ->  (1, 6) [#6H3:1]-[#1H0:6]
   5 (2, 7) [#6H2:2]-[#1H0:7]
   6 (2, 8) [#6H2:2]-[#1H0:8]
   7 (3, 9) [#8H1:3]-[#1H0:9]
####################################
   6    6 match:    6 unmatched:    2
    1 Sj: [*:1]~[#1:2]
    2 Sj: [*:1]~[H0:2]
       (2, 3, 4, 5, 6, 7)
   0 (1, 2) [#6H3:1]-[#6H2:2]
   1 (2, 3) [#6H2:2]-[#8H1:3]
   2  ->  (1, 4) [#6H3:1]-[#1H0:4]
   3  ->  (1, 5) [#6H3:1]-[#1H0:5]
   4  ->  (1, 6) [#6H3:1]-[#1H0:6]
   5  ->  (2, 7) [#6H2:2]-[#1H0:7]
   6  ->  (2, 8) [#6H2:2]-[#1H0:8]
   7  ->  (3, 9) [#8H1:3]-[#1H0:9]
####################################

Here we see there were 8 unique partitions found using the given search
settings. The `Sj` patterns indicate that the splits all produce the same
partitioning and are therefore grouped together. This was done because
`splitter.unique = False`. Then, for each partition, the 8 structures are shown
and an arrow indicates that the structure matches the new split. In order for a
partition to be valid, it must match some, but not all structures.

Notice that even with one molecule and minimal search settings, the output is
somewhat complex. The `BESMARTS` package tries to hide most of this behind the
higher-level functions, such as clustering and force field fitting. This
example shows roughly how such functions work to produce novel SMARTS patterns.

Analytic/Direct search
----------------------

Next the direct split method is shown. As mentioned above, this approach
requires an particular partition, and then the code tries to find a SMARTS that
satisfies the partition.

.. code-block:: python

>>> from besmarts.core import graphs
>>> from besmarts.core import configs
>>> from besmarts.core import topology
>>> from besmarts.core import splits
>>> from besmarts.core import codecs
>>> from besmarts.core import configs
>>> from besmarts.core import compute
>>> 
>>> from besmarts.core.primitives import primitive_key
>>> 
>>> # use the RDKit plugin
>>> from besmarts.codecs import codec_rdkit
>>> 
>>> configs.processors = 1
>>> 
>>> 
>>> prims = (primitive_key.ELEMENT, primitive_key.HYDROGEN), (
>>>     primitive_key.BOND_ORDER,
>>> )
>>> 
>>> # use all available primitives
>>> # prims = (None, None)
>>> 
>>> gcd = codec_rdkit.graph_codec_rdkit(*prims)
>>> 
>>> ###
>>> 
>>> branch_min = 0
>>> branch_limit = 0
>>> branch_depth_min = 0
>>> branch_depth = 0
>>> 
>>> bit_depth_min = 1
>>> bit_depth_max = 1
>>> 
>>> splitter = configs.smarts_splitter_config(
>>>     bit_depth_min,
>>>     bit_depth_max,
>>>     branch_min,
>>>     branch_limit,
>>>     branch_depth_min,
>>>     branch_depth,
>>>     unique=False,
>>>     return_matches=True,
>>>     max_splits=0,
>>>     split_general=True,
>>>     split_specific=True,
>>>     unique_complements=False,
>>>     unique_complements_prefer_min=True,
>>> )
>>> 
>>> # for this to work, we need to extend our graphs to at least the depth of S0
>>> extender = configs.smarts_extender_config(branch_depth, branch_depth, True)
>>> 
>>> spec = configs.smarts_perception_config(
>>>     splitter, extender
>>> )
>>> 
>>> ###
>>> smi = "CCO"
>>> # beg = gcd.smiles_decode(smi)
>>> ###
>>> G = {0: gcd.smiles_decode(smi)}
>>> ic_list = [s for s in graphs.graph_to_structure_bonds(G[0])]
>>> # selections = [(i, x) for i in G for x in graphs.graph_bonds(G[i])]
>>> 
>>> topo = topology.bond
>>> 
>>> ###
>>> matches = (1, 7)
>>> 
>>> for i in range(len(ic_list)):
>>>     if i not in matches:
>>>         print(i, gcd.smarts_encode(ic_list[i]))
>>> for i in matches:
>>>     print(i, "->", gcd.smarts_encode(ic_list[i]))
>>> 
>>> results: splits.split_return_type = splits.split_partition(topo, spec, ic_list, matches, gcd=gcd, maxmoves=0)
0 [#6H3:1](-[#1H0])(-[#1H0])(-[#1H0])-[#6H2:2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
2 [#6H3:1](-[#1H0:4])(-[#1H0])(-[#1H0])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
3 [#6H3:1](-[#1H0])(-[#1H0:5])(-[#1H0])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
4 [#6H3:1](-[#1H0])(-[#1H0])(-[#1H0:6])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
5 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0:7])(-[#1H0])-[#8H1]-[#1H0]
6 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0])(-[#1H0:8])-[#8H1]-[#1H0]
1 -> [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0])(-[#1H0])-[#8H1:3]-[#1H0]
7 -> [#8H1:3](-[#1H0:9])-[#6H2](-[#1H0])(-[#1H0])-[#6H3](-[#1H0])(-[#1H0])-[#1H0]
LUN:  [#1,#6;H0,H2:2]-[#8H1:3]
LHS:  [__:2]-[#8H1:3]
RUN:  [#6;H2,H3:1]-[#1,#6;H0,H2:2]
RHS:  [#6_:1]-[__:2]
LHS_DIFF:  [__:2]_[#8H1:3]
LHS_INVE:  [*:2]-[#8H1:3]
RHS_DIFF:  [__:1]_[__:2]
BESTLHS:  [*:1]-[#8H1:2]

As above, the structures are printed, except the desired partition is indicated
with the arrows. We selected the two structures that have oxygen in the bond, now the goal is to find a SMARTS pattern that matches only these two.
Some informational output is shown, and at the bottom we see BESTLHS is indicated a match was found.

Below is some custom result parsing:

.. code-block:: python

>>> shards = results.value
>>> 
>>> ###
>>> removeA = shards[2]
>>> addA = shards[3]
>>> nummoves = len(removeA) + len(addA)
>>> verbose = True
>>> shard = shards[0]
>>> matches = [x for x in range(len(ic_list)) if x not in removeA and (x in matches or x in addA)]
>>> if shard is not None:
>>>     print(f"Matches only the input with {nummoves} swaps:", gcd.smarts_encode(shard))
>>>     if verbose and (removeA or addA):
>>>         print("RemoveA", removeA)
>>>         print("AddA", addA)
>>>         for i in range(len(ic_list)):
>>>             if i not in matches:
>>>                 print(i, gcd.smarts_encode(ic_list[i]))
>>>         for i in range(len(ic_list)):
>>>             if i in matches:
>>>                 print(i, "->", gcd.smarts_encode(ic_list[i]))
>>> 
>>> shard = shards[1]
>>> if shard is not None:
>>>     print(f"Matches the input complement with {nummoves} swaps:", gcd.smarts_encode(shard))
>>>     if verbose and (removeA or addA):
>>>         print("RemoveA", removeA)
>>>         print("AddA", addA)
>>>         for i in range(len(ic_list)):
>>>             if i in matches:
>>>                 print(i, gcd.smarts_encode(ic_list[i]))
>>>         for i in range(len(ic_list)):
>>>             if i not in matches:
>>>                 print(i, "->", gcd.smarts_encode(ic_list[i]))
Matches only the input with 0 swaps: [*:1]-[#8H1:2]

And so we see that we were able to find a SMARTS pattern that indeed splits the two structures. There are two concepts of interest here. First, we may want an approximate result that satisfies the matches rather than specify an exact partition. In such a case, we can set `maxmoves` to a positive integer. If no SMARTS pattern can be found that matches the exact partition, it tries to find a SMARTS pattern that would match a partition if `maxmoves` structures are included in the original partition. For example, there are 3 indestinguishable CH methyl bonds. If we specify `matches=(2,)`, we get

.. code-block:: python

>>> matches = (2,)
>>> 
>>> for i in range(len(ic_list)):
>>>     if i not in matches:
>>>         print(i, gcd.smarts_encode(ic_list[i]))
>>> for i in matches:
>>>     print(i, "->", gcd.smarts_encode(ic_list[i]))
>>> 
>>> results: splits.split_return_type = splits.split_partition(topo, spec, ic_list, matches, gcd=gcd, maxmoves=0)
0 [#6H3:1](-[#1H0])(-[#1H0])(-[#1H0])-[#6H2:2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
1 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0])(-[#1H0])-[#8H1:3]-[#1H0]
3 [#6H3:1](-[#1H0])(-[#1H0:5])(-[#1H0])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
4 [#6H3:1](-[#1H0])(-[#1H0])(-[#1H0:6])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
5 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0:7])(-[#1H0])-[#8H1]-[#1H0]
6 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0])(-[#1H0:8])-[#8H1]-[#1H0]
7 [#8H1:3](-[#1H0:9])-[#6H2](-[#1H0])(-[#1H0])-[#6H3](-[#1H0])(-[#1H0])-[#1H0]
2 -> [#6H3:1](-[#1H0:4])(-[#1H0])(-[#1H0])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
LUN:  [#6H3:1]-[#1H0:4]
LHS:  [#6H3:1]-[#1H0:4]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#6H3:1]_[#1H0:4]
LHS_INVE:  [#6H3:1]-[#1H0:4]
RHS_DIFF:  [__:1]_[__:2]

a.k.a. nothing. Now if we increase maxmoves to 2 (since we know there are two other identical structures):

.. code-block:: python

>>> matches = (2,)
>>> 
>>> for i in range(len(ic_list)):
>>>     if i not in matches:
>>>         print(i, gcd.smarts_encode(ic_list[i]))
>>> for i in matches:
>>>     print(i, "->", gcd.smarts_encode(ic_list[i]))
>>> 
>>> results: splits.split_return_type = splits.split_partition(topo, spec, ic_list, matches, gcd=gcd, maxmoves=2)
0 [#6H3:1](-[#1H0])(-[#1H0])(-[#1H0])-[#6H2:2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
1 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0])(-[#1H0])-[#8H1:3]-[#1H0]
3 [#6H3:1](-[#1H0])(-[#1H0:5])(-[#1H0])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
4 [#6H3:1](-[#1H0])(-[#1H0])(-[#1H0:6])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
5 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0:7])(-[#1H0])-[#8H1]-[#1H0]
6 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0])(-[#1H0:8])-[#8H1]-[#1H0]
7 [#8H1:3](-[#1H0:9])-[#6H2](-[#1H0])(-[#1H0])-[#6H3](-[#1H0])(-[#1H0])-[#1H0]
2 -> [#6H3:1](-[#1H0:4])(-[#1H0])(-[#1H0])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
LUN:  [#6H3:1]-[#1H0:4]
LHS:  [#6H3:1]-[#1H0:4]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#6H3:1]_[#1H0:4]
LHS_INVE:  [#6H3:1]-[#1H0:4]
RHS_DIFF:  [__:1]_[__:2]
BESTLHS:  [#6H3:1]-[#1H0:2]
Matches only the input with 2 swaps: [#6H3:1]-[#1H0:2]
RemoveA set()
AddA {3, 4}
0 [#6H3:1](-[#1H0])(-[#1H0])(-[#1H0])-[#6H2:2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
1 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0])(-[#1H0])-[#8H1:3]-[#1H0]
5 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0:7])(-[#1H0])-[#8H1]-[#1H0]
6 [#6H2:2](-[#6H3](-[#1H0])(-[#1H0])-[#1H0])(-[#1H0])(-[#1H0:8])-[#8H1]-[#1H0]
7 [#8H1:3](-[#1H0:9])-[#6H2](-[#1H0])(-[#1H0])-[#6H3](-[#1H0])(-[#1H0])-[#1H0]
2 -> [#6H3:1](-[#1H0:4])(-[#1H0])(-[#1H0])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
3 -> [#6H3:1](-[#1H0])(-[#1H0:5])(-[#1H0])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]
4 -> [#6H3:1](-[#1H0])(-[#1H0])(-[#1H0:6])-[#6H2](-[#1H0])(-[#1H0])-[#8H1]-[#1H0]

We now see that `[#6H3:1]-[#1H0:2]` is the desired, expected result. The
algorithm always returns the result with the fewest moves. Note that for a
given partition, there might be two unique results: one that matches the input
partition, and one that matches the complement. This is shown as the "LHS" (the
input partition) and "RHS" (the complement).

Note that, because this method is usually used for smaller datasets with only a
few structures, it has yet to make use of the scaling features in the numerical
search, search as graph compression and distributed compute.

Hybrid search
-------------

There are times when we have a few structures and we want to find their splits,
but using a numerical search becomes too expensive for large search spaces, and
smaller spaces find nothing. Since we also don't have a particular partition in
mind, a direct search is not possible. To overcome this, a function is provided
that will generate all partitions and perform a direct search on each:

.. code-block:: python

>>> from besmarts.core import graphs
>>> from besmarts.core import configs
>>> from besmarts.core import topology
>>> from besmarts.core import splits
>>> from besmarts.core import codecs
>>> from besmarts.core import configs
>>> from besmarts.core import compute
>>> 
>>> from besmarts.core.primitives import primitive_key
>>> 
>>> # use the RDKit plugin
>>> from besmarts.codecs import codec_rdkit
>>> 
>>> configs.processors = 1
>>> 
>>> prims = (primitive_key.ELEMENT, primitive_key.HYDROGEN), (
>>>     primitive_key.BOND_ORDER,
>>> )
>>> 
>>> # use all available primitives
>>> # prims = (None, None)
>>> 
>>> gcd = codec_rdkit.graph_codec_rdkit(*prims)
>>> 
>>> branch_min = 0
>>> branch_limit = 0
>>> branch_depth_min = 0
>>> branch_depth = 0
>>> 
>>> bit_depth_min = 1
>>> bit_depth_max = 1
>>> 
>>> splitter = configs.smarts_splitter_config(
>>>     bit_depth_min,
>>>     bit_depth_max,
>>>     branch_min,
>>>     branch_limit,
>>>     branch_depth_min,
>>>     branch_depth,
>>>     unique=False,
>>>     return_matches=True,
>>>     max_splits=0,
>>>     split_general=True,
>>>     split_specific=True,
>>>     unique_complements=False,
>>>     unique_complements_prefer_min=True,
>>> )
>>> 
>>> # for this to work, we need to extend our graphs to at least the depth of S0
>>> extender = configs.smarts_extender_config(branch_depth, branch_depth, True)
>>> 
>>> spec = configs.smarts_perception_config(
>>>     splitter, extender
>>> )
>>> 
>>> smi = "CCO"
>>> G = {0: gcd.smiles_decode(smi)}
>>> ic_list = [graphs.structure_remove_unselected(s) for s in graphs.graph_to_structure_bonds(G[0])]
>>> 
>>> topo = topology.bond
>>> 
>>> # give a unique label to each for combination generation
>>> labels = [str(i) for i in range(len(ic_list))]
>>> 
>>> # this is k in the nCk partition generation
>>> # will be limited to n//2
>>> spec.splitter.bit_search_limit = 9
>>> results: splits.split_return_type = splits.split_all_partitions(topo, spec, ic_list, labels, gcd=gcd, maxmoves=0)
>>> 
>>> shards = results.value
>>> 
>>> for j, (lhs, rhs, matched, unmatch) in enumerate(shards, 1):
>>>     print(f"###\n{j:2d} Sj: {gcd.smarts_encode(lhs)}")
>>>     for i in range(len(ic_list)):
>>>         if i not in matched:
>>>             print(i, gcd.smarts_encode(ic_list[i]))
>>>         else:
>>>             print(i, "->", gcd.smarts_encode(ic_list[i]))
>>> 
Direct on 1 combo ('0',) depth 0 0
LUN:  [#6H3:1]-[#6H2:2]
LHS:  [#6H3:1]-[#6H2:2]
RUN:  [#6,#8;!H0!H4:2]-[#1,#8;H0,H1:3]
RHS:  [__:2]-[__:3]
LHS_DIFF:  [#6H3:1]_[#6H2:2]
LHS_INVE:  [#6H3:1]-[#6H2:2]
RHS_DIFF:  [__:2]_[__:3]
BESTLHS:  [#6H3:1]-[#6H2:2]
Direct on 1 combo ('1',) depth 0 0
LUN:  [#6H2:2]-[#8H1:3]
LHS:  [#6H2:2]-[#8H1:3]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#6H2:2]_[#8H1:3]
LHS_INVE:  [#6H2:2]-[#8H1:3]
RHS_DIFF:  [__:1]_[__:2]
BESTLHS:  [#6H2:1]-[#8H1:2]
Direct on 1 combo ('2',) depth 0 0
LUN:  [#6H3:1]-[#1H0:4]
LHS:  [#6H3:1]-[#1H0:4]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#6H3:1]_[#1H0:4]
LHS_INVE:  [#6H3:1]-[#1H0:4]
RHS_DIFF:  [__:1]_[__:2]
Direct on 1 combo ('3',) depth 0 0
LUN:  [#6H3:1]-[#1H0:5]
LHS:  [#6H3:1]-[#1H0:5]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#6H3:1]_[#1H0:5]
LHS_INVE:  [#6H3:1]-[#1H0:5]
RHS_DIFF:  [__:1]_[__:2]
Direct on 1 combo ('4',) depth 0 0
LUN:  [#6H3:1]-[#1H0:6]
LHS:  [#6H3:1]-[#1H0:6]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#6H3:1]_[#1H0:6]
LHS_INVE:  [#6H3:1]-[#1H0:6]
RHS_DIFF:  [__:1]_[__:2]
Direct on 1 combo ('5',) depth 0 0
LUN:  [#6H2:2]-[#1H0:7]
LHS:  [#6H2:2]-[#1H0:7]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#6H2:2]_[#1H0:7]
LHS_INVE:  [#6H2:2]-[#1H0:7]
RHS_DIFF:  [__:1]_[__:2]
Direct on 1 combo ('6',) depth 0 0
LUN:  [#6H2:2]-[#1H0:8]
LHS:  [#6H2:2]-[#1H0:8]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#6H2:2]_[#1H0:8]
LHS_INVE:  [#6H2:2]-[#1H0:8]
RHS_DIFF:  [__:1]_[__:2]
Direct on 1 combo ('7',) depth 0 0
LUN:  [#8H1:3]-[#1H0:9]
LHS:  [#8H1:3]-[#1H0:9]
RUN:  [#6,#8;!H0!H4:1]-[#1,#6;H0,H2:2]
RHS:  [__:1]-[__:2]
LHS_DIFF:  [#8H1:3]_[#1H0:9]
LHS_INVE:  [#8H1:3]-[#1H0:9]
RHS_DIFF:  [__:1]_[__:2]
BESTLHS:  [#8H1:1]-[#1H0:2]
Direct on 2 combo ('0', '1') depth 0 0
LUN:  [#6,#8;H1,H3:1]-[#6H2:2]
LHS:  [__:1]-[#6H2:2]
RUN:  [#6,#8;!H0!H4:1]-[#1H0:4]
RHS:  [__:1]-[#1H0:4]
LHS_DIFF:  [__:1]_[#6H2:2]
LHS_INVE:  [*:1]-[#6H2:2]
RHS_DIFF:  [__:1]_[__:4]
... skipped long output...
Direct on 4 combo ('0', '5', '6', '7') depth 0 0
LUN:  [#1,#6;H0,H3:1]-[#6,#8;H1,H2:2]
LHS:  [__:1]-[__:2]
RUN:  [#6;H2,H3:2]-[#1,#8;H0,H1:3]
RHS:  [#6_:2]-[__:3]
LHS_DIFF:  [__:1]_[__:2]
RHS_DIFF:  [#6_:2]_[__:3]
RHS_INVE:  [#6:2]-[*:3]
RHS_INTR:  [#6:2]-[*:3]
###
 1 Sj: [#6H3:1]-[#6H2:2]
0 -> [#6H3:1]-[#6H2:2]
1 [#6H2:2]-[#8H1:3]
2 [#6H3:1]-[#1H0:4]
3 [#6H3:1]-[#1H0:5]
4 [#6H3:1]-[#1H0:6]
5 [#6H2:2]-[#1H0:7]
6 [#6H2:2]-[#1H0:8]
7 [#8H1:3]-[#1H0:9]
###
 2 Sj: [#6H2:1]-[#8H1:2]
0 [#6H3:1]-[#6H2:2]
1 -> [#6H2:2]-[#8H1:3]
2 [#6H3:1]-[#1H0:4]
3 [#6H3:1]-[#1H0:5]
4 [#6H3:1]-[#1H0:6]
5 [#6H2:2]-[#1H0:7]
6 [#6H2:2]-[#1H0:8]
7 [#8H1:3]-[#1H0:9]
###
 3 Sj: [#8H1:1]-[#1H0:2]
0 [#6H3:1]-[#6H2:2]
1 [#6H2:2]-[#8H1:3]
2 [#6H3:1]-[#1H0:4]
3 [#6H3:1]-[#1H0:5]
4 [#6H3:1]-[#1H0:6]
5 [#6H2:2]-[#1H0:7]
6 [#6H2:2]-[#1H0:8]
7 -> [#8H1:3]-[#1H0:9]
###
 4 Sj: [*:1]-[#8H1:2]
0 [#6H3:1]-[#6H2:2]
1 -> [#6H2:2]-[#8H1:3]
2 [#6H3:1]-[#1H0:4]
3 [#6H3:1]-[#1H0:5]
4 [#6H3:1]-[#1H0:6]
5 [#6H2:2]-[#1H0:7]
6 [#6H2:2]-[#1H0:8]
7 -> [#8H1:3]-[#1H0:9]
###
 5 Sj: [#6H2:1]-[#1H0:2]
0 [#6H3:1]-[#6H2:2]
1 [#6H2:2]-[#8H1:3]
2 [#6H3:1]-[#1H0:4]
3 [#6H3:1]-[#1H0:5]
4 [#6H3:1]-[#1H0:6]
5 -> [#6H2:2]-[#1H0:7]
6 -> [#6H2:2]-[#1H0:8]
7 [#8H1:3]-[#1H0:9]
###
 6 Sj: [#6:1]-[#1H0:2]
0 [#6H3:1]-[#6H2:2]
1 [#6H2:2]-[#8H1:3]
2 -> [#6H3:1]-[#1H0:4]
3 -> [#6H3:1]-[#1H0:5]
4 -> [#6H3:1]-[#1H0:6]
5 -> [#6H2:2]-[#1H0:7]
6 -> [#6H2:2]-[#1H0:8]
7 [#8H1:3]-[#1H0:9]
###
 7 Sj: [#6H3:1]-[#1H0:2]
0 [#6H3:1]-[#6H2:2]
1 [#6H2:2]-[#8H1:3]
2 -> [#6H3:1]-[#1H0:4]
3 -> [#6H3:1]-[#1H0:5]
4 -> [#6H3:1]-[#1H0:6]
5 [#6H2:2]-[#1H0:7]
6 [#6H2:2]-[#1H0:8]
7 [#8H1:3]-[#1H0:9]
###
 8 Sj: [*:1]-[#6H2:2]
0 -> [#6H3:1]-[#6H2:2]
1 -> [#6H2:2]-[#8H1:3]
2 [#6H3:1]-[#1H0:4]
3 [#6H3:1]-[#1H0:5]
4 [#6H3:1]-[#1H0:6]
5 -> [#6H2:2]-[#1H0:7]
6 -> [#6H2:2]-[#1H0:8]
7 [#8H1:3]-[#1H0:9]
###
 9 Sj: [#6H3:1]-[*:2]
0 -> [#6H3:1]-[#6H2:2]
1 [#6H2:2]-[#8H1:3]
2 -> [#6H3:1]-[#1H0:4]
3 -> [#6H3:1]-[#1H0:5]
4 -> [#6H3:1]-[#1H0:6]
5 [#6H2:2]-[#1H0:7]
6 [#6H2:2]-[#1H0:8]
7 [#8H1:3]-[#1H0:9]

This hybrid approach took 4.8s to find 9 splits where the numerical took 1.8s
to find 6 splits as shown above. However, the hybrid approach found splits
using a maximum of 4 bits whereas the numerical approach used only 1 bit.
Modifying the numerical search to search up to 4 bits resulted in a runtime of
29.2s and found 23 splits. 

