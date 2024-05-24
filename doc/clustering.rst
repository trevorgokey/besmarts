
Clustering
==========

Here we show how to build data-driven SMARTS hierarchies using data associated
with SMARTS patterns. There are two methods depending on what the format of the
data is since it defines what the objective is. We can either build hierarchies
based on categorical data, such as existing labels, or numerical data, such as
bond lengths. Here we show both methods


Clustering categorical data
--------------------------

In this example, we have a SMILES and we have labeled the bonds with either "a"
or "b". The goal here is to define a SMARTS hierarchy that would assign the
same labels. Here is how to do it:

.. code-block:: python

>>> from besmarts.cluster.cluster_assignment import smiles_assignment_str
>>> from besmarts.core.assignments import smiles_assignment_group_bonds
>>> from besmarts.cluster.cluster_optimization import cluster_classifications
>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.assign.hierarchy_assign_rdkit import smarts_hierarchy_assignment_rdkit
>>> from besmarts.core import hierarchies
>>> 
>>> gcd = graph_codec_rdkit()
>>> labeler = smarts_hierarchy_assignment_rdkit()
>>> 
>>> smi = "[C:1]([H:3])#[C:2][H:4]"
>>> assns = {(1,2): "a", (1,3): "b", (2,4): "b"}
>>> sa  = smiles_assignment_str(smi, assns)
>>> 
>>> sag = smiles_assignment_group_bonds([sa])
>>> cst = cluster_classifications(gcd, labeler, sag)
>>> 
>>> hierarchies.smarts_hierarchy_print(cst.hierarchy)

.. code-block::
    
    Assigning molecule     1/1 at depth 0
    Labels per unique structure that need more depth
    There are 2/3 unique structures at depth 0
    There are 0 problems:
    Max depth is set to 0
    2024-05-23 17:26:26.587125 Labeling subgraphs
    2024-05-23 17:26:26.587688 Checking consistency...
    Optimization strategy is building steps...
    2024-05-23 17:26:26.587762 The optimization strategy has the following iterations:
    ->   1. op= 1 a=[0] b=1->1 d=0->0 n=0->3
         2. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         3. op= 1 a=[0] b=2->2 d=0->0 n=0->3
         4. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         5. op= 1 a=[0] b=3->3 d=0->0 n=0->3
         6. op=-1 a=[0] b=0->0 d=0->0 n=0->0
    Targets for this macro step 1:
    1 p0
    N Targets: 1
    Step tracker for current macro step 1
    
    
    *******************
     2024-05-23 17:26:26.587850 iteration=   1 macro=  1/6 X=        0 params=(1|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'a', 'b'} [*:1]~[*:2]
    =====
    
    2024-05-23 17:26:26.587901 Saving checkpoint to chk.cst.p
    2024-05-23 17:26:26.588329 Collecting SMARTS for p0 N=3/3 and setting to depth=0
     == iteration=   2 macro=  1/6 micro=  1/1 operation=1 params=(1|2) cluster=p0   N= 3 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=3
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    000003 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    2024-05-23 17:26:26.591432 Union merging=3
    2024-05-23 17:26:26.604654 Union merging=2
    2024-05-23 17:26:26.613429 Union merging=1
    2024-05-23 17:26:26.613705 Union is [#6H1X2x0!rA+0:1]!@;-,#[#1,#6;H0,H1;X1,X2;x0;!r;A;+0:3]
    2024-05-23 17:26:26.659793 Generating splits
    2024-05-23 17:26:26.660697 Generating single splits
    2024-05-23 17:26:26.661676 Generated 16 splits
    BIT [_______:1]_;-[_______:2]
    BIT [_______:1]_;!-[_______:2]
    BIT [_______:1]_;#[_______:2]
    BIT [_______:1]_;!#[_______:2]
    BIT [_______:1]_;_[#1______:2]
    BIT [_______:1]_;_[!#1______:2]
    BIT [_______:1]_;_[#6______:2]
    BIT [_______:1]_;_[!#6______:2]
    BIT [_______:1]_;_[_H0_____:2]
    BIT [_______:1]_;_[_!H0_____:2]
    BIT [_______:1]_;_[_H1_____:2]
    BIT [_______:1]_;_[_!H1_____:2]
    BIT [_______:1]_;_[__X1____:2]
    BIT [_______:1]_;_[__!X1____:2]
    BIT [_______:1]_;_[__X2____:2]
    BIT [_______:1]_;_[__!X2____:2]
    2024-05-23 17:26:26.670146 Building tasks
    Started local workspace on ('127.0.0.1', 39777)
    workspace listening on local host. Remote connections prohibited.
    2024-05-23 17:26:28.291644 Searching atoms=2 data=3 bit_depth=1/1 b_j=1/32 hits=0            
    2024-05-23 17:26:28.297140 Searching atoms=2 data=3 bit_depth=1/1 b_j=4/32 hits=2            
    2024-05-23 17:26:28.302437 Searching atoms=2 data=3 bit_depth=1/1 b_j=7/32 hits=2            
    2024-05-23 17:26:28.307670 Searching atoms=2 data=3 bit_depth=1/1 b_j=10/32 hits=2            
    2024-05-23 17:26:28.312911 Searching atoms=2 data=3 bit_depth=1/1 b_j=13/32 hits=2            
    2024-05-23 17:26:28.318257 Searching atoms=2 data=3 bit_depth=1/1 b_j=16/32 hits=2            
    2024-05-23 17:26:28.323462 Searching atoms=2 data=3 bit_depth=1/1 b_j=20/32 hits=2            
    2024-05-23 17:26:28.328734 Searching atoms=2 data=3 bit_depth=1/1 b_j=23/32 hits=2            
    2024-05-23 17:26:28.333925 Searching atoms=2 data=3 bit_depth=1/1 b_j=26/32 hits=2            
    2024-05-23 17:26:28.339113 Searching atoms=2 data=3 bit_depth=1/1 b_j=29/32 hits=2            
    2024-05-23 17:26:28.344373 Searching atoms=2 data=3 bit_depth=1/1 b_j=32/32 hits=2            
    Progress: 100.00%        32/32
    Finished: 100.00%        32/32
    Removing workspace ('127.0.0.1', 39777)
    Closing workspace
    2024-05-23 17:26:28.683668 Calculating partitions for hits=2
    Started local workspace on ('127.0.0.1', 41193)
    workspace listening on local host. Remote connections prohibited.
    Submitting 2 packets of work
    Chunk: 100.00%         2/2
    Finished: 100.00%         2/2
    Removing workspace ('127.0.0.1', 41193)
    Closing workspace
    2024-05-23 17:26:29.189932 Unique hits 1/2
    2024-05-23 17:26:29.189951 Searching atoms done; data=3 hits=1
    2024-05-23 17:26:29.190409 Collecting new candidates
    2024-05-23 17:26:29.190443 Scanning done.
    2024-05-23 17:26:29.190452
    
    
    Generating SMARTS on 1
    2024-05-23 17:26:29.214462 Labeling
    2024-05-23 17:26:29.215407 Rebuilding assignments
    2024-05-23 17:26:29.215503 Rebuilding mappings
    Tree:
    **  0   0 p0   {'a', 'b'} [*:1]~[*:2]
    =====
    
    Scanning 1 candidates for operation=1
    2024-05-23 17:26:29.312445 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 44067)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
    ->    1 Cnd.    1/1 p0      X=   -1.00000 dX=   -1.00000 N=      1 C= Y [*:1]!-[*:2]
    Performing 1 operations
    There are 1 nodes returned
    Operations per parameter for this micro:
    Counter({'p0': 1})
    Micro total: 1 should be 1
    Operations per parameter for this macro:
    Counter({'p0': 1})
    Macro total: 1 should be 1
    Pruned 0 empty nodes; candidates now 1/1
    []
    
    >>>>> New parameter    1/1 p1 parent p0 Objective   -1.00000 Delta   -1.00000 Partition 2|1
     >>>>> (0, None, -1) Local dObj   -1.00000 [*:1]!-[*:2]
    
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    There were 1 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-05-23 17:26:29.676434 Visited {'p1', 'p0'}
    Assignments changed for p1, will retarget
    Restarting optimization search
    Targets for this macro step 1:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 1
    p1 1
    
    
    *******************
     2024-05-23 17:26:29.677039 iteration=   2 macro=  1/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-05-23 17:26:29.677097 Saving checkpoint to chk.cst.p
    2024-05-23 17:26:29.677437 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   3 macro=  1/6 micro=  1/2 operation=1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:26:29.679564 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   4 macro=  1/6 micro=  2/2 operation=1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:26:29.680607 Scanning done.
    2024-05-23 17:26:29.680615
    
    
    Generating SMARTS on 0
    2024-05-23 17:26:29.702862 Labeling
    2024-05-23 17:26:29.703771 Rebuilding assignments
    2024-05-23 17:26:29.703873 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:26:29.749018 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 37805)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-05-23 17:26:30.076852 Visited set()
    Targets for this macro step 2:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 2
    p1 2
    p0 2
    
    
    *******************
     2024-05-23 17:26:30.077510 iteration=   4 macro=  2/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-05-23 17:26:30.077568 Saving checkpoint to chk.cst.p
    2024-05-23 17:26:30.077923 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   5 macro=  2/6 micro=  1/2 operation=-1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:26:30.077994 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   6 macro=  2/6 micro=  2/2 operation=-1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:26:30.078018 Scanning done.
    2024-05-23 17:26:30.078026
    
    
    Generating SMARTS on 1
    2024-05-23 17:26:30.102539 Labeling
    2024-05-23 17:26:30.103404 Rebuilding assignments
    2024-05-23 17:26:30.103495 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:26:30.152613 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 44097)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-05-23 17:26:30.556847 Visited {'p1'}
    Targets for this macro step 3:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 3
    p1 3
    p0 3
    
    
    *******************
     2024-05-23 17:26:30.557473 iteration=   6 macro=  3/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-05-23 17:26:30.557547 Saving checkpoint to chk.cst.p
    2024-05-23 17:26:30.557919 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   7 macro=  3/6 micro=  1/2 operation=1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=2->2 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:26:30.560069 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   8 macro=  3/6 micro=  2/2 operation=1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=2->2 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:26:30.561130 Scanning done.
    2024-05-23 17:26:30.561141
    
    
    Generating SMARTS on 0
    2024-05-23 17:26:30.584125 Labeling
    2024-05-23 17:26:30.585012 Rebuilding assignments
    2024-05-23 17:26:30.585243 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:26:30.629026 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 43481)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-05-23 17:26:30.939930 Visited set()
    Targets for this macro step 4:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 4
    p1 4
    p0 4
    
    
    *******************
     2024-05-23 17:26:30.940536 iteration=   8 macro=  4/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-05-23 17:26:30.940592 Saving checkpoint to chk.cst.p
    2024-05-23 17:26:30.940935 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   9 macro=  4/6 micro=  1/2 operation=-1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:26:30.941005 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  10 macro=  4/6 micro=  2/2 operation=-1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:26:30.941029 Scanning done.
    2024-05-23 17:26:30.941036
    
    
    Generating SMARTS on 1
    2024-05-23 17:26:30.964999 Labeling
    2024-05-23 17:26:30.965867 Rebuilding assignments
    2024-05-23 17:26:30.965958 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:26:31.012444 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 42867)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-05-23 17:26:31.403375 Visited {'p1'}
    Targets for this macro step 5:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 5
    p1 5
    p0 5
    
    
    *******************
     2024-05-23 17:26:31.403986 iteration=  10 macro=  5/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-05-23 17:26:31.404045 Saving checkpoint to chk.cst.p
    2024-05-23 17:26:31.404401 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  11 macro=  5/6 micro=  1/2 operation=1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:26:31.406498 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  12 macro=  5/6 micro=  2/2 operation=1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:26:31.407535 Scanning done.
    2024-05-23 17:26:31.407543
    
    
    Generating SMARTS on 0
    2024-05-23 17:26:31.430259 Labeling
    2024-05-23 17:26:31.431184 Rebuilding assignments
    2024-05-23 17:26:31.431285 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:26:31.475712 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 39675)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-05-23 17:26:31.800129 Visited set()
    Targets for this macro step 6:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 6
    p1 6
    p0 6
    
    
    *******************
     2024-05-23 17:26:31.800932 iteration=  12 macro=  6/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-05-23 17:26:31.800992 Saving checkpoint to chk.cst.p
    2024-05-23 17:26:31.801361 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  13 macro=  6/6 micro=  1/2 operation=-1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:26:31.801432 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  14 macro=  6/6 micro=  2/2 operation=-1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:26:31.801457 Scanning done.
    2024-05-23 17:26:31.801464
    
    
    Generating SMARTS on 1
    2024-05-23 17:26:31.829043 Labeling
    2024-05-23 17:26:31.829926 Rebuilding assignments
    2024-05-23 17:26:31.830018 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:26:31.875788 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 39273)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-05-23 17:26:33.283468 Visited {'p1'}
    Nothing found. Done.
    Start time: 2024-05-23 17:26:26.586739
    End   time: 2024-05-23 17:26:33.285234
    p0 {'b'}
    p1 {'a'}
    ACCURACY: 1.0
    **  0 p0 [*:1]~[*:2]
    **   1 p1 [*:1]!-[*:2]

There is quite a bit going on, but the last output shows the final hierarchy. The solution found was a SMARTS pattern `[*:1]!-[*:2]`.

Clustering numerical data
-------------------------

In this example, we have a SMILES and we have a bond length associated with
each bond. The goal here is to find a hierarchy where a parent and child SMARTS
patterns have a mean bond length difference of greater than the threshhold,
here 0.1 Angstrom. From the data given, we see that the algorithm should find a
hierarchy that separates bond 1-2 from bonds 1-3 and 2-4 since the difference
is 0.2 A and above the 0.1 threshold.

.. code-block:: python
>>> from besmarts.cluster.cluster_assignment import smiles_assignment_float
>>> from besmarts.core.assignments import smiles_assignment_group_bonds
>>> from besmarts.cluster.cluster_optimization import cluster_means
>>> from besmarts.cluster.cluster_objective import clustering_objective_mean_separation
>>> from besmarts.codecs.codec_rdkit import graph_codec_rdkit
>>> from besmarts.assign.hierarchy_assign_rdkit import smarts_hierarchy_assignment_rdkit
>>> from besmarts.core import hierarchies
>>> 
>>> gcd = graph_codec_rdkit()
>>> labeler = smarts_hierarchy_assignment_rdkit()
>>> 
>>> smi = "[C:1]([H:3])#[C:2][H:4]"
>>> assns = {(1,2): [1.1], (1,3): [1.3], (2,4): [1.3]}
>>> sa  = smiles_assignment_float(smi, assns)
>>> 
>>> objective = clustering_objective_mean_separation(split_separation=0.1)
>>> 
>>> sag = smiles_assignment_group_bonds([sa])
>>> cst = cluster_means(gcd, labeler, sag, objective=objective)
>>> 
>>> hierarchies.smarts_hierarchy_print(cst.hierarchy)

.. code-block::

    2024-05-23 17:34:22.988580 Labeling subgraphs
    2024-05-23 17:34:22.989177 Checking consistency...
    Optimization strategy is building steps...
    2024-05-23 17:34:22.989276 The optimization strategy has the following iterations:
    ->   1. op= 1 a=[0] b=1->1 d=0->0 n=0->3
         2. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         3. op= 1 a=[0] b=2->2 d=0->0 n=0->3
         4. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         5. op= 1 a=[0] b=3->3 d=0->0 n=0->3
         6. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         7. op= 1 a=[0] b=1->1 d=1->1 n=0->3
         8. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         9. op= 1 a=[0] b=2->2 d=1->1 n=0->3
        10. op=-1 a=[0] b=0->0 d=0->0 n=0->0
        11. op= 1 a=[0] b=3->3 d=1->1 n=0->3
        12. op=-1 a=[0] b=0->0 d=0->0 n=0->0
    Targets for this macro step 1:
    1 p0
    N Targets: 1
    Step tracker for current macro step 1
    
    
    *******************
     2024-05-23 17:34:22.989398 iteration=   1 macro=  1/12 X=        0 params=(1|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.2333 Var=    0.0089 N=      3 Min=    1.1000 Max=    1.3000 [*:1]~[*:2]
    =====
    
    2024-05-23 17:34:22.989459 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:22.989971 Collecting SMARTS for p0 N=3/3 and setting to depth=0
     == iteration=   2 macro=  1/12 micro=  1/1 operation=1 params=(1|1) cluster=p0   N= 3 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=3
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    000003 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    2024-05-23 17:34:22.993107 Union merging=3
    2024-05-23 17:34:23.006315 Union merging=2
    2024-05-23 17:34:23.015042 Union merging=1
    2024-05-23 17:34:23.015262 Union is [#6H1X2x0!rA+0:1]!@;-,#[#1,#6;H0,H1;X1,X2;x0;!r;A;+0:3]
    2024-05-23 17:34:23.063181 Generating splits
    2024-05-23 17:34:23.064124 Generating single splits
    2024-05-23 17:34:23.065146 Generated 16 splits
    BIT [_______:1]_;-[_______:2]
    BIT [_______:1]_;!-[_______:2]
    BIT [_______:1]_;#[_______:2]
    BIT [_______:1]_;!#[_______:2]
    BIT [_______:1]_;_[#1______:2]
    BIT [_______:1]_;_[!#1______:2]
    BIT [_______:1]_;_[#6______:2]
    BIT [_______:1]_;_[!#6______:2]
    BIT [_______:1]_;_[_H0_____:2]
    BIT [_______:1]_;_[_!H0_____:2]
    BIT [_______:1]_;_[_H1_____:2]
    BIT [_______:1]_;_[_!H1_____:2]
    BIT [_______:1]_;_[__X1____:2]
    BIT [_______:1]_;_[__!X1____:2]
    BIT [_______:1]_;_[__X2____:2]
    BIT [_______:1]_;_[__!X2____:2]
    2024-05-23 17:34:23.073941 Building tasks
    Started local workspace on ('127.0.0.1', 46809)
    workspace listening on local host. Remote connections prohibited.
    2024-05-23 17:34:24.556833 Searching atoms=2 data=3 bit_depth=1/1 b_j=1/32 hits=0            
    2024-05-23 17:34:24.562316 Searching atoms=2 data=3 bit_depth=1/1 b_j=4/32 hits=2            
    2024-05-23 17:34:24.567590 Searching atoms=2 data=3 bit_depth=1/1 b_j=7/32 hits=2            
    2024-05-23 17:34:24.572827 Searching atoms=2 data=3 bit_depth=1/1 b_j=10/32 hits=2            
    2024-05-23 17:34:24.578048 Searching atoms=2 data=3 bit_depth=1/1 b_j=13/32 hits=2            
    2024-05-23 17:34:24.583373 Searching atoms=2 data=3 bit_depth=1/1 b_j=16/32 hits=2            
    2024-05-23 17:34:24.588562 Searching atoms=2 data=3 bit_depth=1/1 b_j=20/32 hits=2            
    2024-05-23 17:34:24.593869 Searching atoms=2 data=3 bit_depth=1/1 b_j=23/32 hits=2            
    2024-05-23 17:34:24.599082 Searching atoms=2 data=3 bit_depth=1/1 b_j=26/32 hits=2            
    2024-05-23 17:34:24.604290 Searching atoms=2 data=3 bit_depth=1/1 b_j=29/32 hits=2            
    2024-05-23 17:34:24.609566 Searching atoms=2 data=3 bit_depth=1/1 b_j=32/32 hits=2            
    Progress: 100.00%        32/32
    Finished: 100.00%        32/32
    Removing workspace ('127.0.0.1', 46809)
    Closing workspace
    2024-05-23 17:34:24.977858 Calculating partitions for hits=2
    Started local workspace on ('127.0.0.1', 46745)
    workspace listening on local host. Remote connections prohibited.
    Submitting 2 packets of work
    Chunk: 100.00%         2/2
    Finished: 100.00%         2/2
    Removing workspace ('127.0.0.1', 46745)
    Closing workspace
    2024-05-23 17:34:25.500784 Unique hits 1/2
    2024-05-23 17:34:25.500811 Searching atoms done; data=3 hits=1
    2024-05-23 17:34:25.501254 Collecting new candidates
    2024-05-23 17:34:25.501293 Scanning done.
    2024-05-23 17:34:25.501302
    
    
    Generating SMARTS on 1
    2024-05-23 17:34:25.526312 Labeling
    2024-05-23 17:34:25.527250 Rebuilding assignments
    2024-05-23 17:34:25.527356 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.2333 Var=    0.0089 N=      3 Min=    1.1000 Max=    1.3000 [*:1]~[*:2]
    =====
    
    Scanning 1 candidates for operation=1
    2024-05-23 17:34:25.635865 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 33345)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
    ->    1 Cnd.    1/1 p0      X=   -0.20000 dX=   -0.20000 N=      1 C= Y [*:1]!-[*:2]
    Performing 1 operations
    There are 1 nodes returned
    Operations per parameter for this micro:
    Counter({'p0': 1})
    Micro total: 1 should be 1
    Operations per parameter for this macro:
    Counter({'p0': 1})
    Macro total: 1 should be 1
    Pruned 0 empty nodes; candidates now 1/1
    []
    
    >>>>> New parameter    1/1 p1 parent p0 Objective   -0.20000 Delta   -0.20000 Partition 2|1
     >>>>> (0, None, -1) Local dObj   -0.20000 [*:1]!-[*:2]
    
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    There were 1 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:27.086791 Visited {'p1', 'p0'}
    Assignments changed for p1, will retarget
    Restarting optimization search
    Targets for this macro step 1:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 1
    p1 1
    
    
    *******************
     2024-05-23 17:34:27.087571 iteration=   2 macro=  1/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:27.087654 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:27.088079 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   3 macro=  1/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:34:27.090249 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   4 macro=  1/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:34:27.091321 Scanning done.
    2024-05-23 17:34:27.091330
    
    
    Generating SMARTS on 0
    2024-05-23 17:34:27.115462 Labeling
    2024-05-23 17:34:27.116401 Rebuilding assignments
    2024-05-23 17:34:27.116507 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:34:27.162630 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 34887)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:27.480271 Visited set()
    Targets for this macro step 2:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 2
    p1 2
    p0 2
    
    
    *******************
     2024-05-23 17:34:27.480978 iteration=   4 macro=  2/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:27.481049 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:27.481450 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   5 macro=  2/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:27.481523 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   6 macro=  2/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:27.481549 Scanning done.
    2024-05-23 17:34:27.481556
    
    
    Generating SMARTS on 1
    2024-05-23 17:34:27.505213 Labeling
    2024-05-23 17:34:27.506082 Rebuilding assignments
    2024-05-23 17:34:27.506177 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:34:27.557830 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 35115)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
          1 Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
    Performing 0 operations
    There are 0 nodes returned
    Operations per parameter for this micro:
    Counter()
    Micro total: 0 should be 0
    Operations per parameter for this macro:
    Counter()
    Macro total: 0 should be 0
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:28.946573 Visited {'p1'}
    Targets for this macro step 3:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 3
    p1 3
    p0 3
    
    
    *******************
     2024-05-23 17:34:28.947405 iteration=   6 macro=  3/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:28.947489 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:28.947924 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   7 macro=  3/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=2->2 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:34:28.950024 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   8 macro=  3/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=2->2 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:34:28.951032 Scanning done.
    2024-05-23 17:34:28.951040
    
    
    Generating SMARTS on 0
    2024-05-23 17:34:28.973661 Labeling
    2024-05-23 17:34:28.974601 Rebuilding assignments
    2024-05-23 17:34:28.974703 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:34:29.022555 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 41995)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:29.360038 Visited set()
    Targets for this macro step 4:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 4
    p1 4
    p0 4
    
    
    *******************
     2024-05-23 17:34:29.360897 iteration=   8 macro=  4/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:29.360995 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:29.361478 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   9 macro=  4/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:29.361577 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  10 macro=  4/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:29.361608 Scanning done.
    2024-05-23 17:34:29.361615
    
    
    Generating SMARTS on 1
    2024-05-23 17:34:29.387180 Labeling
    2024-05-23 17:34:29.388064 Rebuilding assignments
    2024-05-23 17:34:29.388158 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:34:29.435904 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 42075)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
          1 Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
    Performing 0 operations
    There are 0 nodes returned
    Operations per parameter for this micro:
    Counter()
    Micro total: 0 should be 0
    Operations per parameter for this macro:
    Counter()
    Macro total: 0 should be 0
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:30.833520 Visited {'p1'}
    Targets for this macro step 5:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 5
    p1 5
    p0 5
    
    
    *******************
     2024-05-23 17:34:30.834224 iteration=  10 macro=  5/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:30.834296 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:30.834712 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  11 macro=  5/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:34:30.836866 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  12 macro=  5/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:34:30.837929 Scanning done.
    2024-05-23 17:34:30.837937
    
    
    Generating SMARTS on 0
    2024-05-23 17:34:30.861361 Labeling
    2024-05-23 17:34:30.862329 Rebuilding assignments
    2024-05-23 17:34:30.862433 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:34:30.912483 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 39275)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:31.236828 Visited set()
    Targets for this macro step 6:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 6
    p1 6
    p0 6
    
    
    *******************
     2024-05-23 17:34:31.237655 iteration=  12 macro=  6/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:31.237731 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:31.238182 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  13 macro=  6/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:31.238261 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  14 macro=  6/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:31.238292 Scanning done.
    2024-05-23 17:34:31.238300
    
    
    Generating SMARTS on 1
    2024-05-23 17:34:31.265314 Labeling
    2024-05-23 17:34:31.266191 Rebuilding assignments
    2024-05-23 17:34:31.266287 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:34:31.312560 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 39921)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
          1 Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
    Performing 0 operations
    There are 0 nodes returned
    Operations per parameter for this micro:
    Counter()
    Micro total: 0 should be 0
    Operations per parameter for this macro:
    Counter()
    Macro total: 0 should be 0
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:31.687313 Visited {'p1'}
    Targets for this macro step 7:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 7
    p1 7
    p0 7
    
    
    *******************
     2024-05-23 17:34:31.688129 iteration=  14 macro=  7/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:31.688212 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:31.688715 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  15 macro=  7/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=1->1 depth=1->1 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:34:31.692777 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  16 macro=  7/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=1->1 depth=1->1 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1](!@;-[#1H0X1x0!rA+0])!@;#[#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:34:31.695950 Scanning done.
    2024-05-23 17:34:31.695961
    
    
    Generating SMARTS on 0
    2024-05-23 17:34:31.719879 Labeling
    2024-05-23 17:34:31.720772 Rebuilding assignments
    2024-05-23 17:34:31.720871 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:34:31.769163 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 42435)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:32.087089 Visited set()
    Targets for this macro step 8:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 8
    p1 8
    p0 8
    
    
    *******************
     2024-05-23 17:34:32.087841 iteration=  16 macro=  8/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:32.087925 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:32.088387 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  17 macro=  8/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:32.088478 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  18 macro=  8/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:32.088504 Scanning done.
    2024-05-23 17:34:32.088511
    
    
    Generating SMARTS on 1
    2024-05-23 17:34:32.113207 Labeling
    2024-05-23 17:34:32.114077 Rebuilding assignments
    2024-05-23 17:34:32.114171 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:34:32.162481 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 43797)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
          1 Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
    Performing 0 operations
    There are 0 nodes returned
    Operations per parameter for this micro:
    Counter()
    Micro total: 0 should be 0
    Operations per parameter for this macro:
    Counter()
    Macro total: 0 should be 0
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:33.606664 Visited {'p1'}
    Targets for this macro step 9:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 9
    p1 9
    p0 9
    
    
    *******************
     2024-05-23 17:34:33.607350 iteration=  18 macro=  9/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:33.607418 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:33.607869 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  19 macro=  9/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=2->2 depth=1->1 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:34:33.611778 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  20 macro=  9/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=2->2 depth=1->1 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1](!@;-[#1H0X1x0!rA+0])!@;#[#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:34:33.614733 Scanning done.
    2024-05-23 17:34:33.614741
    
    
    Generating SMARTS on 0
    2024-05-23 17:34:33.638134 Labeling
    2024-05-23 17:34:33.639069 Rebuilding assignments
    2024-05-23 17:34:33.639172 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:34:33.689168 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 39949)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:34.010123 Visited set()
    Targets for this macro step 10:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 10
    p1 10
    p0 10
    
    
    *******************
     2024-05-23 17:34:34.010900 iteration=  20 macro= 10/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:34.010976 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:34.011401 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  21 macro= 10/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:34.011479 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  22 macro= 10/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:34.011507 Scanning done.
    2024-05-23 17:34:34.011515
    
    
    Generating SMARTS on 1
    2024-05-23 17:34:34.044524 Labeling
    2024-05-23 17:34:34.045599 Rebuilding assignments
    2024-05-23 17:34:34.045699 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:34:34.102507 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 40569)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
          1 Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
    Performing 0 operations
    There are 0 nodes returned
    Operations per parameter for this micro:
    Counter()
    Micro total: 0 should be 0
    Operations per parameter for this macro:
    Counter()
    Macro total: 0 should be 0
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:34.656675 Visited {'p1'}
    Targets for this macro step 11:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 11
    p1 11
    p0 11
    
    
    *******************
     2024-05-23 17:34:34.657338 iteration=  22 macro= 11/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:34.657408 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:34.657824 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  23 macro= 11/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=3->3 depth=1->1 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-05-23 17:34:34.661917 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  24 macro= 11/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=3->3 depth=1->1 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1](!@;-[#1H0X1x0!rA+0])!@;#[#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0]
    
    Skipping p1 since all graphs are the same
    2024-05-23 17:34:34.665025 Scanning done.
    2024-05-23 17:34:34.665035
    
    
    Generating SMARTS on 0
    2024-05-23 17:34:34.688677 Labeling
    2024-05-23 17:34:34.689611 Rebuilding assignments
    2024-05-23 17:34:34.689714 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-05-23 17:34:34.735826 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 36601)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:35.058439 Visited set()
    Targets for this macro step 12:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 12
    p1 12
    p0 12
    
    
    *******************
     2024-05-23 17:34:35.059726 iteration=  24 macro= 12/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-05-23 17:34:35.060044 Saving checkpoint to chk.cst.p
    2024-05-23 17:34:35.060669 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  25 macro= 12/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:35.060864 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  26 macro= 12/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-05-23 17:34:35.060972 Scanning done.
    2024-05-23 17:34:35.061009
    
    
    Generating SMARTS on 1
    2024-05-23 17:34:35.091613 Labeling
    2024-05-23 17:34:35.092507 Rebuilding assignments
    2024-05-23 17:34:35.092606 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-05-23 17:34:35.139353 workqueue started on ('0.0.0.0', 55555)
    Started local workspace on ('127.0.0.1', 35229)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
          1 Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
    Performing 0 operations
    There are 0 nodes returned
    Operations per parameter for this micro:
    Counter()
    Micro total: 0 should be 0
    Operations per parameter for this macro:
    Counter()
    Macro total: 0 should be 0
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-05-23 17:34:35.542994 Visited {'p1'}
    Nothing found. Done.
    Start time: 2024-05-23 17:34:22.988097
    End   time: 2024-05-23 17:34:35.544876
    **  0 p0 [*:1]~[*:2]
    **   1 p1 [*:1]!-[*:2]

Similar to the categorical case, the same SMARTS pattern was found.
