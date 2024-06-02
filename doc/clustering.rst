
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
>>> from besmarts.core import configs
>>> 
>>> configs.workqueue_port = 54321 # make sure this port is open/unused
>>> configs.remote_compute_enable = False # port is only open to localhost
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
    2024-06-02 10:59:16.684934 Labeling subgraphs
    2024-06-02 10:59:16.685488 Checking consistency...
    Optimization strategy is building steps...
    2024-06-02 10:59:16.685556 The optimization strategy has the following iterations:
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
     2024-06-02 10:59:16.685609 iteration=   1 macro=  1/6 X=        0 params=(1|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'a', 'b'} [*:1]~[*:2]
    =====
    
    2024-06-02 10:59:16.685637 Saving checkpoint to chk.cst.p
    2024-06-02 10:59:16.686071 Collecting SMARTS for p0 N=3/3 and setting to depth=0
     == iteration=   2 macro=  1/6 micro=  1/1 operation=1 params=(1|2) cluster=p0   N= 3 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=3
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    000003 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    2024-06-02 10:59:16.689140 Union merging=3
    2024-06-02 10:59:16.703539 Union merging=2
    2024-06-02 10:59:16.739951 Union merging=1
    2024-06-02 10:59:16.740180 Union is [#6H1X2x0!rA+0:1]!@;-,#[#1,#6;H0,H1;X1,X2;x0;!r;A;+0:3]
    2024-06-02 10:59:16.789138 Generating splits
    2024-06-02 10:59:16.790030 Generating single splits
    2024-06-02 10:59:16.791063 Generated 16 splits
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
    2024-06-02 10:59:16.799974 Building tasks
    workspace listening on local host. Remote connections prohibited.
    Progress: 100.00%        32/32
    Finished: 100.00%        32/32
    Removing workspace ('127.0.0.1', 33721)
    Closing workspace
    2024-06-02 10:59:21.176492 Calculating partitions for hits=2
    workspace listening on local host. Remote connections prohibited.
    Submitting 2 packets of work
    Chunk: 100.00%         2/2
    Finished: 100.00%         2/2
    Removing workspace ('127.0.0.1', 36193)
    Closing workspace
    2024-06-02 10:59:24.460336 Unique hits 1/2
    2024-06-02 10:59:24.460352 Searching atoms done; data=3 hits=1
    2024-06-02 10:59:24.460824 Collecting new candidates
    2024-06-02 10:59:24.460858 Scanning done.
    2024-06-02 10:59:24.460865
    
    
    Generating SMARTS on 1
    2024-06-02 10:59:24.500874 Labeling
    2024-06-02 10:59:24.501874 Rebuilding assignments
    2024-06-02 10:59:24.501982 Rebuilding mappings
    Tree:
    **  0   0 p0   {'a', 'b'} [*:1]~[*:2]
    =====
    
    Scanning 1 candidates for operation=1
    2024-06-02 10:59:24.595275 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=   -1.00000 dX=   -1.00000 N=      1 C= Y [*:1]!-[*:2]
                                                                                
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
    2024-06-02 10:59:27.919517 Visited {'p1', 'p0'}
    Assignments changed for p1, will retarget
    Restarting optimization search
    Targets for this macro step 1:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 1
    p1 1
    
    
    *******************
     2024-06-02 10:59:27.920159 iteration=   2 macro=  1/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-06-02 10:59:27.920196 Saving checkpoint to chk.cst.p
    2024-06-02 10:59:27.920560 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   3 macro=  1/6 micro=  1/2 operation=1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 10:59:27.922583 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   4 macro=  1/6 micro=  2/2 operation=1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-06-02 10:59:27.923556 Scanning done.
    2024-06-02 10:59:27.923561
    
    
    Generating SMARTS on 0
    2024-06-02 10:59:27.947135 Labeling
    2024-06-02 10:59:27.948160 Rebuilding assignments
    2024-06-02 10:59:27.948256 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 10:59:27.995228 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-06-02 10:59:28.310074 Visited set()
    Targets for this macro step 2:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 2
    p1 2
    p0 2
    
    
    *******************
     2024-06-02 10:59:28.310808 iteration=   4 macro=  2/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-06-02 10:59:28.310847 Saving checkpoint to chk.cst.p
    2024-06-02 10:59:28.311265 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   5 macro=  2/6 micro=  1/2 operation=-1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 10:59:28.311331 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   6 macro=  2/6 micro=  2/2 operation=-1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 10:59:28.311351 Scanning done.
    2024-06-02 10:59:28.311355
    
    
    Generating SMARTS on 1
    2024-06-02 10:59:28.338544 Labeling
    2024-06-02 10:59:28.339483 Rebuilding assignments
    2024-06-02 10:59:28.339573 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 10:59:28.388501 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    1.00000 N=      3 C= N [*:1]!-[*:2]
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-06-02 10:59:31.709658 Visited {'p1'}
    Targets for this macro step 3:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 3
    p1 3
    p0 3
    
    
    *******************
     2024-06-02 10:59:31.710212 iteration=   6 macro=  3/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-06-02 10:59:31.710245 Saving checkpoint to chk.cst.p
    2024-06-02 10:59:31.710597 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   7 macro=  3/6 micro=  1/2 operation=1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=2->2 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 10:59:31.712680 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   8 macro=  3/6 micro=  2/2 operation=1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=2->2 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-06-02 10:59:31.713693 Scanning done.
    2024-06-02 10:59:31.713698
    
    
    Generating SMARTS on 0
    2024-06-02 10:59:31.737206 Labeling
    2024-06-02 10:59:31.738138 Rebuilding assignments
    2024-06-02 10:59:31.738233 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 10:59:31.785195 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-06-02 10:59:32.116772 Visited set()
    Targets for this macro step 4:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 4
    p1 4
    p0 4
    
    
    *******************
     2024-06-02 10:59:32.117548 iteration=   8 macro=  4/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-06-02 10:59:32.117592 Saving checkpoint to chk.cst.p
    2024-06-02 10:59:32.118019 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   9 macro=  4/6 micro=  1/2 operation=-1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 10:59:32.118092 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  10 macro=  4/6 micro=  2/2 operation=-1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 10:59:32.118113 Scanning done.
    2024-06-02 10:59:32.118117
    
    
    Generating SMARTS on 1
    2024-06-02 10:59:32.157983 Labeling
    2024-06-02 10:59:32.158872 Rebuilding assignments
    2024-06-02 10:59:32.158971 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 10:59:32.205257 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    1.00000 N=      3 C= N [*:1]!-[*:2]
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-06-02 10:59:35.532537 Visited {'p1'}
    Targets for this macro step 5:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 5
    p1 5
    p0 5
    
    
    *******************
     2024-06-02 10:59:35.533180 iteration=  10 macro=  5/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-06-02 10:59:35.533214 Saving checkpoint to chk.cst.p
    2024-06-02 10:59:35.533577 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  11 macro=  5/6 micro=  1/2 operation=1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 10:59:35.535239 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  12 macro=  5/6 micro=  2/2 operation=1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-06-02 10:59:35.536045 Scanning done.
    2024-06-02 10:59:35.536053
    
    
    Generating SMARTS on 0
    2024-06-02 10:59:35.564996 Labeling
    2024-06-02 10:59:35.565884 Rebuilding assignments
    2024-06-02 10:59:35.565976 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 10:59:35.611708 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-06-02 10:59:35.932799 Visited set()
    Targets for this macro step 6:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 6
    p1 6
    p0 6
    
    
    *******************
     2024-06-02 10:59:35.933382 iteration=  12 macro=  6/6 X=       -1 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    2024-06-02 10:59:35.933415 Saving checkpoint to chk.cst.p
    2024-06-02 10:59:35.933743 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  13 macro=  6/6 micro=  1/2 operation=-1 params=(2|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 10:59:35.933804 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  14 macro=  6/6 micro=  2/2 operation=-1 params=(2|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 10:59:35.933819 Scanning done.
    2024-06-02 10:59:35.933824
    
    
    Generating SMARTS on 1
    2024-06-02 10:59:35.960568 Labeling
    2024-06-02 10:59:35.961440 Rebuilding assignments
    2024-06-02 10:59:35.961532 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 10:59:36.008529 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    1.00000 N=      3 C= N [*:1]!-[*:2]
                                                                                
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]!-[*:2]
    2024-06-02 10:59:39.329489 Visited {'p1'}
    Nothing found. Done.
    Start time: 2024-06-02 10:59:16.684534
    End   time: 2024-06-02 10:59:39.331211
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
>>> from besmarts.core import configs
>>> 
>>> configs.workqueue_port = 54321 # make sure this port is open/unused
>>> configs.remote_compute_enable = False # port is only open to localhost
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

    2024-06-02 11:04:47.408432 Labeling subgraphs
    2024-06-02 11:04:47.409015 Checking consistency...
    Optimization strategy is building steps...
    2024-06-02 11:04:47.409105 The optimization strategy has the following iterations:
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
     2024-06-02 11:04:47.409183 iteration=   1 macro=  1/12 X=        0 params=(1|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.2333 Var=    0.0089 N=      3 Min=    1.1000 Max=    1.3000 [*:1]~[*:2]
    =====
    
    2024-06-02 11:04:47.409223 Saving checkpoint to chk.cst.p
    2024-06-02 11:04:47.409753 Collecting SMARTS for p0 N=3/3 and setting to depth=0
     == iteration=   2 macro=  1/12 micro=  1/1 operation=1 params=(1|1) cluster=p0   N= 3 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=3
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    000003 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    2024-06-02 11:04:47.412902 Union merging=3
    2024-06-02 11:04:47.427718 Union merging=2
    2024-06-02 11:04:47.438833 Union merging=1
    2024-06-02 11:04:47.439097 Union is [#6H1X2x0!rA+0:1]!@;-,#[#1,#6;H0,H1;X1,X2;x0;!r;A;+0:3]
    2024-06-02 11:04:47.485960 Generating splits
    2024-06-02 11:04:47.486797 Generating single splits
    2024-06-02 11:04:47.487659 Generated 16 splits
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
    2024-06-02 11:04:47.496216 Building tasks
    workspace listening on local host. Remote connections prohibited.
    Progress: 100.00%        32/32
    Finished: 100.00%        32/32
    Removing workspace ('127.0.0.1', 43469)
    Closing workspace
    2024-06-02 11:04:51.912694 Calculating partitions for hits=2
    workspace listening on local host. Remote connections prohibited.
    Submitting 2 packets of work
    Chunk: 100.00%         2/2
    Finished: 100.00%         2/2
    Removing workspace ('127.0.0.1', 38483)
    Closing workspace
    2024-06-02 11:04:56.186463 Unique hits 1/2
    2024-06-02 11:04:56.186477 Searching atoms done; data=3 hits=1
    2024-06-02 11:04:56.186912 Collecting new candidates
    2024-06-02 11:04:56.186945 Scanning done.
    2024-06-02 11:04:56.186952
    
    
    Generating SMARTS on 1
    2024-06-02 11:04:56.211579 Labeling
    2024-06-02 11:04:56.212446 Rebuilding assignments
    2024-06-02 11:04:56.212521 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.2333 Var=    0.0089 N=      3 Min=    1.1000 Max=    1.3000 [*:1]~[*:2]
    =====
    
    Scanning 1 candidates for operation=1
    2024-06-02 11:04:56.308235 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=   -0.20000 dX=   -0.20000 N=      1 C= Y [*:1]!-[*:2]
                                                                                
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
    2024-06-02 11:04:59.622629 Visited {'p1', 'p0'}
    Assignments changed for p1, will retarget
    Restarting optimization search
    Targets for this macro step 1:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 1
    p1 1
    
    
    *******************
     2024-06-02 11:04:59.623325 iteration=   2 macro=  1/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:04:59.623368 Saving checkpoint to chk.cst.p
    2024-06-02 11:04:59.623765 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   3 macro=  1/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 11:04:59.625794 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   4 macro=  1/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=1->1 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-06-02 11:04:59.626766 Scanning done.
    2024-06-02 11:04:59.626771
    
    
    Generating SMARTS on 0
    2024-06-02 11:04:59.656996 Labeling
    2024-06-02 11:04:59.657927 Rebuilding assignments
    2024-06-02 11:04:59.658041 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 11:04:59.708607 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-06-02 11:05:00.036512 Visited set()
    Targets for this macro step 2:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 2
    p1 2
    p0 2
    
    
    *******************
     2024-06-02 11:05:00.037438 iteration=   4 macro=  2/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:00.037493 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:00.037989 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   5 macro=  2/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:00.038059 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   6 macro=  2/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:00.038081 Scanning done.
    2024-06-02 11:05:00.038085
    
    
    Generating SMARTS on 1
    2024-06-02 11:05:00.078292 Labeling
    2024-06-02 11:05:00.079384 Rebuilding assignments
    2024-06-02 11:05:00.079491 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 11:05:00.128601 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
                                                                                
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
    2024-06-02 11:05:03.456472 Visited {'p1'}
    Targets for this macro step 3:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 3
    p1 3
    p0 3
    
    
    *******************
     2024-06-02 11:05:03.457515 iteration=   6 macro=  3/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:03.457610 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:03.458065 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   7 macro=  3/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=2->2 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 11:05:03.460152 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   8 macro=  3/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=2->2 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-06-02 11:05:03.461218 Scanning done.
    2024-06-02 11:05:03.461226
    
    
    Generating SMARTS on 0
    2024-06-02 11:05:03.484470 Labeling
    2024-06-02 11:05:03.485383 Rebuilding assignments
    2024-06-02 11:05:03.485477 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 11:05:03.531656 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-06-02 11:05:03.853047 Visited set()
    Targets for this macro step 4:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 4
    p1 4
    p0 4
    
    
    *******************
     2024-06-02 11:05:03.853761 iteration=   8 macro=  4/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:03.853807 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:03.854220 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   9 macro=  4/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:03.854293 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  10 macro=  4/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:03.854310 Scanning done.
    2024-06-02 11:05:03.854314
    
    
    Generating SMARTS on 1
    2024-06-02 11:05:03.886263 Labeling
    2024-06-02 11:05:03.887092 Rebuilding assignments
    2024-06-02 11:05:03.887176 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 11:05:03.931957 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
                                                                                
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
    2024-06-02 11:05:07.255984 Visited {'p1'}
    Targets for this macro step 5:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 5
    p1 5
    p0 5
    
    
    *******************
     2024-06-02 11:05:07.256638 iteration=  10 macro=  5/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:07.256683 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:07.257088 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  11 macro=  5/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 11:05:07.259204 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  12 macro=  5/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-06-02 11:05:07.260226 Scanning done.
    2024-06-02 11:05:07.260231
    
    
    Generating SMARTS on 0
    2024-06-02 11:05:07.285031 Labeling
    2024-06-02 11:05:07.285952 Rebuilding assignments
    2024-06-02 11:05:07.286044 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 11:05:07.335100 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-06-02 11:05:07.649581 Visited set()
    Targets for this macro step 6:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 6
    p1 6
    p0 6
    
    
    *******************
     2024-06-02 11:05:07.650354 iteration=  12 macro=  6/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:07.650402 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:07.650854 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  13 macro=  6/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:07.650937 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  14 macro=  6/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:07.650956 Scanning done.
    2024-06-02 11:05:07.650961
    
    
    Generating SMARTS on 1
    2024-06-02 11:05:07.681759 Labeling
    2024-06-02 11:05:07.682697 Rebuilding assignments
    2024-06-02 11:05:07.682801 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 11:05:07.731889 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
                                                                                
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
    2024-06-02 11:05:11.049060 Visited {'p1'}
    Targets for this macro step 7:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 7
    p1 7
    p0 7
    
    
    *******************
     2024-06-02 11:05:11.049725 iteration=  14 macro=  7/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:11.049769 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:11.050172 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  15 macro=  7/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=1->1 depth=1->1 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 11:05:11.054154 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  16 macro=  7/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=1->1 depth=1->1 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1](!@;-[#1H0X1x0!rA+0])!@;#[#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0]
    
    Skipping p1 since all graphs are the same
    2024-06-02 11:05:11.057292 Scanning done.
    2024-06-02 11:05:11.057299
    
    
    Generating SMARTS on 0
    2024-06-02 11:05:11.080485 Labeling
    2024-06-02 11:05:11.081414 Rebuilding assignments
    2024-06-02 11:05:11.081508 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 11:05:11.128346 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-06-02 11:05:11.442149 Visited set()
    Targets for this macro step 8:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 8
    p1 8
    p0 8
    
    
    *******************
     2024-06-02 11:05:11.442822 iteration=  16 macro=  8/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:11.442864 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:11.443268 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  17 macro=  8/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:11.443329 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  18 macro=  8/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:11.443345 Scanning done.
    2024-06-02 11:05:11.443349
    
    
    Generating SMARTS on 1
    2024-06-02 11:05:11.467130 Labeling
    2024-06-02 11:05:11.467989 Rebuilding assignments
    2024-06-02 11:05:11.468070 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 11:05:11.514927 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
                                                                                
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
    2024-06-02 11:05:14.829189 Visited {'p1'}
    Targets for this macro step 9:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 9
    p1 9
    p0 9
    
    
    *******************
     2024-06-02 11:05:14.829845 iteration=  18 macro=  9/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:14.829889 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:14.830313 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  19 macro=  9/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=2->2 depth=1->1 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 11:05:14.834290 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  20 macro=  9/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=2->2 depth=1->1 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1](!@;-[#1H0X1x0!rA+0])!@;#[#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0]
    
    Skipping p1 since all graphs are the same
    2024-06-02 11:05:14.837327 Scanning done.
    2024-06-02 11:05:14.837332
    
    
    Generating SMARTS on 0
    2024-06-02 11:05:14.861249 Labeling
    2024-06-02 11:05:14.862122 Rebuilding assignments
    2024-06-02 11:05:14.862208 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 11:05:14.911881 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-06-02 11:05:15.226564 Visited set()
    Targets for this macro step 10:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 10
    p1 10
    p0 10
    
    
    *******************
     2024-06-02 11:05:15.227388 iteration=  20 macro= 10/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:15.227439 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:15.227927 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  21 macro= 10/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:15.227994 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  22 macro= 10/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:15.228011 Scanning done.
    2024-06-02 11:05:15.228015
    
    
    Generating SMARTS on 1
    2024-06-02 11:05:15.252244 Labeling
    2024-06-02 11:05:15.253074 Rebuilding assignments
    2024-06-02 11:05:15.253156 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 11:05:15.301668 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
                                                                                
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
    2024-06-02 11:05:18.619994 Visited {'p1'}
    Targets for this macro step 11:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 11
    p1 11
    p0 11
    
    
    *******************
     2024-06-02 11:05:18.620711 iteration=  22 macro= 11/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:18.620759 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:18.621190 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  23 macro= 11/12 micro=  1/2 operation=1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=3->3 depth=1->1 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Var=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2](!@;#[#6H1X2x0!rA+0])!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-06-02 11:05:18.624135 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  24 macro= 11/12 micro=  2/2 operation=1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=3->3 depth=1->1 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]!-[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1](!@;-[#1H0X1x0!rA+0])!@;#[#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0]
    
    Skipping p1 since all graphs are the same
    2024-06-02 11:05:18.626291 Scanning done.
    2024-06-02 11:05:18.626297
    
    
    Generating SMARTS on 0
    2024-06-02 11:05:18.657163 Labeling
    2024-06-02 11:05:18.658140 Rebuilding assignments
    2024-06-02 11:05:18.658227 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    2024-06-02 11:05:18.705096 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    2024-06-02 11:05:19.029217 Visited set()
    Targets for this macro step 12:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 12
    p1 12
    p0 12
    
    
    *******************
     2024-06-02 11:05:19.029938 iteration=  24 macro= 12/12 X=     -0.2 params=(2|1) G=Y S=Y bits=1->3 depth=0->1 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    2024-06-02 11:05:19.029981 Saving checkpoint to chk.cst.p
    2024-06-02 11:05:19.030400 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  25 macro= 12/12 micro=  1/2 operation=-1 params=(2|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:19.030462 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  26 macro= 12/12 micro=  2/2 operation=-1 params=(2|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-06-02 11:05:19.030478 Scanning done.
    2024-06-02 11:05:19.030482
    
    
    Generating SMARTS on 1
    2024-06-02 11:05:19.054996 Labeling
    2024-06-02 11:05:19.055868 Rebuilding assignments
    2024-06-02 11:05:19.055954 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Var=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Var=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]!-[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    2024-06-02 11:05:19.105084 workqueue started on ('127.0.0.1', 54321)
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    
    Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= Y [*:1]!-[*:2]
                                                                                
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
    2024-06-02 11:05:22.430115 Visited {'p1'}
    Nothing found. Done.
    Start time: 2024-06-02 11:04:47.407944
    End   time: 2024-06-02 11:05:22.432080
    **  0 p0 [*:1]~[*:2]
    **   1 p1 [*:1]!-[*:2]

Similar to the categorical case, the same SMARTS pattern was found.
