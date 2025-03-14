
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

    from besmarts.cluster.cluster_assignment import smiles_assignment_str
    from besmarts.core.assignments import smiles_assignment_group_bonds
    from besmarts.cluster.cluster_optimization import cluster_classifications
    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.assign.hierarchy_assign_rdkit import smarts_hierarchy_assignment_rdkit
    from besmarts.core import hierarchies
    from besmarts.core import configs
    
    configs.workqueue_port = 54321 # make sure this port is open/unused
    configs.remote_compute_enable = False # port is only open to localhost
    
    gcd = graph_codec_rdkit()
    labeler = smarts_hierarchy_assignment_rdkit()
    
    smi = "[C:1]([H:3])#[C:2][H:4]"
    assns = {(1,2): "a", (1,3): "b", (2,4): "b"}
    sa  = smiles_assignment_str(smi, assns)
    
    sag = smiles_assignment_group_bonds([sa])
    cst = cluster_classifications(gcd, labeler, sag)
    
    hierarchies.smarts_hierarchy_print(cst.hierarchy)


.. code-block::

    Assigning molecule     1/1 at depth 0
    Labels per unique structure that need more depth
    There are 2/3 unique structures at depth 0
    There are 0 problems:
    Max depth is set to 0
    2024-12-17 10:52:34.477112 Labeling subgraphs
    2024-12-17 10:52:34.477964 Checking consistency...
    Optimization strategy is building steps...
    2024-12-17 10:52:34.478116 The optimization strategy has the following iterations:
    ->   1. op= 1 a=[0] b=1->1 d=0->0 n=0->1
         2. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         3. op= 1 a=[0] b=2->2 d=0->0 n=0->2
         4. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         5. op= 1 a=[0] b=3->3 d=0->0 n=0->3
         6. op=-1 a=[0] b=0->0 d=0->0 n=0->0
    Targets for this macro step 1:
    1 p0
    N Targets: 1
    Step tracker for current macro step 1
    p0 1
    
    
    *******************
     2024-12-17 10:52:34.478233 iteration=   1 macro=  1/6 X=        0 params=(2|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b', 'a'} [*:1]~[*:2]
    =====
    
    2024-12-17 10:52:34.478282 Saving checkpoint to chk.cst.p
    2024-12-17 10:52:34.479079 Collecting SMARTS for p0 N=3/3 and setting to depth=0
     == iteration=   2 macro=  1/6 micro=  1/1 operation=1 params=(2|2) cluster=p0   N= 3 overlap=[0] bits=1->1 depth=0->0 branch=0->1
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=3
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    000003 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    2024-12-17 10:52:34.482381 Union merging=3
    2024-12-17 10:52:34.500785 Union merging=1
    2024-12-17 10:52:34.501043 Union is [#6H1X2x0!rA+0:1]!@;-,#[#1,#6;H0,H1;X1,X2;x0;!r;A;+0:3]
    2024-12-17 10:52:34.501556 Generating splits
    2024-12-17 10:52:34.502301 Generating single splits
    2024-12-17 10:52:34.504652 Generated 16 splits
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
    2024-12-17 10:52:34.513843 Building tasks
    workspace listening on local host. Remote connections prohibited.
    2024-12-17 10:52:34.514065 P:   0.00%    32/32 IQ:    1 OQ:    0 IP:    0 LF:    0 RF:    0 RIQ:    0 ROQ:    0 RIP:    0  ERC:    0.0 
    2024-12-17 10:52:34.614145 P: 100.00%     0/32 IQ:    0 OQ:    0 IP:    0 LF:    0 RF:    0 RIQ:    0 ROQ:    0 RIP:    0  ERC:    0.0 
        1 CND SPLITS=Y  [*:1]!-[*:2]
        2 CND SPLITS=Y  [*:1]-[*:2]
        5 CND SPLITS=Y  [*:1]!#[*:2]
        6 CND SPLITS=Y  [*:1]#[*:2]
        9 CND SPLITS=N  [*:1]~[!#1:2]
       10 CND SPLITS=Y  [*:1]~[#1:2]
       13 CND SPLITS=Y  [*:1]~[!#6:2]
       14 CND SPLITS=N  [*:1]~[#6:2]
       17 CND SPLITS=N  [*:1]~[!H0:2]
       18 CND SPLITS=Y  [*:1]~[H0:2]
       21 CND SPLITS=Y  [*:1]~[!H1:2]
       22 CND SPLITS=N  [*:1]~[H1:2]
       25 CND SPLITS=N  [*:1]~[!X1:2]
       26 CND SPLITS=Y  [*:1]~[X1:2]
       29 CND SPLITS=Y  [*:1]~[!X2:2]
       30 CND SPLITS=N  [*:1]~[X2:2]
    Finished: 100.00%        32/32
    Closing workspace
    2024-12-17 10:52:34.705268 Calculating partitions for hits=10
    workspace listening on local host. Remote connections prohibited.
    Submitting 10 packets of work
    Closing workspace
    2024-12-17 10:52:34.811964 Unique hits 1/10
        7 HIT S0= 2     -> Sj= 1     [*:1]#[*:2]
        1     DUP [*:1]~[!X2:2]
        2     DUP [*:1]~[X1:2]
        3     DUP [*:1]~[!H1:2]
        4     DUP [*:1]~[H0:2]
        5     DUP [*:1]~[!#6:2]
        6     DUP [*:1]~[#1:2]
        8     DUP [*:1]!#[*:2]
        9     DUP [*:1]-[*:2]
       10     DUP [*:1]!-[*:2]
    2024-12-17 10:52:34.812934 Searching atoms done; data=3 hits=1
    2024-12-17 10:52:34.813104 Collecting new candidates
    2024-12-17 10:52:34.813126 Scanning done.
    2024-12-17 10:52:34.813133
    
    
    Generating SMARTS on 1
    2024-12-17 10:52:34.822349 Labeling
    2024-12-17 10:52:34.823310 Rebuilding assignments
    2024-12-17 10:52:34.823394 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b', 'a'} [*:1]~[*:2]
    =====
    
    Scanning 1 candidates for operation=1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    Cnd.    1/1 p0      X=   -1.00000 dX=   -1.00000 N=      1 C= Y [*:1]#[*:2]
                                                                               
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
    ->    1 Cnd.    1/1 p0      X=   -1.00000 dX=   -1.00000 N=      1 C= Y [*:1]#[*:2]
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
     >>>>> (0, None, -1) Local dObj   -1.00000 [*:1]#[*:2]
    
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    There were 1 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    2024-12-17 10:52:34.916340 Visited {'p0', 'p1'}
    Assignments changed for p0, will retarget
    Assignments changed for p1, will retarget
    Restarting optimization search
    Targets for this macro step 1:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 1
    p0 1
    p1 1
    
    
    *******************
     2024-12-17 10:52:34.916927 iteration=   2 macro=  1/6 X=       -1 params=(3|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    2024-12-17 10:52:34.916987 Saving checkpoint to chk.cst.p
    2024-12-17 10:52:34.917330 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   3 macro=  1/6 micro=  1/2 operation=1 params=(3|2) cluster=p0   N= 2 overlap=[0] bits=1->1 depth=0->0 branch=0->1
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-12-17 10:52:34.919388 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   4 macro=  1/6 micro=  2/2 operation=1 params=(3|2) cluster=p1   N= 1 overlap=[0] bits=1->1 depth=0->0 branch=0->1
    
    Attempting to split p1:
    S0: [*:1]#[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-12-17 10:52:34.920404 Scanning done.
    2024-12-17 10:52:34.920411
    
    
    Generating SMARTS on 0
    2024-12-17 10:52:34.926656 Labeling
    2024-12-17 10:52:34.927640 Rebuilding assignments
    2024-12-17 10:52:34.927738 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    2024-12-17 10:52:35.015457 Visited set()
    Targets for this macro step 2:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 2
    p0 2
    p1 2
    
    
    *******************
     2024-12-17 10:52:35.016147 iteration=   4 macro=  2/6 X=       -1 params=(3|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    2024-12-17 10:52:35.016214 Saving checkpoint to chk.cst.p
    2024-12-17 10:52:35.016597 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   5 macro=  2/6 micro=  1/2 operation=-1 params=(3|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:52:35.016677 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   6 macro=  2/6 micro=  2/2 operation=-1 params=(3|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:52:35.016705 Scanning done.
    2024-12-17 10:52:35.016711
    
    
    Generating SMARTS on 1
    2024-12-17 10:52:35.025372 Labeling
    2024-12-17 10:52:35.026315 Rebuilding assignments
    2024-12-17 10:52:35.026420 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    Cnd.    1/1 p0      X=    0.00000 dX=    1.00000 N=      3 C= N [*:1]#[*:2]
                                                                               
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    2024-12-17 10:52:35.115490 Visited {'p1'}
    Targets for this macro step 3:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 3
    p0 3
    p1 3
    
    
    *******************
     2024-12-17 10:52:35.116151 iteration=   6 macro=  3/6 X=       -1 params=(3|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    2024-12-17 10:52:35.116218 Saving checkpoint to chk.cst.p
    2024-12-17 10:52:35.116606 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   7 macro=  3/6 micro=  1/2 operation=1 params=(3|2) cluster=p0   N= 2 overlap=[0] bits=2->2 depth=0->0 branch=0->2
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-12-17 10:52:35.118449 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   8 macro=  3/6 micro=  2/2 operation=1 params=(3|2) cluster=p1   N= 1 overlap=[0] bits=2->2 depth=0->0 branch=0->2
    
    Attempting to split p1:
    S0: [*:1]#[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-12-17 10:52:35.119490 Scanning done.
    2024-12-17 10:52:35.119499
    
    
    Generating SMARTS on 0
    2024-12-17 10:52:35.125876 Labeling
    2024-12-17 10:52:35.126798 Rebuilding assignments
    2024-12-17 10:52:35.126887 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    2024-12-17 10:52:35.215341 Visited set()
    Targets for this macro step 4:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 4
    p0 4
    p1 4
    
    
    *******************
     2024-12-17 10:52:35.215914 iteration=   8 macro=  4/6 X=       -1 params=(3|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    2024-12-17 10:52:35.215972 Saving checkpoint to chk.cst.p
    2024-12-17 10:52:35.216299 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   9 macro=  4/6 micro=  1/2 operation=-1 params=(3|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:52:35.216373 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  10 macro=  4/6 micro=  2/2 operation=-1 params=(3|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:52:35.216399 Scanning done.
    2024-12-17 10:52:35.216406
    
    
    Generating SMARTS on 1
    2024-12-17 10:52:35.229950 Labeling
    2024-12-17 10:52:35.230810 Rebuilding assignments
    2024-12-17 10:52:35.230886 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    Cnd.    1/1 p0      X=    0.00000 dX=    1.00000 N=      3 C= N [*:1]#[*:2]
                                                                               
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    2024-12-17 10:52:35.318658 Visited {'p1'}
    Targets for this macro step 5:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 5
    p0 5
    p1 5
    
    
    *******************
     2024-12-17 10:52:35.319183 iteration=  10 macro=  5/6 X=       -1 params=(3|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    2024-12-17 10:52:35.319237 Saving checkpoint to chk.cst.p
    2024-12-17 10:52:35.319547 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  11 macro=  5/6 micro=  1/2 operation=1 params=(3|2) cluster=p0   N= 2 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))              {'b'} [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))              {'b'} [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-12-17 10:52:35.321560 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  12 macro=  5/6 micro=  2/2 operation=1 params=(3|2) cluster=p1   N= 1 overlap=[0] bits=3->3 depth=0->0 branch=0->3
    
    Attempting to split p1:
    S0: [*:1]#[*:2]
    Matched N=1
    000001 (0, (1, 2))              {'a'} [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-12-17 10:52:35.322564 Scanning done.
    2024-12-17 10:52:35.322571
    
    
    Generating SMARTS on 0
    2024-12-17 10:52:35.328626 Labeling
    2024-12-17 10:52:35.329628 Rebuilding assignments
    2024-12-17 10:52:35.329719 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    2024-12-17 10:52:35.418673 Visited set()
    Targets for this macro step 6:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 6
    p0 6
    p1 6
    
    
    *******************
     2024-12-17 10:52:35.419223 iteration=  12 macro=  6/6 X=       -1 params=(3|2) G=Y S=Y bits=1->3 depth=0->0 branch=0->3
    *******************
    
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    2024-12-17 10:52:35.419276 Saving checkpoint to chk.cst.p
    2024-12-17 10:52:35.419593 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=  13 macro=  6/6 micro=  1/2 operation=-1 params=(3|2) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:52:35.419664 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  14 macro=  6/6 micro=  2/2 operation=-1 params=(3|2) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:52:35.419689 Scanning done.
    2024-12-17 10:52:35.419695
    
    
    Generating SMARTS on 1
    2024-12-17 10:52:35.427191 Labeling
    2024-12-17 10:52:35.428046 Rebuilding assignments
    2024-12-17 10:52:35.428127 Rebuilding mappings
    Tree:
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    Cnd.    1/1 p0      X=    0.00000 dX=    1.00000 N=      3 C= N [*:1]#[*:2]
                                                                               
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0   {'b'} [*:1]~[*:2]
    **  1   1 p1   {'a'} [*:1]#[*:2]
    2024-12-17 10:52:35.515337 Visited {'p1'}
    Nothing found. Done.
    Start time: 2024-12-17 10:52:34.427978
    End   time: 2024-12-17 10:52:35.516856
    p0 {'b'}
    p1 {'a'}
    ACCURACY: 1.0
    **  0 p0 [*:1]~[*:2]
    **   1 p1 [*:1]#[*:2]

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
    from besmarts.cluster.cluster_assignment import smiles_assignment_float
    from besmarts.core.assignments import smiles_assignment_group_bonds
    from besmarts.cluster.cluster_optimization import cluster_means
    from besmarts.cluster.cluster_objective import clustering_objective_mean_separation
    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.assign.hierarchy_assign_rdkit import smarts_hierarchy_assignment_rdkit
    from besmarts.core import hierarchies
    from besmarts.core import configs
    
    configs.workqueue_port = 54321 # make sure this port is open/unused
    configs.remote_compute_enable = False # port is only open to localhost
    
    gcd = graph_codec_rdkit()
    labeler = smarts_hierarchy_assignment_rdkit()
    
    smi = "[C:1]([H:3])#[C:2][H:4]"
    assns = {(1,2): [1.1], (1,3): [1.3], (2,4): [1.3]}
    sa  = smiles_assignment_float(smi, assns)
    
    objective = clustering_objective_mean_separation(split_separation=0.1)
    
    sag = smiles_assignment_group_bonds([sa])
    cst = cluster_means(gcd, labeler, sag, objective=objective)
    
    hierarchies.smarts_hierarchy_print(cst.hierarchy)

.. code-block::

    2024-12-17 10:56:06.570249 Labeling subgraphs
    2024-12-17 10:56:06.571007 Checking consistency...
    Optimization strategy is building steps...
    2024-12-17 10:56:06.571142 The optimization strategy has the following iterations:
    ->   1. op= 1 a=[0] b=1->1 d=0->0 n=0->0
         2. op=-1 a=[0] b=0->0 d=0->0 n=0->0
         3. op= 1 a=[0] b=2->2 d=0->0 n=0->0
         4. op=-1 a=[0] b=0->0 d=0->0 n=0->0
    Targets for this macro step 1:
    1 p0
    N Targets: 1
    Step tracker for current macro step 1
    p0 1
    
    
    *******************
     2024-12-17 10:56:06.571270 iteration=   1 macro=  1/4 X=        0 params=(2|1) G=N S=Y bits=1->2 depth=0->0 branch=0->0
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.2333 Std=    0.0943 N=      3 Min=    1.1000 Max=    1.3000 [*:1]~[*:2]
    =====
    
    2024-12-17 10:56:06.571374 Saving checkpoint to chk.cst.p
    2024-12-17 10:56:06.571891 Collecting SMARTS for p0 N=3/3 and setting to depth=0
     == iteration=   2 macro=  1/4 micro=  1/1 operation=1 params=(2|1) cluster=p0   N= 3 overlap=[0] bits=1->1 depth=0->0 branch=0->0
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=3
    000001 (0, (1, 3))               Mean=    1.3000 Std=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (1, 2))               Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    000003 (0, (2, 4))               Mean=    1.3000 Std=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    2024-12-17 10:56:06.575155 Union merging=3
    2024-12-17 10:56:06.592133 Union merging=1
    2024-12-17 10:56:06.592363 Union is [#6H1X2x0!rA+0:1]!@;-,#[#1,#6;H0,H1;X1,X2;x0;!r;A;+0:3]
    2024-12-17 10:56:06.592841 Generating splits
    2024-12-17 10:56:06.593546 Generating single splits
    2024-12-17 10:56:06.594134 Generated 8 splits
    BIT [_______:1]_;-[_______:2]
    BIT [_______:1]_;#[_______:2]
    BIT [_______:1]_;_[#1______:2]
    BIT [_______:1]_;_[#6______:2]
    BIT [_______:1]_;_[_H0_____:2]
    BIT [_______:1]_;_[_H1_____:2]
    BIT [_______:1]_;_[__X1____:2]
    BIT [_______:1]_;_[__X2____:2]
    2024-12-17 10:56:06.598571 Building tasks
    workspace listening on local host. Remote connections prohibited.
    2024-12-17 10:56:06.598748 P:   0.00%     8/8 IQ:    1 OQ:    0 IP:    0 LF:    0 RF:    0 RIQ:    0 ROQ:    0 RIP:    0  ERC:    0.0 
    2024-12-17 10:56:06.698827 P: 100.00%     0/8 IQ:    0 OQ:    0 IP:    0 LF:    0 RF:    0 RIQ:    0 ROQ:    0 RIP:    0  ERC:    0.0 
        1 CND SPLITS=Y  [*:1]-[*:2]
        2 CND SPLITS=Y  [*:1]#[*:2]
        3 CND SPLITS=Y  [*:1]~[#1:2]
        4 CND SPLITS=N  [*:1]~[#6:2]
        5 CND SPLITS=Y  [*:1]~[H0:2]
        6 CND SPLITS=N  [*:1]~[H1:2]
        7 CND SPLITS=Y  [*:1]~[X1:2]
        8 CND SPLITS=N  [*:1]~[X2:2]
    Finished: 100.00%         8/8
    Closing workspace
    2024-12-17 10:56:06.785257 Calculating partitions for hits=5
    workspace listening on local host. Remote connections prohibited.
    Submitting 5 packets of work
    Closing workspace
    2024-12-17 10:56:06.881913 Unique hits 1/5
        4 HIT S0= 2     -> Sj= 1     [*:1]#[*:2]
        1     DUP [*:1]~[X1:2]
        2     DUP [*:1]~[H0:2]
        3     DUP [*:1]~[#1:2]
        5     DUP [*:1]-[*:2]
    2024-12-17 10:56:06.882460 Searching atoms done; data=3 hits=1
    2024-12-17 10:56:06.882557 Collecting new candidates
    2024-12-17 10:56:06.882584 Scanning done.
    2024-12-17 10:56:06.882590
    
    
    Generating SMARTS on 1
    2024-12-17 10:56:06.890613 Labeling
    2024-12-17 10:56:06.891500 Rebuilding assignments
    2024-12-17 10:56:06.891583 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.2333 Std=    0.0943 N=      3 Min=    1.1000 Max=    1.3000 [*:1]~[*:2]
    =====
    
    Scanning 1 candidates for operation=1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    Cnd.    1/1 p0      X=   -0.20000 dX=   -0.20000 N=      1 C= Y [*:1]#[*:2]
                                                                               
    Nanostep 1: The filtered results of the candidate scan N=1 total=1:
    ->    1 Cnd.    1/1 p0      X=   -0.20000 dX=   -0.20000 N=      1 C= Y [*:1]#[*:2]
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
     >>>>> (0, None, -1) Local dObj   -0.20000 [*:1]#[*:2]
    
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    There were 1 successful operations
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    2024-12-17 10:56:06.979662 Visited {'p0', 'p1'}
    Assignments changed for p0, will retarget
    Assignments changed for p1, will retarget
    Restarting optimization search
    Targets for this macro step 1:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 1
    p0 1
    p1 1
    
    
    *******************
     2024-12-17 10:56:06.980184 iteration=   2 macro=  1/4 X=     -0.2 params=(3|1) G=N S=Y bits=1->2 depth=0->0 branch=0->0
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    =====
    
    2024-12-17 10:56:06.980269 Saving checkpoint to chk.cst.p
    2024-12-17 10:56:06.980596 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   3 macro=  1/4 micro=  1/2 operation=1 params=(3|1) cluster=p0   N= 2 overlap=[0] bits=1->1 depth=0->0 branch=0->0
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Std=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Std=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-12-17 10:56:06.982714 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   4 macro=  1/4 micro=  2/2 operation=1 params=(3|1) cluster=p1   N= 1 overlap=[0] bits=1->1 depth=0->0 branch=0->0
    
    Attempting to split p1:
    S0: [*:1]#[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-12-17 10:56:06.983758 Scanning done.
    2024-12-17 10:56:06.983766
    
    
    Generating SMARTS on 0
    2024-12-17 10:56:06.989963 Labeling
    2024-12-17 10:56:06.990924 Rebuilding assignments
    2024-12-17 10:56:06.991015 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    2024-12-17 10:56:07.078634 Visited set()
    Targets for this macro step 2:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 2
    p0 2
    p1 2
    
    
    *******************
     2024-12-17 10:56:07.079185 iteration=   4 macro=  2/4 X=     -0.2 params=(3|1) G=N S=Y bits=1->2 depth=0->0 branch=0->0
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    =====
    
    2024-12-17 10:56:07.079271 Saving checkpoint to chk.cst.p
    2024-12-17 10:56:07.079586 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   5 macro=  2/4 micro=  1/2 operation=-1 params=(3|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:56:07.079654 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   6 macro=  2/4 micro=  2/2 operation=-1 params=(3|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:56:07.079679 Scanning done.
    2024-12-17 10:56:07.079685
    
    
    Generating SMARTS on 1
    2024-12-17 10:56:07.087497 Labeling
    2024-12-17 10:56:07.088419 Rebuilding assignments
    2024-12-17 10:56:07.088508 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= N [*:1]#[*:2]
                                                                               
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    2024-12-17 10:56:07.178663 Visited {'p1'}
    Targets for this macro step 3:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 3
    p0 3
    p1 3
    
    
    *******************
     2024-12-17 10:56:07.179202 iteration=   6 macro=  3/4 X=     -0.2 params=(3|1) G=N S=Y bits=1->2 depth=0->0 branch=0->0
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    =====
    
    2024-12-17 10:56:07.179308 Saving checkpoint to chk.cst.p
    2024-12-17 10:56:07.179637 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   7 macro=  3/4 micro=  1/2 operation=1 params=(3|1) cluster=p0   N= 2 overlap=[0] bits=2->2 depth=0->0 branch=0->0
    
    Attempting to split p0:
    S0: [*:1]~[*:2]
    Matched N=2
    000001 (0, (1, 3))               Mean=    1.3000 Std=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:1]!@;-[#1H0X1x0!rA+0:3]
    000002 (0, (2, 4))               Mean=    1.3000 Std=    0.0000 N=      1 Min=    1.3000 Max=    1.3000 [#6H1X2x0!rA+0:2]!@;-[#1H0X1x0!rA+0:4]
    
    Skipping p0 since all graphs are the same
    2024-12-17 10:56:07.181858 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=   8 macro=  3/4 micro=  2/2 operation=1 params=(3|1) cluster=p1   N= 1 overlap=[0] bits=2->2 depth=0->0 branch=0->0
    
    Attempting to split p1:
    S0: [*:1]#[*:2]
    Matched N=1
    000001 (0, (1, 2))               Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [#6H1X2x0!rA+0:1]!@;#[#6H1X2x0!rA+0:2]
    
    Skipping p1 since all graphs are the same
    2024-12-17 10:56:07.182951 Scanning done.
    2024-12-17 10:56:07.182959
    
    
    Generating SMARTS on 0
    2024-12-17 10:56:07.189262 Labeling
    2024-12-17 10:56:07.190217 Rebuilding assignments
    2024-12-17 10:56:07.190311 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    =====
    
    Scanning 0 candidates for operation=1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=0 total=0:
    
    Nanostep 1: The filtered results of the candidate scan N=0 total=0:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    2024-12-17 10:56:07.275306 Visited set()
    Targets for this macro step 4:
    1 p0
    2 p1
    N Targets: 2
    Step tracker for current macro step 4
    p0 4
    p1 4
    
    
    *******************
     2024-12-17 10:56:07.275926 iteration=   8 macro=  4/4 X=     -0.2 params=(3|1) G=N S=Y bits=1->2 depth=0->0 branch=0->0
    *******************
    
    Tree:
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    =====
    
    2024-12-17 10:56:07.276021 Saving checkpoint to chk.cst.p
    2024-12-17 10:56:07.276354 Collecting SMARTS for p0 N=2/3 and setting to depth=0
     == iteration=   9 macro=  4/4 micro=  1/2 operation=-1 params=(3|1) cluster=p0   N= 2 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:56:07.276424 Collecting SMARTS for p1 N=1/3 and setting to depth=0
     == iteration=  10 macro=  4/4 micro=  2/2 operation=-1 params=(3|1) cluster=p1   N= 1 overlap=[0] bits=0->0 depth=0->0 branch=0->0
    
    2024-12-17 10:56:07.276447 Scanning done.
    2024-12-17 10:56:07.276453
    
    
    Generating SMARTS on 1
    2024-12-17 10:56:07.285620 Labeling
    2024-12-17 10:56:07.286525 Rebuilding assignments
    2024-12-17 10:56:07.286612 Rebuilding mappings
    Tree:
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    =====
    
    Scanning 1 candidates for operation=-1
    workspace listening on local host. Remote connections prohibited.
    The unfiltered results of the candidate scan N=1 total=1:
    Cnd.    1/1 p0      X=    0.00000 dX=    0.20000 N=      3 C= N [*:1]#[*:2]
                                                                               
    Nanostep 1: The filtered results of the candidate scan N=0 total=1:
    There were 0 successful operations
    **  0   0 p0    Mean=    1.3000 Std=    0.0000 N=      2 Min=    1.3000 Max=    1.3000 [*:1]~[*:2]
    **  1   1 p1    Mean=    1.1000 Std=    0.0000 N=      1 Min=    1.1000 Max=    1.1000 [*:1]#[*:2]
    2024-12-17 10:56:07.378695 Visited {'p1'}
    Nothing found. Done.
    Start time: 2024-12-17 10:56:06.522996
    End   time: 2024-12-17 10:56:07.380206
    **  0 p0 [*:1]~[*:2]
    **   1 p1 [*:1]#[*:2]

Similar to the categorical case, the same SMARTS pattern was found.
