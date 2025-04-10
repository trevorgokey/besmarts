"""
besmarts.core.configs
"""

import os

processors = os.cpu_count()
remote_compute_enable = True
workqueue_port = 55555
compute_verbosity = 0
compute_task_chunksize = 1000

compute_runtime = {
    "is_remote": False,
    "verbosity": 0,
    "task_chunksize": 1000,
}
# for computations that may suffer from error accumulation, round the numbers
# to this precision
precision = 12

match_print_limit = 100

class smiles_perception_config:
    def __init__(
        self,
        allow_unconnected,
        protonate,
        strip_hydrogen,
        aromaticity,
    ):
        self.allow_unconnected: bool = bool(allow_unconnected)
        self.protonate: bool = bool(protonate)
        self.strip_hydrogen: bool = bool(strip_hydrogen)
        self.aromaticity: str = str(aromaticity)


class mapper_config:
    def __init__(self, add_nodes, fill_new_nodes, mode):
        """
        Initialize configuration.

        Parameters
        ----------
        add_nodes : bool
            Whether to add node nodes to the structures if there are unmapped nodes

        fill_new_nodes : bool
            If new nodes are added, determines whether the new node will be filled,
            otherwise it will be empty.
        mode : str, "low" or "high"
            During the mapping, determine whether the score should prefer highest
            overlap or lowest overlap

        """
        self.add_nodes: bool = bool(add_nodes)
        self.fill_new_nodes: bool = bool(fill_new_nodes)
        self.mode: str = str(mode)

    def copy(self):
        return mapper_config(self.add_nodes, self.fill_new_nodes, self.mode)


class smarts_extender_config:
    """Configures how many atoms in the reference data should we consider when
    writing our SMARTS strings

    Determines `structure.select`.
    """
    __slots__ = ("depth_min", "depth_max", "include_hydrogen")

    def __init__(self, depth_min: int, depth_max: int, include_hydrogen: bool):
        self.depth_min: int = depth_min
        self.depth_max: int = depth_max
        self.include_hydrogen: bool = include_hydrogen

    def copy(self):
        return smarts_extender_config(
            self.depth_min, self.depth_max, self.include_hydrogen
        )


class smarts_splitter_config:
    """
    Configuration object for splitting SMARTS patterns.

    Parameters beginning ``bit_search`` choose the number of bits that may be
    considered in the primitives of a SMARTS pattern ($k$ in eq. 11 of `Gokey
    and Mobley 2023`_). This is the number of decisions made in each split. When
    making a pattern more specific, this is the number of bits that should be
    flipped.

    Parameters beginning ``branch`` configure extending the SMARTS pattern by
    adding new wildcard atoms that can then be specialized via bit search. Atoms
    can be added by extending the graph out to ``branch_depth`` from the SMARTS
    pattern's mapped atoms, which index the atoms that form the internal
    coordinate. The mapped atoms have a depth of 0, neighbouring atoms have
    1, and so forth. Once this search space is defined, the actual number of
    atoms to keep is set by ``branch_min`` and ``branch_limit``. The bit search
    is spread across all atoms, so extending the graph by more atoms can dilute
    this search.

    Attributes
    ==========
    bit_search_min
        Minimum number of decisions made when splitting a SMARTS string.
    bit_search_limit
        Maximum number of decisions made when splitting a SMARTS string.
    branch_min
        Minimum number of new atoms to add in each split.
    branch_limit
        Maximum number of new atoms to add in each split.
    self.branch_depth_min
        Minimum ``branch_depth`` of new atoms to consider
    branch_depth_limit
        Maximum ``branch_depth`` of new atoms to consider
    unique
        If ``True``, do not return multiple splits that produce the same
        partitioning of the reference structures; otherwise, do not perform this
        pruning.
    return_matches
        Perform substructure matching on each split, or just return them?
        ``True`` permits an optimization.
    max_splits
        How many SMARTS patterns are returned. Low values allow fast returns,
        but not necessarily the best split.
    split_general
        If ``True``, return splits that tend to produce broad hierarchies by
        turning on bits (using $\hat{b_i}$), otherwise do
        not return these splits.
    split_specific
        If ``True``, return splits that tend to produce deep hierarchies by
        turning off bits (using the partial inverse of $\hat{b_i}$), otherwise
        do not return these splits.
    unique_complements
        If ``true``, do not return multiple splits that produce complementary
        partitions of the reference structures; otherwise, do not perform this
        pruning. Complementary patterns are those that produce the same
        partition by matching the opposite reference structures.
    unique_complements_prefer_min
        If ``True`` and ``unique_complements == True``, keep the complement
        where the SMARTS pattern matches the fewest structures. If ``False``
        and ``unique_complements == True``, keep the complement where the SMARTS
        pattern matches the most structures.
    primitives
        Collection of SMARTS primitives that should be split on.

    .. _Gokey and Mobley 2023:
        https://doi.org/10.26434/chemrxiv-2023-v969f-v3
    """
    __slots__ = (
        "bit_search_min",
        "bit_search_limit",
        "branch_min",
        "branch_limit",
        "branch_depth_min",
        "branch_depth_limit",
        "unique",
        "return_matches",
        "max_splits",
        "split_general",
        "split_specific",
        "unique_complements",
        "unique_complements_prefer_min",
        "primitives",
    )

    def __init__(
        self,
        bit_search_min,
        bit_search_limit,
        branch_min,
        branch_limit,
        branch_depth_min,
        branch_depth_limit,
        unique=True,
        return_matches=True,
        max_splits=None,
        split_general=True,
        split_specific=True,
        unique_complements=True,
        unique_complements_prefer_min=True,
        primitives=None,
    ):
        # Bit search is roughly the number of bits to flip when making a string
        # more specific (ie, number of decisions to make when splitting a SMARTS
        # pattern
        # Minimum number of decisions
        self.bit_search_min: int = int(bit_search_min)
        # Maximum number of decisions
        self.bit_search_limit: None | int = bit_search_limit
        # Branch is about extending the graph
        # Atoms can be added by extending the graph out to branch_depth from
        # the indexed atoms/topology atoms/internal coordinate atoms. The atoms
        # in the internal coordinate have a depth of 0, neighbouring atoms have
        # 1, etc. Once this search space is defined, the actual number of atoms
        # to keep is set by branch_min and branch_limit. The bit search is
        # spread across all atoms, so extending the graph by more atoms can
        # dilute this search.
        # Minimum number of new atoms to add
        self.branch_min: None | int = branch_min
        # Maximum number of new atoms to add
        self.branch_limit: None | int = branch_limit
        # Minimum number of new atoms in each "direction" to consider
        self.branch_depth_min: None | int = branch_depth_min
        # Maximum number of new atoms in each "direction" to consider
        self.branch_depth_limit: None | int = branch_depth_limit
        # Do you want all splits, or just those that provide unique partitions
        self.unique: bool = bool(unique)
        # Perform substructure matching on each split, or just return them?
        # True permits an optimization
        self.return_matches: bool = bool(return_matches)
        # How many SMARTS patterns are returned
        # Low values allow fast returns, but not necessarily the best split
        self.max_splits: int = max_splits
        # Whether to return splits that tend to produce broad hierarchies by
        # turning on bits (using i_hat_j)
        self.split_general: bool = split_general
        # Whether to return splits that tend to produce deep hierarchies by
        # turning off bits (using the inverse of i_hat_j)
        self.split_specific: bool = split_specific
        # Complementary patterns are those that produce the same partition by
        # matching the opposite reference structures. if unique_complements is
        # true, prune complementary partitions; if false, keep them. Independent
        # of ``unique``, which is the equivalent for exactly identical
        # partitions.
        self.unique_complements: bool = bool(unique_complements)
        # True: keep the complement where the SMARTS pattern matches the fewest structures
        # False: keep the complement where the SMARTS pattern matches the most structures
        self.unique_complements_prefer_min: bool = bool(
            unique_complements_prefer_min
        )
        # Which of the SMARTS primitives should be split on?
        self.primitives = primitives

    def copy(self):
        return smarts_splitter_config(
            self.bit_search_min,
            self.bit_search_limit,
            self.branch_limit,
            self.branch_depth_limit,
            self.unique,
            self.return_matches,
            self.split_general,
            self.split_specific,
        )


class smarts_perception_config:
    __slots__ = ("splitter", "extender")

    def __init__(
        self,
        splitter: smarts_splitter_config,
        extender: smarts_extender_config,
    ):
        self.splitter = splitter
        self.extender = extender

    def copy(self):
        return smarts_perception_config_copy(self)


def smarts_perception_config_copy(self):
    split = self.splitter.copy()
    extend = self.extender.copy()
    return smarts_perception_config(split, extend)
