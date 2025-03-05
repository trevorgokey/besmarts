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
        self.bit_search_min: int = int(bit_search_min)
        self.bit_search_limit: None | int = bit_search_limit
        self.branch_min: None | int = branch_min
        self.branch_limit: None | int = branch_limit
        self.branch_depth_min: None | int = branch_depth_min
        self.branch_depth_limit: None | int = branch_depth_limit
        self.unique: bool = bool(unique)
        self.return_matches: bool = bool(return_matches)
        self.max_splits = max_splits
        self.split_general = split_general
        self.split_specific = split_specific
        self.unique_complements: bool = bool(unique_complements)
        self.unique_complements_prefer_min: bool = bool(
            unique_complements_prefer_min
        )
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
