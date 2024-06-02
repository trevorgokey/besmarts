"""
besmarts.core.tree_visitors

Likely a shoe-in for the old offsb tree operations
"""

from typing import Generator
import datetime
import copy

from besmarts.core.chem import bechem
from besmarts.core.primitives import primitive_key, element_tr
from besmarts.core import graphs
from besmarts.core import chem


"""
ideally this creates a task for each node and returns a list
"""

class tree_visitor:
    def on_start(self, g):
        return None

    def on_tree(self):
        return None

    def on_stop(self, result):
        return None

    def on_descend(self):
        return None

    def on_ascend(self):
        return None

    def on_node_enter(self, idx):
        return None

    def on_node(self, idx, atom: bechem, *args, **kwargs):
        return None

    def on_node_edges(self, node, edges):
        return None

    def on_node_exit(self, idx):
        return None

    def on_edge_enter(self, idx):
        return None

    def on_edge(self, idx, bond: bechem, *args, **kwargs):
        return None

    def on_edge_exit(self, idx):
        return None

