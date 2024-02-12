"""
besmarts.core.tree_merge
"""

from besmarts.core import trees
from besmarts.core import tree_iterators


def tree_index_merge(
    hA: trees.tree_index,
    rootA: trees.tree_node,
    hB: trees.tree_index,
    rootB: trees.tree_node,
    index=None,
):
    node = trees.tree_node(0, "", "", rootB.name)
    hA.node_add(rootA.index, node, index=index)
    up = node.index
    mapping = {rootB.index: up}
    for eb in tree_iterators.tree_iter_breadth_first(hB, rootB):
        up = mapping[hB.above[eb.index]]
        node = trees.tree_node(None, "", "", eb.name)
        node = hA.node_add(up, node)
        mapping[eb.index] = node.index

    return hA
