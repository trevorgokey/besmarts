"""
besmarts.core.trees

Graphs with no cycles, a rough concept of directed edges, and have a root nodes
"""

from typing import Dict, List


class tree_node:
    __slots__ = ("index", "category", "type", "name")

    def __init__(self, index, category, type, name):
        self.index: int = index
        self.category: str = category
        self.type: str = type
        self.name: str = name

    def __hash__(self):
        return hash((self.index, self.category, self.type, self.name))

    def copy(self):
        return tree_node_copy(self)


node_id = int


class tree_index:
    """
    A tree data structure that serves as a hierarchical index.
    """

    __slots__ = ("nodes", "above", "below")

    def __init__(self):
        self.nodes: Dict[node_id, tree_node] = {}
        self.above: Dict[node_id, node_id] = {}
        self.below: Dict[node_id, List[node_id]] = {}

    def node_remove(self, idx: int) -> tree_node:
        return tree_index_node_remove(self, idx)

    def node_remove_by_name(self, name: str) -> List[tree_node]:
        return tree_index_node_remove_by_name(self, name)

    def node_add(self, above_idx, tnode: tree_node, index=None) -> tree_node:
        return tree_index_node_add(self, above_idx, tnode, index)

    def node_add_below(self, above_idx, index=None) -> tree_node:
        tnode = tree_node(0, "", "", "")
        return tree_index_node_add(self, above_idx, tnode, index)

    def copy(self):
        return tree_index_copy(self)


def tree_index_node_add(
    tree: tree_index, above_idx, tnode: tree_node, index=None
) -> tree_node:
    """
    Add a entry to the index.

    Parameters
    ----------
    up_idx : int
        The entry index that the incoming entry should be attached to
    hent : hentry
        The entry to add
    index : int
        The position to add the entry of the parent entry has existing children. A
        value of 0 will insert the new entry in front. A value of None will append
        the hentry at the end.

    Returns
    -------
    hentry
        The entry with a new index
    """

    assert isinstance(above_idx, int) or above_idx is None

    if tree.nodes:
        idx = max(tree.nodes) + 1
    else:
        idx = 0

    below = tree.below.get(above_idx, [])

    if above_idx is not None:
        if index is None:
            below.append(idx)
        else:
            below.insert(index, idx)

    tree.below[above_idx] = below

    tree.above[idx] = above_idx
    tnode.index = idx

    tree.nodes[idx] = tnode
    tree.below[idx] = []

    return tnode


def tree_index_node_remove_by_name(tree, name: str) -> List[tree_node]:
    """
    Remove an entry from the index and repair the hierarchy. The children of the
    entry are attached to the entry's parent.

    Parameters
    ----------
    idx : int
        The index of the entry to remove

    Returns
    -------
    hentry
        The removed entry
    """
    removed = []
    to_remove = []
    for idx, node in list(tree.nodes.items()):
        if node.name == name:
            to_remove.append(idx)
    for idx in to_remove:
        node = tree_index_node_remove(tree, idx)
        removed.append(node)

    return removed


def tree_index_node_remove(tree, idx):
    """
    Remove an entry from the index and repair the hierarchy. The children of the
    entry are attached to the entry's parent.

    Parameters
    ----------
    idx : int
        The index of the entry to remove

    Returns
    -------
    hentry
        The removed entry
    """

    tnode = tree.nodes.pop(idx)
    above = tree.above.pop(idx)
    below = tree.below.pop(idx)

    for x in below:
        tree.above[x] = above

    # remove reference from up
    if above is None:
        above_below = None
    else:
        above_below = tree.below[above]
        pos = above_below.index(idx)
        above_below.remove(idx)

        # shift tree into above's children
        above_below = above_below[:pos] + below + above_below[pos:]
        tree.below[above] = above_below

    return tnode


def tree_index_copy(tree: tree_index):
    _nodes = {}
    for k, v in tree.nodes.items():
        _nodes[k] = tree_node(v.index, v.category, v.type, v.name)
    above = tree.above.copy()
    below = {}
    for k, v in tree.below.items():
        below[k] = v.copy()
    t = tree_index()
    t.nodes.update(_nodes)
    t.above.update(above)
    t.below.update(below)
    return t


def tree_index_node_depth(tree: tree_index, node: tree_node):
    """
    Return the depth of the hentry in the hierarchical index.

    Parameters
    ----------
    tree : tree_index
        The input tree

    node : tree_node
        The node to get the depth of

    Returns
    -------
    int
        The node depth.
    """
    n = node.index
    l = 0
    while tree.above[n] is not None:
        l += 1
        n = tree.above[n]
        # print(n, l, tree.below[n])
    return l


def tree_index_roots(t: tree_index):
    roots = [t.nodes[i] for i, x in t.above.items() if x is None]
    return roots


def tree_node_copy(n: tree_node):
    return tree_node(n.index, n.category, n.type, n.name)
