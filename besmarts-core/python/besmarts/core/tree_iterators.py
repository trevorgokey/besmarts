"""
besmarts.core.tree_iterators

Iterators for traversing trees.tree_index objects
"""

from besmarts.core import trees


def yield_single(hidx: trees.tree_index, v: trees.tree_node):
    yield v


def yield_if_single(hidx, v, select, state):
    if select is None or select == v.name:
        if state is None or state == v.state:
            yield from yield_single(hidx, v)


def yield_if(hidx, v, select, state):
    if hasattr(v, "__iter__"):
        for vi in v:
            yield from yield_if(hidx, vi, select, state)
    else:
        yield from yield_if_single(hidx, v, select, state)


def tree_iter_depth_first_single(
    hidx: trees.tree_index, v: trees.tree_node, select=None, state=None
):
    for c in hidx.below[v.index]:
        c = hidx.nodes[c]
        yield from tree_iter_depth_first_single(hidx, c, select, state)
    yield from yield_if(hidx, v, select, state)


def tree_iter_depth_first(
    hidx: trees.tree_index, v: trees.tree_node, select=None, state=None
):
    if hasattr(v, "__iter__"):
        for vi in v:
            yield from tree_iter_depth_first(hidx, vi, select, state)
    else:
        yield from tree_iter_depth_first_single(hidx, v, select, state)


def tree_iter_breadth_first_single(
    hidx: trees.tree_index, v: trees.tree_node, select=None, state=None
):
    if hidx.above.get(v.index) is None:
        yield_single(hidx, v)
    for c in hidx.below[v.index]:
        c = hidx.nodes[c]
        yield from yield_if(hidx, c, select, state)
    for c in hidx.below[v.index]:
        c = hidx.nodes[c]
        yield from tree_iter_breadth_first_single(hidx, c, select, state)


def tree_iter_breadth_first(
    hidx: trees.tree_index, v: trees.tree_node, select=None, state=None
):
    if hasattr(v, "__iter__"):
        for vi in v:
            yield from tree_iter_breadth_first(hidx, vi, select, state)
    else:
        yield from tree_iter_breadth_first_single(hidx, v, select, state)


def tree_iter_dive_single(hidx, v, select=None, state=None):
    yield from yield_if_single(hidx, v, select, state)
    for c in hidx.below[v.index]:
        c = hidx.nodes[c]
        yield from tree_iter_dive_single(hidx, c, select, state)


def tree_iter_dive_single_reverse(hidx, v, select=None, state=None):
    for c in reversed(hidx.below[v.index]):
        c = hidx.nodes[c]
        yield from tree_iter_dive_single_reverse(hidx, c, select, state)
    yield from yield_if_single(hidx, v, select, state)


def tree_iter_dive(hidx, v, select=None, state=None):
    if hasattr(v, "__iter__"):
        for vi in v:
            yield from tree_iter_dive(hidx, vi, select, state)
    else:
        yield from tree_iter_dive_single(hidx, v, select, state)


def tree_iter_dive_reverse(hidx, v, select=None, state=None):
    if hasattr(v, "__iter__"):
        for vi in reversed(v):
            yield from tree_iter_dive_reverse(hidx, vi, select, state)
    else:
        yield from tree_iter_dive_single_reverse(hidx, v, select, state)


def tree_iter_to_root_single(hidx, v, select=None, state=None):
    yield from yield_if_single(hidx, v, select, state)
    if hidx.above[v.index] is not None:
        parent = hidx.nodes[hidx.above[v.index]]
        yield from tree_iter_to_root_single(hidx, parent, select, state)


def tree_iter_to_root(hidx, v, select=None, state=None):
    if hasattr(v, "__iter__"):
        for vi in v:
            yield from tree_iter_to_root(hidx, vi, select, state)
    else:
        yield from tree_iter_to_root_single(hidx, v, select, state)
