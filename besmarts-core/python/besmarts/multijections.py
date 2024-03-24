"""
besmarts.multijections

Multijections are still very experimental
"""

from besmarts.core import graph_bitwise
from besmarts.core import graphs
from besmarts.core import mapper

# TODO: transfer this to use workspaces
# only allow add_nodes=0|1 for now until I can decide what it means to use 2|3
class multijection:
    def __init__(
        self,
        domain: graphs.structure,
        fill=0,
    ):
        self.reference = graphs.structure_copy(domain)
        self.codomain_map = {}

        self.domain = graphs.structure_copy(domain)
        self.codomain = {}
        self.add_nodes = 1
        self.fill = fill

    def reset(self):
        self.domain = graphs.structure_copy(self.reference)
        self.codomain_map.clear()
        self.codomain.clear()

    def refresh(self):
        old = {k: v for k, v in self.codomain_map.items()}
        self.codomain.clear()
        self.codomain_map.clear()
        for oldG in old:
            self.add(oldG)

    def remove(self, G: graphs.structure):
        if G in self.codomain_map:
            g = self.codomain_map.pop(G)
            self.codomain.pop(g)
            self.domain = graphs.structure_copy(self.reference)
            self.refresh()

    def add(self, G: graphs.structure, skip=None):
        g = self.codomain_map.get(G)
        if g is not None:
            return

        T = mapper.map_to(
            self.domain,
            G,
            strict=False,
            equality=False,
            skip=skip,
            mode="high",
            add_nodes=self.add_nodes,
            fill=self.fill,
        )

        if T.map is None:
            self.codomain[G] = None
            self.codomain_map[G] = G
            return

        d, g, m = T.G, T.H, T.map

        # d, g, m = mapper_compose_graphs(
        #     self.domain,
        #     G,
        #     M,
        #     add_nodes=self.add_nodes,
        #     fill_new_nodes=self.fill_new_nodes,
        # )

        # print("added", m, G, g)

        i = 0
        # print("domain", gcd.smarts_encode(self.domain))
        # print("new domain:", gcd.smarts_encode(d))
        while hash(self.domain) != hash(d):

            # print(gcd.smarts_encode(d))
            i += 1
            # print(i, "set to new domain")

            self.domain = d
            old = {k: v for k, v in self.codomain_map.items()}
            for j, (oldG, oldg) in enumerate(old.items()):
                # print(i,j)


                _skip = self.codomain[oldg]
                # print("Visting existing", i,j, gcd.smarts_encode(oldg), "skip is ", _skip)
                if len(_skip) != len(self.domain.nodes):
                    # print("Iter map", i,j, len(_skip), len(self.domain.nodes))
                    _skip = self.codomain.pop(oldg)
                    # if _skip:
                    #     _skip = {k:v for k,v in _skip.items() if k in self.domain.nodes and v in oldg.nodes}
                    if oldG in self.codomain_map:
                        self.codomain_map.pop(oldG)
                    _T = mapper.map_to(
                        self.domain,
                        oldg,
                        strict=False,
                        equality=False,
                        skip=_skip,
                        mode="high",
                        add_nodes=self.add_nodes,
                        fill=self.fill,
                    )

                    # print("Iter compose", i,j)
                    # print("Iter compose", _m)
                    d, _g, _m = _T.G, _T.H, _T.map
                    # d, _g, _m = mapper_compose_graphs(
                    #     self.domain,
                    #     oldg,
                    #     _m,
                    #     add_nodes=self.add_nodes,
                    #     fill_new_nodes=self.fill_new_nodes,
                    # )
                    # print("Iter compose 2", _m)
                    self.codomain[_g] = _m
                    self.codomain_map[oldG] = _g

        self.domain = d
        self.codomain[g] = m
        self.codomain_map[G] = g

def multijection_mcs(A, ref, fill=True):

    M = multijection(ref, fill=fill)
    M.add_nodes = 1
    for s in A:
        M.add(s)

    return M

def multijection_fit_reference(A, ref, fill=True):

    M = multijection(ref, fill=fill)
    M.add_nodes = 3
    for s in A:
        M.add(s)

    return M

def multijection_branch(A, ref, fill=True):

    M = multijection(ref, fill=fill)
    M.add_nodes = 1
    for s in A:
        M.add(s)

    return M

def multijection_union(T: multijection):

    result = graphs.structure_clear(T.domain)
    for codomain, M in T.codomain.items():
        result = graph_bitwise.subgraph_union(result, codomain, M)
    return result


def multijection_intersection_conditional(T: multijection):
    result = graphs.structure_clear(T.domain)

    for codomain, M in T.codomain.items():
        result = intersection_conditional(result, codomain, None, map=M)
    return result


def multijection_intersection(T: multijection):
    result = graphs.structure_clear(T.domain)

    reference = T.domain
    T.add_nodes
    config = configs.mapper_config(T.add_nodes, T.fill, "low")

    return intersection_list_parallel(
        list(T.codomain_map.values()), config, reference=reference
    )

def multijection_save(T: multijection, fname):
    atom_order = {}
    with open(fname, "w") as f:
        g = graphs.structure_remove_unselected(T.domain)
        atom_order = {i: j for i, j in enumerate(g.select)}

        lines = codec_native.graph_save(g, atom_order)
        f.write("\n".join(lines[:3]) + "\n")
        f.write(" ".join(lines[3:]) + "\n")

        order = sorted(atom_order)
        for g, m in T.codomain.items():
            this_order = {i: m[atom_order[i]] for i in sorted(order)}
            lines = codec_native.graph_save(g, this_order)
            f.write(" ".join(lines[3:]) + "\n")


def multijection_load(
    fname,
    mode="high",
    strict=False,
    equality=False,
    add_nodes=False,
    fill=False,
) -> multijection:
    header = []
    domain = None
    domain_order = []
    n_atom_prims = 0
    n_bond_prims = 0
    T = None
    with open(fname) as f:
        for i, line in enumerate(f):
            tokens = line.split()
            if tokens[0] == "#GRAPH":
                header.append(tokens)
            elif tokens[0] == "#ATOM":
                header.append(tokens)
                n_atom_prims = len(tokens[1:])
            elif tokens[0] == "#BOND":
                header.append(tokens)
                n_bond_prims = len(tokens[1:])
            else:
                lines = []
                i = 1
                order = []
                while i < len(tokens):
                    if tokens[i - 1] == tokens[i]:
                        lines.append(tokens[i - 1 : (i + 1 + n_atom_prims)])
                        order.append(abs(int(tokens[i])))
                        i += n_atom_prims + 2
                    else:
                        lines.append(tokens[i - 1 : (i + 1 + n_bond_prims)])
                        i += n_bond_prims + 2
                lines = header + lines
                if domain:
                    g = codec_native.graph_load(lines)
                    m = {i: j for i, j in zip(domain_order, order)}
                    T.codomain[g] = m
                    T.codomain_map[g] = g
                else:
                    domain = codec_native.graph_load(lines)
                    domain_order = order
                    T = multijection(
                        domain,
                        mode=mode,
                        strict=strict,
                        equality=equality,
                        add_nodes=add_nodes,
                        fill=fill,
                    )

    return T
