"""
besmarts.core.graph_multijections

Multijections are still very experimental
"""

from besmarts.core import graph_bitwise

# TODO: transfer this to use workspaces
class multijection:
    def __init__(
        self,
        domain: graphs.structure,
        mode="high",
        strict=False,
        equality=False,
        add_nodes=0,
        fill_new_nodes=0,
    ):
        self.reference = graphs.structure_copy(domain)
        self.codomain_map = {}

        self.domain = graphs.structure_copy(domain)
        self.codomain = {}
        self.strict = strict
        self.equality = equality
        self.mode = mode
        self.add_nodes = add_nodes
        self.fill_new_nodes = fill_new_nodes

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

        T = map_to(
            self.domain,
            G,
            strict=self.strict,
            equality=self.equality,
            skip=skip,
            mode=self.mode,
            add_nodes=self.add_nodes,
            fill=self.fill_new_nodes,
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
        # print("start")
        while hash(self.domain) != hash(d):
            i += 1
            # print(i)
            self.domain = d
            old = {k: v for k, v in self.codomain_map.items()}
            for j, (oldG, oldg) in enumerate(old.items()):
                # print(i,j)

                # print("Iter", i,j)

                _skip = self.codomain[oldg]
                if False or len(_skip) != len(self.domain.nodes):
                    # print("Iter map", i,j, len(_skip), len(self.domain.nodes))
                    _skip = self.codomain.pop(oldg)
                    # if _skip:
                    #     _skip = {k:v for k,v in _skip.items() if k in self.domain.nodes and v in oldg.nodes}
                    if oldG in self.codomain_map:
                        self.codomain_map.pop(oldG)
                    _T = map_to(
                        self.domain,
                        oldg,
                        strict=self.strict,
                        equality=self.equality,
                        skip=_skip,
                        mode=self.mode,
                        add_nodes=self.add_nodes,
                        fill=self.fill_new_nodes,
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

def multijection_union(T: multijection):

    result = graphs.structure_clear(T.domain)
    for codomain, M in T.codomain.items():
        result = graph_bitwise.union(result, codomain, map=M)
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
    config = configs.mapper_config(T.add_nodes, T.fill_new_nodes, "low")

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
    fill_new_nodes=False,
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
                        fill_new_nodes=fill_new_nodes,
                    )

    return T
