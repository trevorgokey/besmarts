"""
besmarts.core.rulesets

Defines the rules for determining whether a given graph is "valid" where the
current state is whether the graph is a valid molecule/SMILES
"""

from besmarts.core.primitives import primitive_key
from besmarts.core.chem import bechem
from besmarts.core.graph_visitors import graph_visitor
from besmarts.core.graphs import graph
from besmarts.core import arrays


class visitor_ruleset(graph_visitor):
    def __init__(self, primitive_codecs):
        self.primitive_codecs = primitive_codecs
        self.beg: graph = None
        self.small_molecule_cutoff = 5
        self.hydrogen_rule = ("n_minus", {None: None})
        self.noncarbon_rule = ("n_max", {None: None})


    def on_start(self, g: graph):
        self.beg = g
        return None

    def on_stop(self, g: graph, result):
        self.beg: graph = None
        return None

    def on_graph(self, nodes, edges):

        codec = self.primitive_codecs[primitive_key.FORMAL_CHARGE]

        if abs(codec.count_charge_smiles(nodes)) > 2:
            return False

        N = len(self.beg.nodes)
        if len(nodes) == N:
            # codec = self.primitive_codecs[primitive_key.ELEMENT]
            # C = codec.count_element_smiles(nodes, n=6)

            name, rule = self.hydrogen_rule
            if name == "n_minus":
                codec = self.primitive_codecs[primitive_key.HYDROGEN]
                H = codec.count_hydrogen_smiles(nodes)

                x = rule.get(N, rule[None])
                if type(x) is float:
                    # print("H", H, "total", N, (N-H)/(N+H))
                    if H/N < x:
                        return False
                elif type(x) is int:
                    if H < N - x:
                        return False

            if self.noncarbon_rule[0] == "n_max":
                codec = self.primitive_codecs[primitive_key.ELEMENT]

                C = codec.count_carbon_smiles(nodes)
                N = len(self.beg.nodes)
                calc = len(nodes) - C
                rule = self.noncarbon_rule[1]

                x = rule.get(N, rule[None])
                # print("Calc:", calc, len(nodes), C, len(self.beg.nodes), x)
                if x is not None and calc > x:
                    return False

        return True

    def on_nodes(self, nodes):


        return True

    def on_node(self, idx, besmarts: bechem, *args, **kwargs):

        valid = True
        name = primitive_key.FORMAL_CHARGE
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        q = codec.decode_int(arr.on_first())

        name = primitive_key.ELEMENT
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        e = codec.decode_int(arr.on_first())

        name = primitive_key.VALENCE
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        valence = codec.decode_int(arr.on_first())

        # name = primitive_key.RING_SMALLEST
        # codec = self.primitive_codecs[name]
        # arr = besmarts.primitives[name]
        # r = codec.decode_int(arr.on_first())

        # name = primitive_key.CONNECTIVITY_RING
        # codec = self.primitive_codecs[name]
        # arr = besmarts.primitives[name]
        # x = codec.decode_int(arr.on_first())

        # name = primitive_key.AROMATIC
        # codec = self.primitive_codecs[name]
        # arr = besmarts.primitives[name]
        # a = codec.decode_int(arr.on_first())

        name = primitive_key.HYDROGEN
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        h = codec.decode_int(arr.on_first())

        if False:
            pass
        # elif a == 1 and (x == 0 or r < 2 or e not in [6,7,8,15]):
        #     valid = False
        elif e == 6:
            if q != 0:
                valid = False
            elif valence != 4:
                valid = False

        elif e == 7:
            if q < 0:
                valid = False
            elif valence not in [3]:
                valid = False

        elif e == 8:
            if q > 0:
                valid = False
            if h > 2:
                valid = False
            elif valence != 2:
                valid = False

        elif e in [1, 9, 17, 35, 53]:
            if q != 0:
                valid = False
            elif valence != 1:
                valid = False
        elif e in [15]:
            # if q != 0:
            #     valid = False
            # if h > 1:
            #     valid = False
            if valence not in [3,5,6]:
                valid = False

        elif e in [16]:
            if q not in [0, 1]:
                valid = False
            # elif h > 1:
            #     valid = False
            elif valence not in [2,4,6]:
                valid = False
        elif e in [5]:
            if q != 0:
                valid = False
            elif h > 0:
                valid = False
            elif valence != 3:
                valid = False
        elif e in [1]:
            # if r != 0:
            #     valid = False
            # elif x != 0:
            #     valid = False
            if h > 0:
                valid = False

        if not valid:
            return False

        # if e == 1:
        #     if valence != 1:
        #         continue
        #     if 0:
        #         if c and not frag.primitives["chirality"][0]:  # chiral
        #             continue
        # elif e == 6:
        #     if 0:
        #         if c and not frag.primitives["chirality"][0]:  # chiral
        #             if h + n_edges != 4 or h > 1:
        #                 continue
        #         else:  # not chiral
        #             if (h == 1 and n_edges == 3) or (n_edges == 4 and h == 0):
        #                 continue
        # elif frag.primitives["element"][7]:
        #     if valence != 3:
        #         continue
        #     if 0:
        #         if c and not frag.primitives["chirality"][0]:  # chiral
        #             if h + n_edges != 3 or h > 1 or (h + n_edges == 3 and q == 1):
        #                 continue
        #         else:  # not chiral
        #             if h <= 1 and n_edges >= 2:
        #                 continue
        # elif frag.primitives["element"][8]:
        #     if valence != 2:
        #         continue
        #     if 0:
        #         if c and not frag.primitives["chirality"][0]:  # chiral
        #             continue

        return valid

    def on_node_edges(self, idx, node, edges):

        bo = 0.0
        naro = 0
        # has_double = False
        has_single = 0
        has_stereo = False

        name = primitive_key.FORMAL_CHARGE
        codec = self.primitive_codecs[name]
        arr = node.primitives[name]
        q = codec.decode_int(arr.on_first())

        name = primitive_key.VALENCE
        codec = self.primitive_codecs[name]
        arr = node.primitives[name]
        valence = codec.decode_int(arr.on_first())

        name = primitive_key.ELEMENT
        codec = self.primitive_codecs[name]
        arr = node.primitives[name]
        sym = codec.decode_int(arr.on_first())

        #         if len(valence):
        #             valence = valence[0]
        #         else:
        #             if e == 6:
        #                 valence = 4
        #             elif e == 7:
        #                 valence = 3
        #             elif e == 8:
        #                 valence = 2
        #             elif e == 1:
        #                 valence = 1

        name = primitive_key.HYDROGEN
        codec = self.primitive_codecs[name]
        arr = node.primitives[name]
        h = codec.decode_int(arr.on_first())

        if sym == 16:
            # breakpoint()
            if len(edges) > 2:
                return False
            if len(edges) == 4 and valence == 2:
                return False
            elif len(edges) <= 2 and valence > 2:
                return False
            elif valence == 2 and len(edges) > 2:
                return False

        if sym == 15:
            # breakpoint()
            if len(edges) < 3:
                return False
            if len(edges) == 4 and valence < 3:
                return False
            elif len(edges) <= 3 and valence > 3:
                return False

        edges = dict(edges)
        # has_double = False
        # for eidx, e in edges.items():
        #     arr = e.primitives[primitive_key.BOND_ORDER]
        #     if arr[2]:
        #         has_double = True
        # for eidx, e in edges.items():
        #     arr = e.primitives[primitive_key.BOND_ORDER]
        #     if arr[1] or h > 0:
        #         has_single = True
        # for eidx, e in edges.items():
        #     arr = e.primitives[primitive_key.BOND_ORDER]
        #     if (arr[6] or arr[7]):
        #         if h == 2:
        #             return False
        #         has_stereo = True

        # if has_stereo and not has_double:
        #     return False
        # elif has_double and has_stereo and not has_single:
        #     return False
        # elif has_double and has_single and not has_stereo:
        #     return False

        for eidx, e in edges.items():

            arr = e.primitives[primitive_key.BOND_ORDER]
            if arr[1]:
                bo += 1
                # if eidx in new_edges:
                # has_single += 1

                # check to see if other bond is set to single already
                # other = eidx[0] if eidx[0] == idx else eidx[1]
                # for otheredge, e2 in edge_set:
                #     if otheredge == eidx:
                #         continue
                #     if e2.primitives['bond_order'][1]:
                #         valid = False

            elif arr[2]:
                bo += 2
                # if eidx not in new_edges:

                # other = eidx[0] if eidx[1] == idx else eidx[0]
                # for otheredge in [edge for edge in edges if other in edge]:
                #     if edges[otheredge].primitives['bond_order'][6]:
                #         if has_stereo:
                #             valid = False
                #             continue
                #         has_stereo = True
                #     if edges[otheredge].primitives['bond_order'][7]:
                #         if has_stereo:
                #             valid = False
                #             # continue
                #         has_stereo = True
            elif arr[3]:
                bo += 3
                if sym == 16 and len(edges) == 1:
                    # no sulfur triple (single) bonds
                    return False
            elif arr[4]:
                bo += 4
            elif arr[5]:
                bo += 1.5
                naro += 1
            elif arr[6]:
                # if not frag.primitives['element'][6]:
                #     continue
                bo += 1
                # if eidx in new_edges:
                if has_stereo:
                    valid = False
                    # continue
                # check to see if other bond is set to stereo already
                # other = eidx[0] if eidx[0] == idx else eidx[1]
                # for otheredge, e2 in edges.items():
                #     if otheredge == eidx:
                #         continue
                #     if e2.primitives[primitive_key.BOND_ORDER][6]:
                #         valid = False
                # has_stereo = True

                # # if eidx in new_edges:
                # other = eidx[0] if eidx[1] == idx else eidx[0]
                # for otheredge in [edge for edge in edges if other in edge]:
                #     if edges[otheredge].primitives['bond_order'][2]:
                #         has_double = True
                #     elif edges[otheredge].primitives['bond_order'][6]:
                #         valid = False
                #     elif edges[otheredge].primitives['bond_order'][7]:
                #         valid = False
                #         # continue
            elif arr[7]:
                bo += 1

                if has_stereo:
                    valid = False
                    # continue

                other = eidx[0] if eidx[0] == idx else eidx[1]
                for otheredge, e2 in edges.items():
                    if otheredge == eidx:
                        continue
                    if e2.primitives[primitive_key.BOND_ORDER][7]:
                        valid = False

                # for otheredge in [edge for edge in edges if other in edge]:
                #     if otheredge == eidx:
                #         continue
                #     if edges[otheredge].primitives['bond_order'][7]:
                #         valid = False
                # if eidx in new_edges:
                # other = eidx[0] if eidx[1] == idx else eidx[0]
                # for otheredge in [edge for edge in edges if other in edge]:
                #     if edges[otheredge].primitives['bond_order'][2]:
                #         has_double = True
                #     elif edges[otheredge].primitives['bond_order'][6]:
                #         valid = False
                #     elif edges[otheredge].primitives['bond_order'][7]:
                #         valid = False
                #         # continue
                has_stereo = True
            else:
                bo += 1

        if int(2 * (valence - h + q)) != int(2 * bo):
            return False
        return True

    def on_node_exit(self, idx):
        return None

    def on_edge_enter(self, idx):
        return None

    def on_edge(self, idx, smarts: bechem, *args, **kwargs):
        return True

    def on_edge_exit(self, idx):
        return None
