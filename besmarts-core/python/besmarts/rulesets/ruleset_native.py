
"""
besmarts.rulesets.ruleset_native

Defines the rules for determining whether a given graph is "valid" where the
current state is whether the graph is a valid molecule/SMILES
"""

from besmarts.core.primitives import primitive_key
from besmarts.core.chem import bechem
from besmarts.core.graph_visitors import graph_visitor
from besmarts.core.graphs import graph


class visitor_ruleset(graph_visitor):
    def __init__(self, primitive_codecs):
        self.primitive_codecs = primitive_codecs
        self.beg: graph = None

    def on_start(self, g: graph):
        self.beg = g
        return None

    def on_stop(self, g: graph, result):
        self.beg: graph = None
        return None

    def on_graph(self, nodes, edges):
        valid = True
        codec = self.primitive_codecs[primitive_key.FORMAL_CHARGE]

        if abs(codec.count_charge_smiles(nodes)) != 0:
            valid = False

        return valid

    def on_nodes(self, nodes):

        codec = self.primitive_codecs[primitive_key.ELEMENT]

        if (
            len(self.beg.nodes) > 3
            and (len(nodes) - codec.count_carbon_smiles(nodes))
            / len(self.beg.nodes)
            > 0.5
        ):
            return False
        return True

    def on_node(self, idx, besmarts: bechem, *args, **kwargs):

        valid = True
        name = primitive_key.FORMAL_CHARGE
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        q = codec.decode_int(arr.on()[0])

        name = primitive_key.ELEMENT
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        e = codec.decode_int(arr.on()[0])

        name = primitive_key.VALENCE
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        valence = codec.decode_int(arr.on()[0])

        name = primitive_key.RING_SMALLEST
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        r = codec.decode_int(arr.on()[0])

        name = primitive_key.CONNECTIVITY_RING
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        x = codec.decode_int(arr.on()[0])

        name = primitive_key.AROMATIC
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        a = codec.decode_int(arr.on()[0])

        name = primitive_key.HYDROGEN
        codec = self.primitive_codecs[name]
        arr = besmarts.primitives[name]
        h = codec.decode_int(arr.on()[0])
        
        if False:
            pass
        # elif a == 1 and (x == 0 or r < 2 or e not in [6,7,8,15]):
        #     valid = False
        elif e in [1]:
            if r != 0:
                valid = False
            elif x != 0:
                valid = False
            elif h > 0:
                valid = False
        elif e in [1, 9]:
            if q != 0:
                valid = False
            elif valence != 1:
                valid = False
        elif e == 6:
            if q != 0:
                valid = False
            elif valence != 4:
                valid = False

        elif e == 7:
            if q < 0:
                valid = False
            elif valence != 3:
                valid = False

        elif e == 8:
            if q > 0:
                valid = False
            elif valence != 2:
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
        q = codec.decode_int(arr.on()[0])

        name = primitive_key.VALENCE
        codec = self.primitive_codecs[name]
        arr = node.primitives[name]
        valence = codec.decode_int(arr.on()[0])

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
        h = codec.decode_int(arr.on()[0])

        for eidx, e in edges:

            arr = e.primitives[primitive_key.BOND_ORDER]
            if arr[1]:
                bo += 1
                # if eidx in new_edges:
                has_single += 1

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

                has_double = True
                # other = eidx[0] if eidx[1] == idx else eidx[0]
                # for otheredge in [edge for edge in edges if other in edge]:
                #     if edges[otheredge].primitives['bond_order'][6]:
                #         if has_stereo:
                #             valid = False
                #             # continue
                #         has_stereo = True
                #     if edges[otheredge].primitives['bond_order'][7]:
                #         if has_stereo:
                #             valid = False
                #             # continue
                #         has_stereo = True
            elif arr[3]:
                bo += 3
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
                for otheredge, e2 in edges:
                    if otheredge == eidx:
                        continue
                    if e2.primitives[primitive_key.BOND_ORDER][6]:
                        valid = False
                has_stereo = True

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
            elif arr[7]:
                bo += 1

                if has_stereo:
                    valid = False
                    # continue

                other = eidx[0] if eidx[0] == idx else eidx[1]
                for otheredge, e2 in edges:
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
