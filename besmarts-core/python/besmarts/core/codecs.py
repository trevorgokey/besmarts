"""
besmarts.core.codecs

Definition of the graph codec, the class responsible for encoding and decoding
SMARTS and SMILES strings to a BESMARTS graph representation.
"""

from typing import Sequence, Dict, List, Tuple
import re

from besmarts.core.configs import smiles_perception_config
from besmarts.core.primitives import primitive_key, primitive_codec, element_tr

from besmarts.core import topology
from besmarts.core import chem
from besmarts.core import arrays
from besmarts.core import graphs
from besmarts.core import graph_visitors

from besmarts.core.arrays import bitvec as bitvec
from besmarts.core.arrays import array_dtype


class graph_codec:
    """A graph codec transforms SMARTS and/or SMILES between string and graph
    representations."""

    def __init__(
        self,
        smiles_config: smiles_perception_config,
        primitive_codecs: Dict[primitive_key, primitive_codec],
        dtype: bitvec,
        atom_primitives: Sequence[primitive_key],
        bond_primitives: Sequence[primitive_key],
    ):
        # how to treat the smiles input for decoding
        self.smiles_config = smiles_config

        self.primitive_codecs = primitive_codecs

        # dtype of the primitives
        self.array = dtype

        # selects the primitives to perceive
        self.atom_primitives: Sequence[primitive_key] = atom_primitives
        self.bond_primitives: Sequence[primitive_key] = bond_primitives

    def smiles_decode(self, smiles: str) -> graphs.graph:
        """
        Transform a SMILES string into a graph

        Parameters
        ----------
        smiles : str
            The SMILES string to decode

        Returns
        -------
        graph
            The graph representation of the SMILES string
        """

        raise NotImplementedError()

    def smarts_decode(self, smarts: str) -> graphs.graph:
        """
        Transform a SMARTS string into a subgraph or structure depending if there are
        mapped atoms present.

        Paramters
        ---------
        smarts : str
            The SMARTS string to decode

        Returns
        -------
        graphs.graph
            The graph representation of the SMILES string. If the SMARTS has mapped
            atoms, then a substructure will be returned instead.
        """
        raise NotImplementedError()

    def smiles_encode(self, g: graphs.graph) -> str:
        """
        Transform a graph into a SMILES string. The graph must be a fragment, i.e. all
        primitives have exactly one value set.

        Parameters
        ----------
        g : graphs.graph
            The graph to encode

        Returns
        -------
        str
            The SMILES representation of the graph
        """

        codecs = {
            k: v
            for k, v in self.primitive_codecs.items()
            if k in self.atom_primitives or self.bond_primitives
        }
        visitor = graph_visitors.smiles_visitor(codecs)

        primary = None
        tag = False
        if hasattr(g, "select"):
            primary = g.select
            tag = True

        smiles = graph_visitors.enter_graph(
            visitor, g, primary=primary, tag=tag
        )

        return smiles

    def smarts_encode(self, g: graphs.graph) -> str:
        """
        Transform a graph into a SMARTS string.

        Parameters
        ----------
        g : graphs.graph
            The graph to encode

        Returns
        -------
        str
            The SMARTS representation of the graph
        """

        codecs = {
            k: v
            for k, v in self.primitive_codecs.items()
            if k in self.atom_primitives or k in self.bond_primitives
        }
        visitor = graph_visitors.smarts_visitor(codecs)

        primary = None
        tag = None

        h = g

        if hasattr(g, "topology"):
            tag = True
            primary = [g.select[i] for i in g.topology.primary]
        elif hasattr(g, "select"):
            # h = graphs.subgraph_to_graph(g)
            h = g
            tag = True
            primary = tuple(h.select)

        smiles: str = graph_visitors.enter_graph(visitor, h, primary, tag=tag)

        return smiles


class intvec_codec:
    """
    A codec that transforms graphs to and from a vector form, allowing for
    compact serialization.
    """

    def __init__(self, primitive_codecs, atom_primitives, bond_primitives):
        self.primitive_codecs = primitive_codecs

        # selects the primitives to perceive
        self.atom_primitives: Sequence[primitive_key] = atom_primitives
        self.bond_primitives: Sequence[primitive_key] = bond_primitives

    def graph_encode(self, g: graphs.graph) -> arrays.intvec:
        return graphs.graph_to_intvec(
            g, self.atom_primitives, self.bond_primitives
        )

    def subgraph_encode(self, g: graphs.subgraph) -> arrays.intvec:
        return graphs.subgraph_to_intvec(
            g, self.atom_primitives, self.bond_primitives
        )

    def structure_encode(self, g: graphs.structure) -> arrays.intvec:
        return graphs.structure_to_intvec(
            g, self.atom_primitives, self.bond_primitives
        )

    def graph_decode(self, intvec: arrays.intvec) -> graphs.graph:
        return intvec_codec_graph_decode(
            intvec, self.atom_primitives, self.bond_primitives
        )

    def subgraph_decode(self, intvec: arrays.intvec) -> graphs.subgraph:
        return intvec_codec_subgraph_decode(
            intvec, self.atom_primitives, self.bond_primitives
        )

    def structure_decode(self, intvec: arrays.intvec) -> graphs.structure:
        return intvec_codec_structure_decode(
            intvec, self.atom_primitives, self.bond_primitives
        )


class primitive_codec_formal_charge(primitive_codec):
    def __init__(self):
        super().__init__("formal_charge", "")

    def decode_int(self, x: int) -> int:
        """
        Transform an encoded integer. This method overrides the super class method by
        mapping non-negative integers to all integers.

        Parameters
        ----------
        x : int
            The integer to decode

        Returns
        -------
        int
            The decoded integer
        """

        if x % 2:
            return -(x // 2 + 1)
        else:
            return x // 2

    def encode_int(self, x: int) -> int:
        """
        Transform a raw integer to an encoded form. This method overrides the super
        class method by handling the domain of charges and maps them to the set of
        non-negative integers.

        Parameters
        ----------
        x : int
            The integer to encode

        Returns
        -------
        int
            The encoded integer
        """

        if x < 0:
            return -x * 2 - 1
        else:
            return x * 2

    def tokenize_smarts(self, smarts: str):
        matches = list(re.finditer(r"(!?)(\+|-)([0-9]*)", smarts))
        return matches

    def decode_smarts(self, dtype, smarts: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        on = dtype()
        off = dtype()

        on.maxbits = 8
        off.maxbits = 8

        if "*" in smarts:
            on[:] = True
            return on

        result = on

        matches = self.tokenize_smarts(smarts)

        # rdkit will turn -1 to -, +1 to +
        val = None
        if matches:
            for ret in matches:
                if ret.group(1) == "!":
                    result = off

                if not ret.group(3):
                    # just a + or - means 1
                    val = 1
                else:
                    val = int(ret.group(3))

                if ret.group(2)[0] == "-":
                    val = -val

                encbit = self.encode_int(val)
                result[encbit] = True

        else:
            matches = re.finditer(r"(!?)(\+\+*|--*)", smarts)
            if matches:
                for ret in matches:
                    if ret.group(1) == "!":
                        result = off
                    val = len(ret.group(2))

                    if ret.group(2)[0] == "-":
                        val = -val

                    encbit = self.encode_int(val)
                    result[encbit] = True

        if on.any():
            result = on - off
        else:
            result = ~off

        return result

    def encode_smarts(self, arr: bitvec) -> str:
        """
        Transform the array to a SMARTS string. This method overrides the super class
        method to handle specifics for this primitive.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMARTS string. Note that an invalid SMARTS string will be returned as
            "_" if the array is null
        """

        if arr.all():
            return ""
        if not arr.any():
            return "_"

        on = list(arr.on())
        off = list(arr.off())

        ret = []
        if len(on) < len(off):
            for i in on:
                i = self.decode_int(i)
                if i == 0:
                    ret.append("+0")
                elif i == 1:
                    ret.append("+")
                elif i == -1:
                    ret.append("-")
                elif i > 1:
                    ret.append("+" + str(i))
                elif i < -1:
                    ret.append(str(i))
            ret = ",".join(ret)

        else:
            for i in off:
                i = self.decode_int(i)
                if i == 0:
                    ret.append("!+0")
                elif i == 1:
                    ret.append("!+")
                elif i == -1:
                    ret.append("!-")
                elif i > 1:
                    ret.append("!+" + str(i))
                elif i < -1:
                    ret.append("!" + str(i))
            ret = "".join(ret)
        return ret

    def encode_smiles(self, arr: bitvec) -> str:
        """
        Transform the array to a SMILES string. This method overrides the super class
        method to handle specifics for this primitive, as most primitives are not used
        in SMILES.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMILES string.
        """
        # assert arr.bits() == 1
        return self.encode_smarts(arr)

    def count_charge_smiles(self, nodes):
        """
        Compute the total charge across a set of nodes. Assumes fragments by  by only
        using the first value set.

        Parameters
        ----------
        nodes : Sequence[bechem]
            The nodes to calculate the total charge from

        Returns
        -------
        int
            The charge of the nodes

        """

        count = 0
        for node in nodes.values():
            arr = node.primitives[self.implements]
            q = self.decode_int(arr.on_first())
            count += q
        return count


class primitive_codec_chirality(primitive_codec):
    def __init__(self):
        super().__init__("chirality", "")

    def encode_smarts(self, arr: bitvec) -> str:
        """
        Transform the array to a SMARTS string. This method overrides the super class
        method to handle specifics for this primitive.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMARTS string. Note that an invalid SMARTS string will be returned as
            "_" if the the array is null
        """

        if arr.all():
            return ""
        if not arr.any():
            return "_"
        smarts = []
        for i in arr.on():
            i = self.decode_int(i)
            if i == 1:
                smarts.append("@@")
            elif i == 2:
                smarts.append("@")
        # print("CHIRALITY:", [i for i in arr.on()])
        if len(smarts) == 0 and arr[0]:
            # rdkit does not like this
            # return "!@!@@"
            return ""
        elif len(smarts) == 2:
            return "@,@@"
        elif smarts:
            return smarts[0]
        else:
            return ""

    def encode_smiles(self, arr: bitvec) -> str:
        """
        Transform the array to a SMILES string. This method overrides the super class
        method to handle specifics for this primitive, as most primitives are not used
        in SMILES.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMILES string.
        """
        assert arr.bits() == 1
        return self.encode_smarts(arr)

    def tokenize_smarts(self, smarts) -> List[str]:
        descriptor = r"(!?)(@?)"
        matches = list(re.finditer(descriptor, smarts))
        return matches

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        descriptor = r"(!?)(@?)"
        matches = list(re.finditer(descriptor, obj))

        on: bitvec = dtype()
        off: bitvec = dtype()

        on.maxbits = 3
        off.maxbits = 3

        result = on

        if matches:
            for ret in matches:
                result = on
                if ret.group(1) == "!":
                    result = off
                if ret.group(2) == "@@":
                    encbit = self.encode_int(0)
                    result[encbit] = True
                elif ret.group(2) == "@":
                    encbit = self.encode_int(1)
                    result[encbit] = True
                else:
                    result[:] = True
        elif "*" in obj:
            on[:] = True
            return on

        if on.any():
            result = on - off
        else:
            result = ~off

        return result


class primitive_codec_element(primitive_codec):
    def __init__(self):
        super().__init__("element", "#")

    def encode_smiles(self, arr: bitvec) -> str:
        """
        Transform the array to a SMILES string. This method overrides the super class
        method to handle specifics for this primitive, as most primitives are not used
        in SMILES.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMILES string.
        """
        # assert arr.bits() == 1
        i = arr.on_first()
        return SMILES_ELEMENT(self.decode_int(i))

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """
        arr = super().decode_smarts(dtype, obj)
        arr.maxbits = 130
        return arr

    def count_carbon_smiles(self, nodes):
        count = 0
        for node in nodes.values():
            arr = node.primitives[self.implements]
            count += self.decode_int(arr.on_first()) == 6
        return count

    def count_element_smiles(self, nodes, n):
        count = 0
        for atom in nodes.values():
            arr = atom.primitives[self.implements]
            count += self.decode_int(arr.on_first()) == n
        return count



class primitive_codec_valence(primitive_codec):
    def __init__(self):
        super().__init__("valence", "v")

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """
        arr = super().decode_smarts(dtype, obj)
        arr.maxbits = 8
        return arr


class primitive_codec_hydrogen(primitive_codec):
    def __init__(self):
        super().__init__("hydrogen", "H")

    def encode_smiles(self, arr: bitvec) -> str:
        """
        Transform the array to a SMILES string. This method overrides the super class
        method to handle specifics for this primitive, as most primitives are not used
        in SMILES.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMILES string.
        """
        # assert arr.bits() == 1
        i = self.decode_int(arr.on_first())
        if i > 0:
            return "H"
        return ""

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """
        arr = super().decode_smarts(dtype, obj)
        arr.maxbits = 5
        return arr

    def count_hydrogen_smiles(self, nodes):
        count = 0
        for chem in nodes.values():
            arr = chem.primitives[self.implements]
            count += self.decode_int(arr.on_first())
        return count


class primitive_codec_connectivity_total(primitive_codec):
    def __init__(self):
        super().__init__("connectivity_total", "X")

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """
        arr = super().decode_smarts(dtype, obj)
        arr.maxbits = 5
        return arr


class primitive_codec_connectivity_ring(primitive_codec):
    def __init__(self):
        super().__init__("connectivity_ring", "x")

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """
        arr = super().decode_smarts(dtype, obj)
        arr.maxbits = 5
        return arr


class primitive_codec_aromatic(primitive_codec):
    def __init__(self):
        super().__init__("aromatic", "")

    def encode_smarts(self, arr: bitvec) -> str:
        """
        Transform the array to a SMARTS string. This method overrides the super class
        method to handle specifics for this primitive.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMARTS string. Note that an invalid SMARTS string will be returned as
            "_" if the array is null
        """

        if arr.all():
            return ""
        if not arr.any():
            return "_"

        zero = arr[self.decode_int(0)]
        one = arr[self.decode_int(1)]

        if zero and not one:
            return "A"
        elif one and not zero:
            return "a"
        else:
            return ""

    def tokenize_smarts(self, smarts: str) -> List[str]:
        descriptor = r"(!?)([Aa])"
        matches = list(re.finditer(descriptor, smarts))
        return matches

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        matches = self.tokenize_smarts(obj)

        on: bitvec = dtype()
        off: bitvec = dtype()

        on.maxbits = 2
        off.maxbits = 2

        result = on

        if matches:
            for ret in matches:
                result = on
                if ret.group(1) == "!":
                    result = off
                if ret.group(2) == "A":
                    encbit = self.encode_int(0)
                    result[encbit] = True
                elif ret.group(2) == "a":
                    encbit = self.encode_int(1)
                    result[encbit] = True
                else:
                    result[:] = True
        elif "*" in obj:
            on[:] = True
            return on

        if on.any():
            result = on - off
        else:
            result = ~off

        return result


class primitive_codec_bond_ring(primitive_codec):
    def __init__(self):
        super().__init__("bond_ring", "")

    def encode_smarts(self, arr: bitvec) -> str:
        """
        Transform the array to a SMARTS string. This method overrides the super class
        method to handle specifics for this primitive.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMARTS string. Note that an invalid SMARTS string will be returned as
            "_" if the array is null
        """

        if arr.all():
            return ""
        if not arr.any():
            return "_"

        zero = arr[self.decode_int(0)]
        one = arr[self.decode_int(1)]

        if zero and not one:
            return "!@"
        elif one and not zero:
            return "@"
        else:
            return ""

    def tokenize_smarts(self, smarts: str) -> List[str]:
        descriptor = r"(!?)(@?)"
        matches = list(re.finditer(descriptor, smarts))
        return matches

    def decode_smarts(self, dtype: array_dtype, smarts: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """
        aA = dtype()
        aA.maxbits = 2
        if "*" in smarts:
            aA[:] = True
            return aA

        if "@" in smarts:
            pat = re.search(r"(!?@)", smarts)
            if pat:
                n_g = len(pat.groups())
                aA[self.encode_int(1)] = (
                    True if n_g == 1 and pat.group(1) == "@" else False
                )
                aA[self.encode_int(0)] = (
                    True if n_g == 1 and pat.group(1) == "!@" else False
                )
            if aA.reduce() == 0:
                aA[:] = True

        else:
            aA[:] = True

        return aA


class primitive_codec_bond_order(primitive_codec):
    def __init__(self):
        super().__init__("bond_order", "")

    def tokenize_smarts(self, smarts) -> List[str]:
        descriptor = r"(!?)(~|-|=|:|#|\$|/|\\)"
        matches = list(re.finditer(descriptor, smarts))
        return matches

    def encode_smarts(self, arr: bitvec) -> str:
        """
        Transform the array to a SMARTS string. This method overrides the super class
        method to handle specifics for this primitive.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMARTS string. Note that an invalid SMARTS string will be returned as
            "_" if the array is null
        """

        if arr.all():
            return "~"
        if not arr.any():
            return "_"

        on = [self.decode_int(i) for i in arr.on()]
        off = [self.decode_int(i) for i in arr.off()]

        decoded = []
        neg = ""
        joiner = ","
        if len(on) > len(off):
            neg = "!"
            on = off
            joiner = ""
        for i in on:
            v = ""
            if i == 1:
                v = "-"
            elif i == 2:
                v = "="
            elif i == 3:
                v = "#"
            elif i == 4:
                v = "$"
            elif i == 5:
                v = ":"
            elif i == 6:
                v = "\\"
            elif i == 7:
                v = r"/"
            if v:
                decoded.append(neg + v)
        ret = joiner.join(decoded)
        if ret == "-,:" or ret == ":,-":
            return ""
        else:
            return ret

    def decode_smarts(self, dtype, smarts: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """
        retall = self.tokenize_smarts(smarts)

        onorder = dtype()
        offorder = dtype()

        if retall:
            for ret in retall:
                sym = ret.group(2)
                if sym == "&":
                    continue
                order = offorder
                if ret.group(1) != "!":
                    order = onorder

                if sym == "~":
                    order[:] = True
                elif sym == "-":
                    order[self.encode_int(1)] = True
                elif sym == "=":
                    order[self.encode_int(2)] = True
                elif sym == "#":
                    order[self.encode_int(3)] = True
                elif sym == "$":
                    order[self.encode_int(4)] = True
                elif sym == ":":
                    order[self.encode_int(5)] = True
                elif sym == "/":
                    order[self.encode_int(6)] = True
                elif sym == r"\\":
                    order[self.encode_int(7)] = True
                else:
                    # This should be for []@[] style bonds
                    order[self.encode_int(1)] = True
                    order[self.encode_int(5)] = True

            if onorder.any():
                order = onorder - offorder
            else:
                order = ~offorder
        else:
            order = onorder
            # not specified means -,:
            order[self.encode_int(1)] = True
            order[self.encode_int(5)] = True

        order.maxbits = 8
        return order

    def encode_smiles(self, arr: bitvec) -> str:
        """
        Transform the array to a SMILES string. This method overrides the super class
        method to handle specifics for this primitive, as most primitives are not used
        in SMILES.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMILES string.
        """

        # assert arr.bits() == 1

        sma = self.encode_smarts(arr)

        if sma in "-:":
            return ""
        else:
            return sma


class primitive_codec_ring_smallest(primitive_codec):
    def __init__(self):
        super().__init__("ring_smallest", "r")

    def decode_int(self, x: int) -> int:
        if x > 0:
            return x + 2
        else:
            return x

    def encode_int(self, x: int) -> int:
        if x > 0:
            return x - 2
        else:
            return x

    def encode_smarts(self, arr: bitvec) -> str:
        """
        Transform the array to a SMARTS string. This method overrides the super class
        method to handle specifics for this primitive.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMARTS string. Note that an invalid SMARTS string will be returned as
            "_" if the array is null
        """

        if arr.all():
            return ""
        if not arr.any():
            return "_"
        zero = self.decode_int(0)

        if not arr[zero]:
            if len(arr) > 1:
                if all(arr[1:]):
                    return "r"
            else:
                return "r"

        elif arr[zero]:
            if len(arr) > 1:
                if not any(arr[1:]):
                    return "!r"
            else:
                return "!r"

        on = list(arr.on())
        off = list(arr.off())

        if len(on) < len(off):
            ret = []
            if arr[zero]:
                ret = ["!r"]
            ret += [
                self.smarts_prefix + str(self.decode_int(i))
                for i in on
                if self.decode_int(i) != zero
            ]
            return ",".join(ret)
        else:
            ret = []
            if arr[zero]:
                # ret = ["r"]
                ret += [
                    "!" + self.smarts_prefix + str(self.decode_int(i))
                    for i in off
                    if i != self.decode_int(0)
                ]
            else:
                ret += [
                    "!" + self.smarts_prefix + str(self.decode_int(i))
                    for i in off
                ]
            # bandaid
            for i, x in enumerate(ret):
                if x == "!r0":
                    ret[i] = "r"

            return "".join(ret)

    def tokenize_smarts(self, smarts: str) -> List[str]:
        descriptor = r"(!?)r([3-9]?[0-9]*)"
        matches = list(re.finditer(descriptor, smarts))
        return matches

    def decode_smarts(self, dtype: array_dtype, smarts: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        on: bitvec = dtype()
        off: bitvec = dtype()

        on.maxbits = 7
        off.maxbits = 7

        if "*" in smarts:
            on[:] = True
            return on

        result = on
        matches = self.tokenize_smarts(smarts)

        if matches:
            for ret in matches:
                result = on
                if ret.group(1) == "!":
                    result = off
                if ret.group(2) != "":
                    encbit = self.encode_int(int(ret.group(2)))
                    result[encbit] = True
                else:
                    result[:] = True
                    result[self.encode_int(0)] = False

        if on.any():
            result = on - off
        else:
            result = ~off

        return result


class primitive_codec_variable(primitive_codec):
    """
    This will refer to a variable primitive, aka a custom primitive
    """

    def __init__(self, name):
        super().__init__(name, "=")

    def tokenize_smarts(self, smarts: str) -> List[str]:
        matches = list(re.finditer(r"(!?)=([A-Za-z]+)([0-9]*)", smarts))
        return matches

    def encode_smarts(self, arr: bitvec) -> str:
        return ""

    def decode_smarts(self, arr: bitvec) -> str:
        return ""

    def encode_supersmarts(self, arr: bitvec) -> str:
        return super().encode_smarts(arr)

    def decode_supersmarts(self, arr: bitvec) -> str:
        return super().decode_smarts(arr)


class primitive_codec_unit(primitive_codec):
    def __init__(self):
        super().__init__("unit", r"@")
        self.units: Dict[str, graphs.graph] = {}

    def tokenize_smarts(self, smarts: str):
        matches = list(re.finditer(r"(!?)$([A-Za-z]+[0-9]*)", smarts))
        return matches

    def encode_supersmarts(self, arr: bitvec) -> str:
        """
        Transform the array to a SMARTS string. This method overrides the super class
        method to handle specifics for this primitive. This will expand the graph primitive to SMARTS

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMARTS string. Note that an invalid SMARTS string will be returned as
            "_" if the the array is null
        """

        if arr.all():
            return ""
        if not arr.any():
            return "_"
        smarts = []
        for i in arr.on():
            i = self.decode_int(i)
            if i == 1:
                smarts.append("@@")
            elif i == 2:
                smarts.append("@")
        if len(smarts) == 0 and arr[0]:
            return ""
        elif len(smarts) == 2:
            return "@,@@"
        elif smarts:
            return smarts[0]
        else:
            return ""

    def decode_smiles(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        return dtype()

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        if type(obj) is not str:
            return bitvec()

        matches = self.tokenize_smarts(obj)

        on: bitvec = dtype()
        off: bitvec = dtype()

        on.maxbits = 3
        off.maxbits = 3

        result = on

        if matches:
            for ret in matches:
                result = on
                if ret.group(1) == "!":
                    result = off
                if ret.group(2) and ret.group(2)[0].isdigit():
                    encbit = self.encode_int(int(ret.group(2)))
                    result[encbit] = True
                elif ret.group(2) and ret.group(2)[0].isalpha():
                    num = self.units[ret.group(2)]
                    encbit = self.encode_int(num)
                    result[encbit] = True
                else:
                    result[:] = True
        elif "*" in obj:
            on[:] = True
            return on

        if on.any():
            result = on - off
        else:
            result = ~off

        return result



class primitive_codec_unit_index(primitive_codec):
    """
    This will refer or extract a subset of the graph primitive
    """

    def __init__(self):
        super().__init__("unit_index", ".")

    def tokenize_smarts(self, smarts: str):
        matches = list(re.finditer(r"(!?)\.([1-9][0-9]*)", smarts))
        return matches

    # def encode_smarts(self, arr: bitvec) -> str:
    #     return ""

    def encode_supersmarts(self, arr: bitvec) -> str:
        return ""

    def decode_supersmiles(self, obj) -> bitvec:
        return bitvec()

    def decode_smiles(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        return dtype()

class primitive_codec_link_src(primitive_codec):
    def __init__(self):
        super().__init__("link_src", "<")

    def tokenize_smarts(self, smarts: str):
        matches = list(re.finditer(r"(!?)<([1-9][0-9]*)", smarts))
        return matches

    def encode_supersmarts(self, arr: bitvec) -> str:
        return ""

    def decode_supersmarts(self, arr: bitvec) -> str:
        return ""

    def decode_smiles(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        return dtype()

class primitive_codec_link_dst(primitive_codec):
    def __init__(self):
        super().__init__("link_dst", ">")

    def tokenize_smarts(self, smarts: str):
        matches = list(re.finditer(r"(!?)>([1-9][0-9]*)", smarts))
        return matches

    def encode_supersmarts(self, arr: bitvec) -> str:
        return ""

    def decode_supersmarts(self, arr: bitvec) -> str:
        return ""

    def decode_smiles(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        return dtype()


class primitive_codec_variable(primitive_codec):
    def __init__(self):
        super().__init__("unit", "Q")
        self.units: Dict[str, graphs.graph] = {}

    def tokenize_smarts(self, smarts: str):
        matches = list(re.finditer(r"(!?)Q([0-9]+)=([^|\[\]]+|)", smarts))
        return matches

    def encode_smarts(self, arr: bitvec) -> str:
        """
        Transform the array to a SMARTS string. This method overrides the super class
        method to handle specifics for this primitive. This will expand the graph primitive to SMARTS

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMARTS string. Note that an invalid SMARTS string will be returned as
            "_" if the the array is null
        """

        if arr.all():
            return ""
        if not arr.any():
            return "_"
        smarts = []
        for i in arr.on():
            i = self.decode_int(i)
            if i == 1:
                smarts.append("@@")
            elif i == 2:
                smarts.append("@")
        if len(smarts) == 0 and arr[0]:
            return ""
        elif len(smarts) == 2:
            return "@,@@"
        elif smarts:
            return smarts[0]
        else:
            return ""

    def encode_smiles(self, arr: bitvec) -> str:
        """
        Transform the array to a SMILES string. This method overrides the super class
        method to handle specifics for this primitive, as most primitives are not used
        in SMILES.

        Parameters
        ----------
        arr : array
            The array to encode

        Returns
        -------
        str
            The SMILES string.
        """
        # assert arr.bits() == 1
        return self.encode_smarts(arr)

    def decode_smiles(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        return self.decode_smarts(dtype, obj)

    def decode_smarts(self, dtype, obj: str) -> bitvec:
        """
        Parse and transform the SMARTS atom to a binary form. This method overrides the
        super class method to handle the parsing this primitive.

        Parameters
        ----------
        dtype : array
            The class to use to build new arrays
        smarts : str
            The SMARTS string of the atom to parse

        Returns
        -------
        array
            The decoded SMARTS string.
        """

        if type(obj) is not str:
            return bitvec()

        matches = self.tokenize_smarts(obj)

        on: bitvec = dtype()
        off: bitvec = dtype()

        on.maxbits = 3
        off.maxbits = 3

        result = on

        if matches:
            for ret in matches:
                result = on
                if ret.group(1) == "!":
                    result = off
                if ret.group(2) and ret.group(2)[0].isdigit():
                    encbit = self.encode_int(int(ret.group(2)))
                    result[encbit] = True
                elif ret.group(2) and ret.group(2)[0].isalpha():
                    num = self.units[ret.group(2)]
                    encbit = self.encode_int(num)
                    result[encbit] = True
                else:
                    result[:] = True
        elif "*" in obj:
            on[:] = True
            return on

        if on.any():
            result = on - off
        else:
            result = ~off

        return result


def SMILES_ELEMENT(i):
    return element_tr[str(i)]


class graph_codec_chain(graph_codec):
    def add_graph_codec(codec):
        pass

    def smiles_decode(self, smiles):
        # TODO: Merge graphs
        for codec in self.codecs:
            g = codec.smiles_decode(smiles)

        return g

    def smarts_decode(self, smarts):
        return NotImplementedError()

    def smiles_encode(self, g: graphs.graph, primary=None, tag=None):
        codecs = {
            k: v
            for k, v in self.primitive_codecs.items()
            if k in self.atom_primitives or self.bond_primitives
        }
        visitor = graph_visitors.smiles_visitor(codecs)

        smiles = graph_visitors.enter_graph(visitor, g, primary, tag=tag)

        return smiles

    def smarts_encode(self, g: graphs.graph, primary=None, tag=None):
        codecs = {
            k: v
            for k, v in self.primitive_codecs.items()
            if k in self.atom_primitives or self.bond_primitives
        }
        visitor = graph_visitors.smarts_visitor(codecs)

        h = g
        if hasattr(g, "topology"):
            tag = True
            primary = [g.select[i] for i in g.topology.primary]
        if hasattr(g, "select"):
            h = graphs.subgraph_to_graph(g)

        smiles = graph_visitors.enter_graph(visitor, h, primary, tag=tag)

        return smiles


def primitive_codecs_get() -> Dict[primitive_key, primitive_codec]:
    codecs = {
        "element": primitive_codec_element(),
        "hydrogen": primitive_codec_hydrogen(),
        "connectivity_total": primitive_codec_connectivity_total(),
        "connectivity_ring": primitive_codec_connectivity_ring(),
        "ring_smallest": primitive_codec_ring_smallest(),
        "aromatic": primitive_codec_aromatic(),
        "chirality": primitive_codec_chirality(),
        "valence": primitive_codec_valence(),
        "formal_charge": primitive_codec_formal_charge(),
        "bond_order": primitive_codec_bond_order(),
        "bond_ring": primitive_codec_bond_ring(),
        "link_src": primitive_codec_link_src(),
        "link_dst": primitive_codec_link_dst(),
        "unit": primitive_codec_unit(),
        "unit_index": primitive_codec_unit_index(),
    }
    return codecs


def intvec_codec_graph_decode_auto(
    intvec: arrays.intvec, atom_primitives, bond_primitives
) -> Tuple[graphs.graph, Sequence[int], int]:
    vec = intvec.v
    n_nodes = vec[0]
    n_nprim = vec[1]
    n_edges = vec[2]
    n_eprim = vec[3]
    graph_t = vec[4]

    offset = 5
    idx = offset
    select = []
    nodes = {}
    s = 0
    e = 0
    for nidx in range(n_nodes):
        nid = vec[idx]
        if nid < 0:
            nid = -nid
            select.append(nid)
        s = idx + 1
        e = s + n_nprim
        node_mask = slice(s, e)
        nodes[nid] = chem.bechem(
            {
                name: arrays.bitvec(x)
                for name, x in zip(atom_primitives, vec[node_mask])
            },
            atom_primitives,
        )
        idx = e

    if select:
        reorder = [s for s in select] + [n for n in nodes if n not in select]
        nodes = {n: nodes[n] for n in select}

    edges = {}
    for eidx in range(n_edges):
        eid = graphs.edge((vec[e], vec[e + 1]))
        s = e + 2
        e = s + n_eprim
        edge_mask = slice(s, e)
        edges[eid] = chem.bechem(
            {
                name: arrays.bitvec(x)
                for name, x in zip(bond_primitives, vec[edge_mask])
            },
            bond_primitives,
        )

    return graphs.graph(nodes, edges), select, graph_t


def intvec_codec_graph_decode(
    intvec: arrays.intvec, atom_primitives, bond_primitives
) -> graphs.graph:
    g, _, _ = intvec_codec_graph_decode_auto(
        intvec, atom_primitives, bond_primitives
    )
    return g


def intvec_codec_subgraph_decode(
    intvec: arrays.intvec, atom_primitives, bond_primitives
) -> graphs.subgraph:
    g, select, _ = intvec_codec_graph_decode_auto(
        intvec, atom_primitives, bond_primitives
    )
    return graphs.graph_to_subgraph(g, select)


def intvec_codec_structure_decode(
    intvec: arrays.intvec, atom_primitives, bond_primitives
) -> graphs.structure:
    g, select, graph_t = intvec_codec_graph_decode_auto(
        intvec, atom_primitives, bond_primitives
    )
    topo = topology.topology_index[graph_t]
    return graphs.graph_to_structure(g, select, topo)

def smiles_decode_list_distributed(smiles: List[str], shm=None) -> graphs.graph:
    gcd: graph_codec = shm.gcd
    icd = intvec_codec(gcd.primitive_codecs, gcd.atom_primitives, gcd.bond_primitives)
    return [icd.graph_encode(gcd.smiles_decode(s)) for s in smiles]

def smiles_decode_distributed(smiles: str, shm=None) -> graphs.graph:
    gcd = shm.gcd
    return gcd.smiles_decode(smiles)

def smarts_decode_distributed(smarts: str, shm=None) -> graphs.graph:
    gcd = shm.gcd
    return gcd.smarts_decode(smarts)
