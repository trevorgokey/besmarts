"""
besmarts.core.primitives

Standard definitions of SMARTS primitives.
"""

from typing import List
import re

from besmarts.core.arrays import bitvec
from besmarts.core import arrays


primitive_key = str
primitive_key_set = set(
    [
        "element",
        "hydrogen",
        "connectivity_total",
        "connectivity_ring",
        "ring_smallest",
        "aromatic",
        "chirality",
        "valence",
        "formal_charge",
        "bond_order",
        "bond_ring",
        "unit",
        "unit_index",
        "link_src",
        "link_dst",
        "variable",
    ]
)


class primitive_key:
    ELEMENT = "element"
    HYDROGEN = "hydrogen"
    CONNECTIVITY_TOTAL = "connectivity_total"
    CONNECTIVITY_RING = "connectivity_ring"
    RING_SMALLEST = "ring_smallest"
    AROMATIC = "aromatic"
    CHIRALITY = "chirality"
    VALENCE = "valence"
    FORMAL_CHARGE = "formal_charge"
    BOND_ORDER = "bond_order"
    BOND_RING = "bond_ring"
    UNIT = "unit"
    UNIT_INDEX = "unit_index"
    LINK_SRC = "link_src"
    LINK_DST = "link_dst"
    VARIABLE = "variable"

    def __init__(self, *args, **kwds):
        pass

    def __repr__(self):
        return f"{self.name}"


class primitive_codec:
    def __init__(self, implements, smarts_prefix):
        self.implements: str = implements
        self.smarts_prefix = smarts_prefix

    def decode_int(self, x: int) -> int:
        """
        Transform a decoded integer into its value for SMARTS/SMILES. The
        default behavior is a passthrough and works for most primitives, such
        as element and hydrogen. Special behavior for certain primitives e.g.
        formal charge, must be handled by subclassing.

        Parameters
        ----------
        x : int
            The integer to decode

        Returns
        -------
        int
            The decoded integer
        """

        return x

    def encode_int(self, x: int) -> int:
        """
        Transform an integer parsed from a SMILES or SMARTS representation to a
        valid value for storing in a bit array. The default behavior is a
        passthrough and works for most primitives, such as element and
        hydrogen. Special behavior for certain primitives e.g. formal charge,
        must be handled by subclassing.

        Parameters
        ----------
        x : int
            The integer to encode

        Returns
        -------
        int
            The encoded integer
        """
        return x

    def tokenize_smarts(self, obj) -> List[str]:
        descriptor = f"(!?){self.smarts_prefix}([0-9][0-9]*)"
        matches = list(re.finditer(descriptor, obj))
        return matches

    def decode_supersmarts(self, dtype, obj) -> bitvec:
        return self.decode_smarts(dtype, obj)

    def decode_smarts(self, dtype, obj) -> bitvec:
        """
        Parse and decode a SMARTS atom or bond for a primitive.

        Parameters
        ----------
        dtype : array
            The class used to store the decoded values

        obj : str
            The atom or bond string to parse and decode

        Returns
        -------
        bitvec
            The decoded SMARTS as an array of bits
        """

        matches = self.tokenize_smarts(obj)
        # descriptor = f"(!?){self.smarts_prefix}([0-9][0-9]*)"
        # matches = list(re.finditer(descriptor, obj))

        on = dtype()
        off = dtype()

        if "*" in obj:
            on[:] = True
            return on

        result = on

        if matches:
            for ret in matches:
                result = on
                if ret.group(1) == "!":
                    result = off
                encbit = self.decode_int(int(ret.group(2)))
                result[encbit] = True
        else:
            on[:]

        if on.any():
            result = on - off
        else:
            result = ~off

        return result

    def encode_supersmarts(self, dtype, obj) -> bitvec:
        return self.encode_smarts(dtype, obj)

    def encode_smarts(self, arr: bitvec) -> str:
        """
        Transform a primitive as a bit array to a SMARTS string.

        Parameters
        ----------
        arr : array
            The primitive in the form of a bit array


        Returns
        -------
        str
            The SMARTS string represented the array. Note that if the primitive
            has no values set, it is not a valid SMARTS string and "_" will be
            returned
        """
        if arrays.bitvec_all(arr):
            return ""
        if not arrays.bitvec_any(arr):
            return "_"
        on = list(arrays.bitvec_on(arr))
        off = list(arrays.bitvec_off(arr))

        if len(on) <= len(off):
            strs = [self.smarts_prefix + str(self.decode_int(i)) for i in on]
            return ",".join(strs)
        else:
            strs = [
                "!" + self.smarts_prefix + str(self.decode_int(i)) for i in off
            ]
            return "".join(strs)

    def encode_smiles(self, arr: bitvec) -> str:
        return ""

    def decode_smiles(self, dtype: bitvec, obj) -> bitvec:
        raise NotImplementedError()


element_tr = {
    "1": "H",
    "2": "He",
    "3": "Li",
    "4": "Be",
    "5": "B",
    "6": "C",
    "7": "N",
    "8": "O",
    "9": "F",
    "10": "Ne",
    "11": "Na",
    "12": "Mg",
    "13": "Al",
    "14": "Si",
    "15": "P",
    "16": "S",
    "17": "Cl",
    "18": "Ar",
    "19": "K",
    "20": "Ca",
    "21": "Sc",
    "22": "Ti",
    "23": "V",
    "24": "Cr",
    "25": "Mn",
    "26": "Fe",
    "27": "Co",
    "28": "Ni",
    "29": "Cu",
    "30": "Zn",
    "31": "Ga",
    "32": "Ge",
    "33": "As",
    "34": "Se",
    "35": "Br",
    "36": "Kr",
    "37": "Rb",
    "38": "Sr",
    "39": "Y",
    "40": "Zr",
    "41": "Nb",
    "42": "Mo",
    "43": "Tc",
    "44": "Ru",
    "45": "Rh",
    "46": "Pd",
    "47": "Ag",
    "48": "Cd",
    "49": "In",
    "50": "Sn",
    "51": "Sb",
    "52": "Te",
    "53": "I",
    "54": "Xe",
    "55": "Cs",
    "56": "Ba",
    "57": "La",
    "58": "Ce",
    "59": "Pr",
    "60": "Nd",
    "61": "Pm",
    "62": "Sm",
    "63": "Eu",
    "64": "Gd",
    "65": "Tb",
    "66": "Dy",
    "67": "Ho",
    "68": "Er",
    "69": "Tm",
    "70": "Yb",
    "71": "Lu",
    "72": "Hf",
    "73": "Ta",
    "74": "W",
    "75": "Re",
    "76": "Os",
    "77": "Ir",
    "78": "Pt",
    "79": "Au",
    "80": "Hg",
    "81": "Tl",
    "82": "Pb",
    "83": "Bi",
    "84": "Po",
    "85": "At",
    "86": "Rn",
    "87": "Fr",
    "88": "Ra",
    "89": "Ac",
    "90": "Th",
    "91": "Pa",
    "92": "U",
    "93": "Np",
    "94": "Pu",
    "95": "Am",
    "96": "Cm",
    "97": "Bk",
    "98": "Cf",
    "99": "Es",
    "100": "Fm",
    "101": "Md",
    "102": "No",
    "103": "Lr",
    "104": "Rf",
    "105": "Db",
    "106": "Sg",
    "107": "Bh",
    "108": "Hs",
    "109": "Mt",
    "110": "Ds",
    "111": "Rg",
    "112": "Cn",
    "113": "Nh",
    "114": "Fl",
    "115": "Mc",
    "116": "Lv",
    "117": "Ts",
    "118": "Og",
}
