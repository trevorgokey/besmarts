"""
besmarts.core.chem

The standard definition of a binary encoded SMARTS component of a graph. All
nodes and edges have a bechem object embedded within. These define the
molecular information contained in the graph.
"""

from typing import Dict, List, Sequence, Generator, Tuple

from besmarts.core import arrays
from besmarts.core.primitives import primitive_key, primitive_key_set

class bechem:
    """
    Represents a group of atom or bond SMARTS primitives in a graph
    """

    __slots__ = ("primitives", "select")

    def __init__(self, primitives, select):
        """
        A base type for atom and bond primitives which use bitfields to
        perform set logic on underlying fields. Fields represent variables
        in chemical space, such as symbol or bond order. These are discrete
        and usually stem from discrete forms of chemical space, e.g. SMILES.

        Parameters
        ----------
        """
        self.primitives: Dict[primitive_key, arrays.bitvec] = primitives
        self.select: Sequence[primitive_key] = select

    def __getitem__(self, name: primitive_key) -> arrays.bitvec:
        return self.primitives[name]

    def __setitem__(self, name: primitive_key, array: arrays.bitvec) -> None:
        self.primitives[name] = array

    def __delitem__(self, k) -> None:
        del self.primitives[k]
        self.select = tuple((x for x in self.select if x != k))

    def __len__(self) -> int:
        return len(self.select)

    def __repr__(self) -> str:
        ret = " ".join(
            [
                "{}: {}".format(name, self.primitives[name])
                for name in self.select
            ]
        )
        return "bechem(" + ret + ")"

    def __iter__(self):
        yield from bechem_iter(self, skip_ones=False)

    def __hash__(self) -> int:
        fields = tuple(
            ((name, hash(self.primitives[name])) for name in self.select)
        )
        return hash(fields)

    def __contains__(self, o) -> bool:
        """
        Return whether the this instance is a super set of the argument by checking if
        all selected primitives are super sets
        """
        bechem_check_sane_compare(self, o)

        fields_a = [f for f in self.select if f in o.select]
        fields_o = [f for f in o.select if f in fields_a]
        ret = True
        for a_field, o_field in zip(fields_a, fields_o):
            a_vec: arrays.bitvec = self.primitives[a_field]
            o_vec: arrays.bitvec = o.primitives[o_field]
            if (o_vec - a_vec).any():
                ret = False
                break
        return ret

    def __and__(self, o) -> "bechem":
        """
        Bitwise AND of selected primitives
        """
        # bitwise and (intersection)
        return bechem_dispatch_op(self, o, arrays.bitvec_and)

    def __or__(self, o) -> "bechem":
        """
        Bitwise OR of selected primitives
        """
        # bitwise or (union)
        return bechem_dispatch_op(self, o, arrays.bitvec_or)

    def __xor__(self, o) -> "bechem":
        """
        Bitwise XOR of selected primitives
        """
        # bitwise xor
        return bechem_dispatch_op(self, o, arrays.bitvec_xor)

    def __invert__(self) -> "bechem":
        """
        Bitwise complement of selected primitives
        """
        # negation (not a)
        return bechem_dispatch_op(self, None, arrays.bitvec_not)

    def __add__(self, o) -> "bechem":
        """
        Bitwise OR of selected primitives
        """
        # a + b is union
        return self | o

    def __sub__(self, o) -> "bechem":
        """
        Bitwise difference of selected primitives
        """
        # a - b is a marginal, note that a - b != b - a
        ans = self ^ o
        return self & ans

    def __eq__(self, o) -> bool:
        """
        We can't use hashes here because hash(-1) == hash(-2), so ~0b0 == ~0b1
        """
        return all([
            o.primitives[name] == self.primitives[name]
            for name in self.select
        ])

    def __lt__(self, o) -> bool:
        a, b = bechem_reduce_longest(self, o)
        return a < b

    def __gt__(self, o) -> bool:
        a, b = bechem_reduce_longest(self, o)
        return a > b

    def __le__(self, o) -> bool:
        a, b = bechem_reduce_longest(self, o)
        return a <= b

    def __ge__(self, o) -> bool:
        a, b = bechem_reduce_longest(self, o)
        return a >= b

    def __ne__(self, o) -> bool:
        """
        We can't use hashes here because hash(-1) == hash(-2), so ~0b0 == ~0b1
        """
        return any([
            o.primitives[name] != self.primitives[name]
            for name in self.select
        ])

    def bits(self, maxbits=False) -> int:
        return bechem_bits(self, maxbits=maxbits)

    def bits_max(self) -> int:
        return bechem_bits_max(self)

    def any(self):
        return bechem_any(self)

    def all(self):
        return bechem_all(self)

    def clear(self):
        bechem_clear(self)

    def fill(self):
        bechem_fill(self)

    def is_valid(self):
        return bechem_is_valid(self)

    def is_null(self):
        return bechem_is_null(self)

    def copy(self):
        return bechem_copy(self)

    def to_fragments(self):
        return bechem_to_fragments(self)

    def disable(self, field):
        return bechem_disable(self, field)

    def enable(self, field):
        return bechem_enable(self, field)


def bechem_get(
    bc: bechem, name: primitive_key, default: arrays.bitvec
) -> arrays.bitvec:
    return bc.primitives.get(name, default)


def bechem_is_fragment(bc: bechem) -> bool:
    """
    Returns whether the SMARTS has exclusively single values per primitive.
    """

    for field in bc.select:
        arr = bc.primitives[field]
        if arr.bits() != 1:
            return False
    return True


def bechem_is_null(bc: bechem) -> bool:
    """
    Returns whether any primitive has no values set.

    Parameters
    ----------

    Returns
    -------
    bool
    """

    for name in bc.select:
        arr = bc.primitives[name]
        if arr.is_null():
            return True

    return False


def bechem_is_valid(bc: bechem) -> bool:
    """
    Returns whether all primitives have at least one value set, i.e. it can be
    encoded into valid SMARTS

    Parameters
    ----------

    Returns
    -------
    bool
    """
    return not bechem_is_null(bc)


def bechem_recurse_fields(
    bc: bechem, fields: dict, pos=None
) -> Generator[Dict, None, None]:

    """
    Generate all fragments given the set primitives
    """

    if len(fields) == 0:
        return

    if pos is None:
        pos = []

    fvals = tuple(fields.values())
    if len(pos) == len(fields):

        ret = [field[i] for field, i in zip(fvals, pos)]


        if all(ret): 
            yield pos.copy()
        return

    l_pos = len(pos)
    pos.append(0)
    # print("field", l+1, "/", len(fields), "has length", len(fields[l]))

    # if there is a maxbits set, do not iterate beyond it
    field_i = fvals[l_pos]

    bits = len(field_i)

    if arrays.bitvec_is_inverted(field_i):
        bits = field_i.maxbits

    for i in range(bits):
        pos[-1] = i
        # print("querying", i, "of pos", len(pos), pos)
        yield from bechem_recurse_fields(bc, fields, pos=pos.copy())


def bechem_clear(bc: bechem) -> None:
    """
    Clear all values in selected primitives
    """
    for name in bc.select:
        bv = bc.primitives[name]
        bv.clear()


def bechem_fill(bc: bechem) -> None:
    """
    Fill all selected primitives
    """

    for name in bc.select:
        arr = bc.primitives[name]
        arr[:] = True


def bechem_all(bc: bechem) -> bool:
    """
    Return whether all primitives are full
    """

    for name in bc.select:
        arr = bc.primitives[name]
        if not arr.all():
            return False

    return True


def bechem_any(bc: bechem) -> bool:
    """
    Return whether any primitive is not null
    """

    for name in bc.select:
        arr = bc.primitives[name]
        if arr.any():
            return True

    return False


def bechem_bits(bc: bechem, maxbits=False, return_all=False) -> int:
    """
    Return the number of bits across all selected primitives
    """

    bits = []
    for i, name in enumerate(bc.select):
        arr: arrays.bitvec = bc.primitives[name]
        bits.append(arr.bits(maxbits=maxbits))
    n_bits = (sum(bits), bits)

    if return_all:
        return n_bits
    else:
        return n_bits[0]

def bechem_bits_max(bc: bechem, return_all=False) -> int:
    """
    Return the number of bits across all selected primitives
    """

    bits = []
    for i, name in enumerate(bc.select):
        arr: arrays.bitvec = bc.primitives[name]
        bits.append(arr.maxbits)
    n_bits = (sum(bits), bits)

    if return_all:
        return n_bits
    else:
        return n_bits[0]


def bechem_align_score(bc: bechem, o: bechem) -> int:
    """
    Return the number of overlapping bits
    """

    return bechem_bits(bc & o, maxbits=True)


def bechem_copy(bc: bechem) -> bechem:
    """
    Return a copy of bc.
    Returns
    -------
    cls: type(cls)
       The new chem_type object
    """

    primitives = bc.primitives.copy()
    for field in primitives:
        primitives[field] = arrays.bitvec_copy(bc.primitives[field])

    return bechem(primitives, tuple(bc.select))


def bechem_to_fragments(bc: bechem) -> Sequence[bechem]:
    """
    Return all fragments of that can be generated from the selected primitives
    """

    terms = []

    for x in bechem_recurse_fields(
        bc,
        {name: bc.primitives[name] for name in bc.select},
        pos=list(),
    ):
        prim = {}
        for k, v in zip(bc.select, x):
            prim[k] = arrays.bitvec(maxbits=bc.primitives[k].maxbits)
            prim[k][v] = True
        chem = bechem(prim, tuple(bc.select))
        terms.append(chem)

    return terms


def bechem_enable(bc: bechem, field: primitive_key) -> bool:
    """
    Select a currently disabled primitive

    Parameters
    ----------
    field : primitive_key
    The primitive to enable

    Returns
    -------
    bool
        Whether the primitive was successfully enabled and the selection has been
        modified
    """

    if field not in bc.primitives:
        status = False
    elif field not in bc.select:
        status = True
        bc.select = tuple(list(bc.select) + [field])
    else:
        status = False

    return status


def bechem_disable(bc: bechem, field: primitive_key) -> bool:
    """
    Deselect a currently enabled primitive

    Parameters
    ----------
    field : primitive_key
    The primitive to disable

    Returns
    -------
    bool
        Whether the operation was successful
    """

    if field in bc.select:
        bc.select = tuple((x for x in bc.select if x != field))

    status = True
    return status


def bechem_iter(
    bc: bechem, skip_ones=False, iter_inverse=False, primitives=None
) -> Generator[bechem, None, None]:
    """
    Generate a copy for every bit across all selected primitives
    """

    blank = bechem_copy(bc)
    bechem_clear(blank)

    if not primitives:
        primitives = bc.select

    # print(f"iterating primitives: {primitives}")
    for field in bc.select:

        if field not in primitives:
            continue

        bv: arrays.bitvec = bc.primitives[field]
        if skip_ones and arrays.bitvec_bits(bv) == 1:
            continue
        cursor = blank.primitives[field]
        for bit in bv:
            cursor += bit
            blank.primitives[field] = cursor
            yield bechem_copy(blank)
            if iter_inverse:
                cursor.clear()
                cursor += ~bit
                blank.primitives[field] = cursor
                yield bechem_copy(blank)
            cursor.clear()
            blank.primitives[field] = cursor


def bechem_reduce(bc: bechem) -> int:
    """
    Get the sum/number of bits across all selected primitives
    """
    return bechem_bits(bc, maxbits=True)


def bechem_reduce_longest(bc: bechem, o: bechem) -> Tuple[int, int]:
    """
    Get the number of bits using another instance to pad if necessary
    """

    suma = 0
    sumb = 0
    for field in bc.select:
        a: arrays.bitvec = bc.primitives[field]
        b: arrays.bitvec = o.primitives[field]
        x, y = bechem_reduce_longest(a, b)
        suma += x
        sumb += y

    return suma, sumb


def bechem_check_sane_compare(bc: bechem, o: bechem) -> bool:
    """
    Check if it is sane to perform an operation against another object

    Raises
    ------
    Exception
    """
    if type(bc) != type(o):
        raise Exception("chem operations must use same type")
    return True


def bechem_dispatch_op(bc: bechem, o: bechem, fn) -> bechem:
    """
    Broadcast an operation across primitives.
    """

    if o:
        bechem_check_sane_compare(bc, o)

    ret = bechem_copy(bc)
    for name_a in bc.select:

        field_a = bc.primitives[name_a]

        if o is None:
            ret.primitives[name_a] = fn(field_a)
        else:
            field_b = o.primitives[name_a]
            ret.primitives[name_a] = fn(field_a, field_b)

    return ret


def bechem_add(a: bechem, b: bechem) -> bechem:
    return a + b


def bechem_iadd(a: bechem, b: bechem) -> bechem:
    a += b
    return a


def bechem_subtract(a: bechem, b: bechem) -> bechem:
    return a - b


def bechem_subtract_conditional(a: bechem, b: bechem) -> bechem:
    c = a - b
    for p in c.primitives:
        if c.primitives[p].is_null():
            c.primitives[p] = a.primitives[p].copy()
    return c


def bechem_isubtract(a: bechem, b: bechem):
    a -= b
    return a


def bechem_xor(a: bechem, b: bechem):
    return a ^ b


def bechem_ixor(a: bechem, b: bechem):
    a ^= b
    return a


def bechem_and(a: bechem, b: bechem):
    return a & b


def bechem_and_conditional(a: bechem, b: bechem):
    c = a & b

    for p in c.primitives:
        if c.primitives[p].is_null():
            c.primitives[p] = a.primitives[p].copy()

    return c


def bechem_iand(a: bechem, b: bechem):
    a &= b
    return a


def bechem_or(a: bechem, b: bechem):
    return a | b


def bechem_ior(a: bechem, b: bechem):
    a |= b
    return a


def bechem_identity(a: bechem, b: bechem = None):
    return a


def bechem_neg(a: bechem, b: bechem = None):
    return ~a
