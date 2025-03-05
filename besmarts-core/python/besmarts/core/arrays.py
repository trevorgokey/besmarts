"""
besmarts.core.arrays

The core data types that store SMARTS primitives.
"""

from typing import List, Tuple
import array
import itertools

INF = -1

class bitvec:

    """
    Holds the bit vector of a single primitive as a single integer.
    """

    __slots__ = "v", "maxbits"

    def __init__(self, val=0, maxbits=64):
        self.v: int = val
        self.maxbits: int = maxbits

    def inv(self):
        return self.v < 0

    def __iter__(self):
        for x in bitvec_on(self):
            b = bitvec(maxbits=self.maxbits)
            b[x] = True
            yield b

    def __getitem__(self, i):
        l = min(len(self), self.maxbits)
        if isinstance(i, slice):
            # [:]
            if i.stop is None and i.start is None and i.step is None:
                return [((self.v >> x) & 1) for x in range(l)]

            start = 0 if i.start is None else i.start
            end = i.stop

            if end is None:
                end = self.maxbits

            if self.maxbits == INF:
                raise IndexError(
                    "Cannot supply an infinite length array "
                    "(input slice had no bounds and maxbits was inf (-1))"
                )

            step = 1 if i.step is None else i.step
            rev = step < 0
            step = abs(step)

            # [i:] and i > maxbits
            if start >= l:
                diff = (end - start + 1) // step
                pad = list([int(self.v < 0)] * diff)
                pad = pad[::-1]
                return pad

            # [i:j:k] and 2**j > v
            if end > l:
                diff = (end - l) // step
                expl = [(self.v >> x) & 1 for x in range(start, l, step)]
                pad = list([int(self.v < 0)] * diff)
                expl += pad
                if rev:
                    expl = expl[::-1]
                return expl

            # [i:j:k] and all within v
            expl = [(self.v >> x) & 1 for x in range(start, end, step)]
            if rev:
                expl = expl[::-1]
            return expl

        # [i]
        elif i >= l:
            return int(self.v < 0)
        else:
            return (self.v >> i) & 1

    def __setitem__(self, i, v: int):
        v: bool = v > 0

        if isinstance(i, slice):
            # [:] = v
            if i.stop is None and i.start is None and i.step is None:
                if v:
                    self.v = -1
                else:
                    self.v = 0
                return

            start = 0 if i.start is None else i.start
            end = int(max(self.maxbits if i.stop is None else i.stop, start))

            if end > self.maxbits:
                end = self.maxbits

            step = 1 if i.step is None else i.step

            # [x:] = v
            if i.stop is None:
                mask = sum(1 << i for i in range(start, end, step))
                if v:
                    self.v |= mask
                else:
                    mask = ~mask
                    self.v &= mask

            else:
                # [x:y:z] = v
                mask = sum(1 << i for i in range(start, end, step))
                if v:
                    self.v |= mask
                else:
                    mask = ~mask
                    self.v &= mask

        # [5] = True
        elif isinstance(i, int):
            mask = 1 << i
            if v:
                self.v |= mask
            else:
                mask = ~mask
                self.v &= mask

        # [[1,0,1,0,1]] = v
        elif isinstance(i, bitvec):
            on = bitvec_on(i)
            mask = sum(1 << x for x in on)
            if v:
                self.v |= mask
            else:
                mask = ~mask
                self.v &= mask

        # [[1,2,3]] = v
        elif hasattr(i, "__iter__"):
            mask = sum(1 << x for x in set(i))
            if v:
                self.v |= mask
            else:
                mask = ~mask
                self.v &= mask

        else:
            raise Exception(
                "Using datatype {} for setitem not supported".format(type(i))
            )

    def __len__(self) -> int:
        if self.v == 0:
            return 0

        l = len(bin(self.v))

        if self.v < 0:
            return l - 3

        return l - 2

    def __repr__(self) -> str:
        v = self.v
        neg = " "

        if v < 0:
            v = -self.v - 1
            neg = "~"

        return neg + "{:>s}".format(bin(v))

    def __hash__(self) -> int:
        return self.v

    def __and__(self, o) -> "bitvec":
        return bitvec_and(self, o)

    def __or__(self, o) -> "bitvec":
        return bitvec_or(self, o)

    def __xor__(self, o) -> "bitvec":
        return bitvec_xor(self, o)

    def __invert__(self) -> "bitvec":
        return bitvec_not(self)

    def __add__(self, o) -> "bitvec":
        return bitvec_or(self, o)

    def __sub__(self, o) -> "bitvec":
        return bitvec_subtract(self, o)

    def __eq__(self, o) -> bool:
        return bitvec_equal(self, o)

    def __ne__(self, o) -> bool:
        return not self == o

    def __contains__(self, o) -> "bool":
        return bitvec_is_null(bitvec_subtract(o, self))

    def any(self) -> bool:
        return bitvec_any(self)

    def all(self) -> bool:
        return bitvec_all(self)

    def is_null(self) -> bool:
        return bitvec_is_null(self)

    def on(self) -> List[int]:
        return bitvec_on(self)

    def on_first(self) -> int:
        return bitvec_on_first(self)

    def off(self) -> List[int]:
        return bitvec_off(self)

    def off_first(self) -> int:
        return bitvec_off_first(self)

    def bits(self, maxbits=False) -> int:
        return bitvec_bits(self, maxbits=maxbits)

    def reduce(self) -> int:
        return bitvec_sum(self)

    def clear(self):
        self.v = 0

    def copy(self):
        return bitvec_copy(self)


class intvec:
    __slots__ = "v"

    def __init__(self):
        self.v = array.array("q")


def bitvec_sum(bv: bitvec) -> int:
    return bv.v


def bitvec_is_inverted(bv: bitvec) -> bool:
    return bv.v < 0


def bitvec_bits(bv: bitvec, maxbits=False) -> int:
    inv = bitvec_is_inverted(bv)
    if inv:
        if maxbits:
            return bv.maxbits
        else:
            return INF

    return len(bitvec_on(bv))


def bitvec_on(bv: bitvec) -> List[int]:
    return [i for i in range(bv.maxbits) if (bv.v >> i) & 1]


def bitvec_on_first(bv: bitvec) -> List[int]:

    for i in range(bv.maxbits):
        if (bv.v >> i) & 1:
            return i
    return None


def bitvec_off(bv: bitvec) -> List[int]:
    return [i for i in range(bv.maxbits) if not (bv.v >> i) & 1]


def bitvec_off_first(bv: bitvec) -> List[int]:

    for i in range(bv.maxbits):
        if not (bv.v >> i) & 1:
            return i
    return None


def bitvec_explicit_flip(bv: bitvec) -> None:
    bv.v = ~bv.v


def bitvec_is_null(bv: bitvec) -> bool:
    return bv.v == 0


def bitvec_clear(bv: bitvec) -> None:
    bv.v = 0


def bitvec_copy(bv: bitvec) -> bitvec:
    return bitvec(bv.v, bv.maxbits)


def bitvec_all(bv: bitvec) -> bool:
    return bv.v == -1


def bitvec_any(bv: bitvec) -> bool:
    return bv.v != 0


def bitvec_reduce(bv: bitvec) -> int:
    return bv.v

def bitvec_not(a: bitvec) -> bitvec:
    return bitvec(~a.v, a.maxbits)


def bitvec_or(a: bitvec, b: bitvec) -> bitvec:
    return bitvec(a.v | b.v, max(a.maxbits, b.maxbits))


def bitvec_and(a: bitvec, b: bitvec) -> bitvec:
    return bitvec(a.v & b.v, min(a.maxbits, b.maxbits))


def bitvec_xor(a: bitvec, b: bitvec) -> bitvec:
    return bitvec(a.v ^ b.v, max(a.maxbits, b.maxbits))


def bitvec_subtract(a: bitvec, b: bitvec) -> bitvec:
    return bitvec(a.v & (a.v ^ b.v), a.maxbits)


def bitvec_equal(a: bitvec, b: bitvec) -> bool:
    return a.v == b.v


def bitvec_subset(a: bitvec, b: bitvec) -> bool:
    return a.v == (a.v & b.v)


def bitvec_superset(a: bitvec, b: bitvec) -> bool:
    return b.v == (a.v & b.v)


array_dtype = type(bitvec)


def batched(iterable, n: int):
    """
    https://docs.python.org/3/library/itertools.html#itertools-recipes
    Batch data into tuples of length n. The last batch may be shorter.

    Parameters
    ----------
    iterable
    n : int

    Returns
    -------
    Generator
    """
    # batched('ABCDEFG', 3) --> ABC DEF G
    it = iter(iterable)
    batch = tuple(itertools.islice(it, n))
    while batch:
        yield batch
        batch = tuple(itertools.islice(it, n))

def argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

def argmin(iterable):
    return min(enumerate(iterable), key=lambda x: x[1])[0]

def argsort(iterable):
    return [x[0] for x in sorted(enumerate(iterable), key=lambda x: x[1])]

def find_unsigned_typecode_min(N: int):
    code = None
    for c in "QLIHB":
        a = array.array(c)
        if N < 2 ** (8 * a.itemsize):
            # print(f"Array max value is {2**(8*a.itemsize) - 1} and itemsize is {a.itemsize}; len(A) is {len(A)}")
            code = c
    assert code, "Number is too large for this machine"
    return code

def flatten_list(l, times=1):
    if times == 0:
        return l
    if len(l) == 0:
        return l
    if times == -1:
        if isinstance(l[0], list):
            ll = [a for b in l if hasattr(b, "__iter__") for a in b]
            return flatten_list(ll, times)
        else:
            return l
    else:
        return flatten_list(
            [a for b in l if hasattr(b, "__iter__") for a in b], times - 1
        )

def array_scale(a, s):
    return [i*s for i in a]

def array_translate(a, s):
    return [i+s for i in a]

def array_sum(a):
    return sum(a)

def array_mean(a):
    return sum(a)/len(N)

def array_add(a, b):
    return [i+j for i,j in zip(a,b)]

def array_abs(a):
    return [abs(i) for i in a]

def array_difference(a, b):
    return [i-j for i,j in zip(a,b)]

def array_multiply(a, b):
    return [i*j for i,j in zip(a,b)]

def array_divide(a, b):
    return [i/j for i,j in zip(a,b)]

def array_inner_product(a, b):
    return sum((i*j for i,j in zip(a,b)))

def array_cross(a, b):
    return (
        (a[1]*b[2] - a[2]*b[1]),
        (a[2]*b[0] - a[0]*b[2]),
        (a[0]*b[1] - a[1]*b[0])
    )

def array_unit(a, b):
    """
    unit vector from a to b
    """
    return array_scale(array_difference(b, a), 1/array_distance(b, a))

def array_basis(a, b):
    """
    unit vector from a to b and its projection (magnitude)
    """
    r = array_distance(a, b)
    return array_scale(array_difference(b, a), 1/r), r

def array_magnitude(a) -> float:
    return sum([x*x for x in a])**.5


def array_distance(a,b) -> float:
    return sum([x*x for x in array_difference(b,a)])**.5

def array_round(a, b) -> List[float]:
    return [round(x, b) for x in a]

def measure_distance(xyz1, xyz2):

    result = [[array_distance(a,b)] for a, b in zip(xyz1, xyz2)]

    return result
