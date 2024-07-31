"""
besmarts.cluster.cluster_objective

Objective functions to use for optimizing cluster hierarchies.
"""

from typing import Sequence
import itertools
from besmarts.core import arrays
from besmarts.core.clusters import clustering_objective

def sum_obj(A):

    obj = []
    for a in A:
        if not obj:
            obj.extend(a)
        else:
            obj = arrays.array_add(obj, a)
    return obj

def format_array(a, s: str):
    return "".join(map(s.format, a))


class clustering_objective_mean_separation(clustering_objective):
    def __init__(self, split_separation=0.1, merge_separation=0.1):
        self.split_separation = split_separation
        self.merge_separation = merge_separation

    def is_discrete(self):
        return False

    def split(
        self, A: Sequence[Sequence[float]], B: Sequence[Sequence[float]], overlap=0.0
    ) -> float:

        # A: Sequence[float] = list(itertools.chain.from_iterable(A))
        # B: Sequence[float] = list(itertools.chain.from_iterable(B))
        A: Sequence[float] = list(A)
        B: Sequence[float] = list(B)
        if len(A) == 0 or len(B) == 0:
            return overlap

        abar = arrays.array_scale(sum_obj(A), 1/len(A))
        bbar = arrays.array_scale(sum_obj(B), 1/len(B))
        d = max(map(abs, arrays.array_difference(abar, bbar)))

        if d < self.split_separation:
            return d
        else:
            return -d

    def merge(
        self, A: Sequence[Sequence[float]], B: Sequence[Sequence[float]], overlap=0.0
    ) -> float:

        # A: Sequence[float] = list(itertools.chain.from_iterable(A))
        # B: Sequence[float] = list(itertools.chain.from_iterable(B))
        A: Sequence[float] = list(A)
        B: Sequence[float] = list(B)
        if len(A) == 0:
            abar = [0.0]
        else:
            abar = arrays.array_scale(sum_obj(A), 1/len(A))
            # abar = sum(A) / len(A)

        if len(B) == 0:
            bbar = [0.0]
        else:
            bbar = arrays.array_scale(sum_obj(B), 1/len(B))
            # bbar = sum(B) / len(B)

        d = max(map(abs, arrays.array_difference(abar, bbar)))
        # d = abs(abar - bbar)

        if d < self.merge_separation:
            return -d
        else:
            return d

    def report(self, A: Sequence[Sequence[float]]) -> str:

        # A: Sequence[float] = list(itertools.chain.from_iterable(A))
        A: Sequence[float] = list(A)
        if len(A) == 0:
            abar = [0]
            avar = [0]
            amin = [0]
            amax = [0]
        else:
            M = len(A[0])
            abar = arrays.array_scale(sum_obj(A), 1/len(A))
            avar = []
            for i in range(M):
                x = arrays.array_translate([a[i] for a in A], -abar[i])
                x = arrays.array_inner_product(x, x) / len(A)
                avar.append(x**.5)
            # avar = [arrays.array_sum(arrays.array_translate([a[i] for a in A], -abar[i]/len(A))) for i in range(M)]
            # abar = sum(A) / len(A)
            # avar = sum([(x-abar)**2 for x in A])/len(A) 
            amin = [min([a[i] for a in A]) for i in range(M)]
            amax = [max([a[i] for a in A]) for i in range(M)]
            # amax = max(A)
        s = "{:9.4f}"
        return (
            f" Mean= {format_array(abar, s)}"
            + f" Std= {format_array(avar, s)}"
            + f" N= {len(A):6d}"
            + f" Min= {format_array(amin, s)}"
            + f" Max= {format_array(amax, s)}"
        )

    def single(self, A: Sequence[Sequence[float]], overlap=0.0) -> float:
        A: Sequence[float] = list(itertools.chain.from_iterable(A))
        if len(A) == 0:
            return 0.0
        abar = sum(A) / len(A)
        return sum([x * x for x in A]) / len(A) - abar * abar

class clustering_objective_variance_separation(clustering_objective):
    def __init__(self, split_separation=1.0, merge_separation=1.0):
        self.split_separation = split_separation
        self.merge_separation = merge_separation

    def is_discrete(self):
        return False

    def split(
        self, A: Sequence[Sequence[float]], B: Sequence[Sequence[float]], overlap=0.0
    ) -> float:

        A: Sequence[float] = list(itertools.chain.from_iterable(A))
        B: Sequence[float] = list(itertools.chain.from_iterable(B))
        if len(A) == 0 or len(B) == 0:
            return overlap
        abar = sum(A) / len(A)
        avar = sum([(ai - abar)**2 for ai in A]) / len(A)
        bbar = sum(B) / len(B)
        bvar = sum([(bi - bbar)**2 for bi in B]) / len(B)
        d = abs(avar - bvar)

        if d < self.split_separation:
            return d
        else:
            return -d

    def merge(
        self, A: Sequence[Sequence[float]], B: Sequence[Sequence[float]], overlap=0.0
    ) -> float:

        A: Sequence[float] = list(itertools.chain.from_iterable(A))
        B: Sequence[float] = list(itertools.chain.from_iterable(B))
        if len(A) == 0:
            avar = 0.0
        else:
            abar = sum(A) / len(A)
            avar = sum([(ai - abar)**2 for ai in A]) / len(A)

        if len(B) == 0:
            bvar = 0.0
        else:
            bbar = sum(B) / len(B)
            bvar = sum([(bi - bbar)**2 for bi in B]) / len(B)

        d = abs(avar - bvar)

        if d < self.merge_separation:
            return -d
        else:
            return d

    def report(self, A: Sequence[Sequence[float]]) -> str:

        A: Sequence[float] = list(itertools.chain.from_iterable(A))
        if len(A) == 0:
            abar = 0
            avar = 0
            amin = 0
            amax = 0
        else:
            abar = sum(A) / len(A)
            avar = sum([(x-abar)**2 for x in A])/len(A) 
            amin = min(A)
            amax = max(A)
        return (
            f" Mean= {abar:9.4f}"
            + f" Var= {avar:9.4f}"
            + f" N= {len(A):6d}"
            + f" Min= {amin:9.4f}"
            + f" Max= {amax:9.4f}"
        )

    def single(self, A: Sequence[Sequence[float]], overlap=0.0) -> float:
        A: Sequence[float] = list(itertools.chain.from_iterable(A))
        if len(A) == 0:
            return 0.0
        abar = sum(A) / len(A)
        avar = sum([(ai - abar)**2 for ai in A]) / len(A)
        return avar

class clustering_objective_classification(clustering_objective):

    def __init__(self):
        pass

    def is_discrete(self):
        return True

    def split(self, A: Sequence[str], B: Sequence[str], overlap=0) -> float:
        a = set(A)
        b = set(B)
        i = len(a.intersection(b))
        # return (len(a)-overlap)**2 + (len(b)-overlap)**2 - len(a.union(b))**2
        return i - overlap - 1.0 # len(a.symmetric_difference(b))

    def merge(self, A: Sequence[str], B: Sequence[str], overlap=0) -> float:

        a = set(A)
        b = set(B)
        i = len(a.intersection(b))
        # print(f"\nA: {a} B: {b} i: {i}: obj: {overlap - i}")
        # return self.split(A, B, overlap=overlap)
        return -(i - overlap) # + len(a.symmetric_difference(b))

    def overlap(self, A: Sequence[str], B: Sequence[str]) -> float:

        a = set(A)
        b = set(B)
        i = len(a.intersection(b))
        return i

    def single(self, A: Sequence[str]) -> float:

        return (len(set(A)) - 1)**2

    def report(self, A: Sequence[str]) -> str:
        return str(set(A))
