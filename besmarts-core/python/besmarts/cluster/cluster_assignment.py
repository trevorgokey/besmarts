"""
besmarts.cluster.cluster_assignment

Data structures for associating data to SMILES.
"""

from typing import Dict, Sequence, List
from besmarts.core import assignments


class smiles_assignment_str(assignments.smiles_assignment):

    __slots__ = "smiles", "selections"

    def __init__(self, smiles, selections):
        self.smiles: str = smiles
        self.selections: Dict[Sequence[int], str] = selections

    def copy(self):
        return smiles_assignment_str_copy(self)

    def compute(self, idx) -> str:
        return smiles_assignment_str_compute(self, idx)

def smiles_assignment_str_diff(A: smiles_assignment_str, B: smiles_assignment_str):
    assert A.smiles == B.smiles
    c = {}
    e = ""

    for sa, lbl in A.selections.items():
        if lbl != B.selections.get(sa, ""):
            c[sa] = lbl
    
    return smiles_assignment_str(A.smiles, c)

def smiles_assignment_str_modified(A: List[smiles_assignment_str], B: List[smiles_assignment_str]):

    lbls = set()

    for a, b in zip(A, B):
        c = smiles_assignment_str_diff(a, b)
        lbls.update(c.selections.values())

    return list(lbls)


class smiles_assignment_float(assignments.smiles_assignment):

    __slots__ = "smiles", "selections"

    def __init__(self, smiles, selections):
        self.smiles: str = smiles
        self.selections: Dict[Sequence[int], List[float]] = selections

    def copy(self):
        return smiles_assignment_float_copy(self)

    def compute(self, idx) -> List[float]:
        return smiles_assignment_float_compute(self, idx)

def smiles_assignment_str_compute(
    sa: smiles_assignment_str, idx: Sequence[int]
) -> str:
    return sa.selections.get(idx, "")


def smiles_assignment_float_compute(
    sa: smiles_assignment_float, idx: Sequence[int]
) -> List[float]:
    return sa.selections.get(idx, [])


def smiles_assignment_float_copy(sa) -> smiles_assignment_float:
    return smiles_assignment_float(sa.smiles, sa.selections.copy())


def smiles_assignment_str_copy(sa) -> smiles_assignment_str:
    return smiles_assignment_str(sa.smiles, sa.selections.copy())
