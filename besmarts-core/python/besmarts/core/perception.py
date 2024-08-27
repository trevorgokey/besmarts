
"""
besmarts.core.perception
"""

from besmarts.core import codecs
from besmarts.core import assignments

class perception_model:
    def __init__(self, gcd, labeler):
        self.gcd: codecs.graph_codec = gcd
        self.labeler: assignments.smarts_hierarchy_assignment = labeler
        self.icd: codecs.intvec_codec = codecs.intvec_codec(
            gcd.primitive_codecs, gcd.atom_primitives, gcd.bond_primitives
        )
