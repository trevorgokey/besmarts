
"""
besmarts.core.perception
"""

from besmarts.core import codecs
from besmarts.core import assignments

class perception_model:
    """Functions that can be used to load, label and serialize graphs

    Attributes
    ----------
    gcd: codecs.graph_codec
        A graph codec, responsible for loading graphs from SMILES/SMARTS
        strings. Initialized from the ``gcd`` parameter.
    labeler: assignments.smarts_hierarchy_assignment
        A labeler, responsible for assigning labels to graphs by matching them
        against a SMARTS hierarchy. Initialized from the ``labeler`` parameter.
    icd: codecs.intvec_codec
        A vector-of-integers codec, responsible for serializing graphs into a
        compact format. Initialized automatically from the primitives specified
        by the ``gcd`` parameter.
    """
    def __init__(self, gcd: codecs.graph_codec, labeler: assignments.smarts_hierarchy_assignment):
        """Functions that can be used to load, label and serialize graphs

        Parameters
        ----------
        gcd
            A graph codec, responsible for loading graphs from SMILES/SMARTS
            strings.
        labeler
            A labeler, responsible for assigning labels to graphs by matching them
            against a SMARTS hierarchy.
        """
        self.gcd: codecs.graph_codec = gcd
        self.labeler: assignments.smarts_hierarchy_assignment = labeler
        self.icd: codecs.intvec_codec = codecs.intvec_codec(
            gcd.primitive_codecs, gcd.atom_primitives, gcd.bond_primitives
        )
