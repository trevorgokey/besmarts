
from besmarts.core import perception
from besmarts.codecs import codec_rdkit
from besmarts.assign import hierarchy_assign_rdkit


def perception_model_rdkit():
    gcd = codec_rdkit.graph_codec_rdkit()
    labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
    return perception.perception_model(gcd, labeler)
