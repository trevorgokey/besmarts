
"""
besmarts.mechanics.smirnoff_xml
"""

from typing import List, Dict
from xml.etree import ElementTree, Element

def smirnoff_xml_read(fname: str) -> Dict:

    xml = {}

    has_parameters = [
        "Bonds",
        "Angles",
        "ProperTorsions",
        "ImproperTorsions",
        "vdW",
        "LibraryCharges",
        "Constraints"
        "ChargeIncrementModel"
    ]

    active = ""

    for event, element in ElementTree.iterparse(fname, ["start", "end"]):
        # print(element.tag, event, active, element.attrib)

        if active and event == "start":
            if element.tag in has_parameters:
                xml[active]["options"] = dict(element.attrib)
            else:
                xml[active]["parameters"].append(dict(element.attrib))
        elif element.tag in has_parameters:
            if event == "start":
                active = element.tag
                xml[active] = {}
                xml[active]["options"] = dict(element.attrib)
                xml[active]["parameters"] = []
            elif event == "end":
                active = ""
        elif not active:
            xml[element.tag] = {}
            xml[element.tag]["options"] = dict(element.attrib)
            xml[element.tag]["parameters"] = []

    return xml

def smirnoff_xml_write_version_0p3(xml: Dict, fname: str):
    """
    xml = {
        Bonds: tag: {
            options: {k:v}, # store in attribs
            parameters: [ # these are subelements
                Bond: tag: {k:v} # terms are attribs
            ]
    }

    """
    model_to_parameter = {
        "Bonds": "Bond",
        "Angles": "Angle",
        "ProperTorsions": "ProperTorsion",
        "ImproperTorsions": "ImproperTorsion",
        "vdW": "Atom",
        "LibraryCharges": "LibraryCharge",
        "Constraints": "Constraint",
        "ChargeIncrementModel": "ChargeIncrement"
    }
    model_keys = {
        "Bonds": ["smirks", "id", "k", "length"],
        "Angles": ["smirks", "id", "k", "angle"],
        "ProperTorsions": ["smirks", "id", "periodicity", "k", "phase"],
        "ImproperTorsions": ["smirks", "id", "periodicity", "k", "phase"],
        "vdW": ["smirks", "id", "epsilon", "sigma"],
        "LibraryCharges": ["smirks", "id", "charge"],
        "Constraints": ["smirks", "id", "c"],
        "ChargeIncrementModel": ["smirks", "id", "charge"]
    }
    model_multi_parameter = {
        "Bonds": False,
        "Angles": False,
        "ProperTorsions": True,
        "ImproperTorsions": True,
        "vdW": False,
        "LibraryCharges": True,
        "Constraints": True,
        "ChargeIncrementModel": True
    }

    root_attrib = {
        "version"="0.3" "aromaticity_model":"MDL"
    }

    root = Element("SMIRNOFF", attrib=root_attrib)
    tree = ElementTree(root)

    for toptag, topvals in xml.items():
        node = Element(toptag)
        if toptag in has_parameters:
            param_root = Element(model_to_parameter[toptag])
            node.attrib = topvals.get("options", {})
            for ptag, pvals in topvals["parameters"]:
                for x in ["smirks", "id"]:
                    assert x in pvals
                param = Element(ptag, attrib=pvals)
                param_root.append(param)
        else:
            node.text = topvals.pop("parameters", "")
            attrib = topvals.pop("options", {})
            root.append(node, attrib=attrib)

    tree.write(fname)
    return tree
