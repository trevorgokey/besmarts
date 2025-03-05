
"""
besmarts.mechanics.smirnoff_xml
"""

from typing import Dict
from xml.etree.ElementTree import ElementTree, Element, iterparse, indent

def smirnoff_xml_read(fname: str) -> Dict:

    xml = {}

    has_parameters = {
        "Bonds": "Bond",
        "Angles": "Angle",
        "ProperTorsions": "Proper",
        "ImproperTorsions": "Improper",
        "vdW": "Atom",
        "LibraryCharges": "LibraryCharge",
        "Constraints": "Constraint",
        "ChargeIncrementModel": "ChargeIncrement"
    }

    active = ""

    for event, element in iterparse(fname, ["start", "end"]):
        # print(element.tag, event, active, element.attrib)

        if active and event == "start":
            if element.tag in has_parameters:
                xml[active]["options"] = dict(element.attrib)
            else:
                p = {has_parameters[active]: dict(element.attrib)}
                xml[active]["parameters"].append(p)
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


def smirnoff_xml_write(xml: Dict, fname: str):
    """
    """

    root_tag = "SMIRNOFF"
    root_attrib = xml.pop(root_tag)["options"]
    root = Element(root_tag, attrib=root_attrib)
    tree = ElementTree(root)

    for name, p in xml.items():
        proot = Element(name, attrib=p["options"])
        root.append(proot)
        for param in p["parameters"]:
            for pname, terms in param.items():
                node = Element(pname, attrib=terms)
                proot.append(node)

    indent(tree, space="  ", level=0)
    tree.write(fname, method="xml")
    return tree
