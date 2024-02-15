
"""
besmarts.mechanics.smirnoff_xml
"""

from typing import List, Dict
from xml.etree import ElementTree

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
