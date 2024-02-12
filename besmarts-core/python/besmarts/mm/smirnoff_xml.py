
"""
besmarts.mm.smirnoff_xml
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

# import pprint
# pprint.pprint(smirnoff_read_xml("/home/tgokey/pol_split.offxml"))

# this should possibly just be load/save?
# later we can open it open to mess with it

def smirnoff_xml_write(fname: str, d: dict) -> bool:
    pass

# def forcefield_smirnoff_save_xml(hidx: forcefields.forcefield, fname: str) -> bool:

#     out = ""
#     root = ElementTree.Element("SMIRNOFF", attrib=hidx.db["SMIRNOFF"])
#     for force_tag, idx in hidx.force_entry_idx.items():
#         if idx is None:
#             continue
#         hent = hidx[idx]
#         key = hent.key
#         force = hidx.forces.get(key)
#         if force is None:
#             print("Warning: force is missing for", key)
#             continue
#         force_root = ElementTree.SubElement(root, force_tag, hidx.db[key])
#         params = hindex_iter_dive(hidx, hent)
#         _ = next(params)
#         for param in params:

#             lbl = param.key
#             tag = hidx.db[lbl]["tag"]
#             param = hidx.db[lbl]["parameter"]
#             attribs = {"id": lbl, "smirks": param.smarts}
#             for dim, unit in force.units.items():
#                 if param.term_type[dim] == term_type.SCALAR:
#                     term_str = str(param.terms[dim][0])
#                     if unit:
#                         term_str = term_str + " * " + unit
#                     attribs[dim] = term_str
#                 elif param.term_type[dim] == term_type.ARRAY:
#                     for i, term in enumerate(param.terms[dim], 1):
#                         term_str = str(param.terms[dim][i - 1])
#                         if unit:
#                             term_str = term_str + " * " + unit
#                         attribs[dim + str(i)] = term_str
#             node = ElementTree.SubElement(force_root, tag, attribs)
#     for idx in hidx.root.down:
#         if idx in hidx.force_entry_idx.values():
#             continue
#         hent = hidx[idx]
#         attrib = dict(hidx.db[hent.key])
#         attrib.pop("tag")
#         node = ElementTree.SubElement(root, hent.key, hidx.db[hent.key])

#     # electr = ElementTree.SubElement(root, "Electrostatics",
#     #     dict(
#     #         version="0.3",
#     #         scale12="0.0",
#     #         scale13="0.0",
#     #         scale14=str(1/1.2),
#     #         scale15="1.0",
#     #         cutoff="9.0 * angstrom",
#     #         switch_width="0.0 * angstrom",
#     #         method="PME",
#     #     )
#     # )
#     # am1bcc = ElementTree.SubElement(root, "ToolkitAM1BCC", dict(version="0.3"))

#     out = "\n".join(
#         ElementTree.tostringlist(root, encoding="unicode", short_empty_elements=False)
#     )
#     # ElementTree.indent(tree)
#     with open(fname, "w") as f:
#         f.write(out)
#     return True
    # tree.write(fname, encoding='utf')

