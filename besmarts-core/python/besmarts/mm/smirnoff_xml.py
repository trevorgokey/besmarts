
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

# def forcefield_smirnoff_load_xml(fname: str, gcd: graph_codec) -> forcefields.forcefield:

#     f = forces.force()
#     hidx = forcefields.forcefield()
#     db = dict()
#     hidx.db = db
#     xml = ElementTree.iterparse(fname, ["start", "end"])
#     active = ""
#     up = None
#     entries = []
#     topology = None
#     for event, element in xml:
#         print(element.tag, event, active)
#         if element.tag == "SMIRNOFF":
#             hidx.db[element.tag] = dict(element.attrib)
#             continue
#         elif element.tag in ["Author", "Date"]:
#             continue
#         elif event == "end" and element.tag in hidx.forces:
#             active = ""
#             up.key = element.tag
#             if entries:
#                 hent = hidx.hentry_add(up.index, hentry())
#                 lbl = entries.pop(0)
#                 hent.key = lbl
#                 hidx = hindex_organize_smarts(hidx, hent, entries)
#                 up = None
#             # for idx, e in hidx.entries.items():
#             #     print(idx, e, hidx.node_depth(e))
#             # breakpoint()
#             continue

#         elif element.tag in hidx.forces:
#             active = element.tag
#             topology = hidx.forces[active].topology
#             up = hidx.hentry_add(0, hentry())
#             up.key = element.tag
#             hidx.force_entry_idx[element.tag] = up.index
#             db[element.tag] = dict(element.attrib)
#             db[element.tag]["tag"] = element.tag
#             entries.clear()

#         elif not active and event == "start":
#             hent = hidx.hentry_add(0, hentry())
#             hent.key = element.tag
#             db[element.tag] = dict(element.attrib)
#             db[element.tag]["tag"] = element.tag


#         elif active and event == "start":

#             smarts = element.attrib.get("smirks")
#             lbl = element.attrib.get("id")
#             if lbl is None:
#                 lbl = element.attrib.get("name")
#             if lbl is None:
#                 lbl = smarts
#             if lbl is None:
#                 continue
#             assert lbl not in db, f"Duplicate parameter {lbl}"
#             db[lbl] = {}
#             # g = None
#             # try:
#             #     if "$" not in smarts:
#             #         pass
#             #         # g: subgraph = gcd.smarts_decode(smarts)
#             #         # g: structure = structure(g.nodes, g.edges, g.select, topology)
#             #     else:
#             #         print(f"WARNING: Recursive SMARTS for {lbl}: {smarts}")
#             # except Exception as e:
#             #     print("Could not parse", e)
#             # db[lbl]["graph"] = g

#             entries.append(lbl)
#             # db[lbl] = dict(element.attrib)
#             terms = dict(element.attrib)
#             key = lbl
#             db[lbl]["parameter"] = parameters.topology_parameter(hidx.forces[active], key, smarts, terms)
#             db[lbl]["tag"] = element.tag
#             # print("Loaded:", lbl)
#             # print(db[lbl])
#     return hidx

# def forcefield_smirnoff_load_xml(fname: str, gcd: graph_codec) -> forcefields.forcefield:

#     f = forces.force()
#     hidx = forcefields.forcefield()
#     db = dict()
#     hidx.db = db
#     xml = ElementTree.iterparse(fname, ["start", "end"])
#     active = ""
#     up = None
#     entries = []
#     topology = None
#     force = []
#     for event, element in xml:
#         print(element.tag, event, active)
#         if element.tag == "SMIRNOFF":
#             hidx.db[element.tag] = dict(element.attrib)
#             continue
#         elif element.tag in ["Author", "Date"]:
#             continue
#         elif event == "end" and element.tag in hidx.forces:
#             active = ""
#             up.key = element.tag
#             if entries:
#                 hent = hidx.hentry_add(up.index, hentry())
#                 lbl = entries.pop(0)
#                 hent.key = lbl
#                 hidx = hindex_organize_smarts(hidx, hent, entries)
#                 up = None
#             # for idx, e in hidx.entries.items():
#             #     print(idx, e, hidx.node_depth(e))
#             # breakpoint()
#             continue

#         elif element.tag in hidx.forces:
#             active = element.tag
#             topology = hidx.forces[active].topology
#             up = hidx.hentry_add(0, hentry())
#             up.key = element.tag
#             hidx.force_entry_idx[element.tag] = up.index
#             db[element.tag] = dict(element.attrib)
#             db[element.tag]["tag"] = element.tag
#             entries.clear()

#         elif not active and event == "start":
#             if element.tag == "Bonds":
#                 f = smirnoff_xml_bonds_to_force(xml)
#                 force[f.name] = f  
#             hent = hidx.hentry_add(0, hentry())
#             hent.key = element.tag
#             db[element.tag] = dict(element.attrib)
#             db[element.tag]["tag"] = element.tag


#         elif active and event == "start":

#             smarts = element.attrib.get("smirks")
#             lbl = element.attrib.get("id")
#             if lbl is None:
#                 lbl = element.attrib.get("name")
#             if lbl is None:
#                 lbl = smarts
#             if lbl is None:
#                 continue
#             assert lbl not in db, f"Duplicate parameter {lbl}"
#             db[lbl] = {}
#             # g = None
#             # try:
#             #     if "$" not in smarts:
#             #         pass
#             #         # g: subgraph = gcd.smarts_decode(smarts)
#             #         # g: structure = structure(g.nodes, g.edges, g.select, topology)
#             #     else:
#             #         print(f"WARNING: Recursive SMARTS for {lbl}: {smarts}")
#             # except Exception as e:
#             #     print("Could not parse", e)
#             # db[lbl]["graph"] = g

#             entries.append(lbl)
#             # db[lbl] = dict(element.attrib)
#             terms = dict(element.attrib)
#             key = lbl
#             db[lbl]["parameter"] = parameters.topology_parameter(hidx.forces[active], key, smarts, terms)
#             db[lbl]["tag"] = element.tag
#             # print("Loaded:", lbl)
#             # print(db[lbl])
#     return hidx

# def make_forcefield_smirnoff_empty() -> forcefields.forcefield:
#     models = {}
#     models["Bonds"] = [forces.chemical_model_force_bond_harmonic()]
#     if False:
#         models["Angles"] = [ forces.chemical_model_force_angle_harmonic() ]
#         models["Torsions"] = [ forces.chemical_model_force_torsion_periodic_2term() ]
#         models["Impropers"] = [ forces.chemical_model_force_outofplane_periodic_2term() ]
#         models["vdW"] = [forces.chemical_model_force_vdw_lennard_jones_cutoff()]

#         # need to handle electrostatics somehow
#         # thinking I can make besmarts-offtk and then charge and save assignments
#         # this means I need a method to save assignments
#         # ASSIGNMENT FLOAT ATOM ; x y z charge 
#         # SMILES ASDASDAS
#         # 1 1.0 2.0 3.0 0.0
#         models["Electrostatics"] = [forces.chemical_model_force_coulomb_cutoff()]

#     ff = forcefields.forcefield(models)
#     # ff.model.forces["LibraryCharges"] = [forces.chemical_model_coulomb_cutoff()]
#     # ff.model.forces["ChargeIncrementModel"] = forces.chemical_model_coulomb_cutoff()
#     return ff


#def chemical_model_bond_harmonic_smarts_smirnoff_xml(gcd, labeler, xml) -> forces.chemical_model:
#    """
#    unit_hierarchy
#        unit 1 0 bond  [#6:1]~[#6:2]            ; ETHANE
#        unit 2 0 angle [#6:1]~[#6:2]~[#6:3]     ; PROPANE
#        unit 3 0 bond  [@1:1]<2>1~[@1:2]        ; BUTANE
#        unit 4 0 angle [#6X4:1]([#1:2])([#1:3]) ; CH2

#    force bonds bond_harmonic
#    topology_term l length kJ/mol float
#    topology_term k stiffness kJ/mol/A float
#    smarts_hierarchy k l
#        smarts 1 0 [*:1]~[*:2] 1 1.5
#    index_hierarchy u
#        unit single [@1:1] k l
#            unit_parameter 1-2 1 1.5
#        unit single [@2:1] k l
#            unit_parameter 1-2 1 1.5
#            unit_parameter 2-3 1 1.5
#        unit single [@3:1] k l ; this means we expand to the CG
#            unit_parameter 1-2 1 1.5
#        unit pair [@1:1]~[@1:2] k l ; pairs and beyond do not expand
#            unit_parameter 1-2 1 1.5
#            unit_parameter 2-1 1 1.5
        
#            hierarchy a bond k l
#            parameter 1 1 1.5 [*:1]~[*:2] ; comment
#    """

#    cm = chemical_model("B", "Bonds", topology.bond)

#    cm.energy_function = energy_function_spring 
#    cm.force_function = force_function_spring 
#    cm.internal_function = assignments.smiles_assignment_geometry_bonds

#    # define the terms of this model
#    cm.topology_terms = {
#        "k": topology_term("k", "stiffness", "float", "kJ/mol/A/A", {}, "", {}),
#        "l": topology_term("l", "length", "float", "kJ/mol/A", {}, "", {}),
#    }

#    ############################################################################
#    # the terms are determined by a smarts matching
#    proc = chemical_model_procedure_smarts_assignment(gcd, labeler, cm.topology_terms)
#    # define the perception model for this force
#    # the main hierarchy is first searched, then the inner hierarchy is searched

#    proc.unit_hierarchy = None

#    # we descend into the unit hierarchy when the SMARTS matches at least once,
#    # and I think we only search superSMARTS here and so disable all atomic
#    # primitives
#    # we may want to provide an explicit mapping here

#    # proc.unit_hierarchy.smarts[0] = "[*:1]"
#    uindex = 0

#    # this is the hierarchy within the first main hierarchy node
#    proc.smarts_hierarchies = {
#        uindex: hierarchies.structure_hierarchy(
#            trees.tree_index(), {}, {}, topology.bond
#        )
#    }

#    # # create a default param
#    for bondline in bondlines:
#        i = proc.smarts_hierarchies[uindex].index.node_add_below(None)
#        i.name = "b1"
#        proc.smarts_hierarchies[uindex].smarts[i.index] = "[*:1]~[*:2]"

#        # define the parameter order for this force
#        proc.topology_parameters[(uindex, i.name)] = {"k": i.name, "l": i.name}

#        # add the term values
#        cm.topology_terms["k"].values[i.name] = [500.0]
#        cm.topology_terms["l"].values[i.name] = [1.5]

#    cm.procedures.append(proc)

#    return cm

