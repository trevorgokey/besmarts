
"""
besmarts.data.smiles_graphs
"""

from besmarts.codecs import codec_native
from besmarts.core import primitives

def graph_water():
    s = [x.split() for x in [
        "#GRAPH",
        "#ATOM element",
        "#BOND bond_order",
        "1 1 256", # O
        "2 2 2",   # H
        "3 3 2",   # H
        "1 2 1",   # single bond
        "1 3 1"    # single bond
    ]]
    g = codec_native.graph_load(s)

    # codecs = codec_native.primitive_codecs_get()
    # codecs = {k:codecs[k] for k in ("element", "bond_order")}
    # gcd = codec_native.graph_codec_native(codecs, ["element"], ["bond_order"])
    # print(g.nodes[1].primitives)
    # print(gcd.smarts_encode(g))
    # print(gcd.smiles_encode(g))
    
    return g

def graph_dinitrogen():
    smi = "[N:1]#[N:2]"
    sel = {(1,): [[1.0, 0.0, 0.0]], (2,): [[0.0, 0.0, 0.0]]}
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos

def pos_water():
    smi = "[H:1]-[O:2]-[H:3]"
    sel = {(1,): [[1.0, 0.0, 0.0], [1.2, 0.0, 0.0]], (2,): [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]], (3,): [[0.0, 1.0, 0.0], [0.0, 1.1, 0.0]]}
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos

def ethane_pos():

    # unit: kJ/mol
    # Frame,NonbondedForce,PeriodicTorsionForce,HarmonicAngleForce,HarmonicBondForce,TotalEnergy
    # 0,4.107396125793457,0.15501700341701508,40.698150634765625,0.273992121219635,45.23455588519573
    # LJ is 0.3149093985557556
    # QQ is 3.7924864292144775
    # using the chem sys below. should be sage 2.0 with oe am1bcc

    smi = "[C:1]([H:3])([H:4])([H:5])[C:2]([H:6])([H:7])([H:8])"
    sel = {
        (1,): [[10*-0.18969499, 10*-0.3937415 , 10*-1.1148261  ]],
        (2,): [[10*-0.05805168, 10*-0.31429192, 10*-1.1157967  ]],
        (3,): [[10*-0.27382693, 10*-0.32386214, 10*-1.1170473  ]],
        (4,): [[10*-0.2049986 , 10*-0.45749822, 10*-1.0272645  ]],
        (5,): [[10*-0.1928178 , 10*-0.45454127, 10*-1.2057095  ]],
        (6,): [[10* 0.0315621 , 10*-0.3762089 , 10*-1.1258872  ]],
        (7,): [[10*-0.04893475, 10*-0.25069237, 10*-1.0272638  ]],
        (8,): [[10*-0.05737103, 10*-0.25367138, 10*-1.206851   ]],
    }
    pos = assignments.smiles_assignment_float(smi, sel)
    return pos
