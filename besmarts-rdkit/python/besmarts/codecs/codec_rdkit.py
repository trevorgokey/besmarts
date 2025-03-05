"""
besmarts.besmarts_rdkit.codecs

BESMARTS graph encoding using the RDKit perception model
"""

from typing import Dict, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from besmarts.core import configs
from besmarts.core.arrays import bitvec as array
from besmarts.core.arrays import array_dtype
from besmarts.core.chem import bechem
from besmarts.core import graphs
from besmarts.core import assignments
from besmarts.core.graphs import graph
from besmarts.core.primitives import primitive_key, primitive_codec
from besmarts.core.codecs import graph_codec
from besmarts.core.codecs import (
    primitive_codec_element,
    primitive_codec_hydrogen,
    primitive_codec_connectivity_total,
    primitive_codec_connectivity_ring,
    primitive_codec_ring_smallest,
    primitive_codec_aromatic,
    primitive_codec_chirality,
    primitive_codec_valence,
    primitive_codec_formal_charge,
    primitive_codec_bond_order,
    primitive_codec_bond_ring,
)



class graph_codec_rdkit(graph_codec):
    def __init__(self, atom_primitives=None, bond_primitives=None):

        smiles_config = configs.smiles_perception_config(
            False,
            True,
            False,
            "MDL",
        )

        primitive_codecs: Dict[
            primitive_key, primitive_codec
        ] = primitive_codecs_get()

        if atom_primitives is None:

            atom_primitives = tuple(
                (
                    "element",
                    "hydrogen",
                    "connectivity_total",
                    "connectivity_ring",
                    "ring_smallest",
                    "aromatic",
                    # primitive_key.CHIRALITY,
                    "formal_charge",
                )
            )

        if bond_primitives is None:
            bond_primitives = tuple(
                (
                    "bond_ring",
                    "bond_order",
                )
            )

        super().__init__(
            smiles_config,
            primitive_codecs,
            array,
            atom_primitives,
            bond_primitives,
        )

    def smiles_decode(self, smi) -> graphs.graph:
        return rdkit_smiles_decode(
            self.smiles_config,
            self.primitive_codecs,
            self.array,
            self.atom_primitives,
            self.bond_primitives,
            smi,
        )

    def smarts_decode(self, sma) -> graphs.graph:
        return rdkit_smarts_decode(
            self.primitive_codecs,
            self.array,
            self.atom_primitives,
            self.bond_primitives,
            sma,
        )
    def sdf_decode(self, sdf) -> assignments.graph_assignment:
        """
        """
        sa, extras = rdkit_sdf_to_smiles_assignment(sdf)
        g = rdkit_smiles_decode(
            self.smiles_config,
            self.primitive_codecs,
            self.array,
            self.atom_primitives,
            self.bond_primitives,
            sa.smiles,
        )
        g = graphs.subgraph_as_graph(g)
        smiles = sa.smiles
        # smiles = self.smiles_encode(g)
        return assignments.graph_assignment(smiles, sa.selections, g), extras

    def xyz_decode(self, smi, xyz) -> assignments.graph_assignment:
        """
        """
        sa = rdkit_xyz_to_smiles_assignment(smi, xyz)
        g = rdkit_smiles_decode(
            self.smiles_config,
            self.primitive_codecs,
            self.array,
            self.atom_primitives,
            self.bond_primitives,
            sa.smiles,
        )
        g = graphs.subgraph_as_graph(g)
        smiles = sa.smiles
        # smiles = self.smiles_encode(g)
        return assignments.graph_assignment(smiles, sa.selections, g)

    def rdmol_decode(self, mol: Chem.Mol) -> graphs.graph:
        """
        Build a graph directly from an RDKit molecule. This sidesteps all
        sanitization and manipulation; the graph is built by directly querying
        the molecule for primmitives. Make sure the aromaticity model is
        compatible!

        Parameters
        ----------
        mol: Chem.Mol
            The RDKit molecule to parse

        Returns
        -------
        graphs.graph
        """
        return rdkit_mol_decode(
            self.primitive_codecs,
            self.array,
            self.atom_primitives,
            self.bond_primitives,
            mol
        )

    @staticmethod
    def list_implemented_atom_primitives() -> Sequence[primitive_key]:
        return tuple(list_atom_primitives())

    @staticmethod
    def list_implemented_bond_primitives() -> Sequence[primitive_key]:
        return tuple(list_bond_primitives())


def rdkit_smarts_decode(
    codecs, arr: array, atom_primitives, bond_primitives, sma
) -> graphs.graph:

    if r"$" in sma:
        # print(f"Warning, recursive SMARTS {sma} detected, skipping")
        return sma
    mol = Chem.MolFromSmarts(sma)

    # https://sourceforge.net/p/rdkit/mailman/message/29261087/
    # JP,

    # On Mon, May 14, 2012 at 12:57 PM, JP <jeanp...@in...> wrote:
    # >
    # >
    # > I create a molecule without sanitization (red flag) - because, oh well, I
    # > downloaded this sd file from the web so it must be perfectly curated and
    # > what not.
    # >
    # > When I try Chem.AddHs, I get
    # >
    # > <rdkit.Chem.rdchem.Mol object at 0x140a3d0>
    # > [11:47:00]
    # >
    # > ****
    # > Pre-condition Violation
    # > getNumImplicitHs() called without preceding call to calcImplicitValence()
    # > Violation occurred on line 167 in file
    # > /opt/RDKit_trunk/Code/GraphMol/Atom.cpp
    # > Failed Expression: d_implicitValence>-1
    # > ****
    # >
    # > Traceback (most recent call last):
    # >   File "./test.py", line 107, in <module>
    # >     Chem.AddHs(m, addCoords=True)
    # > RuntimeError: Pre-condition Violation
    # >
    # > I understand that this may be related to having sanitization switched off
    # > (in fact if I turn it on it works), but my question is - is this the correct
    # > error message?
    # > Perhaps calcImplicitValence() should be called regardless of sanitization?

    # In order for the code to add coordinates for the added Hs, it needs to
    # calculate the valence at each atom. You don't need to sanitize the
    # molecules to work around this, the following snippet should work just
    # fine:
    #   m = Chem.MolFromMolBlock(mol_block, sanitize=False)
    #   m.UpdatePropertyCache(strict=False)
    #   mh=Chem.AddHs(m, addCoords=True)

    # Best,
    # -greg

    mol.UpdatePropertyCache(strict=False)

    nodes = {}

    chem_codecs = {name: codecs[name] for name in atom_primitives}
    idx = 1
    selection = get_tags(mol)
    indices = get_indices(mol)
    for atom in mol.GetAtoms():
        idx = indices[atom.GetIdx()]
        primitives = parse_atom(chem_codecs, arr, atom)
        nodes[idx] = bechem(primitives, atom_primitives)

    edges = {}
    chem_codecs = {name: codecs[name] for name in bond_primitives}

    for bond in mol.GetBonds():
        idx_i = indices[bond.GetBeginAtom().GetIdx()]
        idx_j = indices[bond.GetEndAtom().GetIdx()]

        primitives = parse_bond(chem_codecs, arr, bond)
        edges[graphs.edge((idx_i, idx_j))] = bechem(primitives, bond_primitives)

    if selection:
        select = tuple(
            (
                *(selection[i] for i in sorted(selection)),
                *(i for i in sorted(nodes) if i not in selection.values()),
            )
        )
        return graphs.subgraph(nodes, edges, select)
    else:
        return graphs.graph(nodes, edges)

def rdkit_xyz_to_smiles_assignment(smiles, xyz) -> Tuple[assignments.smiles_assignment, Dict]:

    molsmi = Chem.MolFromSmiles(smiles)
    charge = Chem.GetFormalCharge(molsmi)

    mol = Chem.Mol(Chem.MolFromXYZBlock(xyz))
    rdDetermineBonds.DetermineBonds(mol, charge=charge, useAtomMap=True)
    
    indices = get_indices(mol)

    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(indices[atom.GetIdx()])

    smi = Chem.MolToSmiles(mol)

    sel = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        i = indices[idx],
        if i not in sel:
            sel[i] = []
        for conf in mol.GetConformers():
            xyz = conf.GetAtomPosition(idx)
            sel[i].append(list(xyz))

    # extras = mol.GetPropsAsDict()

    return assignments.smiles_assignment_float(smi, sel)

def rdkit_sdf_to_smiles_assignment(sdf) -> Tuple[assignments.smiles_assignment, Dict]:

    mol = next(Chem.SDMolSupplier(sdf, sanitize=False))

    indices = get_indices(mol)

    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(indices[atom.GetIdx()])

    smi = Chem.MolToSmiles(mol)

    sel = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        i = indices[idx],
        if i not in sel:
            sel[i] = []
        for conf in mol.GetConformers():
            xyz = conf.GetAtomPosition(idx)
            sel[i].append(list(xyz))

    extras = mol.GetPropsAsDict()

    return assignments.smiles_assignment_float(smi, sel), extras


def rdkit_smiles_decode(
    pcp, codecs, arr: array, atom_primitives, bond_primitives, smi
) -> graphs.graph:

    global aromaticity_incompatible_warning

    mol = Chem.MolFromSmiles(smi, sanitize=False)

    flags = (
        Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
    )

    if not pcp.protonate:
        flags ^= Chem.SanitizeFlags.SANITIZE_ADJUSTHS

    Chem.SanitizeMol(mol, flags)

    if pcp.protonate:
        mol = Chem.AddHs(mol)

    lut = {
        "OEAroModel_MDL": Chem.AromaticityModel.AROMATICITY_MDL,
        "MDL": Chem.AromaticityModel.AROMATICITY_MDL,
        "RDKIT_MDL": Chem.AromaticityModel.AROMATICITY_MDL,
        "RDKIT": Chem.AromaticityModel.AROMATICITY_RDKIT,
        "RDKIT_SIMPLE": Chem.AromaticityModel.AROMATICITY_SIMPLE,
        "RDKIT_DEFAULT": Chem.AromaticityModel.AROMATICITY_DEFAULT,
    }
    lut.update(Chem.AromaticityModel.names)

    warn = aromaticity_incompatible_warning
    if not warn and pcp.aromaticity == "OEAroModel_MDL":
        print("Warning, aromaticity set to OEAroModel_MDL with RDKit.", end="")
        print(" Model set to RDKit MDL. Expect differences.")
        aromaticity_incompatible_warning = True

    elif pcp.aromaticity not in lut:
        print(f"Could not set aromaticity to {pcp.aromaticity}")
        print("Choose from:")
        for n in lut:
            print(n)
        assert pcp.aromaticity in lut, "Unknown aromaticity model"

    Chem.SetAromaticity(mol, lut[pcp.aromaticity])

    if pcp.strip_hydrogen:
        mol = Chem.RemoveHs(mol)

    return rdkit_mol_decode(
        pcp,
        codecs,
        arr,
        atom_primitives,
        bond_primitives,
        mol
    )


def rdkit_mol_decode(
    pcp, codecs, arr: array, atom_primitives, bond_primitives, mol
) -> graphs.graph:

    nodes = {}

    selection = get_tags(mol)
    indices = get_indices(mol)
    for atom in mol.GetAtoms():
        primitives = {}

        for name in atom_primitives:
            codec: primitive_codec = codecs[name]
            primitives[name] = codec.decode_smiles(arr, atom)
        chem = bechem(primitives, atom_primitives)

        idx = indices[atom.GetIdx()]
        nodes[idx] = chem

    edges = {}

    for bond in mol.GetBonds():
        idx_i = indices[bond.GetBeginAtom().GetIdx()]
        idx_j = indices[bond.GetEndAtom().GetIdx()]
        primitives = {}
        for name in bond_primitives:
            codec: primitive_codec = codecs[name]
            primitives[name] = codec.decode_smiles(arr, bond)

        chem = bechem(primitives, bond_primitives)
        edges[graphs.edge((idx_i, idx_j))] = chem

    if selection:
        select = tuple(
            (
                *(selection[i] for i in sorted(selection)),
                *(i for i in sorted(nodes) if i not in selection.values()),
            )
        )
        return graphs.subgraph(nodes, edges, select)
    else:
        return graphs.graph(nodes, edges)

def list_atom_primitives() -> Sequence[primitive_key]:
    atom_primitives = tuple(
        (
            primitive_key.ELEMENT,
            primitive_key.HYDROGEN,
            primitive_key.CONNECTIVITY_TOTAL,
            primitive_key.CONNECTIVITY_RING,
            primitive_key.RING_SMALLEST,
            primitive_key.AROMATIC,
            primitive_key.CHIRALITY,
            primitive_key.VALENCE,
            primitive_key.FORMAL_CHARGE,
        )
    )
    return atom_primitives


def list_bond_primitives() -> Sequence[primitive_key]:
    bond_primitives = tuple(
        (
            primitive_key.BOND_RING,
            primitive_key.BOND_ORDER,
        )
    )
    return bond_primitives


def parse_chem(codecs, arr: array, chem) -> Dict[primitive_key, array]:

    string = chem.GetSmarts()

    prims = {}
    codec: primitive_codec
    for prim, codec in codecs.items():
        chem = string.replace(";", "").replace("&", "")
        prims[prim] = codec.decode_smarts(arr, chem)

    return prims


def parse_bond(codecs, arr: array, bond) -> Dict[primitive_key, array]:

    return parse_chem(codecs, arr, bond)


def parse_atom(codecs, arr, atom) -> Dict[primitive_key, array]:
    """"""
    return parse_chem(codecs, arr, atom)


def get_tags(mol):
    tag_map = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        tag = atom.GetAtomMapNum()
        if tag > 0:
            tag_map[idx] = tag
    return tag_map


def get_indices(mol):

    tag_map = get_tags(mol)

    nidx = 1
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in tag_map:
            while nidx in tag_map.values():
                nidx += 1
            tag_map[idx] = nidx

    return tag_map


class primitive_codec_element_rdkit(primitive_codec_element):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        x = int(obj.GetAtomicNum())
        array = dtype()
        array[self.encode_int(x)] = True
        array.maxbits = 118
        return array


class primitive_codec_hydrogen_rdkit(primitive_codec_hydrogen):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:

        array = dtype()
        h1 = obj.GetTotalNumHs(includeNeighbors=True)
        h2 = obj.GetNumExplicitHs()

        # h1 seems to work for SMILES, but h2 works for SMARTS.
        # If it doesn't work, then it gives 0, so it seems reasonable
        # to accept a nonzero answer if it is given
        if h1 > h2:
            array[self.encode_int(h1)] = True
        else:
            array[self.encode_int(h2)] = True
        array.maxbits = 5

        return array


class primitive_codec_connectivity_total_rdkit(
    primitive_codec_connectivity_total
):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        x = obj.GetTotalDegree()
        array = dtype()
        array[self.encode_int(x)] = True
        array.maxbits = 5
        return array


class primitive_codec_connectivity_ring_rdkit(
    primitive_codec_connectivity_ring
):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        x = sum([0] + [1 for b in obj.GetBonds() if b.IsInRing()])
        array = dtype()
        array[self.encode_int(x)] = True
        array.maxbits = 5
        return array


class primitive_codec_ring_smallest_rdkit(primitive_codec_ring_smallest):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        ring = 0
        if obj.IsInRing():
            for i in range(3, 103):
                if obj.IsInRingSize(i):
                    ring = i
                    break
        array = dtype()
        array[self.encode_int(ring)] = True
        array.maxbits = 99
        return array


class primitive_codec_aromatic_rdkit(primitive_codec_aromatic):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        x = int(obj.GetIsAromatic())
        array = dtype()
        array[self.encode_int(x)] = True
        array.maxbits = 2
        return array


class primitive_codec_chirality_rdkit(primitive_codec_chirality):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        x = int(obj.GetChiralTag())
        array = dtype()
        array[self.encode_int(x)] = True
        array.maxbits = 3
        return array


class primitive_codec_valence_rdkit(primitive_codec_valence):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        x = int(obj.GetTotalValence())
        array = dtype()
        array[self.encode_int(x)] = True
        array.maxbits = 8
        return array


class primitive_codec_formal_charge_rdkit(primitive_codec_formal_charge):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        x = int(obj.GetFormalCharge())
        array = dtype()
        array[self.encode_int(x)] = True
        array.maxbits = 8
        return array


class primitive_codec_bond_ring_rdkit(primitive_codec_bond_ring):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        x = int(obj.IsInRing())
        array = dtype()
        array[self.encode_int(x)] = True
        array.maxbits = 2
        return array


class primitive_codec_bond_order_rdkit(primitive_codec_bond_order):
    def decode_smiles(self, dtype: array_dtype, obj) -> array:
        order_map = {
            0.0: 0,
            1.0: 1,
            2.0: 2,
            3.0: 3,
            4.0: 4,
            1.5: 5,
        }
        stereomap = {
            2: 6,
            3: 7,
        }
        bo = order_map[obj.GetBondTypeAsDouble()]
        chiral = int(obj.GetStereo())
        if chiral in stereomap:
            bo = stereomap[chiral]

        array = dtype()
        array[self.encode_int(bo)] = True
        array.maxbits = 8
        return array


def primitive_codecs_get() -> Dict[primitive_key, primitive_codec]:
    codecs = {
        "element": primitive_codec_element_rdkit(),
        "hydrogen": primitive_codec_hydrogen_rdkit(),
        "connectivity_total": primitive_codec_connectivity_total_rdkit(),
        "connectivity_ring": primitive_codec_connectivity_ring_rdkit(),
        "ring_smallest": primitive_codec_ring_smallest_rdkit(),
        "aromatic": primitive_codec_aromatic_rdkit(),
        "chirality": primitive_codec_chirality_rdkit(),
        "valence": primitive_codec_valence_rdkit(),
        "formal_charge": primitive_codec_formal_charge_rdkit(),
        "bond_order": primitive_codec_bond_order_rdkit(),
        "bond_ring": primitive_codec_bond_ring_rdkit(),
    }
    return codecs

aromaticity_incompatible_warning = False
