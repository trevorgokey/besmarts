"""
examples/init_bonds_angles.py
"""

import os
import tempfile
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import smirnoff_models
from besmarts.perception import perception_rdkit

pcp = perception_rdkit.perception_model_rdkit()

FF = """
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds
        version="0.4"
        potential="harmonic"
        fractional_bondorder_method="AM1-Wiberg"
        fractional_bondorder_interpolation="linear"
    >
        <Bond
            smirks="[*:1]~[*:2]"
            id="b1"
            length="1.5 * angstrom"
            k="500.0 * angstrom**-2 * mole**-1 * kilocalorie"
        ></Bond>
        <Bond
            smirks="[*:1]~[#1:2]"
            id="b2"
            length="1.0 * angstrom"
            k="700.0 * angstrom**-2 * mole**-1 * kilocalorie"
        ></Bond>
    </Bonds>
    <Angles version="0.3" potential="harmonic">
        <Angle
            smirks="[*:1]~[X4:2]~[*:3]"
            angle="109.5 * degree"
            k="60.0 * mole**-1 * radian**-2 * kilocalorie"
            id="a1"
        ></Angle>
        <Angle
            smirks="[*:1]~[X3:2]~[*:3]"
            angle="120.0 * degree"
            k="80.0 * mole**-1 * radian**-2 * kilocalorie"
            id="a2"
        ></Angle>
        <Angle
            smirks="[*:1]~[X2:2]~[*:3]"
            angle="180.0 * degree"
            k="10.0 * mole**-1 * radian**-2 * kilocalorie"
            id="a3"
        ></Angle>
    </Angles>
    <ProperTorsions
        version="0.4"
        potential="k*(1+cos(periodicity*theta-phase))"
        default_idivf="auto"
        fractional_bondorder_method="AM1-Wiberg"
        fractional_bondorder_interpolation="linear"
    >
        <Proper
            smirks="[*:1]~[*:2]~[*:3]~[*:4]"
            id="t1"
            periodicity1="3"
            phase1="0.0 * degree"
            k1="0.15 * mole**-1 * kilocalorie"
            idivf1="1.0"
        ></Proper>
    </ProperTorsions>
    <ImproperTorsions 
        version="0.3"
        potential="k*(1+cos(periodicity*theta-phase))"
        default_idivf="auto"
    > </ImproperTorsions>
    <vdW
        version="0.3"
        potential="Lennard-Jones-12-6"
        combining_rules="Lorentz-Berthelot"
        scale12="0.0"
        scale13="0.0"
        scale14="0.5"
        scale15="1.0"
        cutoff="9.0 * angstrom"
        switch_width="1.0 * angstrom"
        method="cutoff"
    >
        <Atom
            smirks="[*:1]"
            id="n1"
            epsilon="0.01 * mole**-1 * kilocalorie"
            rmin_half="0.75 * angstrom"
        ></Atom>
    </vdW>
    <Electrostatics
        version="0.3"
        scale12="0.0"
        scale13="0.0"
        scale14="0.8333333333"
        scale15="1.0" cutoff="9.0 * angstrom"
        switch_width="0.0 * angstrom"
        method="PME"
    ></Electrostatics>
    <LibraryCharges version="0.3"></LibraryCharges>
    <ToolkitAM1BCC version="0.3"></ToolkitAM1BCC>
</SMIRNOFF>"""

fd, ff_fname = tempfile.mkstemp(suffix=".offxml", text=True)
with os.fdopen(fd, 'w') as f:
    f.write(FF)
csys = smirnoff_models.smirnoff_load(ff_fname, pcp)
os.remove(ff_fname)

sdf = """
  -OEChem-04032417453D

 20 20  0     0  0  0  0  0  0999 V2000
    2.9658   -1.6819    0.0797 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1414   -1.8277    1.4584 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0599   -0.6392    1.8273 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8571   -1.0249   -0.4689 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9149   -0.5111    0.4452 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1828   -1.3026    2.3212 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7114   -0.8926   -1.9660 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3862    0.7898    1.1682 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0681    1.5755   -0.8685 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3432   -0.6385   -0.9442 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5277    0.3426   -0.2210 S   0  0  0  0  0  0  0  0  0  0  0  0
    2.4023   -1.4901    4.2151 Br  0  0  0  0  0  0  0  0  0  0  0  0
    3.7154   -2.0925   -0.5917 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.0124   -2.3426    1.8506 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3115   -0.2291    2.4936 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6653    0.1586   -2.2648 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5580   -1.3612   -2.4724 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7902   -1.3683   -2.3144 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4655    1.8053    1.2001 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3037    0.3467    1.1518 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  4  2  0  0  0  0
  1 13  1  0  0  0  0
  2  6  2  0  0  0  0
  2 14  1  0  0  0  0
  3  5  2  0  0  0  0
  3  6  1  0  0  0  0
  3 15  1  0  0  0  0
  4  5  1  0  0  0  0
  4  7  1  0  0  0  0
  5 11  1  0  0  0  0
  6 12  1  0  0  0  0
  7 16  1  0  0  0  0
  7 17  1  0  0  0  0
  7 18  1  0  0  0  0
  8 11  1  0  0  0  0
  8 19  1  0  0  0  0
  8 20  1  0  0  0  0
  9 11  2  0  0  0  0
 10 11  2  0  0  0  0
M  END
> <qm_energy (kcal/mol)>
-2164061.7164133266

> <qcarchive_id>
20909050

$$$$"""

fd, sdf_fname = tempfile.mkstemp(suffix=".sdf", text=True)
with os.fdopen(fd, 'w') as f:
    f.write(sdf)

pos, extras = pcp.gcd.sdf_decode(sdf_fname)
os.remove(sdf_fname)

psys = mm.chemical_system_to_physical_system(csys, [pos])

angles = mm.chemical_system_reset_angles(csys, [psys])
bonds = mm.chemical_system_reset_bond_lengths(csys, [psys])

show = set((k[2] for k in list(bonds) + list(angles)))
# Angles are in radians
mm.chemical_system_print(csys, show_parameters=show)
