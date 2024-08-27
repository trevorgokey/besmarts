"""
examples/10_vibfreq.py

Calculate the MM vibrational frequencies of ethane. Also shows how to calculate
energies. Compares between the SciPy and OpenMM backends for minimization.
"""

import os
import tempfile
import time

# using scipy so numpy is present
import numpy as np
from besmarts.core import graphs

from besmarts.mechanics import vibration
from besmarts.mechanics import objectives
from besmarts.mechanics import smirnoff_models
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import optimizers_scipy
from besmarts.mechanics import optimizers_openmm

from besmarts.core import assignments
from besmarts.perception import perception_rdkit

# Load the RDKit perception model
pcp = perception_rdkit.perception_model_rdkit()

# The molecule definition
smi = "[C:1]([H:3])([H:4])([H:5])[C:2]([H:6])([H:7])([H:8])"
ethane_xyz = """8

C   0.93354  0.05210 -0.06101
C   2.46698  0.05208 -0.06093
H   0.53566 -0.55247 -0.88337
H   0.53569  1.06658 -0.17346
H   0.53558 -0.35772  0.87374
H   2.86494  0.46192 -0.99569
H   2.86483 -0.96239  0.05151
H   2.86486  0.65667  0.76142"""

# The force field
xml = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="AROMATICITY_MDL">
    <Bonds version="0.4" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
        <Bond smirks="[#6X4:1]-[#6X4:2]" id="b1" length="1.527940216866 * angstrom" k="419.9869268191 * angstrom**-2 * mole**-1 * kilocalorie"></Bond>
        <Bond smirks="[#6X4:1]-[#1:2]" id="b84" length="1.090139506109 * angstrom" k="719.6424928981 * angstrom**-2 * mole**-1 * kilocalorie"></Bond>
    </Bonds>
    <Angles version="0.3" potential="harmonic">
        <Angle smirks="[*:1]~[#6X4:2]-[*:3]" angle="110.0631999136 * degree" k="121.1883270155 * mole**-1 * radian**-2 * kilocalorie" id="a1"></Angle>
        <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="108.5839257083 * degree" k="75.08254435747 * mole**-1 * radian**-2 * kilocalorie" id="a2"></Angle>
    </Angles>
    <ProperTorsions version="0.4" potential="k*(1+cos(periodicity*theta-phase))" default_idivf="auto" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
        <Proper smirks="[#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]" periodicity1="3" phase1="0.0 * degree" id="t3" k1="0.2516073078789 * mole**-1 * kilocalorie" idivf1="1.0"></Proper>
    </ProperTorsions>
    <ImproperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))" default_idivf="auto">
    </ImproperTorsions>
    <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1.0" cutoff="9.0 * angstrom" switch_width="1.0 * angstrom" method="cutoff">
        <Atom smirks="[#1:1]-[#6X4]" epsilon="0.01577948280971 * mole**-1 * kilocalorie" id="n2" rmin_half="1.48419980825 * angstrom"></Atom>
        <Atom smirks="[#6X4:1]" epsilon="0.1088406109251 * mole**-1 * kilocalorie" id="n16" rmin_half="1.896698071741 * angstrom"></Atom>
    </vdW>
    <Electrostatics version="0.3" scale12="0.0" scale13="0.0" scale14="0.8333333333" scale15="1.0" cutoff="9.0 * angstrom" switch_width="0.0 * angstrom" method="PME"></Electrostatics>
    <LibraryCharges version="0.3">
    </LibraryCharges>
    <ToolkitAM1BCC version="0.3"></ToolkitAM1BCC>
</SMIRNOFF>"""

# Normally, this is already a file and it can be loaded directly
fd, fname = tempfile.mkstemp(suffix=".offxml")
with open(fd, 'w') as f:
    f.write(xml)
csys = smirnoff_models.smirnoff_load(fname, pcp)
os.remove(fname)

# read the molecule
pos = assignments.xyz_to_graph_assignment(pcp.gcd, smi, ethane_xyz)

# parameterize
psys = mm.chemical_system_to_physical_system(csys, [pos])

# optimize
t = time.perf_counter_ns()
# Set this to 1 to use OpenMM
if 0:
    print("Minimizing ethane with SciPy")
    minpos = optimizers_scipy.optimize_positions_scipy(csys, psys)
    energy = objectives.physical_system_energy(psys, csys)
else:
    print("Minimizing ethane with OpenMM")
    minpos = optimizers_openmm.optimize_positions_openmm(csys, psys)
    energy = optimizers_openmm.physical_system_energy_openmm(psys, csys)
t = time.perf_counter_ns() - t

print(f"Energy: {energy:.8f} kJ/mol. Elapsed: {t*1e-9:.4f} sec")

# A quick way to update the psys with the new positions.
# Reuse everything so that all existing parameters are kept and not reprocessed
psys = mm.chemical_system_to_physical_system(
    csys,
    [minpos],
    ref=psys,
    reuse=list(psys.models)
)


def freq(hess, pos):

    hess = np.array(hess) / 4.184
    syms = graphs.graph_symbols(pos.graph)
    xyzs = np.vstack([*pos.selections.values()])
    mass = [3*[vibration.mass_table[syms[s]]] for s in syms]

    # set to verbose=True to get all vibrational modes saved to xyz files
    # The zero is a number to make the files unique
    f, vib = vibration.hessian_modes(hess, syms, xyzs, mass, 0, verbose=False)

    return f


# Hessian using pure python to evaluate energy
t = time.perf_counter_ns()
hess_bes = objectives.physical_system_hessian(psys, csys, h=1e-4)
t = time.perf_counter_ns() - t
print(f"Hessian in pure Python time: {t*1e-9:.4f} sec")
fb = freq(hess_bes, minpos)

# Hessian using OpenMM to evalulate energy
t = time.perf_counter_ns()
hess_omm = optimizers_openmm.physical_system_hessian_openmm(psys, csys, h=1e-4)
t = time.perf_counter_ns() - t
print(f"Hessian in OpenMM time:      {t*1e-9:.4f} sec")
fo = freq(hess_omm, minpos)

print("\n")
print("       Frequencies (cm-1)       ")
print("================================")
print("Vib.  besmarts   openmm    delta")
print("--------------------------------")
for i, (a, b) in enumerate(zip(fb, fo), 1):
    prefix = "V"
    if i < 4:
        prefix = "T"
    elif i < 7:
        prefix = "R"
    else:
        i -= 6

    lhs = f"{prefix + str(i):>4s}  "

    cols = [a, b, a-b]
    rhs = map("{:8.2f}".format, cols)

    print(lhs + " ".join(rhs))
print("================================")
