
Examples
========

__1.  SMILES IO (01_smiles_io.py)__

A simple example of getting a graph representation from a SMILES and how to use
the native BESMARTS format.

__2. Match a SMARTS from a molecule with to parameter SMARTS (02_matching.py)__

This example loads a target SMARTS and performs a few queries on it using
a native BESMARTS function.

__3. Bond union (03_bond_union.py)__

Example of combining the bonds of propane into a single SMARTS.

__4. Initializing FF parameters (04_init_bonds_angles.py)__

Resetting FF bonds and angles according to a molecule geometry

__5. Build a force field for ethane (05_forcefield_build.py)__

This example builds a force field with parameters specific to ethane.

__6. Cluster bond lengths (06_smarts_cluster_bond_lengths.py)__

Use automated chemical perception in BESMARTS to determine the SMARTS patterns
needed to cluster bond lengths that differ 0.1 Angstroms using SMARTS patterns
as the cluster labels.

__7. Automated chemical perception for force field design (07_besmarts_fit.py)__

This performs automated chemical perception starting from OpenFF Sage 2.1. For
expedience, only b4 is targeted and only the lengths are fit. The resulting
force field was fit on positions and gradients of a single molecule.

__8. Find all splits using numerical SMARTS search (08_smarts_split.py)__

Enumerate all SMARTS patterns that can split (partition) the bonds of CCO into
two groups.

__9. Minimize a molecules energy by optimizing positions (09_minimize.py)__

Perform a standard geometry minimization.

__10. Vibrational frequencies of ethane with OpenMM (10_vibfreq.py)__

Calculate the MM vibrational frequencies of ethane. Also shows how to calculate
energies using the pure python implementation versus OpenMM. Compares between
the SciPy and OpenMM backends for minimization.

__11. Ab-Initio (valence) chemical perception (AICP) (11_aicp.py)__

Ab-Initio (valence) chemical perception (AICP)

This is a full valence fit of parameters starting from a base, empty set of
parameters. The initial parameters are first rapidly expanded by clustering
bonds, angles, and torsions based on Hessian projection, followed by a BESMARTS
parameter search which will refine the parameters to fit the objective. Here
we fit on the geometry, gradient, and frequencies of a single molecule. The
molecule is chemically diverse and has 25 atoms, so it presents a typical
challenge in force field design.

Estimated time to complete is 1 hour.
