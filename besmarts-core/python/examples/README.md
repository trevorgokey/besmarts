
Examples
========

1. SMILES IO (1_smiles_io.py)
-----------------------------
A simple example of getting a graph representation from a SMILES and how to
use the native BESMARTS format.

2. Match a SMARTS from a molecule with to parameter SMARTS (4_matching.py)
--------------------------------------------------------------------------
This example loads a target SMARTS and performs a few queries on it using
a native BESMARTS function.

3. Bond union (3_bond_union.py)
-------------------------------
Example of combining the bonds of propane into a single SMARTS.

4. Initializing FF parameters (4_init_bonds_angles.py)
------------------------------------------------------
Resetting FF bonds and angles according to a molecule geometry

5. Perform a force field fit on ethane (5_forcefield_fit.py) 
------------------------------------------------------------
This example builds a force field and fits the parameters to ethane. No
parameter search or chemical perception is performed.

6. Cluster bond lengths (6_smarts_cluster_bond_lengths.py)
----------------------------------------------------------
Use automated chemical perception in BESMARTS to determine the SMARTS patterns
needed to cluster bond lengths that differ 0.1 Angstroms using SMARTS patterns
as the cluster labels.

7. Automated chemical perception for force field design (7_besmarts_fit.py)
---------------------------------------------------------------------------
This performs automated chemical perception starting from OpenFF Sage 2.1. For
expedience, only b4 is targeted and only the lengths are fit. The resulting
force field was fit on positions and gradients of a single molecule.
