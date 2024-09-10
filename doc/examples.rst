Examples
========

These examples are also available on github.

**1.  SMILES IO (01_smiles_io.py)**

A simple example of getting a graph representation from a SMILES and how to use
the native BESMARTS format.

.. literalinclude:: examples/01_smiles_io.py
   :language: python

**2. Match a SMARTS from a molecule with to parameter SMARTS (02_matching.py)**

This example loads a target SMARTS and performs a few queries on it using
a native BESMARTS function.

.. literalinclude:: examples/02_matching.py
   :language: python

**3. Bond union (03_bond_union.py)**

Example of combining the bonds of propane into a single SMARTS.

.. literalinclude:: examples/03_bond_union.py
   :language: python


**4. Initializing FF parameters (04_init_bonds_angles.py)**

Resetting FF bonds and angles according to a molecule geometry

.. literalinclude:: examples/04_init_bonds_angles.py
   :language: python

**5. Build a force field for ethane (05_forcefield_build.py)**

This example builds a force field with parameters specific to ethane.

.. literalinclude:: examples/05_forcefield_build.py
   :language: python

**6. Cluster bond lengths (06_smarts_cluster_bond_lengths.py)**

Use automated chemical perception in BESMARTS to determine the SMARTS patterns
needed to cluster bond lengths that differ 0.1 Angstroms using SMARTS patterns
as the cluster labels.

.. literalinclude:: examples/06_smarts_cluster_bond_lengths.py
   :language: python

**7. Automated chemical perception for force field design (07_besmarts_fit.py)**

This performs automated chemical perception starting from OpenFF Sage 2.1. For
expedience, only b4 is targeted and only the lengths are fit. The resulting
force field was fit on positions and gradients of a single molecule.

.. literalinclude:: examples/07_besmarts_fit.py
   :language: python

**8. Find all splits using numerical SMARTS search (08_smarts_split.py)**

Enumerate all SMARTS patterns that can split (partition) the bonds of CCO into
two groups.

.. literalinclude:: examples/08_smarts_split.py
   :language: python


**9. Minimize a molecules energy by optimizing positions (09_minimize.py)**

Perform a standard geometry minimization.

.. literalinclude:: examples/09_minimize.py
   :language: python

**10. Vibrational frequencies of ethane with OpenMM (10_vibfreq.py)**

Calculate the MM vibrational frequencies of ethane. Also shows how to calculate
energies using the pure python implementation versus OpenMM. Compares between
the SciPy and OpenMM backends for minimization.

.. literalinclude:: examples/10_vibfreq.py
   :language: python

**11. Ab-Initio (valence) chemical perception (AICP) (11_aicp.py)**

Ab-Initio (valence) chemical perception (AICP)

.. literalinclude:: examples/11_aicp.py
   :language: python

