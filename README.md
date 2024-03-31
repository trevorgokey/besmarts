
# BESMARTS

*A toolkit for data-driven force field design based on binary-encoded SMARTS*

Highlights of this package are:

* Map and perform bitwise operations between two molecular substructures of
  arbitrary size
* Search/iterate a substructure at the SMARTS-primitive level, using both
  _numerical_ and _analytic_ approaches
* Cluster molecular data by SMARTS using a SMARTS hierarchy
* Calculate energy and gradients using a classical force field based on (a
  basic implementation of) the SMIRNOFF format
* Geometry optimization
* Force field parameter optimization (under development)

See the ChemRxiv [preprint](https://doi.org/10.26434/chemrxiv-2023-v969f-v3)
for the theoretical unpinnings on which this package is based.

# Installation

Currently, the best way to install is to clone and then install with pip.

For environment users (e.g. venv or conda), one should probably create an empty
environment first:
```
conda create -n besmarts python
conda activate besmarts
```
or
```
python -m venv besmarts
. besmarts/bin/activate

```

followed by the actual install:

```
git clone https://github.com/trevorgokey/besmarts besmarts-git
cd besmarts-git/besmarts-core/python
python -m pip install .
cd ../../besmarts-rdkit/python
python -m pip install .
```

RDKit is needed to decode SMILES into graphs and offers a faster implementation
of SMARTS matching when labeling from a SMARTS hierarchy.

Geometry optimization uses the SciPy minimizer and can be installed using
using a similar process as above with `besmarts-scipy`.

Energy and gradient evaluations are implemented, but require partial charges. By
default, `besmarts` with try to charge molecules using the `sqm` program from
`ambertools` suite. Consequently, make sure `sqm` is in your path by installing
via `conda` or by other means.

# Documentation

Documentation in this repository is hosted on [RTD](https://besmarts.readthedocs.io)

