
# BESMARTS

*A toolkit for data-driven force field design based on binary-encoded SMARTS*

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



