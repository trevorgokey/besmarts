
# BESMARTS

*A toolkit for force field design based on binary-encoded SMARTS*

# Installation

Currently, the best way to install is to clone and then run

```
git clone https://github.com/trevorgokey/besmarts
cd besmarts/besmarts-core/python
python -m pip install .
cd ../../besmarts-rdkit/python
python -m pip install .
```

RDKit is needed to decode SMILES into graphs and offers a faster implementation
of SMARTS matching when labeling from a SMARTS hierarchy.


