
Installation
============

Currently, the best way to install is to clone and then install with pip.

For environment users (e.g. venv or conda), one should probably create an empty
environment first:

.. code-block:: bash

    conda create -n besmarts python
    conda activate besmarts

or (if not installing ambertools)

.. code-block:: bash

    python -m venv besmarts
    . besmarts/bin/activate


followed by the actual install:

.. code-block:: bash

    git clone https://github.com/trevorgokey/besmarts besmarts-git
    cd besmarts-git/besmarts-core/python
    python -m pip install .
    cd ../../besmarts-rdkit/python
    python -m pip install .

RDKit is needed to decode SMILES into graphs and offers a faster implementation
of SMARTS matching when labeling from a SMARTS hierarchy.

Geometry optimization in addition to force field optimization uses the SciPy
minimizer and can be installed using using a similar process as above with
`besmarts-scipy`:

.. code-block:: bash

    cd besmarts-git/besmarts-scipy/python
    python -m pip install .

Molecular mechanics energy and gradient evaluations are implemented, but
require partial charges. By default, `besmarts` will try to charge molecules
with `am1bcc` using the `sqm` program from `ambertools` suite. Consequently,
make sure `sqm` is in your `PATH` by installing via `conda` (`conda install
-c conda-forge ambertools`) or by other means.


