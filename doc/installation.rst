Installation
============

BESMARTS is on PyPI, but certain dependencies are not and only exist on conda
unless they are built from source. Currently, the best way to install all
functionality is to use conda to create a python environment with Ambertools:

.. code-block:: bash

     conda create -n besmarts -c conda-forge python ambertools 
     conda activate besmarts

and then install BESMARTS via pip:

.. code-block:: bash

    pip install 'besmarts[rdkit,scipy,openmm]'

For a more custom installation, one can install from git. If ambertools is not
needed, one may set up an environment via

.. code-block:: bash

    python -m venv besmarts
    . besmarts/bin/activate

followed by the actual install in develop mode:

.. code-block:: bash

    git clone https://github.com/trevorgokey/besmarts besmarts-git
    cd besmarts-git/besmarts-core/python
    python -m pip install -e .
    cd ../../besmarts-rdkit/python
    python -m pip install -e .
    cd ../../besmarts-mechanics/python
    python -m pip install -e .
    cd ../../besmarts-openmm/python
    python -m pip install -e .
    cd ../../besmarts-scipy/python
    python -m pip install -e .

RDKit is needed to decode SMILES into graphs and offers a faster implementation
of SMARTS matching when labeling from a SMARTS hierarchy.

Geometry optimization uses the SciPy minimizer and can be installed using
a similar process as above with `besmarts-scipy`. There is also an
interface to OpenMM and its molecule energy minimizer can be used instead after installing
`besmarts-openmm`. The OpenMM plugin is quite a bit faster and is recommended
if large or heavily numerical computations are needed. The native interface
to calculating energies and gradients is useful if novel functional forms are
needed and not in standard packages (e.g. OpenMM). We recommend using OpenMM if
it is supported and available on your system. The energies, hessians, and
gradients compared between the native and OpenMM implementations are nearly
exact; for energy/gradient down to 12 decimal places for linear terms, around 6
places for torsions, and between 4-12 places for Hessians. Included in the native
implementation is a very fast analytic MM hessian method that can produce the
entire matrix in the time it takes OpenMM to evalulate a single energy.

Molecular mechanics energy and gradient evaluations are implemented, but
require partial charges. By default, `besmarts` will try to charge molecules
with `am1bcc` using the `sqm` program from `ambertools` suite. Consequently,
make sure `sqm` is in your `PATH` by installing via `conda` or by other means.
