.. besmarts documentation master file, created by
   sphinx-quickstart on Sat Dec  3 18:54:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to besmarts's documentation!
====================================


Description
===========

This package was written to handle binary-encoded SMARTS (BESMARTS) with an
emphasis on handling aspects which are needed for molecular mechanics (MM)
force field (FF) parameterization. Standard MM force fields apply forces to
atoms, bonds, angles, and dihedrals (both "torsions" and "out-of-planes"), termed
*structures*. BESMARTS considers all SMARTS, SMILES, and everywhere in between
as graphs, and *structures* are *subgraphs* of some underlying *graph*. Much
effort has gone into manipulating SMARTS patterns which describe *structures*.
The SMIRNOFF force field specification uses SMARTS structures to assign
parameters to molecules, and a major difficulty using this specification is how
to programmatically create new parameters. This work stemmed from the need to
manipulate SMARTS patterns in fine detail, primitive at a time. We hope that
you find this work useful!

Getting Started
===============


.. toctree::
   :maxdepth: 2

   installation.rst
   usage.rst
   examples.rst
   contents.rst

Important considerations
========================

This package only implements and processes a subset of SMARTS for encoding.
Recursive SMARTS are not handled, and all SMARTS primitives are assumed to be
independent. This means that it is not possible to encode something like 
``[#6H3,#7H2]``, but it is possible to encode ``[#6,#7;H2,H3]``.
The former case is harder to reason with using functionality such as taking the
union and other mapping tasks, and has been put to the side for now.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
