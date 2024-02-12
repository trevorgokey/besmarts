Description
===========

Goal: Publish the binary encoded SMARTS (besmarts) package

Use cases:

Main functionality
==================

Perception using RDKIT
----------------------
	- Name: besmarts_rdkit

	1. Perceive the SMARTS information to build the core model 

Data definition
---------------
	- Name: besmarts_core

	1. Generate fragments from molecules (atoms/bonds/etc)
		- atoms
		- bonds
		- angles
		- dihedrals
		- impropers
		- chains
		- rings

	2. Aligning two structures
		- Equality
		- Subset/Superset
		- Best using overlap
		- Align using reference
		- Perform mapping with bit overlap
		- Enumerate valid mappings

	2. Perform operations on two environments depending on mapper
		- union (or; |)
		- intersection (and; &)
		- xor (^)
		- difference (-)

	3. Generate a unique SMARTS for each environment (uniquify_fragments)
		- mapper_extend_uniquify
		- mapper_extend_extend


Specialized functionality
=========================


Splitting and partitioning
--------------------------
	- Name: besmarts-splitter

	1. Find the SMARTS partition between two groups (partition)
		- extends subgraph as needed
		* find closest parition with number of max moves

	2. Find the partitions from iterating shards
		- can automatically extend branching

Enumerate SMILES from a SMARTS pattern
--------------------------------------
	- Name: besmarts-resolve
	
	1. Use rulesets to determine how to build molecules

Hierarchy
---------
	- Name: besmarts-harare

	1. Tree data structure

	1. Build a hierarchy of SMARTS

	2. Combine 2 hierarchies


Getting Started
===============

Installation
------------
