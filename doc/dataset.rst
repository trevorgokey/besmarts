
Datasets
========


Expanding one SMARTS to many SMILES
-----------------------------------

The last useful task that this package performs is related to dataset design.
Using a set of rules, one can take a SMARTS pattern and filter all satisfactory
graphs. This ultimately produces a list of SMILES from a single graph:

.. code-block:: python

    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core.codecs import primitive_key as p
    from besmarts.core import graphs
    
    from besmarts.resolve.enumeration import enumerate_smiles
    from besmarts.resolve.rulesets import visitor_ruleset
    
    atom = (p.ELEMENT, p.HYDROGEN, p.FORMAL_CHARGE, p.VALENCE)
    bond = (p.BOND_ORDER,)
    gcd = graph_codec_rdkit(atom, bond)
    
    gcd.smiles_config.strip_hydrogen = True
    
    rulesets = [visitor_ruleset(gcd.primitive_codecs)]
    
    for l, smi in enumerate(("C1CCCCC1",)):
    
        b = gcd.smiles_decode(smi)
        print(gcd.smiles_encode(b))
    
        graphs.graph_fill(b)
    
        for chem in b.nodes.values():
            chem.primitives[p.ELEMENT][:] = False
            chem.primitives[p.ELEMENT][6] = True
            # chem.primitives[p.ELEMENT][7] = True
            # chem.primitives[p.ELEMENT][8] = True
            # chem.primitives[p.ELEMENT][9] = True
            chem.primitives[p.FORMAL_CHARGE][:] = False
            chem.primitives[p.FORMAL_CHARGE][0] = True
            chem.primitives[p.VALENCE][:] = False
            # chem.primitives[p.VALENCE][1:5] = True
            chem.primitives[p.VALENCE][4] = True
        for chem in b.edges.values():
            chem.primitives[p.BOND_ORDER][:] = False
            chem.primitives[p.BOND_ORDER][1:4] = True
    
    
        mols = []
        smis = set()
        for i, P in enumerate(enumerate_smiles(b, rulesets)):
            if P not in mols:
                mols.append(P)
                smi = gcd.smiles_encode(P)
                if smi not in smis:
                    smis.add(smi)
                    print(f"{len(smis):4d} {smi}")


Here we see that we start from the graph of benzene and then fill the entire 
graph with wildcards. We then set the graph such that only 0-charge carbon with
a valence of 4 are allowed, and then permit bond orders from 1 to 3 (single, 
double, and triple bonds). The ruleset with then ensure that all atoms and bonds
are valid and will produce a meaningful SMILES pattern.
