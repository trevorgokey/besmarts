Maximum common substructure (MCS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python


    from besmarts.core import graphs
    from besmarts.core import mapper
    from besmarts.core import configs
    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    from besmarts.core import primitives
    
    
    def MCS(g, h):
    
        g_atoms = graphs.graph_to_structure_atoms(g)
        h_atoms = graphs.graph_to_structure_atoms(h)
    
        intersect = configs.mapper_config(False, False, "low")
    
        mcs = {}
        pruned = {gi: [] for gi,_ in enumerate(g_atoms)}
    
        added = True
        i = -1
        while added:

            added = False
            i += 1
            config = configs.smarts_extender_config(i, i, False)
            modified = mapper.mapper_smarts_extend(config, [*g_atoms, *h_atoms])

            if not modified and i > 0:
                break

            for gi, ga in enumerate(g_atoms):
                for hi, ha in enumerate(h_atoms):

                    if hi in pruned[gi]:
                        continue
    
                    r = mapper.intersection(ga, ha, intersect)
                    r = graphs.structure_prune_null_nodes(r)
    
                    if r is None:
                        pruned[gi].append(hi)
                        continue
    
                    n = len(r.nodes)
                    if n not in mcs:
                        mcs[n] = []

                    mcs[n].append((gi, hi, r))
                    added = True
        
        return mcs
    
    gcd = graph_codec_rdkit()
    gcd.atom_primitives = tuple(x for x in gcd.atom_primitives if x != primitives.primitive_key.HYDROGEN)
    gcd.smiles_config.strip_hydrogen = True
    
    g = gcd.smiles_decode("C[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O)N") # trialanine
    h = gcd.smiles_decode("C(C(=O)NCC(=O)NCC(=O)O)N") #triglycine
    
    mcs = MCS(g, h)
    if mcs:
        for (gi, hi, r) in mcs[max(mcs)]:
            print(gcd.smarts_encode(r))
            break
