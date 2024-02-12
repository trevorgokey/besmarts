Codecs
------

Parsing SMILES to build a complete description is complicated. Many variables
go into perceiving how a molecule should be constructed using SMILES, as most
but not all information needed is given. Because of this, a SMILES string must
be *decoded*, and most implementations agree but there are many cases where
they do not. Rather than make another SMILES decoder, we instead require that
an existing implementation is used, such as RDKit. BESMARTS considers RDKit to
be a SMILES and SMARTS decoder which in turn produces BESMARTS graphs.
Once a pattern has been decoded, it no longer requires any interaction with the
codec. This decouples BESMARTS from the way graphs are built, and makes it easy
to implement additional codecs.

SMILES and SMARTS both represent a graph, and each piece of information about
the node/edge is described by primitives. These primitives also need their own
codecs. Thus, a graph codec is responsible for generating the graphs, but also
decoded the primitives that the pattern contains. As an example, a simple 
passthrough type of codec can be used to represent the atomic number, where 
'#6' is encoded as 6. However, for data that is not positive integers such
as formal charge, a custom codec is needed. In this case, a scheme to enumerate
all integers to the positive integers is used. The decoders of these primitives
are then used to determine how to represent the data in SMILES or SMARTS format.

