# sphinx-apidoc --implicit-namespaces -f -e -o . ../besmarts ../besmarts/backends/io ../besmarts/backends/graph
# sphinx-apidoc --implicit-namespaces -f -e -o . ../besmarts/backends/graph/besmarts_networkx/
# sphinx-apidoc --implicit-namespaces -f -e -o . ../besmarts/backends/io/besmarts_rdkit

sphinx-apidoc --implicit-namespaces -e -f  -o besmarts ~/projects/besmarts/besmarts-core/python/besmarts/
# sphinx-apidoc --implicit-namespaces -e -f  -o besmarts-hierarchy ~/projects/besmarts/besmarts-hierarchy/besmarts/
sphinx-apidoc --implicit-namespaces -e -f  -o besmarts-rdkit ~/projects/besmarts/besmarts-rdkit/besmarts/
# sphinx-apidoc --implicit-namespaces -e -f  -o besmarts-resolve ~/projects/besmarts/besmarts-resolve/besmarts/
# sphinx-apidoc --implicit-namespaces -e -f  -o besmarts-splitter ~/projects/besmarts/besmarts-splitter/besmarts/

