opts="-M --ext-autodoc --implicit-namespaces -d 3 --remove-old -e -f"
sphinx-apidoc $opts -o besmarts        ../besmarts-core/python/besmarts
sphinx-apidoc $opts -o besmarts-rdkit  ../besmarts-rdkit/python/besmarts
sphinx-apidoc $opts -o besmarts-scipy  ../besmarts-scipy/python/besmarts
sphinx-apidoc $opts -o besmarts-openmm ../besmarts-scipy/python/besmarts
