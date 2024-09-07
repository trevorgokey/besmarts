"""
besmarts.tests.test_import

Import all files. Mostly designed as a test to list every object in the codebase
but can also detect circular imports and check if pypy3 can run everything
"""

import unittest

class test_import(unittest.TestCase):
    def test_import_core(self):
        from besmarts.core import (
            arrays,
            assignments,
            chem,
            clusters,
            codecs,
            compute,
            configs,
            geometry,
            graph_visitors,
            graphs,
            hierarchies,
            hierarchy_merge,
            mapper,
            optimization,
            primitives,
            returns,
            splits,
            topology,
            tree_iterators,
            trees,
        )

    def test_import_codecs(self):
        from besmarts.codecs import (
            codec_native,
        )

    def test_import_cluster(self):
        from besmarts.cluster import (
            cluster_assignment,
            cluster_objective,
            cluster_optimization,
        )

    def test_import_assign(self):
        from besmarts.assign import (
            hierarchy_assign_native,
        )

    def test_import_mechanics(self):
        from besmarts.mechanics import (
            force_harmonic,
            force_pairwise,
            force_periodic,
            masses,
            smirnoff_models,
            smirnoff_xml,
        )

