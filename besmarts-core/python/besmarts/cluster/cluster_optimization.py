"""
besmarts.cluster.cluster_optimization

Optimization strategies to form cluster hierarchies.
"""

from typing import List, Callable

from besmarts.core import (
    clusters,
    configs,
    codecs,
    assignments,
    optimization,
    tree_iterators,
)

from besmarts.cluster import cluster_objective


class optimization_iteration_default(optimization.optimization_iteration):
    def __init__(self, steps):
        self.cursor = 0
        self.steps: List[optimization.optimization_step] = steps
        self.repeat = False

    def is_done(self):
        return self.cursor >= len(self.steps)

    def next(self):
        return optimization.optimization_iteration_is_done(self)

    def repeat_step(self):
        return optimization.optimization_iteration_repeat_step(self)


class optimization_strategy_default(optimization.optimization_strategy):

    """
    determines how to step the optimization forward, choosing which hyperparameter
    to try next
    """

    def __init__(
        self,
        bounds: configs.smarts_perception_config = None,
        tree_iterator: Callable = None,
    ):


        if bounds is None:
            splitter = configs.smarts_splitter_config(
                1, 3, 0, 3, 0, 2, True, True, 100, True, True
            )
            extender = configs.smarts_extender_config(
                0, 2, True
            )
            bounds = configs.smarts_perception_config(splitter, extender)
        super().__init__(bounds)

        if tree_iterator is None:
            self.tree_iterator = tree_iterators.tree_iter_dive

        self.bounds = bounds
        self.cursor = -1
        self.maxedits_limit = 0
        self.repeat = False
        self.steps = []
        self.tree_iterator: Callable = tree_iterators.tree_iter_dive


def cluster_classifications(
    gcd: codecs.graph_codec,
    labeler: assignments.smarts_hierarchy_assignment,
    sag: assignments.smiles_assignment_group,
    objective: clusters.clustering_objective = None,
    optimization: optimization.optimization_strategy = None,
    initial_conditions: clusters.smarts_clustering = None,
):

    if objective is None:
        objective = cluster_objective.clustering_objective_classification()

    if optimization is None:
        max_depth = clusters.smarts_clustering_find_max_depth(
            assignments.smiles_assignment_group_to_structure_assignment_group(
                sag, gcd
            ),
            8,
        )
        splitter = configs.smarts_splitter_config(
            1, 3, 0, 3, 0, max_depth, True, True, 100, True, True
        )
        extender = configs.smarts_extender_config(
            0, max_depth, True
        )
        cfg = configs.smarts_perception_config(splitter, extender)
        optimization = optimization_strategy_default(cfg)

    if initial_conditions is None:
        initial_conditions = clusters.clustering_initial_conditions(gcd, sag)

    clst = clusters.smarts_clustering_optimize(
        gcd, labeler, sag, objective, optimization, initial_conditions
    )

    groups = clusters.clustering_build_ordinal_mappings(clst, sag)
    problems = []
    correct = 0
    params = set()
    for albl, blbls in groups.items():
        x = set(tuple(blbls))
        print(albl, x)
        if len(x) > 1:
            problems.append((albl, list(x)))
        else:
            correct += 1
    if problems:
        print("PROBLEMS:", problems)
    print("ACCURACY:", correct/len(groups))

    return clst


def cluster_means(
    gcd: codecs.graph_codec,
    labeler: assignments.smarts_hierarchy_assignment,
    sag: assignments.smiles_assignment_group,
    objective: clusters.clustering_objective = None,
    optimization: optimization.optimization_strategy = None,
    initial_conditions: clusters.smarts_clustering = None,
) -> clusters.smarts_clustering:

    if objective is None:
        objective = cluster_objective.clustering_objective_mean_separation()

    if optimization is None:
        splitter = configs.smarts_splitter_config(
            1, 2, 0, 0, 0, 0, True, True, 0, False, True, True, True
        )
        extender = configs.smarts_extender_config(
            0, 0, True
        )
        cfg = configs.smarts_perception_config(splitter, extender)
        optimization = optimization_strategy_default(cfg)

    if initial_conditions is None:
        initial_conditions = clusters.clustering_initial_conditions(gcd, sag)

    clst = clusters.smarts_clustering_optimize(
        gcd, labeler, sag, objective, optimization, initial_conditions
    )

    return clst
