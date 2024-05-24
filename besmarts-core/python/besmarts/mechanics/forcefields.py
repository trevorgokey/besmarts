
"""
besmarts.mechanics.forcefields
"""

from besmarts.mechanics import fits

def forcefield_optimize(
    gcd: codecs.graph_codec,
    labeler: assignments.smarts_hierarchy_assignment,
    gdb: assignments.graph_db,
    tiers: List[objective_tier],
    strategies: List[optimization.optimization_strategy],
    initial_conditions: chemical_system,
) -> chemical_system:
    """
    the macro step will find the nodes
        need a way to iterate over the models
        try with a strategy for each... will need a api layer on top

    now we will generate all splits and then pass them through the filters
    then we accept the operations
    
    for tier in tiers:
        scores.append(tier.compute(gdb, csys) for csys in candidates)

    """

    candidates = []
    keys = [gen(csys, strategy) for csys in candidates]
    psysref = []
    wq = None
    for tier, strategy in zip(tiers, strategies):
        results = fits.objective_tier_take(
            tier,
            gdb,
            keys,
            candidates,
            strategy,
            psysref=psysref,
            wq=wq
        )
        # asdads
        candidates = results.values()
