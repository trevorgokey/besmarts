
def residual_squared_error(x1, x0):

    result = {}
    for k, r1 in x1.selections.items():
        r0 = x0.selections[k]
        result[k] = (r1 - r0)**2
    return result

def measure_rse(new_sys: system, ref_sys: system, measure):
    smi = ref_sys.positions.smiles
    graph = ref_sys.positions.graph

    x0 = measure(smi, graph, ref_sys.positions.selections)
    x1 = measure(smi, graph, new_sys.positions.selections)

    result = residual_squared_error(x1, x0)

    return assignments.graph_assignment(smi, result, g)

def bond_rse(new_sys: system, ref_sys: system):
    result = measure_rse(new_sys, ref_sys, assignments.graph_assignment_geometry_bonds)
    return result

def default_chem_objective(cm):
    pass

class chemical_objective:
    def __init__(self, cm: chemical_model):
        self.total = 0
        self.results = Dict[str, Dict[str, float]]

    def create_compute_tasks(cm):
        return default_chem_objective

class physical_objective:
    def __init__(self, indices):
        self.indices: List[int] = indices
        self.total = 0.0
        self.results: Dict[str, List[assignments.structure_assignment_float]] = {}

    def create_ref_tasks(self, cm, states):
        # assume reference is constant
        # in this case the state is the result, don't send!
        return []

    def create_fit_tasks(self, cm, states):
        return []

    def create_compute_tasks(self, new_sys, ref_sys):
        return []

po = physical_objective([0])
tasks = po.create_ref_tasks(states)


class physical_objective_qm_optgeo(physical_objective):
    """
    This will take the reference and do nothing to it
    """
    def create_fit_tasks(cm, state):
        return [minimization, (cm, state)]

    def create_compute_tasks(self, new_sys, ref_sys):
        return [bond_rse]


class physical_objective_qm_energy(physical_objective):

    def create_fit_tasks(self, cm, state):
        return [(esp_calc, {})]

    def create_compute_tasks(self, new_sys, ref_sys):
        return [(esp_rse, {})]

class physical_objective_qm_force(physical_objective):

    def create_fit_tasks(self, cm, state):
        return [(esp_calc, {})]

    def create_compute_tasks(self, new_sys, ref_sys):
        return [(esp_rse, {})]
    

class qm_physical_objective_const_resp(physical_objective):
    """
    This will take the reference and calc esp, assumes esp will never change
    """
    def create_fit_tasks(self, cm, state):
        return [(esp_calc, {})]

    def create_compute_tasks(self, new_sys, ref_sys):
        return [(esp_rse, {})]

def esp_calc(pm, state) -> systems.state:
    pass

class qm_physical_objective_dynam_resp(physical_objective):
    """
    This will take the reference and calc esp using extra points defined
    by the chemical model. Useful if we want to fit virtual sites and want
    a reference charge in the qm where the points are decided by the FF
    """
    def create_fit_tasks(self, cm, states):
        return [(esp_calc, {})]

    def create_compute_tasks(self, new_sys, ref_sys):
        return [(esp_rse, {})]

def forcefield_optimize(
    ref_sys: List[system],
    cobjective: chemical_objective,
    pobjective: physical_objective,
    strategy: optimization.optimization_strategy,
    initial_conditions: forcefields.forcefield,
) -> forcefields.forcefield:
    """
    Run the ref_sys with the initial chemical model to get a set of new sys,
    then give the two sys to the objective and compute

    then proceed "as normal"

    Note.. we can do this with bond lengths because we can take the positions,
    calc the bonds and return the average per label
    """
