
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



# we need an objective from a pm that produces a sag
# lets say the sag is 3d, then we need:
# objective(psys, pos) -> sag
# the sig would then be sag which are pos
# we then also have an initial csys
# then we would go csys, sag0 -> psys then to obj(psys) -> value
"""
sag0 could be a state with pos and grads and etc
call it ref0
then we would need a calculator for each sag, or in other words a sag provider
csys2pys(ref, csys) -> pm -> obj(pm) -> sag

ref is positions(3d), energy (mol) and grads (3d) for example

csys is the initial ff

so then we would get the psys and then calc pos, grad, and energy

then the objective would take the ref and the psys and return a sag

have a psys, which has a form for the analytic energies and gradients

def modsem(pos, hes) -> sag:
    pass
    (then label, then mean)

    # may want some initializers/computers such as bond length and fc

# refs should be a state?
# smiles_state = Dict[str, structure_assignment]

sst.smiles
sst.graph
sst.assignments[POSITIONS].selections[(1,)] = [[0,0,0]]
sst.assignments[GRADIENTS].selections[(1,)] = [[0,0,0]]
sst.assignments[HESSIANS].selections[(1,1)] = [[[0,0,0],[0,0,0],[0,0,0]]]

# for elecs
sst.assignments[GRID].selections[None] = [[0,0,0,0,0,0,0,0,0,0,0,0]]
sst.assignments[ESP].selections[None] = [[0,0,0,0,0,0,0,0,0,0,0,0]]
# sst.assignments[RADII].selections[(1,)] = [[1.2]]


# sst.assignments["bonds"].selections[(1,)] = [[0]]
# sst.assignments["pairs"].selections[(1,)] = [[0,0,0]]
# sst.assignments["angles"].selections[(1,)] = [[0,0,0]]
# sst.assignments["torsions"].selections[(1,)] = [[0,0,0]]
# sst.assignments["outofplanes"].selections[(1,)] = [[0,0,0]]

def opt(csys, ref_states, objectives, config):


    list_of_energies: list[task] = objective_energy_calcs(refs, compute_config)
    list_of_forces = objective_force_calcs(refs, compute_config)

    task.refs: list of reference data
    task.config
    task.config.procedure = "single"
    task.config.energy = True
    task.config.force = True
    task.config.minimize = False
    task.config.constraints = None
    task.config.restraints = None

    jobs = set(loe, lof)
    compute(psys, config):

    psys = t(csys, pos) # now psys will have the complete desc for energy and force
    states: List[smiles_state]

    crds = states[0].assignments[POSITIONS]
    crds.topology
    crds.selections[(1,)]

    pos = smiles_assignment(states[0].smiles, crds)

    submit_jobs(jobs)

    objective_energy_calc(jobs)

    # gen the the p_models, done
    psys = [csys(ref_i) for ref_i in csys]

    # get the observables, pos, grad, ene
    # one task per psys?
    tgts = [compute(psys_i) for psys_i in psys]

    # then for each 
    o = [obj(psys_i, ref_i) for obj in csys]

    for sag in o:
        total += sum(sum(x) for x in sag.selections.values())

    we will then need to go through the optimization ops and proceed as usual

can i do a concerted objective?

def run(csys, pos, objectives):
    psys
    calc
    obj
    return



def objective(pm) -> sag:

"""

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

def ff_opt(
    csys: chemical_system,
    ref_ssg: List[assignments.smiles_state],
    pobjective: physical_objective,
    cobjective: chemical_objective,
    strategy: optimization.optimization_strategy,
) -> chemical_system:

    """
    need a method to generate the candidates and thats it
    otherwise, i need a method that will go over all models and genreate splits
    it will then reoptimize and generate the new ssg

    optimization strategy will need ot be modified to select procedures

    we will obviously

    psys = csys(ref_ssg)
    which labels everything

    then 
    tgt_ssg = calcs(psys)

    then
    X2 = objective(ref_ssg, tgt_ssg)
    """

def objective_hessian_create_tasks(p, config):
    return (compute_freq_modes_openmm, config, p)

def objective_energy(p, config):
    config.energy = True
    return (compute_single_point, config, p)

def objective_force(p, config):
    config.assignments.add(FORCES)
    return (compute_single_point, config, p)

result = fn(confg, p)
forces = result.assignments[FORCES]



def objective_hessian(reference: psys, predicted: psys):

    h0 = reference.assignments["hessian"]
    
    predicted.assignments["hessian"] = calc_hess(psys)
    h1 = predicted.assignments["hessian"]

    

