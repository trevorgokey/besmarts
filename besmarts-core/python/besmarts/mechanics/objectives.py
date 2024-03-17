
"""
besmarts.mechanics.objectives
"""

import copy
from besmarts.core import assignments
from besmarts.mechanics import molecular_models as mm

# def residual_squared_error(x1, x0):

#     result = {}
#     for k, r1 in x1.selections.items():
#         r0 = x0.selections[k]
#         result[k] = (r1 - r0)**2
#     return result

# def measure_rse(new_sys: system, ref_sys: system, measure):
#     smi = ref_sys.positions.smiles
#     graph = ref_sys.positions.graph

#     x0 = measure(smi, graph, ref_sys.positions.selections)
#     x1 = measure(smi, graph, new_sys.positions.selections)

#     result = residual_squared_error(x1, x0)

#     return assignments.graph_assignment(smi, result, g)

# def bond_rse(new_sys: system, ref_sys: system):
#     result = measure_rse(new_sys, ref_sys, assignments.graph_assignment_geometry_bonds)
#     return result

# class chemical_objective:
#     def __init__(self, cm: chemical_model):
#         self.total = 0
#         self.results = Dict[str, Dict[str, float]]

#     def create_compute_tasks(cm):
#         return default_chem_objective



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
# smiles_state = Dict[str, topology_assignment]

sst.smiles
sst.graph
sst.assignments[POSITIONS].selections[(1,)] = [[0,0,0]]
sst.assignments[GRADIENTS].selections[(1,)] = [[0,0,0]]
sst.assignments[HESSIANS].selections[(1,1)] = [[[0,0,0],[0,0,0],[0,0,0]]]

[smiles, assignments, selections]
(0, POSITIONS, (1,), 0): [x, y, z]

gid = db_graph_assignment_add_graph(db, g)
aid = db_graph_assignment_add_(g)

gid, aid, sid, vid (eid)

topology_assignment_{str, float, unit, comment}

xyz = smiles_dataset_get_values((0, POSITIONS, (1,), 0))

db_graph_assignment:
    graphs: list[graph_assignment]
    assignments = {str: topology_assignment_terms}

db_graph_assignment.graphs[0]
db_graph_assignment.assignments[POSITIONS][0][(1,)] = [[x,y,z]]

db_graph_topology_assignment
    topology
    selections[gid][sid]

CLASSES:
    db_graph_assignment
    db_graph_topology_assignment


POSITIONS
0, 1, 0 x y z
0, 1, 1 x y z
0, 1, 1 x y z
GRADIENTS
0, 1, 0 x y z
0, 1, 1 x y z
0, 1, 1 x y z
HESSIANS
0 1 1 0 xx xy xz yx yy yz zx zy zz

TORSIONS
0 1 2 3 4 0 a
0 1 2 3 5 0 a

ANGLES
0 1 2 3 0 a
0 1 2 4 0 a

# 
ELEMENT
0 1 8

select x, y, z from positions where id=0 and index=0 selection=1

topology_assignment_terms
    name
    units
    selections = (0, (1,)): [x,y,z]
    comment

graph_assignment
    graph
    dict[str: topology_assignment]

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

# class physical_objective:
#     def __init__(self, indices):
#         self.indices: List[int] = indices
#         self.total = 0.0
#         self.results: Dict[str, List[assignments.structure_assignment_float]] = {}

#     def create_ref_tasks(self, cm, states):
#         # assume reference is constant
#         # in this case the state is the result, don't send!
#         return []

#     def create_fit_tasks(self, cm, states):
#         return []

#     def create_compute_tasks(self, new_sys, ref_sys):
#         return []


# class physical_objective_qm_optgeo(physical_objective):
#     """
#     This will take the reference and do nothing to it
#     """
#     def create_fit_tasks(cm, state):
#         return [minimization, (cm, state)]

#     def create_compute_tasks(self, new_sys, ref_sys):
#         return [bond_rse]


# class physical_objective_qm_energy(physical_objective):

#     def create_fit_tasks(self, cm, state):
#         return [(esp_calc, {})]

#     def create_compute_tasks(self, new_sys, ref_sys):
#         return [(esp_rse, {})]

# class physical_objective_qm_force(physical_objective):

#     def create_fit_tasks(self, cm, state):
#         return [(esp_calc, {})]

#     def create_compute_tasks(self, new_sys, ref_sys):
#         return [(esp_rse, {})]
    

# class qm_physical_objective_const_resp(physical_objective):
#     """
#     This will take the reference and calc esp, assumes esp will never change
#     """
#     def create_fit_tasks(self, cm, state):
#         return [(esp_calc, {})]

#     def create_compute_tasks(self, new_sys, ref_sys):
#         return [(esp_rse, {})]


# class qm_physical_objective_dynam_resp(physical_objective):
#     """
#     This will take the reference and calc esp using extra points defined
#     by the chemical model. Useful if we want to fit virtual sites and want
#     a reference charge in the qm where the points are decided by the FF
#     """
#     def create_fit_tasks(self, cm, states):
#         return [(esp_calc, {})]

#     def create_compute_tasks(self, new_sys, ref_sys):
#         return [(esp_rse, {})]

# def forcefield_optimize(
#     ref_sys: List[system],
#     cobjective: chemical_objective,
#     pobjective: physical_objective,
#     strategy: optimization.optimization_strategy,
#     initial_conditions: forcefields.forcefield,
# ) -> forcefields.forcefield:
#     """
#     Run the ref_sys with the initial chemical model to get a set of new sys,
#     then give the two sys to the objective and compute

#     then proceed "as normal"

#     Note.. we can do this with bond lengths because we can take the positions,
#     calc the bonds and return the average per label
#     """

# def ff_opt(
#     csys: chemical_system,
#     ref_ssg: List[assignments.smiles_state],
#     pobjective: physical_objective,
#     cobjective: chemical_objective,
#     strategy: optimization.optimization_strategy,
# ) -> chemical_system:

#     """
#     need a method to generate the candidates and thats it
#     otherwise, i need a method that will go over all models and genreate splits
#     it will then reoptimize and generate the new ssg

#     optimization strategy will need ot be modified to select procedures

#     we will obviously

#     psys = csys(ref_ssg)
#     which labels everything

#     then 
#     tgt_ssg = calcs(psys)

#     then
#     X2 = objective(ref_ssg, tgt_ssg)
#     """

# def objective_hessian_create_tasks(p, config):
#     return (compute_freq_modes_openmm, config, p)

# def objective_energy(p, config):
#     config.energy = True
#     return (compute_single_point, config, p)

# def objective_force(p, config):
#     config.assignments.add(FORCES)
#     return (compute_single_point, config, p)

# result = fn(confg, p)
# forces = result.assignments[FORCES]


# def objective_hessian(reference: psys, predicted: psys):
#     h0 = reference.assignments["hessian"]
#     predicted.assignments["hessian"] = calc_hess(psys)
#     h1 = predicted.assignments["hessian"]

    
# class objective_function_total_energy:
#     def __init__(self, csys, states, objectives):
#         """
#         I would need to iterate all keys, then filter only those that we applied.
#         Optionally, accept an additional filter that pulls even more.

#         """
#         self.csys = csys

#         # the reference coordinates, topo is atom
#         self.pos: assignments.graph_topology_db_table = None

#         # the reference energy, topo is undef
#         self.reference: assignments.graph_topology_db_table = None

#         # compute will need all energy tables
#         # should these be tables, and a lookup will map compute to tables
#         self.tasks = ["bonds", "angles", "torsions"]


#####
# vectorized forms
#####
def array_flatten_assignment(selections):
    lst = []
    keys = []
    n_confs = len(list(selections.values())[0])
    for c in range(n_confs):
        for n, data in selections.items():
            lst.extend(data[c])
            keys.extend(((c, n, i) for i in range(len(data[c]))))
    return lst, keys

def array_geom_energy(args, keys, csys, psys: mm.physical_system):
    energy = 0

    for m, pm in enumerate(psys.models):
        refpos = psys.models[m].positions[0]
        pos = copy.deepcopy(refpos)
        # pos = assignments.graph_assignment_float(
        #         refpos.graph,
        #         {k: v.copy() for k, v in refpos.selections.items()}
        # )

        n_confs = len(list(refpos.selections.values())[0])
        i = 0
        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v


        ic = csys.models[m].internal_function(pos)

        # build system terms
        system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}
        params = pm.values

        # build topology terms
        ene = mm.smiles_assignment_function(csys.models[m].energy_function, system_terms, params, ic)
        energy += sum([x for y in ene.values() for z in y for x in z])

    return energy

def array_geom_gradient(args, keys, csys, psys: mm.physical_system):
    return [-x for x in array_geom_force(args, keys, csys, psys)]

def array_geom_force(args, keys, csys, psys: mm.physical_system):
    """
    csys is for the reference functions and system params
    psys is for the positions and param values
    """
    force = list([0.0]*len(args))

    for m, pm in enumerate(psys.models):

        refpos = pm.positions[0]
        pos = assignments.graph_assignment_float(
                refpos.graph,
                {k: v.copy() for k, v in refpos.selections.items()}
        )

        n_confs = len(list(refpos.selections.values())[0])

        i = 0

        for (c, n, i), v in zip(keys, args):
            pos.selections[n][c][i] = v


        if csys.models[m].derivative_function:

            icq = csys.models[m].internal_function(pos)

            system_terms = {k: v.values for k, v in csys.models[m].system_terms.items()}

            params = pm.values
            f = mm.smiles_assignment_function(csys.models[m].force_function, system_terms, params, icq)

            jac = csys.models[m].derivative_function(pos)

            for ic in f:
                confs_dq = jac.selections[ic]
                nic = len(confs_dq)-1
                for idx, dq in enumerate(confs_dq):
                    fq = f[ic][idx][0]
                    for j, i in enumerate(ic):
                        force[nic*idx + (i-1)*3 + 0] += fq*dq[j][0]*4.184
                        force[nic*idx + (i-1)*3 + 1] += fq*dq[j][1]*4.184
                        force[nic*idx + (i-1)*3 + 2] += fq*dq[j][2]*4.184

    return force
