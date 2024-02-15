"""
besmarts.mm.openmm
"""

from typing import Dict
from openmm.app.topology import Topology

from openmm.app.topology import Topology
from openmm.app.element import Element

from besmarts.core import mm


def physical_system_to_openmm_system(csys, psys: mm.physical_system):

    # make the topology
    otop = Topology()
    for n in psys.positions.graph.nodes
        otop.addAtom(str(n), 
    # then make the forces

    # then make the system
    

class physical_model_openmm(forces.physical_model):
    
    __slots__ = "masses", "forces", "energies", "positions", "internals", "opts"


    def __init__(masses, forces, energies):

        self.masses: List[forces.physical_model] = []
        self.forces: List[assignments.graph_assignment] = []
        self.energies: Dict[str, List[assignments.graph_assignment]] = []
        self.positions: List[assignments.graph_assignment] = None

        self.opts = {"get_force": False, "get_energy": False, "get_pos": True, "mapping": {}}
        self.systems = List[openmm.app.simulation.Simulation]

        # selections is the actual graph indices in the system
        # self.forces["Bonds"][0].selections[(1,2)] = {"function": force_function_bond_harmonic, "k": ((1, 1, 1), None), "l": ((1, 1, 1), None)}
        # self.energies["Bonds"][0].selections[(1,2)] = {"function": energy_function_bond_harmonic, "k": ((1, 1, 1), None), "l": ((1, 1, 1), None)}

        # cm.evaluate(self.energies["Bonds"][0].selections[(1,2)], pos)

    def build(self):

        integrator = openmm.VerletIntegrator(0.1 * simtk.unit.femtoseconds)
        sim = openmm.app.simulation.Simulation(top, system, integrator)
        sim.context.setPositions(xyz * const.angstrom2nm)
        xml = openmm.XmlSerializer.serialize(system)
        
    def evaluate_energy(self, pos) -> List[assignments.graph_assignments]:
        pass

    def evaluate_force(self, pos) -> List[assignments.graph_assignments]:
        pass

    def evaluate_momentum(self, vel) -> List[assignments.graph_assignments]:
        pass

    def minimize(self, pos) -> List[assignments.graph_assignments]:
        return minimization_openmm(self, pos)

    def integrate(self, pos, dt, steps) -> List[assignments.graph_assignments]:
        pass

def energy_openmm(cm_list: List[forces.chemical_model], pos: List[assignments.smiles_assignment_group]) -> float:
    # create the system
    # add the forces
    # get the energy
    # return
    pass

def energy_smirnoff(fname, pos: List[assignments.smiles_assignment_group]) -> float:
    pass

def minimization_openmm(pm: forces.physical_model) -> systems.system_state:
    # create the system, minimize, return new state

    # assume vdW section
    vdw_f = pm.forces['vdW']
    q_f = pm.forces['partial_charge']

    for vdw_g, pc_g in zip(vdw_f, q_f):
        for idx, params in vdw_g.selections.items():

            q = pc_g.selections[idx].partial_charge
            eps = params.eps
            sigma = params.sigma
            mm.create_ptl(idx, eps, sigma, q)

    mm.minimize()

    # need an openmm topology too?!?

    if self.minimize:
        sim.minimizeEnergy()
    state = sim.context.getState(getEnergy=True, getPositions=getPositions)
    energy = state.getPotentialEnergy().in_units_of(
        simtk.unit.kilocalories_per_mole
    )
    pos = None
    if self.minimize:
        pos = state.getPositions(asNumpy=True) / simtk.unit.angstroms


def physical_model_procedure_init_openmm_system(psys, config) -> Dict[str, topology_assignment]:
    """
    creates a system and stores it as a null topology assignment
    """
    pass

def physical_model_procedure_single_point_openmm(psys, config):
    """
    will grab the precomputed openmm system and then calculate the next props
    """
    pass
