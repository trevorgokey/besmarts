
"""
besmarts.mm.systems
"""

from typing import List

from besmarts.mm import forcefields, forces
from besmarts.core import assignments
from besmarts.core import topology


"""
chemical model defines how to parameterize
the physical model contains the parameterized system
the system contains both the chemical_model and physical model

pm = pm(cm)
graph =      transcode_graph(transcoders, graph)
assignment = transcode_positions(transcoders, assignment)

"""
class system_state:
    def __init__(self, smiles, graph):
        # these are the actual numbers, not the functions
        smiles: List[str] = smiles
        graph: List[graphs.graph] = graph
        self.positions: Dict[int, List[assignments.structure_assignment_float]] = {}
        self.velocities: Dict[int, List[assignments.structure_assignment_float]] = {}
        self.forces: Dict[str, List[assignments.structure_assignment_float]] = {}
        self.energies: Dict[str, List[assignments.structure_assignment_float]] = {}

        # single values as a function of entire system
        self.values = {}

        # single values as a function of vector
        self.points = {}

        # basis vectors
        self.spaces = {}

        # vector values as a function of vector
        self.fields = {}

class system_component:
    """
    this is a full description of the potential on a graph
    """

    def __init__(self,
        cm: forcefields.chemical_model,
        pm: forcefields.physical_model,
        state: system_state
    ):
        """
        """
        pass

    def evaluate_energy(self, pos) -> List[assignments.structure_assignment_group]: pass

    def evaluate_force(self, pos) -> List[assignments.structure_assignment_group]: pass

    def evaluate_momentum(self, vel) -> List[assignments.structure_assignment_group]: pass

class system:
    """

    """
    def __init__(self,
        components: Dict[str, system_component],
    ):
        # this is the FF that builds a physical model, essentially the FF
        self.chemical_model = chemical_model

        # this is the parameterized graph, essentially after applying the FF
        # references to the cm for the values
        self.physical_model = physical_model

        self.state = state



class system_nve(system):
    def __init__(self):
        self.box = None

class system_nvt(system):
    def __init__(self):
        self.box = None
        self.thermostat = ""

class system_npt(system):
    def __init__(self):
        self.box = None
        self.thermostat = ""
        self.barostat = ""



