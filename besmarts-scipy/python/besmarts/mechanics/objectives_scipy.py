
from besmarts.mechanics import objectives 
from besmarts.mechanics import optimizers_scipy

class physical_objective_bond_deviation_scipy(
        objectives.physical_objective_bond_deviation):
    def __init__(self):
        self.scale = 1.0
        self.keys = []
        self.get_task_physical_system = None

    def get_task_compute(self, csys, psys):
        return optimizers_scipy.optimize_positions_scipy
