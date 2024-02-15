"""
besmarts.mechanics.molecular_models
"""

from typing import Dict, List, Any
import datetime

from besmarts.core import assignments
from besmarts.core import perception


class topology_term:
    def __init__(self, symbol, name, unit, cast, values, comment, value_comments):
        self.symbol: str = symbol
        self.name: str = name
        self.unit: str = unit
        self.cast: str = cast
        self.values: Dict[int, Any] = values
        self.comment: str = comment
        self.value_comments: Dict[int, str] = value_comments

    def copy():
        return topology_term_copy(self)

def topology_term_copy(t: topology_term):
    return topology_term(
        t.symbol,
        t.name,
        t.unit,
        t.cast,
        t.values.copy(),
        t.comment, 
        t.value_coments.copy()
    )

class system_term:
    def __init__(self, name, symbol, unit, cast, values, comment=""):
        self.name: str = name
        self.symbol: str = symbol
        self.unit: str = unit
        self.cast: str = cast
        self.values: List = values
        self.comment: str = comment

class physical_model:
    """
    The functional form that can be evaluated as a function of positions.
    Also tracks the positions as a cache mechanism to avoid recomputing
    Sort of assumes energies/forces are a function of internal coordinates
    and masses and positions are cartesian
    """

    __slots__ =  "labels", "values", "positions", "topology"


    def __init__(self, positions, labels, values):

        # this will have the current graph and ic positions
        self.positions: List[assignments.graph_assignment_float] = positions

        # the actual parameter values
        # keeps track of the values per procedure
        self.labels: List[Dict] = labels
        self.values: List[Dict] = values

class chemical_model_procedure:
    def __init__(self, name, topo):
        self.name = name
        self.topology = topo
        self.procedure_parameters: Dict[str, int] = {}
    
    def assign(self, pm: physical_model) -> physical_model:
        assert False

    def get_term_labels(self, key) -> Dict:
        assert False

    def get_term_values(self, key) -> Dict:
        assert False


class physical_model_procedure:
    """
    calculates one or more physical properties of a system
    """
    def __init__(self, name, topo):
        self.name = name # HESSIANS
        self.topology = topo
        self.procedure_parameters: Dict[str, int] = {}
    
    def assign(self, pm: physical_model):
        """
        this will return the compute function and config (i.e. the task)
        """

    def get_term_labels(self, key) -> Dict:
        assert False

    def get_term_values(self, key) -> Dict:
        assert False

class chemical_model:
    def __init__(self, symbol, name, topo):
        self.symbol = symbol
        self.name = name
        self.topology = topo

        self.topology_terms: Dict[str, topology_term] = {}
        self.system_terms: Dict[str, system_term] = {}

        self.procedures: List[chemical_model_procedure] = []

        # # the interaction functions
        self.force_func: Callable = None
        self.energy_func: Callable = None


class physical_system:
    def __init__(self, models: List[physical_model]):
        self.models = models

class chemical_system:
    def __init__(self, pcp_model: perception.perception_model, models: List[chemical_model]):
        self.perception = pcp_model
        self.models = models

class chemical_model_procedure_smarts_assignment(chemical_model_procedure):
    def __init__(self, pcp_model: perception.perception_model, topology_terms):
        self.name = ""
        self.perception = pcp_model

        self.topology_parameters: Dict[int, Dict[str, int|str]] = {}
        self.system_parameters: Dict[str, int] = {}

        self.topology_terms = topology_terms
        self.smarts_hierarchies: Dict[int, hierarchies.smarts_assignment_hierarchy] = {}

    def assign(self, cm, pm: physical_model) -> physical_model:
        """
        this will return, for each selection, the reference
        """

        # print(self.topology_terms)
        smiles = [x.smiles for x in pm.positions]
        topo = cm.topology
        unit_i = 0

        lbls = self.perception.labeler.assign(
            self.smarts_hierarchies[unit_i],
            self.perception.gcd,
            smiles,
            self.smarts_hierarchies[unit_i].topology
        )
        
        assn = []
        vals = []
        for x in lbls.assignments:
            p = {}
            v = {}
            for ic, lbl in x.selections.items():
                if lbl is None:
                    continue
                names = self.get_term_labels((unit_i, lbl))
                values = self.get_term_values((unit_i, lbl))
                p[ic] = {term: l for (term, l), x  in zip(names.items(), values.values())}
                v[ic] = {term: x for (term, l), x  in zip(names.items(), values.values())}
            assn.append(p)
            vals.append(v)

        pm.labels.extend(assn)
        pm.values.extend(vals)

        return pm
        

    def get_term_labels(self, k):
        unit_i, smarts_i = k
        terms = self.topology_parameters[(unit_i, smarts_i)]

        return terms

    def get_term_values(self, k):
        values = {}
        terms = self.get_term_labels(k)

        for term_name, term_lbl in terms.items():
            values[term_name] = self.topology_terms[term_name].values[term_lbl]

        return values

class forcefield_metadata:
    def __init__(self):
        self.authors: str = ""
        self.date: str = ""
        self.name: str = ""
        self.version: str = ""
        self.doi: str = ""
        self.description: str = ""
        self.coverage: str = ""
        self.training_molecules: str = ""
        self.training_methods: str = ""
        self.aux = {}

class forcefield:

    __slots__ = ("metadata", "models", "perception")

    def __init__(self, models: Dict[str,chemical_model], pcp_model):

        self.metadata: forcefield_metadata = forcefield_metadata()
        self.models: Dict[str, mm.chemical_model] = None
        self.perception: perception.perception_model = pcp_model

def chemical_system_to_physical_system(cs, pos: assignments.smiles_assignment) -> physical_model:
    ps = physical_system([])

    for cm in cs.models:
        pm = physical_model(pos, [], [])
        print(f"{datetime.datetime.now()} Processing", cm.name)
        for proc in cm.procedures:
            print(f"{datetime.datetime.now()}     Procedure", proc.name)
            procedure: chemical_model_procedure
            pm = proc.assign(cm, pm)
        ps.models.append(pm)
    return ps

def smiles_assignment_function(fn, sys_params, top_params, pos):
    result = {}
    for ic, x in pos.selections.items():
        ic_params = sys_params

        for t_params in top_params:
            p = t_params.get(ic, {})
            ic_params.update(p)

        if ic_params:
            try:
                result[ic] = fn(**ic_params, x=x)
            except TypeError as e:
                print("Partial parameterization: skipping. Error was:")
                print(e)

    return result
