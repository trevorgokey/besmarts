"""
besmarts.mechanics.molecular_models
"""

from typing import Dict, List, Any
import datetime
import copy


from besmarts.core import graphs
from besmarts.core import trees
from besmarts.core import tree_iterators
from besmarts.core import hierarchies
from besmarts.core import assignments
from besmarts.core import perception
from besmarts.core import db


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
    """The functional form that can be evaluated as a function of positions.
    Also tracks the positions as a cache mechanism to avoid recomputing
    Sort of assumes energies/forces are a function of internal coordinates
    and masses and positions are cartesian"""

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
    returns a bunch tables

    """
    def __init__(self, name, topo):
        self.name = name # HESSIANS
        self.topology = topo
        self.procedure_parameters: Dict[str, int] = {}
    
    def assign(self, pm: physical_model):
        """
        this will return the compute function and config (i.e. the task)
        the task will return a assignments.graph_topology_db_table
        """
        return

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

        self.energy_function = None 
        self.force_function = None 
        self.internal_function = None
        self.derivative_function = None


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
        self.smarts_hierarchies: Dict[int, hierarchies.structure_hierarchy] = {}

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

def chemical_model_iter_smarts_hierarchies_nodes(cm):

    for proc in cm.procedures:
        if hasattr(proc, "smarts_hierarchies"):
            proc: chemical_model_procedure_smarts_assignment
            for hidx in proc.smarts_hierarchies.values():
                for root in trees.tree_index_roots(hidx.index):
                    yield from tree_iterators.tree_iter_dive(hidx.index, root)

def chemical_system_iter_smarts_hierarchies_nodes(csys):

    for cm in csys.models:
        yield from chemical_model_iter_smarts_hierarchies_nodes(cm)

def chemical_model_iter_smarts_hierarchies(cm: chemical_model):
    for proc in cm.procedures:
        if hasattr(proc, "smarts_hierarchies"):
            proc: chemical_model_procedure_smarts_assignment
            yield from proc.smarts_hierarchies.values()

def chemical_system_iter_smarts_hierarchies(csys):

    for cm in csys.models:
        yield from chemical_model_iter_smarts_hierarchies(cm)

def chemical_system_get_smarts_node(csys, S):
   return csys.models[int(S.model)].procedures[int(S.type)].smarts_hierarchies[int(S.category)].index


def chemical_system_get_node_hierarchy(csys, node):
    if node is None:
        return None
    m = int(node.category[0])
    cm = csys.models[m]
    for hidx in chemical_model_iter_smarts_hierarchies(cm):
        existing = hidx.index.nodes.get(node.index, None)
        if existing is None:
            continue
        if node.name == existing.name:
            return hidx

def chemical_system_get_node_model(csys, node):
    if node is None:
        return None
    m = int(node.category[0])
    cm = csys.models[m]
    return cm

def chemical_system_smarts_hierarchy_get_node_keys(cm, cid, pid, uid, node):
    kv = {}

    l = node.name
    for t, tv in cm.topology_terms.items():
        lval = tv.values.get(l)
        if lval is None:
            continue
        for i, v in enumerate(lval):
            kv[(cid, t, l, i)] = v

    return kv

def chemical_model_smarts_hierarchy_remove_node(cm: chemical_model, cid, pid, uid, node):

    proc: chemical_model_procedure = cm.procedures[pid]

    h: hierarchies.structure_hierarchy = proc.smarts_hierarchies[uid]
    h.index.node_remove(node.index)

    if node.index in h.smarts:
        h.smarts.pop(node.index)

    if node.index in h.subgraphs:
        h.subgraphs.pop(node.index)

    pkey = (uid, node.name)

    for tname in list(proc.topology_parameters[(uid, node.name)]):
        if node.name in cm.topology_terms[tname].values:
            cm.topology_terms[tname].values.pop(node.name)
    proc.topology_parameters.pop(pkey)
    return

def chemical_model_smarts_hierarchy_copy_node(cm: chemical_model, cid, pid, uid, parent, name):
    proc: chemical_model_procedure = cm.procedures[pid]

    h: hierarchies.structure_hierarchy = proc.smarts_hierarchies[uid]
    node = h.index.node_add_below(
        parent.index, index=0
    )
    if name is None:
        name = f"{cm.symbol}{node.index}"

    node.name = str(name)
    node.category = tuple(parent.category)
    node.type = str(parent.type)
    h.smarts[node.index] = str(h.smarts[parent.index])
    assert h.subgraphs[parent.index].select
    h.subgraphs[node.index] = graphs.subgraph_copy(h.subgraphs[parent.index])

    pkey = (uid, node.name)
    assert pkey not in proc.topology_parameters

    newparms = {}
    for k, v in proc.topology_parameters[(uid, parent.name)].items():
        if v == parent.name:
            newparms[k] = node.name
        else:
            newparms[k] = v
    proc.topology_parameters[pkey] = newparms
    
    
    for tname in list(proc.topology_parameters[(uid, parent.name)]):
        cm.topology_terms[tname].values[node.name] = copy.deepcopy(
            cm.topology_terms[tname].values[parent.name]
        )

    return node


def chemical_model_smarts_hierarchy_add_node(cm, cid, pid, uid, parentid, node_ref, smarts, vals):

    proc = cm.procedures[pid]

    h = proc.smarts_hierarchies[uid]
    node = h.index.node_add_below(
        parentid
    )
    node.name = str(node_ref.name)
    node.category = str(node_ref.category)
    node.type = str(cid)
    h.smarts[node.index] = smarts

    pkey = (uid, node.name)
    assert pkey not in proc.topology_parameters

    proc.topology_parameters[pkey] = {}
    for tname, tvals in vals.items():

        # Make sure that the term is recognized
        term = cm.topology_terms.get(tname)
        assert term
        # Store the values
        term.values[node.name] = tvals.copy()

        # Inform the cm/proc that this node can assign this term name
        proc.topology_parameters[pkey][tname] = node.name

    return node


def chemical_system_smarts_hierarchy_add_node(csys, cid, pid, uid, node_ref, smarts, vals: Dict[str, List]):
    cm = csys.models[cid]
    return chemical_model_smarts_hierarchy_add_node(cm, cid, pid, uid, node_ref, smarts, vals)
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

def chemical_system_iter_keys(csys):
    kv = {}
    for m, cm in enumerate(csys.models):
        for t in cm.system_terms:
            for i, v in enumerate(cm.system_terms[t].values):
                kv[(m, t, i)] = v
        for t in cm.topology_terms:
            for l, vl in cm.topology_terms[t].values.items():
                for i, v in enumerate(vl):
                    kv[(m, t, l, i)] = v
    return kv

def chemical_system_get_value_list(csys, key):
    if len(key) == 3:
        m, t, l = key
        return csys.models[m].topology_terms[t].values[l]
    elif len(key) == 2:
        m, t = key
        return csys.models[m].system_terms[t].values

def physical_system_iter_keys(psys_list: physical_system, csys: chemical_system):

    """
    Generate a flat mapping of keys and values of only the parameters that were
    applied to the physical systems
    """
    kv = {}

    for m, cm in enumerate(csys.models):
        for t in cm.system_terms:
            for i, v in enumerate(cm.system_terms[t].values):
                kv[(m, t, i)] = v

    for psys in psys_list:
        for m, pm in enumerate(psys.models):
            for proc in pm.labels:
                for p in proc.values():
                    for t, l in p.items():
                        for i, v in enumerate(csys.models[m].topology_terms[t].values[l]):
                            kv[(m, t, l, i)] = v
    return kv

def chemical_system_set_value_list(csys, key, values):

    if len(key) == 3:
        m, t, l = key
        csys.models[m].topology_terms[t].values[l] = values
    elif len(key) == 2:
        m, t = key
        csys.models[m].system_terms[t].values = values

def chemical_system_set_value(csys, key, value):

    if len(key) == 4:
        m, t, l, i = key
        csys.models[m].topology_terms[t].values[l][i] = value

    elif len(key) == 3:
        m, t, i = key
        csys.models[m].system_terms[t].values[i] = value

def physical_system_set_value(psys: physical_system, key, value):

    if len(key) == 4:
        m, t, l, i = key
        psys.models[m].topology_terms[t].values[l][i] = value

    elif len(key) == 3:
        m, t, i = key
        psys.models[m].system_terms[t].values[i] = value

    unit_i = 0
    lbl = key[2]
    lbls = pm.get_term_labels((unit_i, lbl))

    pm = psys.models[key[0]]

    if 0:
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

def chemical_system_get_value(csys, key):
    if len(key) == 4:
        m, t, l, i = key
        return csys.models[m].topology_terms[t].values[l][i]
    elif len(key) == 3:
        m, t, i = key
        return csys.models[m].system_terms[t].values[i]


def chemical_system_to_physical_system(cs, pos: assignments.graph_assignment, ref=None, reuse=None) -> physical_model:
    ps = physical_system([])

    for ci, cm in enumerate(cs.models):
        if (ref is not None and reuse is not None) and ci in reuse:
            ps.models.append(ref.models[ci])
        else:
            pm = physical_model(pos, [], [])
            # print(f"{datetime.datetime.now()} Processing", cm.name)
            for proc in cm.procedures:
                #print(f"{datetime.datetime.now()}     Procedure", proc.name)
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
                raise e

    return result

def chemical_system_smarts_complexity(csys: chemical_system):
    """
    get all smarts and calculate the smarts complexity
    """
    C0 = 0
    for ei, hidx in enumerate(
        chemical_system_iter_smarts_hierarchies(csys)
    ):
        for root in trees.tree_index_roots(hidx.index):
            for node in tree_iterators.tree_iter_dive(hidx.index, root):
                g = hidx.subgraphs.get(node.index)
                if g is None:
                    s = hidx.smarts.get(node.index)
                    if s is not None and s:
                        g = csys.perception.gcd.smarts_decode(s)
                        if type(g) is str:
                            # this means we could not parse the str,
                            # e.g. recursive smarts
                            continue
                        hidx.subgraphs[node.index] = g
                        C0 += graphs.graph_bits(g) / len(g.nodes) /1000 + len(g.nodes)
                elif type(g) is not str:
                    C0 += graphs.graph_bits(g) / len(g.nodes)/1000 + len(g.nodes)

                
    return 0.01*C0
