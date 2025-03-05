"""
besmarts.mechanics.molecular_models
"""

import copy
from typing import Dict, List, Any

from besmarts.core import graphs
from besmarts.core import trees
from besmarts.core import tree_iterators
from besmarts.core import hierarchies
from besmarts.core import assignments
from besmarts.core import perception
from besmarts.core import topology


class topology_term:
    def __init__(
        self,
        symbol,
        name,
        unit,
        cast,
        values,
        comment,
        value_comments
    ):
        self.symbol: str = symbol
        self.name: str = name
        self.unit: str = unit
        self.cast: str = cast
        self.values: Dict[int, Any] = values
        self.comment: str = comment
        self.value_comments: Dict[int, str] = value_comments

    def copy(self):
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

    __slots__ = "labels", "values", "positions", "topology"

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
        self.name = name
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

        self.enable = 1


class physical_system:
    def __init__(self, models: List[physical_model]):
        self.models = models


class chemical_system:
    def __init__(
        self,
        pcp_model: perception.perception_model,
        models: List[chemical_model]
    ):
        self.perception = pcp_model
        self.models = models


class chemical_model_procedure_smarts_assignment(chemical_model_procedure):
    def __init__(self, pcp_model: perception.perception_model, topology_terms):
        self.name = ""
        self.perception = pcp_model

        self.topology_parameters: Dict[int, Dict[str, int | str]] = {}
        self.system_parameters: Dict[str, int] = {}

        self.topology_terms = topology_terms
        self.smarts_hierarchies: Dict[int, hierarchies.structure_hierarchy] = {}

        # if we don't find a match, use this instead
        self.default_parameter = None

    def assign(self, cm, pm: physical_model, overrides=None) -> physical_model:
        """
        this will return, for each selection, the reference
        """

        if overrides is None:
            overrides = {}
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
        for xi, x in enumerate(lbls.assignments):
            p = {}
            v = {}
            for ic, lbl in x.selections.items():
                if lbl is None:
                    if self.default_parameter:
                        lbl = self.default_parameter
                    else:
                        if self.smarts_hierarchies[unit_i].topology == topology.pair:
                            print(ic, lbl)
                        continue
                names = self.get_term_labels((unit_i, lbl))
                values = self.get_term_values((unit_i, lbl)).copy()
                for l, lv in overrides.items():
                    if l[1] == names.get(l[0]):
                        # print(f"Override: {l}: {values[l[0]][l[2]]} -> {lv}")
                        values[l[0]][l[2]] = lv
                if len(ic) == 1:
                    ic = xi, ic[0]
                else:
                    ic = tuple([(xi, ici) for ici in ic])
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
    m = int(S.model)
    p = int(S.type)
    u = int(S.category)
    return csys.models[m].procedures[p].smarts_hierarchies[u].index


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


def chemical_model_smarts_hierarchy_remove_node(
    cm: chemical_model,
    cid,
    pid,
    uid,
    node
):

    proc: chemical_model_procedure = cm.procedures[pid]

    h: hierarchies.structure_hierarchy = proc.smarts_hierarchies[uid]
    nodes = [x for x in h.index.nodes.values() if x.name == node.name]
    if len(nodes) > 1:
        print("Multiple nodes have the same name:")
        for n in nodes:
            print(n.index, n.name)
    assert len(nodes) == 1, "Multiple nodes have the same name"
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


def chemical_model_smarts_hierarchy_copy_node(
    cm: chemical_model,
    pid,
    uid,
    parent,
    name
):
    proc: chemical_model_procedure = cm.procedures[pid]

    h: hierarchies.structure_hierarchy = proc.smarts_hierarchies[uid]
    node = h.index.node_add_below(
        parent.index, index=0
    )
    assert (uid, name) not in proc.topology_parameters, f"{name} already exists"

    if name is None:
        i = max(h.index.nodes) + 1
        name = f"{cm.symbol}{i}"
        while (uid, name) in proc.topology_parameters:
            i += 1
            name = f"{cm.symbol}{i}"

    node.name = str(name)
    node.category = tuple(parent.category)
    node.type = str(parent.type)

    nodes = [x.name for x in h.index.nodes.values() if x.name == name]
    if len(nodes) > 1:
        print("Duplicate names:")
        for n in nodes:
            print(n.index, n.name)
    assert len(nodes) == 1, "Duplicate names"

    h.smarts[node.index] = str(h.smarts[parent.index])
    assert h.subgraphs[parent.index].select
    h.subgraphs[node.index] = graphs.subgraph_copy(h.subgraphs[parent.index])

    pkey = (uid, node.name)
    assert pkey not in proc.topology_parameters, f"{pkey} already present"

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


def chemical_model_smarts_hierarchy_add_node(
    cm,
    cid,
    pid,
    uid,
    parentid,
    node_ref,
    smarts,
    vals
):

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


def chemical_system_smarts_hierarchy_add_node(
    csys,
    cid,
    pid,
    uid,
    node_ref,
    smarts,
    vals: Dict[str, List]
):
    cm = csys.models[cid]
    return chemical_model_smarts_hierarchy_add_node(
        cm,
        cid,
        pid,
        uid,
        node_ref,
        smarts,
        vals
    )


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

    def __init__(self, models: Dict[str, chemical_model], pcp_model):

        self.metadata: forcefield_metadata = forcefield_metadata()
        self.models: Dict[str, chemical_model] = None
        self.perception: perception.perception_model = pcp_model


def chemical_system_iter_keys(csys):
    kv = {}
    for m, cm in enumerate(csys.models):
        for t in cm.system_terms:
            for i, v in enumerate(cm.system_terms[t].values):
                kv[(m, t, i)] = v
        for t in cm.topology_terms:
            for lbl, vl in cm.topology_terms[t].values.items():
                for i, v in enumerate(vl):
                    kv[(m, t, lbl, i)] = v
    return kv


def chemical_system_get_value_list(csys, key):
    if len(key) == 3:
        m, t, lbl = key
        return csys.models[m].topology_terms[t].values[lbl]
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
        m, t, lbl = key
        csys.models[m].topology_terms[t].values[lbl].clear()
        csys.models[m].topology_terms[t].values[lbl].extend(values)
    elif len(key) == 2:
        m, t = key
        csys.models[m].system_terms[t].values.clear()
        csys.models[m].system_terms[t].values.extend(values)


def chemical_system_set_value(csys, key, value):

    if len(key) == 4:
        m, t, l, i = key
        if t not in csys.models[m].topology_terms:
            raise IndexError(f"chemical_system_set_value: {t} not present")
        if l not in csys.models[m].topology_terms[t].values:
            raise IndexError(f"chemical_system_set_value: {l} not present")
        N = len(csys.models[m].topology_terms[t].values[l])
        if i < N:
            csys.models[m].topology_terms[t].values[l][i] = value
        else:
            print(f"Warning: chemical_system_set_value: index {i} not present. Adding")
            if csys.models[m].topology_terms[t].values[l] is None:
                csys.models[m].topology_terms[t].values[l] = []
            to_add = list([None]*(i-N)) + [value]
            csys.models[m].topology_terms[t].values[l].extend(to_add)

    elif len(key) == 3:
        m, t, i = key
        csys.models[m].system_terms[t].values[i] = value


def physical_system_set_value(psys: physical_system, key, value):

    if len(key) == 4:
        m, t, l, i = key
        for lbls, values in zip(psys.models[m].labels, psys.models[m].values):
            for ic, ic_lbls in lbls.items():
                if t not in ic_lbls:
                    continue
                term_lbl = ic_lbls[t]
                if term_lbl == l:
                    # print(f"PSYS SETTING {ic}:{t}:{i} = {value} from {values[ic][t][i]}")
                    values[ic][t][i] = value
    elif len(key) == 3:
        # TODO lol
        assert False


def chemical_system_get_value(csys, key, missing=None):
    if len(key) == 4:
        m, t, l, i = key
        try:
            v = csys.models[m].topology_terms[t].values[l][i]
            return v
        except IndexError:
            return missing
    elif len(key) == 3:
        m, t, i = key
        return csys.models[m].system_terms[t].values[i]


def physical_model_values_copy(pm):
    values = []
    for vals in pm.values:
        ic_vals = dict.fromkeys(vals)
        for ic, terms in vals.items():
            # ic_vals[ic] = {t: [*val_array] for t, val_array in terms.items()}
            ic_vals[ic] = {t: val_array.copy() for t, val_array in terms.items()}
        values.append(ic_vals)
    return values


warn_linear = True


def chemical_system_groupby_names(
    csys,
    m,
    psystems,
    selections,
    names=None
) -> dict:
    """
    from chemical_model m, group the assn by the labels in physical_model m
    """
    kv = {k[2]: [] for k in chemical_system_iter_keys(csys) if k[0] == m}
    global warn_linear
    warned = False

    for i, (psys, measure) in enumerate(zip(psystems, selections), 1):
        pm: physical_model = psys.models[m]
        pos = pm.positions[0]
        for ic, ic_terms in pm.labels[0].items():
            lbl = ic_terms['k']
            if names and lbl not in names:
                continue
            if ic not in measure:
                if warn_linear:
                    print(f"Warning, key {ic} did not have data (linear torsion?). Skipping.")
                    warned = True
            else:
                x = measure[ic][0]
                if lbl not in kv:
                    kv[lbl] = []
                kv[lbl].extend(x)
    if warned:
        warn_linear = False

    return kv


def chemical_system_get_ic_measure(csys, psystems, m, fn, names=None) -> dict:
    """
    """
    kv = {
        (k[0], 'l', k[2], None): []
        for k in chemical_system_iter_keys(csys)
        if k[0] == m and k[1] == 'l'
    }
    for psys in psystems:
        pm: physical_model = psys.models[m]
        measure = fn(pm.positions)
        for pi, pos in enumerate(pm.positions):
            for ic, ic_terms in pm.labels[pi].items():
                lbl = ic_terms['l']
                if names and lbl not in names:
                    continue
                x = measure.selections[ic][0]
                key = (m, 'l', lbl, None)
                if key not in kv:
                    kv[key] = []
                kv[key].extend(x)
    return kv


def chemical_system_get_bond_lengths(csys, psystems, names=None) -> dict:
    """
    """
    m = 0
    fn = assignments.graph_assignment_geometry_bond_matrix
    return chemical_system_get_ic_measure(csys, psystems, m, fn, names=names)


def chemical_system_get_angles(csys, psystems, names=None) -> dict:
    """
    """
    m = 1
    fn = assignments.graph_assignment_geometry_angle_matrix
    return chemical_system_get_ic_measure(csys, psystems, m, fn, names=names)


def chemical_system_get_ic_measure_means(csys, m, kv) -> dict:

    assert all((k[0] == m for k in kv))
    means = {}
    for k, v in chemical_system_iter_keys(csys).items():
        if k[0] == m and k[1] == 'l':
            r = kv[(m, 'l', k[2], None)]
            if r:
                means[k] = sum(r)/len(r)
    return means


def chemical_system_get_bond_length_means(csys, psystems, names=None) -> dict:
    kv = chemical_system_get_bond_lengths(csys, psystems, names=names)
    return chemical_system_get_ic_measure_means(csys, 0, kv)


def chemical_system_get_angle_means(csys, psystems, names=None) -> dict:
    kv = chemical_system_get_angles(csys, psystems, names=names)
    return chemical_system_get_ic_measure_means(csys, 1, kv)


def chemical_system_reset_angles(csys, psystems, names=None, skip=None) -> dict:
    kv = chemical_system_get_angle_means(csys, psystems, names=names)
    if skip is None:
        skip = []
    for k, v in kv.items():
        if k[2] not in skip:
            chemical_system_set_value(csys, k, v)
    return kv


def chemical_system_reset_bond_lengths(csys, psystems, names=None, skip=None) -> dict:
    kv = chemical_system_get_bond_length_means(csys, psystems, names=names)
    if skip is None:
        skip = []
    for k, v in kv.items():
        if k[2] not in skip:
            chemical_system_set_value(csys, k, v)
    return kv


def chemical_system_to_physical_system(
    cs,
    pos: List[assignments.graph_assignment],
    ref=None,
    reuse=None
) -> physical_model:
    ps = physical_system([])

    for ci, cm in enumerate(cs.models):
        if (ref is not None and reuse is not None) and ci in reuse:
            values = physical_model_values_copy(ref.models[ci])
            pm = physical_model(pos, ref.models[ci].labels.copy(), values)
            # ps.models.append(ref.models[ci])
        else:
            if cm.enable:
                pm = physical_model(pos, [], [])
                # print(f"{datetime.datetime.now()} Processing", cm.name)
                for proc in cm.procedures:
                    #print(f"{datetime.datetime.now()}     Procedure", proc.name)
                    procedure: chemical_model_procedure
                    pm = proc.assign(cm, pm)
            else:
                pm = physical_model(pos, [{}], [{}])
        ps.models.append(pm)
    return ps


def smiles_assignment_function(fn, sys_params, top_params, pos):
    result = {}
    out = []
    for ic, x in pos.selections.items():
        ic_params = dict(sys_params)

        for t_params in top_params:
            if len(ic[0]) == 3:
                p = dict(t_params.get(tuple((x[:2] for x in ic)), {}))
                p.update(dict(t_params.get(tuple((x[:3] for x in ic)), {})))
                if 's' not in p or len(set(x[2] for x in ic)) > 1:
                    p['s'] = [1.0]
                ic_params.update(p)
            elif len(ic[0]) == 2:
                p = t_params.get(ic, {})
            ic_params.update(p)

        # if 's' in ic_params:
        #     print(ic, ic_params)
        if ic_params:
            try:
                result[ic] = fn(**ic_params, x=x)
                out.append(f"{ic} {ic_params} {x} {result[ic]}")
                # print(ic, ic_params, x, result[ic])
            except TypeError as e:
                print("\n".join(out))
                breakpoint()
                print("Partial parameterization: skipping. Error was:")
                print(e)
                raise e

    # print("\n".join(out))
    return result

def smiles_assignment_matrix_function(fn, sys_params, top_params, posmat):
    result = {}

    for ic, x in posmat.selections.items():
        ic_params = dict(sys_params)

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


def chemical_system_smarts_complexity(csys: chemical_system, B=1.0, C=1.0):
    """
    get all smarts and calculate the smarts complexity
    """
    C0 = []
    atoms = 0
    parameters = 0
    terms = 0
    for ei, hidx in enumerate(
        chemical_system_iter_smarts_hierarchies(csys)
    ):
        for root in trees.tree_index_roots(hidx.index):
            M = len(hidx.topology.primary)
            for node in tree_iterators.tree_iter_dive(hidx.index, root):

                if node.type != 'parameter':
                    continue
                cm = chemical_system_get_node_model(csys, node)
                if cm.symbol in "IT":
                    t = sum([2 for x in cm.topology_terms['n'].values[node.name] if x < 4])
                    t += sum([x*x for x in cm.topology_terms['n'].values[node.name] if x > 3])
                elif cm.symbol in "ABN":
                    t = 2
                elif cm.symbol in "Q":
                    t = 1
                else:
                    t = 2

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
                        c = graphs.graph_complexity(g, scale=.01)
                        C0.append(c)
                        atoms += len(g.nodes)
                        parameters += 1
                        terms += t
                elif type(g) is not str:
                    c = graphs.graph_complexity(g, scale=.01)
                    C0.append(c)
                    # atoms += len(g.nodes) - M
                    atoms += len(g.nodes)
                    parameters += 1
                    terms += t

    # this is average (average bits per atom) of
    BX = sum(C0)/len(C0)
    # this is the scaled number of all atoms
    # CY = atoms*C

    c = terms * BX / 1000 * C
    # print(f"{c=} {terms=} {BX=} {parameters=}")
    return c


def chemical_system_print(csys, show_parameters=None):
    print("Model:")
    for ei, hidx in enumerate(
        chemical_system_iter_smarts_hierarchies(csys)
    ):
        print("Tree:")
        for root in trees.tree_index_roots(hidx.index):
            for e in tree_iterators.tree_iter_dive(hidx.index, root):
                s = trees.tree_index_node_depth(hidx.index, e)
                w = " "*s
                obj_repo = ""
                if e.type != 'parameter' or (show_parameters is None) or e.name in show_parameters:
                    sma = hidx.smarts.get(e.index, "")
                    if sma is None:
                        sma = ""

                    cm: chemical_model = chemical_system_get_node_model(csys, e)
                    params = []
                    for term_sym, term in cm.topology_terms.items():
                        param_vals = term.values.get(e.name)
                        if param_vals is not None:
                            params.append(f"{term_sym}: {param_vals}")
                    sma = hidx.smarts.get(e.index)
                    if sma is None:
                        sma = ""
                    print(
                        f"{s:2d} {int(e.category[0]):3d} {w}{e.name:4s}", sma, ' '.join(params)
                    )
