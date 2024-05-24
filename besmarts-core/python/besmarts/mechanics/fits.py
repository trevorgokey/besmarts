"""
besmarts.mechanics.fits
"""

import copy
import heapq
import datetime
import pickle
import multiprocessing
import collections
import sys
from typing import List, Dict, Tuple
from besmarts.core import configs
from besmarts.core import arrays
from besmarts.core import assignments
from besmarts.core import topology
from besmarts.core import graphs
from besmarts.core import codecs
from besmarts.core import compute
from besmarts.core import clusters
from besmarts.core import mapper
from besmarts.core import splits
from besmarts.core import trees
from besmarts.core import tree_iterators
from besmarts.core import optimization
from besmarts.mechanics import optimizers_scipy 
from besmarts.mechanics import objectives
from besmarts.mechanics import molecular_models as mm


class calculator_config_minimization:
    def __init__(self):
        self.minimize = True
        self.minimize_force_convergence = 0.10
        self.minimize_energy_convergence = 0.10
        self.minimize_position_convergence = 0.10
        self.step_limit = 100
        self.bias_force = {}
        self.constraints = {}


class graph_db_address:
    def __init__(self, *, eid=None, aid=None, gid=None, rid=None, sid=None, xid=None):
        self.eid: list = [] if not eid else eid
        self.aid: list = [] if not aid else aid
        self.gid: list = [] if not gid else gid
        self.rid: list = [] if not rid else rid
        self.sid: list = [] if not sid else sid
        self.xid: list = [] if not xid else xid

class compute_config:

    """
    Default is single point energy
    """

    def __init__(self, addr):
        self.addr = addr 
        self.keys = {}

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_table]]:
        # print("Starting calculation", self)
        csys = self.csys
        GDB = self.GDB
        tbl_idx = assignments.ENERGY
        tid = assignments.POSITIONS
        

        results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []

        for k, v in self.keys.items():
            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.6g} to {v:.6g} d={v-mm.chemical_system_get_value(csys, k):.6g}")
            mm.chemical_system_set_value(csys, k, v)
            # mm.physical_system_set_value(self.psys, k, v)

        for eid, gdb in GDB.entries.items():
            tbl = assignments.graph_db_table(topology.null)
            for gid in gdb.graphs:
                # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, tbl_idx, gid)
                pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid)

                # this means fold all graphs into a single graph
                # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, eid, tid)

                psys = mm.chemical_system_to_physical_system(csys, [pos0], ref=self.psys, reuse=self.reuse)
                ene = objectives.physical_system_energy(psys, csys)
                
                # ene_ga = assignments.graph_assignment("", {(0,): [[energy]]}, gdb.graphs[gid])
                # gdg = assignments.graph_assignment_to_graph_db_graph(ene_ga, topology.null)
                tbl[gid] = gdg
            r = {tbl_idx: tbl}
            # print("Done calculating", self)
            # print("Result is\n", r)
            results.append(r)
        return results


class compute_config_energy(compute_config):

    """
    Default is single point energy
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        # print("Starting calculation", self)
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.ENERGY
        tid = assignments.POSITIONS

        for k, v in self.keys.items():
            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.6g} to {v:.6g} d={v-mm.chemical_system_get_value(csys, k):.6g}")
            mm.chemical_system_set_value(csys, k, v)

        results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []
        for eid, gde in gdb.entries.items():
            tbl = assignments.graph_db_table(topology.null)
            rids = assignments.graph_db_table_get_row_ids(gde[tid])
            for rid in rids:
                system = []
                for gid in gdb.graphs:
                # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, assignments.POSITIONS, gid)

                    pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, rid, gid)
                    system.append(pos0)


                # create a ga but the indices are (gid, rid, sid)
                # just keep an index? 
                # system = graphs.graph_assignment_system(system)

                # psys now has positions in each model indexed by graph
                # this means I will need to make indices for graphs and confs
                # such as (gid, sid, rid)

                psys = mm.chemical_system_to_physical_system(csys, system, ref=self.psys, reuse=self.reuse)

                # this must put everything in the psys into a single system
                # this must mean that we flatten
                # pos = optimizers_scipy.optimize_positions_scipy(csys, psys)
                # psys = mm.chemical_system_to_physical_system(csys, [pos])
                ene = objectives.physical_system_energy(psys, csys)

                tbl.values.append(ene)
                # print("Calculated energy is", ene)
                # ene = objectives.physical_system_energy(psys, csys)
                # ene_ga = assignments.graph_assignment(
                #     "",
                #     {(0,): [[ene]]},
                #     graphs.subgraph_as_graph(gdb.graphs[gid])
                # )
                # # struct = assignments.graph_assignment_to_graph_db_struct(ene_ga, topology.null)
                # gdg = assignments.graph_assignment_to_graph_db_graph(ene_ga, topology.null)
                # tbl.graphs[gid] = gdg
            r = {tbl_idx: tbl}
            # print("Result is\n", r)
            results.append(r)
        # print("Done calculating", self)
        return results

class compute_config_gradient(compute_config):

    """
    Default is single point energy
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        # print("Starting calculation", self)
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.GRADIENTS
        tid = assignments.POSITIONS

        for k, v in self.keys.items():
            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.6g} to {v:.6g} d={v-mm.chemical_system_get_value(csys, k):.6g}")
            mm.chemical_system_set_value(csys, k, v)

        results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []
        for eid, gde in gdb.entries.items():
            tbl = assignments.graph_db_table(topology.atom)
            rids = assignments.graph_db_table_get_row_ids(gde[tid])
            for rid in rids:
                system = []
                for gid in gdb.graphs:
                # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, assignments.POSITIONS, gid)

                    pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, rid, gid)
                    system.append(pos0)


                # create a ga but the indices are (gid, rid, sid)
                # just keep an index? 
                # system = graphs.graph_assignment_system(system)

                # psys now has positions in each model indexed by graph
                # this means I will need to make indices for graphs and confs
                # such as (gid, sid, rid)

                psys = mm.chemical_system_to_physical_system(csys, system, ref=self.psys, reuse=self.reuse)

                # this must put everything in the psys into a single system
                # this must mean that we flatten
                # pos = optimizers_scipy.optimize_positions_scipy(csys, psys)
                # psys = mm.chemical_system_to_physical_system(csys, [pos])
                gx = objectives.physical_system_gradient(psys, csys)
                tbl.values.extend(gx)
                # print("Calculated energy is", ene)
                # ene = objectives.physical_system_energy(psys, csys)
                # ene_ga = assignments.graph_assignment(
                #     "",
                #     {(0,): [[ene]]},
                #     graphs.subgraph_as_graph(gdb.graphs[gid])
                # )
                # # struct = assignments.graph_assignment_to_graph_db_struct(ene_ga, topology.null)
                # gdg = assignments.graph_assignment_to_graph_db_graph(ene_ga, topology.null)
                # tbl.graphs[gid] = gdg
            r = {tbl_idx: tbl}
            # print("Result is\n", r)
            results.append(r)
        # print("Done calculating", self)
        return results

class compute_config_position(compute_config):

    """
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_table]]:
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.POSITIONS
        tid = assignments.POSITIONS

        for k, v in self.keys.items():
            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.6g} to {v:.6g} d={v-mm.chemical_system_get_value(csys, k):.6g}")
            mm.chemical_system_set_value(csys, k, v)
        # need to make this work for multi conformations
        results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []
        for eid, gde in gdb.entries.items():
            tbl = assignments.graph_db_table(topology.atom)
            for gid in gdb.graphs:
                for rid in gde[tid][gid].rows:
                    pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, rid, gid)
                    psys = mm.chemical_system_to_physical_system(csys, [pos0], ref=self.psys, reuse=self.reuse)
                    pos = optimizers_scipy.optimize_positions_scipy(csys, psys)
                    # pos = copy.deepcopy(pos0)
                    # struct = assignments.graph_assignment_to_graph_db_struct(pos, topology.atom)
                    gdg = assignments.graph_assignment_to_graph_db_graph(pos, topology.atom)
                    tbl.graphs[gid] = gdg
            r = {tbl_idx: tbl}
            results.append(r)
        # print("Done calculating", self)
        return results


class objective_config:

    def __init__(self, addr, scale=1):
        self.addr = addr
        self.scale = scale
        self.step_limit = None

    def get_task(self, GDB, csys, keys=None, psys=None, reuse=None) -> compute_config:
        cc = compute_config(self.addr)

        cc.GDB = {i: GDB[i] for i in self.addr.eid}
        cc.csys = csys
        if keys:
            cc.keys.update(keys)
        cc.psys = psys
        cc.reuse = reuse
        return cc

    def compute_gradient_2pt(self, X0, E, D: List[Dict[int, assignments.graph_db_table]], h):
        
        """
        D is the results for some objective idx
        so it is a list of tables (with their own tid)
                                   
        which now means that each is a table with structs which have indices
        """
        dx = []
        for etbls, dtbls  in zip(E, D):
            for tid, dtbl in dtbls.items():
                etbl = etbls[tid]
                for gid, dgdg in dtbl.graphs.items():
                    egdg = etbl.graphs[gid]
                    for rid, dgdr in dgdg.rows.items():
                        egdr = egdg.rows[rid]
                        for cid, dgdc in dgdr.columns.items():
                            egdc = egdr.columns[cid]
                            for v, v0 in zip(egdc.selections.values(), dgdc.selections.values()):
                                # vstr = "["+",".join([f"{vi:10.5f}" for vi in v]) + "]"
                                # v0str = "["+",".join([f"{vi:10.5f}" for vi in v0]) + "]"
                                x = arrays.array_difference(v0, v)
                                dx.extend(x)
                                # print(f"{rid:3d} Position SSE: {sse:10.5f} {vstr} {vstr}")
        # print(f"Total Position SSE: {obj:10.5f} A^2")
        dx = arrays.array_scale(dx, 1.0/(2*h))
        DX = arrays.array_inner_product(X0, dx)

        return DX

    def compute_diff(self, GDB: assignments.graph_db, D: List[Dict[int, assignments.graph_db_table]]):
        """
        D is the results for some objective idx
        so it is a list of tables (with their own tid)
                                   
        which now means that each is a table with structs which have indices
        """
        obj = []
        for gde, tbls  in zip(GDB.entries.values(), D):
            for tid, tbl in tbls.items():
                tbl0 = gde.tables[tid]
                for gid, gdg in tbl.graphs.items():
                    gdg0 = tbl0.graphs[gid]
                    for rid, gdr in gdg.rows.items():
                        gdr0 = gdg0.rows[rid]
                        for cid, gdc in gdr.columns.items():
                            gdc0 = gdr0.columns[cid]
                            for v, v0 in zip(gdc.selections.values(), gdc0.selections.values()):
                                vstr = "["+",".join([f"{vi:10.5f}" for vi in v]) + "]"
                                v0str = "["+",".join([f"{vi:10.5f}" for vi in v0]) + "]"
                                x = arrays.array_difference(v, v0)
                                obj.extend(x)
                                # print(f"{rid:3d} Position SSE: {sse:10.5f} {vstr} {vstr}")
        # print(f"Total Position SSE: {obj:10.5f} A^2")
        return obj

    def compute_total(self, GDB: assignments.graph_db, D: List[Dict[int, assignments.graph_db_table]]):
        """
        D is the results for some objective idx
        so it is a list of tables (with their own tid)
                                   
        which now means that each is a table with structs which have indices
        """
        obj = 0
        for gde, tbls  in zip(GDB.entries.values(), D):
            for tid, tbl in tbls.items():
                tbl0 = gde.tables[tid]
                for gid, gdg in tbl.graphs.items():
                    gdg0 = tbl0.graphs[gid]
                    for rid, gdr in gdg.rows.items():
                        gdr0 = gdg0.rows[rid]
                        for cid, gdc in gdr.columns.items():
                            gdc0 = gdr0.columns[cid]
                            for v, v0 in zip(gdc.selections.values(), gdc0.selections.values()):
                                vstr = "["+",".join([f"{vi:10.5f}" for vi in v]) + "]"
                                v0str = "["+",".join([f"{vi:10.5f}" for vi in v0]) + "]"
                                x = arrays.array_difference(v, v0)
                                sse = arrays.array_inner_product(x, x)
                                obj += sse
                                # print(f"{rid:3d} Position SSE: {sse:10.5f} {vstr} {vstr}")
        # print(f"Total Position SSE: {obj:10.5f} A^2")
        return obj

class objective_config_gradient(objective_config):

    def get_task(self, GDB: assignments.graph_db, csys, keys=None, psys=None, reuse=None) -> compute_config:
        cc = compute_config_gradient(self.addr)
        cc.GDB = assignments.graph_db_get_entries(GDB, self.addr.eid)
        cc.csys = csys
        cc.psys = psys
        if keys:
            cc.keys.update(keys)
        cc.reuse = reuse
        return cc

    def compute_gradient_2pt(self, X0, E0, E1, h):

        dx = []
        for etbls, dtbls  in zip(E0, E1):
            for tid, dtbl in dtbls.items():
                etbl = etbls[tid]
                x = arrays.array_difference(dtbl.values, etbl.values)
                dx.extend(x)
        dx = arrays.array_scale(dx, 1.0/(2*h))
        DX = arrays.array_inner_product(X0, dx)

        return DX

    def compute_diff(self, GDB: assignments.graph_db, D: List[Dict[int, assignments.graph_db_table]]):
        """
        """
        X = []
        for gde, tbls  in zip(GDB.entries.values(), D):
            for tid, tbl in tbls.items():
                tbl0 = gde.tables[tid]
                x = arrays.array_difference(tbl.values, tbl0.values)
                X.extend(x)

        return X

    def compute_total(self, GDB: assignments.graph_db, D: List[Dict[int, assignments.graph_db_table]]):
        """
        D is the results for some objective idx
        so it is a list of tables (with their own tid)
                                   
        which now means that each is a table with structs which have indices
        """
        obj = 0
        g = [0, 0, 0]
        g0 = [0, 0, 0]
        for gde, tbls  in zip(GDB.entries.values(), D):
            for tid, tbl in tbls.items():
                tbl0 = gde.tables[tid]
            
                for i in range(0,len(tbl.values),3):
                    v = tbl.values[i:i+3]
                    v0 = tbl0.values[i:i+3]
                    estr = "["+",".join([f"{vi:10.5f}" for vi in v]) + "]"
                    e0str = "["+",".join([f"{vi:10.5f}" for vi in v0]) + "]"
                    x = arrays.array_difference(v, v0)
                    g[0] += v[0]
                    g[1] += v[1]
                    g[2] += v[2]
                    g0[0] += v0[0]
                    g0[1] += v0[1]
                    g0[2] += v0[2]

                    sse = arrays.array_inner_product(x, x)
                    # print(f"Gradient SSE: {sse:10.5f} kJ/mol/A {estr} {e0str}")
                    obj += sse
        estr = "["+",".join([f"{vi:10.5f}" for vi in g]) + "]"
        e0str = "["+",".join([f"{vi:10.5f}" for vi in g0]) + "]"
        x = arrays.array_difference(g, g0)
        sse = arrays.array_inner_product(x, x)
        # print(f"G summed SSE: {sse:10.5f} kJ/mol/A {estr} {e0str}")
                
        # print(f"Total gradient SSE: {obj:10.5f} kJ/mol/A")
        return obj


class objective_config_energy(objective_config):

    def get_task(self, GDB: assignments.graph_db, csys, keys=None, psys=None, reuse=None) -> compute_config:
        cc = compute_config_energy(self.addr)
        cc.GDB = assignments.graph_db_get_entries(GDB, self.addr.eid)
        cc.csys = csys
        if keys:
            cc.keys.update(keys)
        cc.psys = psys
        cc.reuse = reuse
        return cc

    def compute_total(self, GDB: assignments.graph_db, D: List[Dict[int, assignments.graph_db_table]]):
        """
        D is the results for some objective idx
        so it is a list of tables (with their own tid)
                                   
        which now means that each is a table with structs which have indices
        """
        obj = 0
        for gde, tbls  in zip(GDB.entries.values(), D):
            for tid, tbl in tbls.items():
                tbl0 = gde.tables[tid]
                
                estr = "["+",".join([f"{vi:10.5f}" for vi in tbl.values]) + "]"
                e0str = "["+",".join([f"{vi:10.5f}" for vi in tbl0.values]) + "]"
                x = arrays.array_difference(tbl.values, tbl0.values)
                sse = arrays.array_inner_product(x, x)
                # print(f"Energy SSE: {sse:10.5f} kJ/mol {estr} {e0str}")
                obj += sse
                
        # print(f"Total Energy SSE: {obj:10.5f}")
        return obj

class objective_config_position(objective_config):

    def get_task(self, GDB, csys, keys=None, psys=None, reuse=None) -> compute_config:
        cc = compute_config_position(self.addr)
        cc.GDB = assignments.graph_db_get_entries(GDB, self.addr.eid)
        cc.csys = csys
        if keys:
            cc.keys.update(keys)
        cc.psys = psys
        cc.reuse = reuse
        return cc


class objective_tier:
    def __init__(self):
        self.objectives: Dict[int, objective_config] = {}
        self.bounds = {
            "s": (0,None), # scale
            "c": (0,None), # cutoff
            "r": (0,None), # sigma
            "e": (0,None), # epsilon
            "k": (0,None), # epsilon
            None: (None, None)
        }
        self.key_filter = lambda x: True if x[0] in [0,1,2,3] else False
        self.step_limit = None
        self.accept: int = 0 # 0 keeps all (essentially disables)

def objective_tier_get_keys(ot, csys):
    keys = *filter(ot.key_filter, mm.chemical_system_iter_keys(csys)),
    return keys
    
def gdb_to_physical_systems(gdb, csys):
    psysref = {}
    for eid, gde in gdb.entries.items():
        tid = assignments.POSITIONS
        gid = 0
        rid = 0
        pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
        psysref[eid] = mm.chemical_system_to_physical_system(csys, [pos])
    return psysref

def objective_tier_run(
    ot, gdb: assignments.graph_db, csys: mm.chemical_system, keys,
    oid=None, psysref=None, reuse=None, wq=None, verbose=False
):

    # csys = copy.deepcopy(csys)
    # build the dataset and input ff
    ws = None
    if wq:
        ws = compute.workqueue_new_workspace(wq, shm={})

    

    # build the initial psys (for charges)
    if verbose:
        print(datetime.datetime.now(), "Building physical systems")
    if psysref is None:
        psysref = gdb_to_physical_systems(gdb, csys)

    objectives = ot.objectives
    if oid:
        objectives = {i: objectives[i] for i in oid}
    
    # build the keys to fit
    # kv = mm.chemical_system_iter_keys(csys)
    # print("## All parameter terms")
    # for k,v in kv.items():
    #     print(k, v)
    # model_ids = list(range(len(csys.models)))
    # fit_ids = [0,1,2,3,4,5]
    # kv = {k:v for k, v in kv.items() if k[0] in fit_ids and k[1] in "s"}
    keys = [k for k in keys if ot.key_filter(k)]
    assigned = [(i, k, v) for psys in psysref.values() for i,m in enumerate(psys.models) for a in m.labels[0].values() for k, v in a.items()]
    keys = [k for k in keys if tuple(k[:3]) in assigned or k[1] in "s"]
    # reuse = []
    # # reuse = [i for i, m in enumerate(csys.models) if m.name in ["vdW", "Electrostatics"]]
    # # reuse = [i for i, m in enumerate(csys.models) if m.name in ["Electrostatics"]]
    # # reuse.extend([i for i in model_ids if i not in fit_ids])
    x0 = [mm.chemical_system_get_value(csys, k) for k in keys]

    # keys = list(kv.keys())
    # x0 = list(kv.values())
    # del kv
    # print("## Fitting parameter terms")
    # for k, v0 in zip(keys, x0):
    #     print(k, v0)
    
    args = (keys, csys, gdb, objectives, psysref, reuse, ws, verbose)
    # y0 = optimizers_scipy.fit_gdb(x0, *args) 

    bounds = []
    for k in keys:
        bounds.append(ot.bounds.get(k[1], (None, None)))

    result, y0, y1, gx = optimizers_scipy.optimize_forcefield_gdb_scipy(x0, args, bounds=bounds, step_limit=ot.step_limit)


    if verbose:
        print(f">>> Initial Objective {y0:10.5g}")
        print(f">>> Final Objective   {y1:10.5g}")
        print(f">>> Percent change    {(y1-y0)/y0*100:10.5g}%")

    if ws:
        ws.close()

    kv = {k: v for k,v in zip(keys, result)}

    return kv, y0, y1, gx

def objective_tier_take(
    ot: objective_tier,
    gdb,
    candidate_keys: List[Tuple],
    candidates: List[mm.chemical_system],
    strategy,
    psysref=None,
    wq=None
):
    """
    evaluate each candidate with the given objective_tier
    what should call this?
    something that just generated all candidates should call this b
    for the strat, we will continuously add params until the limit is hit, then
    return those candidates... probably just want to cap this to 1 per macro for now
    """

    results = []
    params = []

    if take_n is None:
        take_n = len(params)

    if wq is None:
        tasks = {i: csys for i, csys in enumerate(candidates)}
        for i, csys in tasks.items() :
            # need a method to get the keys from the csys
            keys = candidate_keys[i]
            p, x = objective_tier_run(ot, gdb, csys, keys, psysref, wq=wq)
            
            heapq.heappush(results, (x, i))
            params.append({k: v for k, v in zip(keys, p)})
        results = heapq.nsmallest(take_n, results)
        results = {i: (x, params[i]) for x, i in results}
    else:
        ws = compute.workqueue_new_workspace(wq)
        
        tasks = {}
        for i, cys in enumerate(candidates):
            keys = candidate_keys[i]
            tasks[i] = ((ot, gdb, csys, keys), {"psysref": psysref, "wq": wq})
        ret = compute.workspace_submit_and_flush(ws, objective_tier_run, tasks)
        results = ((x, i) for i, (_, x) in ret.items())
        results = heapq.heapify(results)
        results = heapq.nsmallest(take_n, results)
        params = zip(tasks[i][0][3], ret[i][0])
        results = {i: (x, {k: v for k, v in params}) for x, i in results}
        ws.close()

    return results

def objective_run_distributed(obj, shm=None):
    return obj.run()


# we pass the dataset because we need to grab the compute settings
# a deduplicator should combine 0,2 and 1,3 since it has the same everything
# except aid
#objectives = {
#    0 :objective_config(
#        graph_db_address(
#            xid=[0],
#            aid=[assignments.POSITIONS]
#        ),
#        min=False
#    ),
#    1: objective_config(
#        graph_db_address(
#            xid=[0],
#            aid=[assignments.POSITIONS]
#        ),
#        min=True
#    ),
#    2: objective_config(
#        graph_db_address(
#            xid=[0],
#            aid=[assignments.GRADIENTS]
#        ),
#        min=False
#    ),
#    3: objective_config(
#        graph_db_address(
#            xid=[0],
#            aid=[assignments.GRADIENTS]
#        ),
#        min=True
#    )
#}
#
#D: list[graph_topology_db] = [graph_topology_db()]

def forcefield_optimization_strategy_default(csys, models=None):

    """
    determines how to step the optimization forward, choosing which hyperparameter
    to try next
    """

    bounds = {}
    for m, cm in enumerate(csys.models):
        if models is not None and m not in models:
            continue
        for h, hidx in enumerate(mm.chemical_model_iter_smarts_hierarchies(cm)):
            splitter = configs.smarts_splitter_config(
                1, 1, 0, 0, 0, 0, True, True, 0, True, True, True, True
            )
            extender = configs.smarts_extender_config(
                0, 0, True
            )
            bounds[m] = configs.smarts_perception_config(splitter, extender)
            break
    
    return forcefield_optimization_strategy(bounds)


class forcefield_optimization_strategy(optimization.optimization_strategy):
    def __init__(self, bounds):
        # build the 
        self.bounds = bounds

        self.objective_accept_total = [0]

        # Only consider the top N clusters of each objective state
        self.objective_accept_clusters = [0]

        # Update objective on each evaluation. Some objectives change if new
        # clusters are added. This option determines whether accepting causes a refresh
        self.objective_update_on_each_accept = True

        self.cursor = -1
        self.maxedits_limit = 0
        self.repeat = False

        self.direct_enable = False
        self.direct_limit = 10

        self.iterative_enable = True

        self.enable_merge = True
        self.enable_split = True
        self.enable_modify = False

        self.steps: List[optimization.optimization_iteration] = []
        self.tree_iterator: Callable = mm.chemical_system_iter_smarts_hierarchies_nodes

        self.step_tracker = {}

        # Number of operations to accept per macro step
        # Relabeling is done here
        # self.accept_max = 1 will give best performance, but will cost the most
        # self.accept_max = 0 is no max
        self.macro_accept_max_total: int = 1

        # Number of operations to accept per micro step
        # We do not relabel here but instead just keep this many
        # self.accept_max = 1 will give best performance, but will cost the most
        # self.accept_max = 0 is no max
        self.micro_accept_max_total: int = 1

        # Number of operations to accept per step per cluster
        # self.accept_max = 1 will give best performance, but will cost the most
        # self.accept_max = 0 is no max
        self.macro_accept_max_per_cluster: int = 1

        # Number of operations to accept per step per cluster
        # self.accept_max = 1 will give best performance, but will cost the most
        # self.accept_max = 0 is no max
        self.micro_accept_max_per_cluster: int = 1

        # If we accept too many operations, some can match nothing due to
        # unexpected occlusion. With this enabled, we shortcut merging and 
        # prevent zero-matching SMARTS from being added.
        self.prune_empty = True


        # Do not merge these
        self.merge_protect_list = []

        # Only operate on these
        self.target_list = []

        # This removes candidates which have an estimated objective diff above
        # this value
        # None disables
        # 0.0 will prune anything that is deemed useless
        self.filter_above: float = 0.0

        self.keep_below: float = 0.0
    def macro_iteration(
        self, clusters: List[trees.tree_node]
    ) -> optimization.optimization_iteration:
        """
        Return a list of iterations that form a macro iteration, where
        we may want to analyze a group of candidates before proceeding
        to the next level of searching

        Parameters
        ----------
        clusters: List[trees.tree_node]
            The nodes of a trees.tree_index to consider in the step

        Returns
        -------
        optimization_step
        """
        if not self.steps:
            self.build_steps()
        return optimization.optimization_strategy_iteration_next(self, clusters)

    def build_steps(self):
        self.steps.extend(
            forcefield_optimization_strategy_build_macro_iterations(self)
        )

def forcefield_optimization_strategy_build_macro_iterations(strat: forcefield_optimization_strategy):

    macro_iters = []
    boundlist = [x.splitter for x in strat.bounds.values()]
    bounds = boundlist[0]
    bounds.branch_depth_min = min((x.branch_depth_min for x in boundlist))
    bounds.branch_depth_limit = max((x.branch_depth_limit for x in boundlist))
    bounds.branch_limit = max((x.branch_limit for x in boundlist))
    bounds.branch_min = min((x.branch_min for x in boundlist))
    bounds.bit_search_min = min((x.bit_search_min for x in boundlist))
    bounds.bit_search_limit = min((x.bit_search_limit for x in boundlist))

    search_cursor = -1

    for branch_d in range(
        bounds.branch_depth_min, bounds.branch_depth_limit + 1
    ):
        for branches in range(bounds.branch_limit, bounds.branch_limit + 1):
            bits = bounds.bit_search_min - 1
            while bits < bounds.bit_search_limit:
                bits += 1

                search_cursor += 1
                if search_cursor < strat.cursor:
                    continue

                if strat.enable_split:
                    steps = []
                    for m, mbounds in strat.bounds.items():

                        mbounds: configs.smarts_perception_config

                        # go through the model bounds and add if current settings
                        # are a subset of bounds
                        if mbounds.splitter.bit_search_limit < bits:
                            continue
                        if mbounds.splitter.bit_search_min > bits:
                            continue
                        if mbounds.splitter.branch_limit < branches:
                            continue
                        if mbounds.splitter.branch_min > branches:
                            continue
                        if mbounds.splitter.branch_depth_limit < branch_d:
                            continue
                        if mbounds.splitter.branch_depth_min > branch_d:
                            continue

                        s = optimization.optimization_step()
                        s.index = len(steps)
                        s.cluster = None
                        s.overlap = [0]
                        s.models.append(m)

                        s.direct_enable = strat.direct_enable
                        s.direct_limit = strat.direct_limit
                        s.iterative_enable = strat.iterative_enable

                        s.operation = strat.SPLIT

                        splitter = configs.smarts_splitter_config(
                            bits,
                            bits,
                            0,
                            branches,
                            branch_d,
                            branch_d,
                            mbounds.splitter.unique,
                            False,
                            0,
                            mbounds.splitter.split_general,
                            mbounds.splitter.split_specific,
                            mbounds.splitter.unique_complements,
                            mbounds.splitter.unique_complements_prefer_min,
                        )
                        extender = configs.smarts_extender_config(
                            branches, branches, mbounds.extender.include_hydrogen
                        )
                        config = configs.smarts_perception_config(
                            splitter, extender
                        )
                        s.pcp = config

                        steps.append(s)
                    macro_iters.append(optimization.optimization_iteration(steps))
                if strat.enable_merge:
                    steps = []
                    for m, mbounds in strat.bounds.items():

                        mbounds: configs.smarts_perception_config

                        # go through the model bounds and add if current settings
                        # are a subset of bounds
                        if mbounds.splitter.bit_search_limit > bits:
                            continue
                        if mbounds.splitter.bit_search_min < bits:
                            continue
                        if mbounds.splitter.branch_limit > branches:
                            continue
                        if mbounds.splitter.branch_min < branches:
                            continue
                        if mbounds.splitter.branch_depth_limit > branch_d:
                            continue
                        if mbounds.splitter.branch_depth_min < branch_d:
                            continue
                        s = optimization.optimization_step()
                        s.index = len(steps)
                        s.models.append(m)
                        s.cluster = None
                        s.overlap = [0]
                        s.direct_enable = strat.direct_enable
                        s.direct_limit = strat.direct_limit
                        s.operation = strat.MERGE
                        s.iterative_enable = strat.iterative_enable

                        splitter = configs.smarts_splitter_config(
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            True,
                            True,
                            0,
                            mbounds.splitter.split_general,
                            mbounds.splitter.split_specific,
                            mbounds.splitter.unique_complements,
                            mbounds.splitter.unique_complements_prefer_min,
                        )
                        extender = configs.smarts_extender_config(0, 0, True)
                        config = configs.smarts_perception_config(
                            splitter, extender
                        )
                        s.pcp = config
                        steps.append(s)

                    macro_iters.append(optimization.optimization_iteration(steps))

                if strat.enable_modify:
                    steps = []
                    for m, mbounds in strat.bounds.items():

                        mbounds: configs.smarts_perception_config

                        # go through the model bounds and add if current settings
                        # are a subset of bounds
                        if mbounds.splitter.bit_search_limit > bits:
                            continue
                        if mbounds.splitter.bit_search_min < bits:
                            continue
                        if mbounds.splitter.branch_limit > branches:
                            continue
                        if mbounds.splitter.branch_min < branches:
                            continue
                        if mbounds.splitter.branch_depth_limit > branch_d:
                            continue
                        if mbounds.splitter.branch_depth_min < branch_d:
                            continue
                        s = optimization.optimization_step()
                        s.index = len(steps)
                        s.models.append(m)
                        s.cluster = None
                        s.overlap = [0]
                        s.direct_enable = False
                        s.direct_limit = False
                        s.operation = strat.MODIFY
                        s.iterative_enable = False

                        splitter = configs.smarts_splitter_config(
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            True,
                            True,
                            0,
                            mbounds.splitter.split_general,
                            mbounds.splitter.split_specific,
                            mbounds.splitter.unique_complements,
                            mbounds.splitter.unique_complements_prefer_min,
                        )
                        extender = configs.smarts_extender_config(0, 0, True)
                        config = configs.smarts_perception_config(
                            splitter, extender
                        )
                        s.pcp = config
                        steps.append(s)

                    macro_iters.append(optimization.optimization_iteration(steps))

                # print("MACRO MERGE")

    strat.cursor = 0

    return macro_iters

def chemical_system_to_optimization_strategy(csys) -> forcefield_optimization_strategy:
    """
    build a strat for each, which is setting the bounds
    just have a list of bounds for each?

    """

# create the objective
# the objectives will define the required task
# and in the future we can condense tasks
# {oid: o}
# obj = {i: objective_config(addr, min=False) for i, addr in enumerate(addrs)} 

# these are the tasks that need to be run
# tasks = {i: sp.get_task(D) for i, sp in objectives.items()}

# this will return a gdb for each
# {oid: {aid: table}}
# results = {i: task.run() for i, task in tasks.items()}

# now we give the result to the objective and it computes the total
# X = sum((x.compute_total([D[xi] for xi in x.addr.xid], results[i]) for i, x in obj.items()))

eid_t = int
def ff_optimize(
    csys0: mm.chemical_system,
    gdb: assignments.graph_db,
    psystems: Dict[eid_t, mm.physical_system],
    strategy: forcefield_optimization_strategy,
    chemical_objective: clusters.clustering_objective,
    initial_objective: objective_tier,
    tiers: List[objective_tier],
    final_objective: objective_tier,
) -> mm.chemical_system:

    """
    initial and final is for the "gold standard" fit
    Dive into a strategy and it will... return a list of force field candidates
    right now it only we would only return 1...
    we ideally want a system where we consider ffs with multiple modifications

    So..:
        1. generate all splits for each model
        2. for tier in tiers: score and filter
        3. accept all remaining
        4. repeat until convergence

    Previously, we would calculate the objective based on assn and the 
    objective split function. This should now just be based on the objective function.
    Split will obviously mean X1 - X0. In this case I will want:
        1. total single+grad
    This means that I will just need a single point. If I do a split, then it will
    be for chemical objective. Then it means I should pack the data with the chemical
    objective. In this case I can give it the two parameters and it can compute the total objective
    """

    started = datetime.datetime.now()
    max_line = 0

    # smiles = [a.smiles for a in sag.assignments]

    # topo = sag.topology

    # group_prefix_str = initial_conditions.group_prefix_str

    # would need to get all clusters... then get the macro for it
    # each node will reference the model

    # hidx = initial_conditions.hierarchy.copy()

    # this now refers to the chemical objective
    # will need to refactor this
    # groups = clustering_build_ordinal_mappings(initial_conditions, sag)
    # print(groups)

    # match = clustering_build_assignment_mappings(initial_conditions, initial_conditions.group)
    # match = initial_conditions.mappings
    # assn = get_assns(sag.assignments, topo)

    gcd = csys0.perception.gcd
    labeler = csys0.perception.labeler
    icd = codecs.intvec_codec(gcd.primitive_codecs, gcd.atom_primitives, gcd.bond_primitives)

    csys = copy.deepcopy(csys0)

    # check_lbls_data_selections_equal(initial_conditions.group, sag)


    """
    nomenclature

    match is node_name:data_idx (mol, ic)
    assn is data_idx:data
    group is node_name:data
    mapping node_name:data_idx
    """

    # group_number = max(
    #     [
    #         int(x.name[1:])
    #         for x in hidx.index.nodes.values()
    #         if x.name[0] == group_prefix_str
    #     ]
    # )

    # group_number += 1

    # gc.collect()
    # lets assume we already have graphs built
    # G0 = gdb.graphs
    G0 = {i: icd.graph_encode(g) for i, g in gdb.graphs.items()}

    # G0 = {i: icd.graph_encode(gcd.smiles_decode(a.smiles)) for i, a in enumerate(sag.assignments)}
    # n_ics = sum((len(s.selections) for s in sag.assignments))


    N = 1
    n_ics = 1
    # N = len(hidx.index.nodes)
    # try:
    #     N = len(set(assn.values()))
    # except Exception:
    #     pass

    repeat = set()
    visited = set()
    iteration = 1
    N_str = "{:" + str(len(str(n_ics))) + "d}"
    success = False

    # roots = trees.tree_index_roots(hidx.index)

    # print(f"{datetime.datetime.now()} Labeling subgraphs")
    # assignments = labeler.assign(hidx, gcd, smiles, topo)
    # csys = csys0
    # cst = smarts_clustering(
    #     hidx,
    #     assignments,
    #     clustering_build_assignment_mappings(hidx, assignments),
    # )
    # print(f"{datetime.datetime.now()} Checking consistency...")
    # check_lbls_data_selections_equal(cst.group, sag)
    wq = compute.workqueue_local('0.0.0.0', configs.workqueue_port)

    # this is now chemical objective
    print(f"{datetime.datetime.now()} Computing chemical objective")
    # TODO: Compute chemical and physical objectives
    # and the CO would 

    C0 = chemical_objective(csys)
    C = C0

    # _, X0 = get_objective(
    #     cyss, assn, objective.split, strategy.overlaps[0], splitting=True
    # )

    print(f"{datetime.datetime.now()} Computing physical objective")
    # need a physical objective



    fitkeys = objective_tier_get_keys(initial_objective, csys)
    fitkeys = [k for k in fitkeys if k[1] in "skeler"]
    fitting_models = set((k[0] for k in fitkeys))
    reuse0=[k for k,_ in enumerate(csys0.models) if k not in fitting_models]
    kv, P00, P0, gp0 = objective_tier_run(
        initial_objective,
        gdb,
        csys,
        fitkeys,
        psysref=psystems,
        reuse=reuse0,
        wq=wq,
        verbose=True
    )
    X0 = P0 + C0
    P = P0
    print(datetime.datetime.now(), f"Initial objective: {P0:13.6g} C={C0:13.6g}")
    for k, v in kv.items():
        v0 = mm.chemical_system_get_value(csys0, k)
        mm.chemical_system_set_value(csys, k, v)
        print(f"{str(k):20s} | New: {v:12.6g} Ref {v0:12.6g} Diff {v-v0:12.6g}")
    # for i, tier in enumerate(tiers):
    #     _, P00, P0, gp0 = objective_tier_run(
    #         tier,
    #         gdb,
    #         csys0,
    #         keys,
    #         psysref=psystems,
    #         reuse=[k for k,_ in enumerate(csys0.models) if k not in fitting_models],
    #         wq=wq
    #     )
    #     print(datetime.datetime.now(), f"Tier {i} objective: P={P0:13.6g} C={C0:13.6g}")
    # here is where I would go through each strategy and return the candidates

    # print(f"Total objective C={C0} P={P0}")
    # print(f"Parameters:")
    # for k,v in kv.items():
    #     print(k, "Val:", v, "Diff:", v - mm.chemical_system_get_value(csys0, k))
    # return C0 + P0

    if not strategy.steps:
        print("Optimization strategy is building steps...")
        strategy.build_steps()
    print(f"{datetime.datetime.now()} The optimization strategy has the following iterations:")
    for ma_i, macro in enumerate(strategy.steps, 1):
        cur = "  "
        if ma_i == strategy.cursor + 1:
            cur = "->"
        for mi_i, micro in enumerate(macro.steps):
            s = micro.pcp.splitter
            a = micro.overlap
            m = micro.models
            b0 = s.bit_search_min
            b1 = s.bit_search_limit
            d0 = s.branch_depth_min
            d1 = s.branch_depth_limit
            n0 = s.branch_min
            n1 = s.branch_limit
            # probably print the models
            print(
                f"{cur} {ma_i:3d}:{micro.index:02d}. op={micro.operation:2d} m={m} a={a} b={b0}->{b1} d={d0}->{d1} n={n0}->{n1}"
            )

    step_tracker = strategy.step_tracker

    while True:
        if success:
            print("Restarting optimization search")
            strategy = optimization.optimization_strategy_restart(strategy)
            success = False

        elif optimization.optimization_strategy_is_done(strategy):
            print("Nothing found. Done.")
            break

        # groups = clustering_build_ordinal_mappings(cst, sag)

        # roots = trees.tree_index_roots(cst.hierarchy.index)
        nodes = [
            x
            for x in strategy.tree_iterator(csys)
            if x.type == "parameter"
            # if strategy.cursor == -1
            # or strategy.cursor >= step_tracker.get((x.category, x.name), 0)
            # and x.type == "parameter"
        ]
        # remove any that are not in the models
        
        print(f"Targets for this macro step {strategy.cursor+1}:")
        for nidx, n in enumerate(nodes, 1):
            print(nidx, n.category, n.name)
        print(f"N Targets: {len(nodes)}")

        print(f"Step tracker for current macro step {strategy.cursor+1}")
        for n, v in step_tracker.items():
            print(n, v + 1)

        # has all nodes
        macro: optimization.optimization_iteration = strategy.macro_iteration(
            nodes
        )

        candidates = {}
        pq = []
        n_added = 0
        n_macro = len(strategy.steps)

        t = datetime.datetime.now()
        for micro in macro.steps:
            config = micro.pcp

            spg = "Y" if config.splitter.split_general else "N"
            sps = "Y" if config.splitter.split_specific else "N"
            print(
                f"\n*******************\n {t}"
                f" iteration={iteration:4d}"
                f" macro={strategy.cursor:3d}/{n_macro}"
                f" X={X0:9.5g}"
                # f" params=({len(cst.mappings)}|{N})"
                f" G={spg}"
                f" S={sps}"
                f" bits={config.splitter.bit_search_min}->{config.splitter.bit_search_limit}"
                f" depth={config.splitter.branch_depth_min}->{config.splitter.branch_depth_limit}"
                f" branch={config.splitter.branch_min}->{config.splitter.branch_limit}"
            )
            print(f"*******************")
            print()
        print("Tree:")
        for ei, e in enumerate(
            mm.chemical_system_iter_smarts_hierarchies_nodes(csys)
        ):
            # s = trees.tree_index_node_depth(cst.hierarchy.index, e)
            s = 0
            # obj_repo = ""
            # if groups[e.name]:
            # obj_repo = objective.report(groups[e.name])
            obj_repo = ""
            print(
                f"** {s:2d} {ei:3d} {e.category} {e.name:4s}",
                # obj_repo,
                # cst.hierarchy.smarts.get(e.index),
            )
        print("=====\n")

        # compute the current psystems

        reuse = reuse0
        psys = {}
        for eid, gde in gdb.entries.items():
            tid = assignments.POSITIONS
            gid = 0
            rid = 0
            pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
            psys[eid] = mm.chemical_system_to_physical_system(csys, [pos], ref=psystems[eid], reuse=reuse)
        psystems = psys


        print(f"{datetime.datetime.now()} Saving checkpoint to chk.cst.p")
        pickle.dump([gdb, csys, strategy, psystems], open("chk.cst.p", "wb"))

        
        step = None

        while not optimization.optimization_iteration_is_done(macro):
            # t = datetime.datetime.now()
            # print(f"{t} Initializing new loop on macro {strategy.cursor}")

            step: optimization.optimization_step = (
                optimization.optimization_iteration_next(macro)
            )

            n_micro = len(macro.steps)
            config: configs.smarts_perception_config = step.pcp
            S = step.cluster
            tkey = (S.category, S.name)
            hidx = mm.chemical_system_get_node_hierarchy(csys, S)
            topo = hidx.topology


            if S.type == 'parameter' and S.index not in hidx.subgraphs:
                hidx.subgraphs[S.index] = gcd.smarts_decode(hidx.smarts[S.index])
            if type(hidx.subgraphs[S.index]) is str:
                # could not parse smarts (e.g. recursive) so we skip
                step_tracker[tkey] = strategy.cursor
                continue
            S0 = graphs.subgraph_to_structure(hidx.subgraphs[S.index], topo)
            # S0 = graphs.subgraph_to_structure(
            #     cst.hierarchy.subgraphs[S.index], topo
            # )

            cfg = config.extender.copy()

            S0_depth = graphs.structure_max_depth(S0)
            d = max(S0_depth, config.splitter.branch_depth_limit)
            # print("Depth is", d)
            cfg.depth_max = d
            cfg.depth_min = S0_depth

            t = datetime.datetime.now()
            # need to get the graphs that match S0
            print(
                f"{t} Collecting SMARTS for {S.name} and setting to depth={S0_depth}"
            )

            # find only unique graphs!

            # now I need to reparameterize and find all matching subgraphs/indices
            selected_ics = []
            selected_graphs = set()
            for eid, ps in psystems.items():
                for gid, lbls in enumerate(ps.models[int(S.category[0])].labels):
                    for ic, term_lbls in lbls.items():
                        if S.name in term_lbls.values():
                            selected_ics.append((eid, ic))
                            selected_graphs.add(eid)
            # aa = cst.mappings[S.name]
            # selected_graphs = set((x[0] for x in aa))
            aa = selected_ics
            G = {k:v for k,v in G0.items() if k in selected_graphs}
            del selected_graphs


            assn_s = aa
            # assn_s = {i: assn[i] for i in cst.mappings[S.name]}

            iteration += 1

            t = datetime.datetime.now()
            print(
                f" =="
                f" iteration={iteration:4d}"
                f" macro={strategy.cursor:3d}/{n_macro}"
                f" micro={macro.cursor:3d}/{n_micro}"
                # f" overlap={strategy.overlaps:3d}"
                f" operation={step.operation}"
                # f" params=({len(cst.mappings)}|{N})"
                f" cluster={S.name:4s}"
                f" N= " + N_str.format(len(aa)) + ""
                f" overlap={step.overlap}"
                f" bits={config.splitter.bit_search_min}->{config.splitter.bit_search_limit}"
                f" depth={config.splitter.branch_depth_min}->{config.splitter.branch_depth_limit}"
                f" branch={config.splitter.branch_min}->{config.splitter.branch_limit}"
            )
            print()

            if step.operation == strategy.SPLIT:
                new_pq = []
                new_candidates = {}
                new_candidates_direct = {}
                direct_success = False

                print(f"Attempting to split {S.name}:")
                print("S0:", gcd.smarts_encode(S0))



                if not aa:
                    print("No matches.")
                    step_tracker[tkey] = strategy.cursor
                    continue

                print(f"Matched N={len(aa)}")
                seen = set()
                # depth = graphs.structure_max_depth(S0)
                extend_config = config.extender.copy()
                extend_config.depth_max = config.splitter.branch_depth_limit
                extend_config.depth_min = config.splitter.branch_depth_min

                # For each node, I present just.. the chemical objective
                # until I can make a case for IC objective accounting

                for seen_i, i in enumerate(aa, 1):
                    g = graphs.graph_to_structure(icd.graph_decode(G[i[0]]), i[1], topo) 
                    graphs.structure_extend(extend_config, [g])
                    g = graphs.structure_remove_unselected(g)
                    # if g not in seen:
                    # seen_g.add(g)
                    if seen_i < 100:
                        print(
                            f"{seen_i:06d} {str(i):24s}",
                            # objective.report([x]),
                            gcd.smarts_encode(g),
                        )
                        seen.add(g)
                print()
                if len(seen) < 2 and len(assn_s) < 100:
                    print(f"Skipping {S.name} since all graphs are the same")
                    step_tracker[tkey] = strategy.cursor
                    continue

                # if len(set(map(tuple, groups[S.name]))) == 1:
                #     print(f"Skipping {S.name} since all data are the same")
                #     step_tracker[tkey.name] = strategy.cursor
                #     continue


                if len(seen) < 2 and len(assn_s) < 100:
                    print(f"Skipping {S.name} since all graphs are the same")
                    step_tracker[tkey] = strategy.cursor
                    continue

                # if objective.single(assn_s.values()) == 0.0:
                #     print(f"Skipping {S.name} due to no objective")
                #     step_tracker[tkey] = strategy.cursor
                #     continue

                if graphs.structure_max_depth(S0) > config.splitter.branch_depth_min:
                    print("This parameter exceeds current depth. Skipping")
                    step_tracker[tkey] = strategy.cursor
                    continue


                # this is where I generate all candidates
                if step.direct_enable and config.splitter.bit_search_limit < len(assn_s):
                    assn_i = []
                    if len(assn_s) < step.direct_limit:
                        # or form matches based on unique smarts
                        a = []
                        extend_config = config.extender.copy()
                        extend_config.depth_max = config.splitter.branch_depth_limit
                        extend_config.depth_min = config.splitter.branch_depth_min
                        for seen_i, (i, x) in enumerate(assn_s.items(), 1):
                            g = graphs.graph_to_structure(icd.graph_decode(G[i[0]]), i[1], topo) 
                            graphs.structure_extend(extend_config, [g])
                            g = graphs.structure_remove_unselected(g)
                            a.append(g)
                        # if len(a) < step.direct_limit:
                        #     lbls = dict()
                            # for (i, idx), sg in a.items():
                            #     lbl_i = lbls.get(a, len(lbls))
                            #     lbls[a] = lbl_i
                            #     assn_i.append(lbl_i)

                        if objective.is_discrete():
                            assn_i.extend(groups[S.name])
                        else:
                            assn_i = list(range(len(a)))

                        pcp = step.pcp.copy()
                        pcp.extender = cfg
                        print("Direct splitting....")
                        # pcp.extender.depth_max = pcp.splitter
                        # pcp.extender.depth_min = strategy.bounds.extender.depth_min
                        ret = splits.split_all_partitions(
                            topo,
                            pcp,
                            a,
                            assn_i,
                            gcd=gcd,
                            maxmoves=0,
                        )

                        for p_j, (Sj, Sj0, matches, unmatches) in enumerate(ret.value, 0):
                            print(f"Found {p_j+1} {gcd.smarts_encode(Sj)}")

                            edits = 0
                            matches = [
                                y
                                for x, y in enumerate(aa)
                                if x in matches
                            ]
                            unmatches = [
                                y
                                for x, y in enumerate(aa)
                                if x in unmatches
                            ]
                            matched_assn = tuple((assn[i] for i in matches))
                            unmatch_assn = tuple((assn[i] for i in unmatches))

                            new_candidates_direct[(step.overlap[0], None, p_j)] = (
                                S,
                                graphs.subgraph_as_structure(Sj, topo),
                                step,
                                matches,
                                unmatches,
                                matched_assn,
                                unmatch_assn,
                            )
                        if len(new_candidates_direct):
                            direct_success = True
                        else:
                            print("Direct found nothing")


                if step.iterative_enable and not direct_success:
                    # Q = mapper.union_list_parallel(
                    #     list(a.values()),
                    #     reference=S0,
                    #     max_depth=graphs.structure_max_depth(S0),
                    #     icd=icd if len(a) > 100000 else None
                    # )

                    # here I need to get the graphs from the gdb
                    # where aa refers to the indices
                    Q = mapper.union_list_parallel(
                        G, aa, topo,
                        # list(a.values()),
                        reference=S0,
                        max_depth=graphs.structure_max_depth(S0),
                        icd=icd
                        # icd=icd if len(a) > 100000 else None
                    )

                    t = datetime.datetime.now()
                    print(f"{t} Union is {gcd.smarts_encode(Q)}")

                    return_matches = config.splitter.return_matches
                    config.splitter.return_matches = True

                    ret = splits.split_structures_distributed(
                        config.splitter,
                        S0,
                        G,
                        aa,
                        wq,
                        icd,
                        Q=Q,
                    )
                    config.splitter.return_matches = return_matches

                    # backmap = {i: j for i, j in enumerate(cst.mappings[S.name])}
                    print(
                        f"{datetime.datetime.now()} Collecting new candidates"
                    )
                    new_candidates = clusters.clustering_collect_split_candidates_serial(
                        S, ret, step
                    )
                    for k in new_candidates:
                        v = list(new_candidates[k])
                        v[1] = graphs.subgraph_as_structure(new_candidates[k][1], topo)
                        new_candidates[k] = tuple(v)

                p_j_max = -1
                if candidates:
                    p_j_max = max(x[2] for x in candidates) + 1
                for k, v in new_candidates.items():
                    k = (k[0], k[1], k[2]+p_j_max)
                    candidates[k] = v
                new_candidates.clear()

                p_j_max = -1
                if candidates:
                    p_j_max = max(x[2] for x in candidates) + 1
                for k, v in new_candidates_direct.items():
                    k = (k[0], k[1], k[2]+p_j_max)
                    candidates[k] = v
                new_candidates_direct.clear()

            elif step.operation == strategy.MERGE:

                for p_j, jidx in enumerate(hidx.index.below[S.index]):
                    J = hidx.index.nodes[jidx]
                    if J.type != "parameter":
                        continue

                    for overlap in step.overlap:
                        key = (overlap, macro.cursor, p_j)
                        cnd = (S, J, step, None, None, None, None)
                        candidates[key] = cnd

        if step is None:
            print(
                f"{datetime.datetime.now()} Warning, this macro step had no micro steps"
            )
            continue

        print(f"{datetime.datetime.now()} Scanning done.")

        print(datetime.datetime.now())
        print(f"\n\nGenerating SMARTS on {len(candidates)}")
        Sj_sma = []

        with multiprocessing.pool.Pool(configs.processors) as pool:
            if step.operation == strategy.SPLIT:
                # Sj_lst = [candidates[x[1]][1] for x in pq]
                # Sj_lst = [graphs.subgraph_as_structure(x[1], topo) for x in candidates.values()]
                Sj_lst = [x[1] for x in candidates.values()]
            elif step.operation == strategy.MERGE:
                Sj_lst = [
                    graphs.subgraph_as_structure(mm.chemical_system_get_node_hierarchy(csys, x[1]).subgraphs[x[1].index], mm.chemical_system_get_node_hierarchy(csys, x[1]).topology)
                    for x in candidates.values()
                ]
                # Sj_lst = [
                #     graphs.subgraph_as_structure(cst.hierarchy.subgraphs[x[1].index], topo)
                #     for x in candidates.values()
                # ]
            Sj_sma = pool.map_async(gcd.smarts_encode, Sj_lst).get()
            del Sj_lst

        # at this point, I need to iterate over candidates, create the new
        # csys, rebuild the psys, and calculate the tier score



        # print(f"{datetime.datetime.now()} Labeling")
        # cur_assignments = labeler.assign(cst.hierarchy, gcd, smiles, topo)
        # print(f"{datetime.datetime.now()} Rebuilding assignments")
        # cur_mappings = clustering_build_assignment_mappings(
        #     cst.hierarchy, cur_assignments
        # )
        # cur_cst = smarts_clustering(
        #     cst.hierarchy.copy(), cur_assignments, cur_mappings
        # )
        # print(f"{datetime.datetime.now()} Rebuilding mappings")
        # groups = clustering_build_ordinal_mappings(cur_cst, sag)
        # check_lbls_data_selections_equal(cst.group, sag)

        cnd_n = len(candidates)

        t = datetime.datetime.now()

        print_chemical_system(csys)

        visited.clear()
        repeat.clear()


        pq_idx = 0
        procs = (
            os.cpu_count() if configs.processors is None else configs.processors
        )
        
        n_keep = None

        print(f"Scoring and filtering {len(candidates)} candidates for operation={step.operation}")
        for t, tier in enumerate(tiers):
            print(f"Tier {t}: Scoring and filtering {len(candidates)} candidates for operation={step.operation}")
            cnd_keys = {i: k for i, k in enumerate(candidates, 1)}

            fitkeys = objective_tier_get_keys(tier, csys)
            fitkeys = [k for k in fitkeys if k[1] in "skeler" and tier.key_filter(k)]
            fitting_models = set([x[0] for x in fitkeys])
            reuse=[k for k,_ in enumerate(csys.models) if k not in fitting_models]

            # reuse = [x for x in range(len(csys.models)) if x != cid]
            tier_psystems = {
                i: mm.chemical_system_to_physical_system(
                    csys,
                    psystems[i].models[0].positions,
                    ref=psystems[i],
                    reuse=reuse
                ) for i in psystems
            }
            reuse = [x for x in range(len(csys.models))]
            shm = compute.shm_local(1, data={
                "objective": tier,
                "csys": csys,
                "gdb": gdb,
                "reuse": reuse,
                "psysref": tier_psystems,
            }) 

            j = tuple(tier.objectives)
            iterable = {
                # (i, j): ((S, Sj, step.operation, edits, [j]), {})
                (i, j): ((S, Sj, step.operation, edits, j), {})
                for i, (
                    (edits, _, p_j),
                    (S, Sj, step, _, _, _, _),
                ) in enumerate(candidates.items(), 1)
                # ) in enumerate(candidates.items(), 1) for j in tiers[0].objectives
            }
            print(datetime.datetime.now(), f"Generated {len(candidates)} x {len(tiers[0].objectives)//len(j)} = {len(iterable)} candidate evalulation tasks")

            chunksize = 10

            if n_ics > 100000000:
                procs = max(1, procs // 10)
            elif n_ics > 50000000:
                procs = max(1, procs // 5)
            elif n_ics > 10000000:
                procs = max(1, procs // 3)
            elif n_ics > 5000000:
                procs = max(1, procs // 2)
            if n_ics > len(candidates)*10:
                shm.procs_per_task = 0
                chunksize = 1

            addr = ("", 0)
            if len(iterable)*(shm.procs_per_task or procs) <= procs:
                addr = ('127.0.0.1', 0)
                procs=len(iterable)

            if configs.processors == 1:
                work = {}
                for k, v in iterable.items():
                    r = calc_tier_distributed(*v[0], **v[1], shm=shm)
                    work[k] = r
            else:
                ws = compute.workqueue_new_workspace(wq, address=addr, nproc=procs, shm=shm)
                # # this modifies the csys, relabels and computes objective
                work = compute.workspace_submit_and_flush(
                    ws,
                    calc_tier_distributed,
                    iterable,
                    chunksize,
                    1.0,
                    len(iterable),
                )
                ws.close()
                ws = None
            # now just sum over the jobs
            # return keep, X, obj, match_len
            work_new = {}
            for i, _ in enumerate(candidates, 1):
                if i not in work_new:
                    work_new[i] = [0, 0, 0, 0]
                for ij, j in work:
                    if i == ij:
                        line = work[(i,j)]
                        work_new[i][0] |= int(line[0])
                        work_new[i][1] += line[1]
                        work_new[i][2] += line[2]
                        work_new[i][3] += line[3]
                    
            work_full = work
            work = work_new

            tier_scores = []
            max_line = 0
            for j, cnd_i in enumerate(sorted(work), 1):
                (keep, X, obj, match_len) = work[cnd_i]
                # cnd_i, key, unit = unit
                (S, Sj, step, _, _, _, _) = candidates[cnd_keys[cnd_i]]

                cX = X + obj
                dcX = X + obj - X0
                if keep:
                    heapq.heappush(tier_scores, (cX, cnd_i))
            if tier.accept:
                accepted_keys = [
                    x[1] for x in heapq.nsmallest(tier.accept, tier_scores)
                ]
                print(f"Accepted {len(accepted_keys)} candidates from tier:")
                for j, cnd_i in enumerate(sorted(work), 1):
                    if cnd_i not in accepted_keys:
                        continue
                    (keep, X, obj, match_len) = work[cnd_i]
                    # cnd_i, key, unit = unit
                    (S, Sj, step, _, _, _, _) = candidates[cnd_keys[cnd_i]]
                    cX = X + obj
                    dcX = X + obj - X0
                    K = "Y" if keep else "N"
                    cout_line = (
                        f"Cnd. {cnd_i:4d}/{len(work)}"
                        f" {S.name:6s}  " 
                        f" X= {cX:10.5f}"
                        f" X0= {X0:10.5f}"
                        f" dX= {dcX:10.5f} N= {match_len:6d} K= {K} {Sj_sma[cnd_i-1]}"
                    )
                    max_line = max(len(cout_line), max_line)
                    # print(datetime.datetime.now())
                    print(cout_line)
                    sys.stdout.flush()

                Sj_sma = [Sj_sma[k-1] for k in accepted_keys]
                candidates = {
                    cnd_keys[k]: candidates[cnd_keys[k]]
                        for k in accepted_keys
                }
                cnd_keys = {i: k for i, k in enumerate(candidates, 1)}
        print(f"Scanning {len(candidates)} candidates for operation={step.operation}")

        macroamt = strategy.macro_accept_max_total
        macroampc = strategy.macro_accept_max_per_cluster
        microamt = strategy.micro_accept_max_total
        microampc = strategy.micro_accept_max_per_cluster

        cnd_n = len(candidates)
        cnd_keys = {i: k for i, k in enumerate(candidates, 1)}
        n_added = 0
        added = True
        kept = set()
        macro_count = collections.Counter()
        ignore = set()
        reuse_cnd = {}
        # wq = compute.workqueue_local("", configs.workqueue_port)
        # print(f"{datetime.datetime.now()} workqueue started on {wq.mgr.address}")
        n_nano = 0
        while added:

            case1 = macroamt == 0 or n_added < macroamt
            case2 = macroampc == 0 or all([x < macroampc for x in macro_count.values()])
            if not (case1 and case2):
                break

            n_nano += 1

            added = False
            best = {}

            cout = {}
            cout_sorted_keys = []


            shm = compute.shm_local(1, data={
                "objective": initial_objective,
                "csys": csys,
                "gdb": gdb,
                "reuse": reuse0,
                "psysref": psystems,
            }) 

            j = tuple(initial_objective.objectives)
            iterable = {
                # (i, j): ((S, Sj, step.operation, edits, [j]), {})
                (i, j): ((S, Sj, step.operation, edits, j), {})
                for i, (
                    (edits, _, p_j),
                    (S, Sj, step, _, _, _, _),
                ) in enumerate(candidates.items(), 1)
                # ) in enumerate(candidates.items(), 1) for j in tiers[0].objectives
            }
            print(datetime.datetime.now(), f"Generated {len(candidates)} x {len(tiers[0].objectives)//len(j)} = {len(iterable)} candidate evalulation tasks")

            chunksize = 10

            if n_ics > 100000000:
                procs = max(1, procs // 10)
            elif n_ics > 50000000:
                procs = max(1, procs // 5)
            elif n_ics > 10000000:
                procs = max(1, procs // 3)
            elif n_ics > 5000000:
                procs = max(1, procs // 2)
            if n_ics > len(candidates)*10:
                shm.procs_per_task = 0
                chunksize = 1

            addr = ("", 0)
            if len(iterable)*(shm.procs_per_task or procs) <= procs:
                addr = ('127.0.0.1', 0)
                procs=len(iterable)



            for k in kept:
                for iterkey in list(iterable):
                    if k == iterkey[0]:
                        iterable.pop(iterkey)

            for k in ignore:
                for iterkey in list(iterable):
                    if k == iterkey[0]:
                        iterable.pop(iterkey)

            for k in reuse_cnd:
                for iterkey in list(iterable):
                    if k == iterkey[0]:
                        iterable.pop(iterkey)

            # for k, v in list(candidates.items()):
            #     S = v[0]
            #     if S.name in cur_cst.mappings:
            #         if objective.single([assn[i] for i in cur_cst.mappings[S.name]]) == 0.0:
            #             if k in iterable:
            #                 iterable.pop(k)
            #     elif k in iterable:
            #         iterable.pop(k)

            if step.operation != strategy.MERGE and (macroamt + macroampc + microamt + microampc == 0):
                # use X0 so later dX will be 0 and kept
                # if we do this for merges, every merge will be taken..
                work = {i: (1, X0, 0.0, 1) for i in cnd_keys}

            else:
                if configs.processors == 1:
                    work = {}
                    for k, v in iterable.items():
                        r = calc_tier_distributed(*v[0], **v[1], shm=shm)
                        work[k] = r
                else:
                    ws = compute.workqueue_new_workspace(wq, address=addr, nproc=procs, shm=shm)
                    # # this modifies the csys, relabels and computes objective
                    work = compute.workspace_submit_and_flush(
                        ws,
                        calc_tier_distributed,
                        iterable,
                        chunksize,
                        1.0,
                        len(iterable),
                    )
                    ws.close()
                    ws = None
                # now just sum over the jobs
                # return keep, X, obj, match_len
                work_new = {}
                for i, _ in enumerate(candidates, 1):
                    if i not in work_new:
                        work_new[i] = [0, 0, 0, 0]
                    for ij, j in work:
                        if i == ij:
                            line = work[(i,j)]
                            work_new[i][0] |= int(line[0])
                            work_new[i][1] += line[1]
                            work_new[i][2] += line[2]
                            work_new[i][3] += line[3]
                        
                work_full = work
                work = work_new
            
            print(f"The unfiltered results of the candidate scan N={len(work)} total={len(iterable)}:")


            best_reuse = None
            if reuse_cnd:
                # just append the best to work and let the loop figure it out
                best_reuse = sorted(reuse_cnd.items(), key=lambda y: (-y[1][0], y[1][1], y[1][2], y[1][3]))[0]
                work[best_reuse[0]] = best_reuse[1]
            
            for j, cnd_i in enumerate(sorted(work), 1):
                (keep, X, obj, match_len) = work[cnd_i]
                # cnd_i, key, unit = unit
                (S, Sj, step, _, _, _, _) = candidates[cnd_keys[cnd_i]]

                dX = X + obj - X0
                keep = keep and dX <= 0.0

                if step.operation == strategy.SPLIT:
                    visited.add(S.name)
                elif step.operation == strategy.MERGE:
                    visited.add(Sj.name)


                reused_line = ""
                if best_reuse is not None and cnd_i == best_reuse[0]:
                    reused_line="*"
                K = "Y" if keep else "N"
                cout_line = (
                    f"Cnd. {cnd_i:4d}/{len(work)}"
                    f" {S.name:6s} {reused_line}" 
                    f" X= {X:10.5f}"
                    f" X0= {X0:10.5f}"
                    f" dX= {dX:10.5f} N= {match_len:6d} K= {K} {Sj_sma[cnd_i-1]}"
                )
                max_line = max(len(cout_line), max_line)
                # print(datetime.datetime.now())
                print('\r' + cout_line, end=" " * (max_line - len(cout_line)))
                sys.stdout.flush()

                if match_len == 0:
                    if step.operation == strategy.SPLIT:
                        keep = False
                        ignore.add(cnd_i)
                        continue

                if not keep:
                    ignore.add(cnd_i)
                    continue


                if cnd_i in kept:
                    ignore.add(cnd_i)
                    continue


                # We prefer to add in this order
                cout_key = None

                # print sorted at the end but only for new
                # this is to speed things up
                cout_key = (-int(keep), X, match_len, cnd_i, S.name)
                cout[cout_key] = cout_line

                # use these below to determine the best ones to keep
                heapq.heappush(cout_sorted_keys, cout_key)

            print("\r" + " " * max_line)

            # if best:
            #     for name, v in best.items():
            #         kept.add(v[0])
            #     added = True
            # else:
            #     break

            # print sorted at the end
            print(f"Nanostep {n_nano}: The filtered results of the candidate scan N={len(cout)} total={len(iterable)}:")
            if len(cout) == 0:
                continue
            ck_i = 1

            cnd_keep = []
            best_params = [x[0] for x in best.values()]
            macroamt = strategy.macro_accept_max_total
            macroampc = strategy.macro_accept_max_per_cluster
            microamt = strategy.micro_accept_max_total
            microampc = strategy.micro_accept_max_per_cluster
            micro_added = 0
            micro_count = collections.Counter()

            while len(cout_sorted_keys):
                ck = heapq.heappop(cout_sorted_keys)

                keeping = "  "

                dX = ck[1] - X0
                case0 = not (strategy.filter_above is not None and strategy.filter_above <= dX)

                if case0:
                    ignore.add(ck[0])

                case1 = macroamt == 0 or n_added < macroamt
                case2 = microamt == 0 or micro_added < microamt
                if case0 and case1 and case2:
                    sname = ck[4]
                    case3 = macroampc == 0 or macro_count[sname] < macroampc
                    case4 = microampc == 0 or micro_count[sname] < microampc
                    if case3 and case4:
                        cnd_keep.append(ck)
                        micro_count[sname] += 1
                        macro_count[sname] += 1
                        micro_added += 1
                        n_added += 1
                # if ck[3] in best_params:
                        keeping = "->"
                        kept.add(ck[0])
                print(f"{keeping} {ck_i:4d}", cout[ck])
                ck_i += 1
            ck = None

            # keys = {x[0]: cnd_keys[x[0]] for x in best.values()}
            keys = {x[3]: cnd_keys[x[3]] for x in cnd_keep}

            # group_number = max(cur_cst.hierarchy.index.nodes)+1
            print(f"Performing {len(keys)} operations")

            csys_ref = copy.deepcopy(csys)
            csys, nodes = perform_operations(
                csys,
                candidates,
                keys,
                Sj_sma,
            )

            print(f"There are {len(nodes)} nodes returned")

            print("Operations per parameter for this micro:")
            print(micro_count)
            print(f"Micro total: {sum(micro_count.values())} should be {micro_added}")

            print("Operations per parameter for this macro:")
            print(macro_count)
            print(f"Macro total: {sum(macro_count.values())} should be {n_added}")

            if len(nodes) == 0:
                success = False
                added = False
                csys = csys_ref
                continue

            # new_assignments = labeler.assign(hidx, gcd, smiles, topo)
            # # print(datetime.datetime.now(), '*** 5')
            # new_match = clustering_build_assignment_mappings(hidx, new_assignments)

            # cst = smarts_clustering(hidx, new_assignments, new_match)

            # groups = clustering_build_ordinal_mappings(cst, sag)

            # if strategy.prune_empty:
            #     prune_count = 0
            #     pruned = []
            #     new_lbls = [n.name for n in nodes.values()]
            #     for n  in list(nodes):
            #         # only consider new nodes since we have other state like
            #         # the tracker that needs to be managed and would be more
            #         # complicated to handle
            #         n = nodes[n]
            #         lbl = n.name

            #         # also consider the case where we occlude the old pattern
            #         m = cst.hierarchy.index.above.get(n.index, None)
            #         if m is None:
            #             continue
            #         m = cst.hierarchy.index.nodes[m]
                    
            #         if (len(groups.get(m.name, [])) == 0 or len(groups[lbl]) == 0) and lbl in new_lbls:
            #             n_list = trees.tree_index_node_remove_by_name(cst.hierarchy.index, lbl)
            #             for n in n_list:
            #                 cst.hierarchy.subgraphs.pop(n.index)
            #                 cst.hierarchy.smarts.pop(n.index)
            #                 prune_count += 1
            #                 pruned.append(n.index)
            #                 micro_count[m.name] -= 1
            #                 micro_count[m.name] = max(micro_count[m.name], 0)
            #                 macro_count[m.name] -= 1
            #                 macro_count[m.name] = max(macro_count[m.name], 0)
            #                 print(f"Pruned {n.name} parent len {len(groups.get(m.name, []))} sj len {len(groups[lbl])}")
            #             # del groups[lbl]
            #             # if lbl in cst.mappings:
            #             #     del cst.mappings[lbl]
                        
            #     for cnd_i in list(nodes):
            #         hent = nodes[cnd_i]
            #         if hent.index in pruned:
            #             del nodes[cnd_i]
            #             ignore.add(cnd_i)
            #             del keys[cnd_i]
            #     if prune_count:
            #     #     groups = clustering_build_ordinal_mappings(cst, sag)
            #         new_assignments = labeler.assign(hidx, gcd, smiles, topo)
            #         # print(datetime.datetime.now(), '*** 5')
            #         new_match = clustering_build_assignment_mappings(hidx, new_assignments)

            #         cst = smarts_clustering(hidx, new_assignments, new_match)

            #         groups = clustering_build_ordinal_mappings(cst, sag)

            #     print(f"Pruned {prune_count} empty nodes; candidates now {len(keys)}/{len(cnd_keep)}")
            #     print(pruned)
            #     del pruned
            #     del prune_count
            #     del new_lbls

            


            success = False
            added = False
            
            fitkeys = objective_tier_get_keys(initial_objective, csys)
            fitkeys = [k for k in fitkeys if k[1] in "skeler" and initial_objective.key_filter(k)]
            fitting_models = set((n.category[0] for n in nodes.values()))
            reuse=[k for k,_ in enumerate(csys.models) if k not in fitting_models]

            # reuse = [x for x in range(len(csys.models)) if x != cid]
            psystems = {
                i: mm.chemical_system_to_physical_system(
                    csys,
                    psystems[i].models[0].positions,
                    ref=psystems[i],
                    reuse=reuse
                ) for i in psystems
            }
            reuse = [x for x in range(len(csys.models))]
            kv0 = {k: mm.chemical_system_get_value(csys, k) for k in fitkeys}
            for k, v in kv0.items():
                print(f"{str(k):20s} | v0 {v:12.6g}")
            kv, P00, P, gp = objective_tier_run(
                initial_objective,
                gdb,
                csys,
                fitkeys,
                psysref=psystems,
                reuse=reuse,
                wq=wq,
                verbose=True
            )
            C = chemical_objective(csys)
            X = P + C
            dX = X - X0
            print(datetime.datetime.now(), f"Accepting objective: {P:13.6g} C={C:13.6g} DX={P+C-X0:13.6g}")
            if dX > 0:
                print(datetime.datetime.now(), f"Objective raised. Skipping")
                success = False
                added = False
                csys = csys_ref
                for c in cnd_keep:
                    ignore.add(c[3])
                    if c[3] in kept:
                        kept.remove(c[3])
                    sname = c[4]
                    micro_count[sname] -= 1
                    macro_count[sname] -= 1
                    micro_added -= 1
                    n_added -= 1
                continue

            for k, v in kv.items():
                v0 = kv0[k]
                # v0 = mm.chemical_system_get_value(csys, k)
                mm.chemical_system_set_value(csys, k, v)
                print(f"{str(k):20s} | New: {v:12.6g} Ref {v0:12.6g} Diff {v-v0:12.6g}")


            recalc = False
            for cnd_i, hent in nodes.items():
                key = keys[cnd_i]

                (S, Sj, step, _, _, _, _) = candidates[key]
                repeat.add(S.name)
                sma = Sj_sma[cnd_i-1]
                kept.add(cnd_i)
                # cnd_i = best[S.name][0] 
                # hent = Sj
                visited.add(hent.name)
                edits = step.overlap[0]
                
                hidx = mm.chemical_system_get_node_hierarchy(csys, S)

                if step.operation == strategy.SPLIT:

                    # obj = edits
                    # if groups[S.name] and groups[hent.name]:
                    #     obj = objective.split(
                    #         groups[S.name], groups[hent.name], overlap=edits
                    #     )

                    # if obj >= 0.0:
                    #     # we get here if we lazily add many params, so now 
                    #     # some will no longer satisfy the constraint
                    #     kept.remove(cnd_i)
                    #     repeat.remove(S.name)
                        
                    #     n_list = trees.tree_index_node_remove_by_name(cst.hierarchy.index, hent.name)
                    #     for n in n_list:
                    #         recalc = True
                    #         cst.hierarchy.subgraphs.pop(n.index)
                    #         cst.hierarchy.smarts.pop(n.index)
                    #         micro_count[m.name] -= 1
                    #         micro_count[m.name] = max(micro_count[m.name], 0)
                    #         macro_count[m.name] -= 1
                    #         macro_count[m.name] = max(macro_count[m.name], 0)
                    #     print(
                    #         f"\n>>>>> Skipping parameter {cnd_i:4d}/{cnd_n}",
                    #         hent.name,
                    #         "parent",
                    #         S.name,
                    #         "Objective",
                    #         f"{X:10.5f}",
                    #         "Delta",
                    #         f"{dX:10.5f}",
                    #         f"Partition {len(cst.mappings[S.name])}|{len(cst.mappings[hent.name])}",
                    #     )
                    #     print(" >>>>>", key, f"Local dObj {obj:10.5f}", sma, end="\n\n")
                    # else:
                    success = True
                    added = True
                    print(
                        f"\n>>>>> New parameter {cnd_i:4d}/{cnd_n}",
                        hent.name,
                        "parent",
                        S.name,
                        "Objective",
                        f"{X:10.5f}",
                        "Delta",
                        f"{dX:10.5f}",
                        # f"Partition {len(cst.mappings[S.name])}|{len(cst.mappings[hent.name])}",
                    )
                    print(" >>>>>", key, f"Local dObj {obj:10.5f}", sma, end="\n\n")

                    repeat.add(hent.name)
                    step_tracker[(hent.category, hent.name)] = 0


                elif step.operation == strategy.MERGE:

                    if (hent.category, hent.name) in step_tracker:
                        step_tracker.pop((hent.category, hent.name))
                    else:
                        print("WARNING", hent.name, "missing from the tracker")

                    visited.add((S.category, S.name))
                    if (hent.category, hent.name) in visited:
                        visited.remove((hent.category, hent.name))

                    above = hidx.index.above.get(S.index)
                    if above is not None:
                        repeat.add((hidx.index.nodes[above].category, hidx.index.nodes[above].name))

                    success = True
                    added = True
                    print(
                        f">>>>> Delete parameter {cnd_i:4d}/{cnd_n}",
                        hent.name,
                        "parent",
                        S.name,
                        "Objective",
                        f"{X:10.5f}",
                        "Delta",
                        f"{dX:10.5f}",
                    )
                    print(" >>>>>", key, f"Local dObj {obj:10.5f}", sma, end="\n\n")

            if False and recalc:
                print("Detected change in result")
                print("Operations per parameter for this micro:")
                print(micro_count)
                print(f"Micro total: {sum(micro_count.values())}")

                print("Operations per parameter for this macro:")
                print(macro_count)
                print(f"Macro total: {sum(macro_count.values())}")

                new_assignments = labeler.assign(hidx, gcd, smiles, topo)
                # print(datetime.datetime.now(), '*** 5')
                new_match = clustering_build_assignment_mappings(hidx, new_assignments)

                cst = smarts_clustering(hidx, new_assignments, new_match)

                groups = clustering_build_ordinal_mappings(cst, sag)


            # print the tree
            print_chemical_system(csys)
            # mod_lbls = cluster_assignment.smiles_assignment_str_modified(
            #     cur_cst.group.assignments, cst.group.assignments
            # )
            # repeat.update(mod_lbls)
            # cur_cst = cst
            # cst = None
            X0 = X
            P0 = P
            C0 = C

        # wq.close()
        # wq = None
        

        if strategy.macro_accept_max_total > 0 and n_added > 0:
            strategy.repeat_step()

        print(f"There were {n_added} successful operations")
        # csys = cur_csys
        # for ei, e in enumerate(
        #     tree_iterators.tree_iter_dive(
        #         csys.hierarchy.index,
        #         trees.tree_index_roots(cur_cst.hierarchy.index),
        #     )
        # ):
        #     s = trees.tree_index_node_depth(cur_cst.hierarchy.index, e)
        #     obj_repo = ""
        #     # if groups[e.name]:
        #     obj_repo = objective.report(groups[e.name])
        #     print(
        #         f"** {s:2d} {ei:3d} {e.name:4s}",
        #         obj_repo,
        #         cur_cst.hierarchy.smarts.get(e.index),
        #     )

        print(f"{datetime.datetime.now()} Visited", visited)
        # for name in ((node.category, node.name) for node in hidx.index.nodes.values()):
        #     if name not in step_tracker:
        #         continue

        #     if name not in repeat:
        #         step_tracker[name] = max(strategy.cursor, step_tracker[name])
        #     else:
        #         print(f"Assignments changed for {name}, will retarget")
        #         step_tracker[name] = 0

        pickle.dump([gdb, csys, strategy, psystems], open("chk.cst.p", "wb"))


    # new_assignments = labeler.assign(cst.hierarchy, gcd, smiles, topo)
    # mappings = clustering_build_assignment_mappings(
    #     cst.hierarchy, new_assignments
    # )
    # cst = smarts_clustering(cst.hierarchy, new_assignments, mappings)
    # pickle.dump([gdb, csys, strategy, psys], open("chk.cst.p", "wb"))

    ended = datetime.datetime.now()

    print(f"Start time: {started}")
    print(f"End   time: {ended}")

    # gc.enable()
    # return csys
    wq.close()

    print_chemical_system(csys)

    return csys, P, C

def calc_tier_distributed(S, Sj, operation, edits, oid, shm=None):
    csys = copy.deepcopy(shm.csys)

    hidx = mm.chemical_system_get_node_hierarchy(csys, S)
    cm = mm.chemical_system_get_node_model(csys, S)
    cid = S.category[0]
    pid = S.category[1]
    uid = S.category[2]

    # node_ref = trees.tree_node_copy(S)
    parentid = S.index

    labeler = csys.perception.labeler
    gcd = csys.perception.gcd
    objective = shm.objective

    reuse = shm.reuse
    gdb = assignments.graph_db_get_entries(shm.gdb, [*set((e for o in oid for e in objective.objectives[o].addr.eid))])
    
    psysref = None
    if shm.psysref:
        psysref = {eid: shm.psysref[eid] for eid in gdb.entries}
    else:
        psysref = gdb_to_physical_systems(gdb, csys)
    
    keep = True
    kv, y0, X, gx, C = {}, 0, 0, [], 0
    # need to perform the operation and then add to keys
    # would also need to add the node to the FF
    if operation == optimization.optimization_strategy.SPLIT:
        node = mm.chemical_model_smarts_hierarchy_copy_node(cm, cid, pid, uid, S, None)
        hidx.subgraphs[node.index] = Sj
        sma = gcd.smarts_encode(Sj)
        hidx.smarts[node.index] = sma


    elif operation == optimization.optimization_strategy.MERGE:
        node = Sj
        sma = hidx.smarts[node.index]
        mm.chemical_model_smarts_hierarchy_remove_node(cm, cid, pid, uid, Sj)

    reuse = [x for x in range(len(csys.models)) if x != cid]
    psysref = {
        i: mm.chemical_system_to_physical_system(
            csys,
            psysref[i].models[0].positions,
            ref=psysref[i],
            reuse=reuse
        ) for i in psysref
    }
    reuse = [x for x in range(len(csys.models))]

    keys = mm.chemical_system_iter_keys(csys)
    # print_chemical_system(csys)
    keys = [k for k in keys if objective.key_filter(k)]
    assigned = [(i, k, v) for psys in psysref.values() for i,m in enumerate(psys.models) for a in m.labels[0].values() for k, v in a.items()]
    keys = [k for k in keys if tuple(k[:3]) in assigned or k[1] in "s"]
    # print("Fitting keys are:")
    # print(keys)

    if operation == optimization.optimization_strategy.SPLIT:

        match_len = 0
        old_match = 0
        for psys in psysref.values():
            pm = psys.models[cid]
            for ic, lbls in pm.labels[0].items():
                lbls = set(lbls.values())
                if node.name in lbls:
                    match_len += 1
                if S.name in lbls:
                    old_match += 1

        new_keys = mm.chemical_system_smarts_hierarchy_get_node_keys(cm, cid, pid, uid, node)
        if old_match == 0 or match_len == 0:   
            keep = False
    elif operation == optimization.optimization_strategy.MERGE: 

        match_len = 0
        for psys in psysref.values():
            pm = psys.models[cid]
            for ic, lbls in pm.labels[0].items():
                lbls = set(lbls.values())
                if S.name in lbls:
                    match_len += 1

    # print("Matches", match_len, "Old matches", old_match)
    # C = graphs.graph_bits(Sj) / len(Sj.nodes) /1000 + len(Sj.nodes)
    if keep:
        C = mm.chemical_system_smarts_complexity(csys)
        kv, y0, X, gx = objective_tier_run(objective, gdb, csys, keys, oid=oid, psysref=psysref, reuse=reuse, wq=None, verbose=False)
    print(f"{S.name}->{sma:40s} OID={oid} {keep} {X} {C} {match_len}")

    return keep, X, C, match_len 

def print_chemical_system(csys):
    print("Model:")
    for ei, hidx in enumerate(
        mm.chemical_system_iter_smarts_hierarchies(csys)
    ):
        print("Tree:")
        for root in trees.tree_index_roots(hidx.index):
            for e in tree_iterators.tree_iter_dive(hidx.index, root):
                s = trees.tree_index_node_depth(hidx.index, e)
                obj_repo = ""
                # if groups[e.name]:
                # obj_repo = objective.report(groups[e.name])
                print(
                    f"** {s:2d} {ei:3d} {e.name:4s}",
                    hidx.smarts.get(e.index, ""),
                )


def chemical_objective(csys):
    return mm.chemical_system_smarts_complexity(csys)

def perform_operations(
        csys: mm.chemical_system,
        candidates,
        keys,
        Sj_sma,
    ):


    nodes = {}
    for cnd_i, key in keys.items():
        (S, Sj, step, _, _, _, _) = candidates[key]
        (edits, _, p_j) = key
        # param_name = "p."
        sma = ""
        added = False
        cid, pid, uid = S.category
        cm = csys.models[cid]

        hidx = mm.chemical_system_get_node_hierarchy(csys, S)
        topo = hidx.topology

        if step.operation == optimization.optimization_strategy.SPLIT:
            # param_name = "p" + str(group_number)
            node = mm.chemical_model_smarts_hierarchy_copy_node(cm, cid, pid, uid, S, None)

            # print(datetime.datetime.now(), '*** 2')
            assert Sj.select
            Sj = graphs.subgraph_relabel_nodes(
                Sj, {n: i for i, n in enumerate(Sj.select, 1)}
            )
            # Sj = graphs.subgraph_to_structure(Sj, topo)

            hidx.subgraphs[node.index] = graphs.subgraph_copy(Sj)
            hidx.smarts[node.index] = Sj_sma[cnd_i-1]

            nodes[cnd_i] = node

        elif step.operation == optimization.optimization_strategy.MERGE:

            mm.chemical_model_smarts_hierarchy_remove_node(cm, cid, pid, uid, Sj)
            nodes[cnd_i] = Sj

    return csys, nodes
