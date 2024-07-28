"""
besmarts.mechanics.fits
"""

import os
import copy
import heapq
import datetime
import pickle
import multiprocessing
import collections
import sys
import math
import itertools
import numpy as np
from typing import List, Dict, Tuple, Callable
from besmarts.core import configs
from besmarts.core import arrays
from besmarts.core import assignments
from besmarts.core import topology
from besmarts.core import graphs
from besmarts.core import codecs
from besmarts.core import geometry
from besmarts.core import compute
from besmarts.core import clusters
from besmarts.core import mapper
from besmarts.core import splits
from besmarts.core import trees
from besmarts.core import tree_iterators
from besmarts.core import optimization
from besmarts.core import primitives
from besmarts.core import logs
from besmarts.core import returns
from besmarts.cluster import cluster_objective
from besmarts.cluster import cluster_optimization
from besmarts.core import array_numpy
from besmarts.mechanics import optimizers_scipy
from besmarts.mechanics import optimizers_openmm
from besmarts.mechanics import objectives
from besmarts.mechanics import hessians
from besmarts.mechanics import vibration
from besmarts.mechanics import molecular_models as mm


class compute_config:

    """
    Default is single point energy
    """

    def __init__(self, addr):
        self.addr = addr 
        self.keys = []

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_table]]:
        # print("Starting calculation", self)
        csys = self.csys
        GDB = self.GDB
        tbl_idx = assignments.ENERGY
        tid = assignments.POSITIONS
        

        all_results = []
        for keyset in self.keys:
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []

            # for k, v in self.keys.items():
                # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.6g} to {v:.6g} d={v-mm.chemical_system_get_value(csys, k):.6g}")
                # mm.chemical_system_set_value(csys, k, v)
                # mm.physical_system_set_value(self.psys, k, v)

            for eid, gdb in GDB.entries.items():
                tbl = assignments.graph_db_table(topology.null)
                for gid in gdb.graphs:
                    # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, tbl_idx, gid)
                    rid = 0
                    pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)

                    # this means fold all graphs into a single graph
                    # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, eid, tid)

                    psys = mm.chemical_system_to_physical_system(csys, [pos0], ref=self.psys, reuse=self.reuse)
                    reapply = set()
                    for k, v in keyset.items():
                        # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.10g} to {v:.10g} d={v-mm.chemical_system_get_value(csys, k):.10g}")
                        # mm.chemical_system_set_value(csys, k, v)
                        # mm.chemical_system_set_value(csys, k, v)
                        mm.physical_system_set_value(psys, k, v)
                        reapply.add(k[0])
                    for m in reapply:
                        procs = csys.models[m].procedures
                        if len(procs) > 1:
                            for _ in range(1, len(psys.models[m].values)):
                                psys.models[m].values.pop()
                                psys.models[m].labels.pop()
                            for proc in procs[1:]:
                                psys.models[m] = proc.assign(csys.models[m], psys.models[m], overrides={k[1:]: v for k, v in keyset.items() if k[0] == m})
                    ene = objectives.physical_system_energy(psys, csys)
                    
                    # ene_ga = assignments.graph_assignment("", {(0,): [[energy]]}, gdb.graphs[gid])
                    # gdg = assignments.graph_assignment_to_graph_db_graph(ene_ga, topology.null)
                    tbl[gid] = gdg
                r = {tbl_idx: tbl}
                # print("Done calculating", self)
                # print("Result is\n", r)
                results.append(r)
            all_results.append(results)
        return all_results


class compute_config_energy(compute_config):

    """
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        # print("Starting calculation", self)
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.ENERGY
        tid = assignments.POSITIONS

        # for k, v in self.keys.items():
            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.6g} to {v:.6g} d={v-mm.chemical_system_get_value(csys, k):.6g}")
            # mm.chemical_system_set_value(csys, k, v)

        all_results = []
        for keyset in self.keys:
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []
            for eid, gde in gdb.entries.items():
                tbl = assignments.graph_db_table(topology.null)
                rids = assignments.graph_db_table_get_row_ids(gde[tid])
                for rid in rids:
                    system = []
                    for gid in gdb.graphs:
                    # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, assignments.POSITIONS, gid)

                        pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
                        system.append(pos0)


                    # create a ga but the indices are (gid, rid, sid)
                    # just keep an index? 
                    # system = graphs.graph_assignment_system(system)

                    # psys now has positions in each model indexed by graph
                    # this means I will need to make indices for graphs and confs
                    # such as (gid, sid, rid)

                    psys = mm.chemical_system_to_physical_system(csys, system, ref=self.psys, reuse=self.reuse)
                    reapply = set()
                    for k, v in keyset.items():
                        # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.10g} to {v:.10g} d={v-mm.chemical_system_get_value(csys, k):.10g}")
                        # mm.chemical_system_set_value(csys, k, v)
                        mm.physical_system_set_value(psys, k, v)
                        reapply.add(k[0])
                    for m in reapply:
                        procs = csys.models[m].procedures
                        if len(procs) > 1:
                            for _ in range(1, len(psys.models[m].values)):
                                psys.models[m].values.pop()
                                psys.models[m].labels.pop()
                            for proc in procs[1:]:
                                psys.models[m] = proc.assign(csys.models[m], psys.models[m], overrides={k[1:]: v for k, v in keyset.items() if k[0] == m})

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
            all_results.append(results)
        # print("Done calculating", self)
        return all_results

class compute_config_gradient(compute_config):

    """
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        """
        compute_config_gradient::run
        """
        # print("Starting calculation", self)
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.GRADIENTS
        tid = assignments.POSITIONS

        verbose=True

        # for k, v in self.keys.items():
            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.6g} to {v:.6g} d={v-mm.chemical_system_get_value(csys, k):.6g}")
            # mm.chemical_system_set_value(csys, k, v)

        all_results = []
        for keyset in self.keys:
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []
            for eid, gde in gdb.entries.items():
                tbl = assignments.graph_db_table(topology.atom)
                rids = assignments.graph_db_table_get_row_ids(gde[tid])
                for rid in rids:
                    system = []
                    for gid in gde.tables[tid].graphs:
                    # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, assignments.POSITIONS, gid)

                        pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
                        # psys = mm.chemical_system_to_physical_system(csys, [pos0], ref=self.psys, reuse=self.reuse)
                        system.append(pos0)


                    # create a ga but the indices are (gid, rid, sid)
                    # just keep an index? 
                    # system = graphs.graph_assignment_system(system)

                    # psys now has positions in each model indexed by graph
                    # this means I will need to make indices for graphs and confs
                    # such as (gid, sid, rid)

                    psys = mm.chemical_system_to_physical_system(csys, system, ref=self.psys, reuse=self.reuse)
                    reapply = set()
                    for k, v in keyset.items():
                        # dprint(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.10g} to {v:.10g} d={v-mm.chemical_system_get_value(csys, k):.10g}", on=verbose)
                        # mm.chemical_system_set_value(csys, k, v)
                        mm.physical_system_set_value(psys, k, v)
                        reapply.add(k[0])
                    for m in reapply:
                        procs = csys.models[m].procedures
                        if len(procs) > 1:
                            for _ in range(1, len(psys.models[m].values)):
                                psys.models[m].values.pop()
                                psys.models[m].labels.pop()
                            for proc in procs[1:]:
                                psys.models[m] = proc.assign(csys.models[m], psys.models[m], overrides={k[1:]: v for k, v in keyset.items() if k[0] == m})

                    # this must put everything in the psys into a single system
                    # this must mean that we flatten
                    # pos = optimizers_scipy.optimize_positions_scipy(csys, psys)
                    # psys = mm.chemical_system_to_physical_system(csys, [pos])

                    #gx = objectives.physical_system_gradient(psys, csys)
                    gx = optimizers_openmm.physical_system_gradient_openmm(csys, psys)
                    gx = arrays.array_round(gx, 12)
                    tbl.values.extend(gx)
                    # if not keyset:
                        # ene = objectives.physical_system_energy(psys, csys)
                        # print(f"Calculated energy for EID {eid} is ", ene)
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
            all_results.append(results)
        # print("Done calculating", self)
        return all_results

class compute_config_hessian(compute_config):

    """
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        """
        compute_config_hessian::run

        This should compute the gradient and hessian
        Also calculate the B2 matrix if necessary
        Do not transform, since the point is to just calculate what is needed
        Measurements are done in the objective
        """
        # print("Starting calculation", self)
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.HESSIANS
        tid = assignments.POSITIONS

        verbose=True

        # for k, v in self.keys.items():
            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.6g} to {v:.6g} d={v-mm.chemical_system_get_value(csys, k):.6g}")
            # mm.chemical_system_set_value(csys, k, v)

        all_results = []
        for ki, keyset in enumerate(self.keys, 1):
            # print(f"Hessian keyset {ki} of {len(self.keys)}")
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []
            for eid, gde in gdb.entries.items():
                tbl_hess = assignments.graph_db_table(topology.null)
                tbl_grad = assignments.graph_db_table(topology.null)
                rids = assignments.graph_db_table_get_row_ids(gde[tid])
                for rid in rids:
                    system = []
                    for gid in gde.tables[tid].graphs:
                    # pos0 = assignments.graph_db_entry_to_graph_assignment(gdb, assignments.POSITIONS, gid)

                        pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
                        # psys = mm.chemical_system_to_physical_system(csys, [pos0], ref=self.psys, reuse=self.reuse)
                        system.append(pos0)


                    # create a ga but the indices are (gid, rid, sid)
                    # just keep an index? 
                    # system = graphs.graph_assignment_system(system)

                    # psys now has positions in each model indexed by graph
                    # this means I will need to make indices for graphs and confs
                    # such as (gid, sid, rid)

                    psys = mm.chemical_system_to_physical_system(csys, system, ref=self.psys, reuse=self.reuse)
                    reapply = set()
                    for k, v in keyset.items():
                        # dprint(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.10g} to {v:.10g} d={v-mm.chemical_system_get_value(csys, k):.10g}", on=verbose)
                        # mm.chemical_system_set_value(csys, k, v)
                        mm.physical_system_set_value(psys, k, v)
                        reapply.add(k[0])
                    for m in reapply:
                        procs = csys.models[m].procedures
                        if len(procs) > 1:
                            for _ in range(1, len(psys.models[m].values)):
                                psys.models[m].values.pop()
                                psys.models[m].labels.pop()
                            for proc in procs[1:]:
                                psys.models[m] = proc.assign(csys.models[m], psys.models[m], overrides={k[1:]: v for k, v in keyset.items() if k[0] == m})

                    pos = psys.models[0].positions[0]
                    # this must put everything in the psys into a single system
                    # this must mean that we flatten
                    # pos = optimizers_scipy.optimize_positions_scipy(csys, psys)
                    # psys = mm.chemical_system_to_physical_system(csys, [pos])
                    hess_qm = gdb.entries[eid].tables[assignments.HESSIANS].values

                    # posmm = optimizers_scipy.optimize_positions_scipy(csys, psys, tol=1e-10) 
                    # posmm = optimizers_openmm.optimize_positions_openmm(csys, psys) 

                    # xyz = np.vstack([x[0] for x in posmm.selections.values()], dtype=float)
                    # xyz = xyz.round(12)
                    args, keys = objectives.array_flatten_assignment(pos.selections)
                    # hess_mm = objectives.array_geom_hessian(args, keys, csys, psys, h=1e-7)

                    hess_mm = optimizers_openmm.physical_system_hessian_openmm(csys, psys, h=1e-4)

                    # hess_mm = objectives.physical_system_hessian(psys, csys)
                    hx = []
                    for row in hess_mm:
                        hx.append(arrays.array_scale(row, 1/4.184))
                    tbl_hess.values.extend(hx)

                    grad_mm = optimizers_openmm.physical_system_gradient_openmm(csys, psys)
                    # grad_mm = objectives.array_geom_gradient(args, keys, csys, psys)
                    gx = arrays.array_scale(grad_mm, 1/4.184)
                    tbl_grad.values.extend(gx)

                r = {
                    assignments.HESSIANS: tbl_hess,
                    assignments.GRADIENTS: tbl_grad
                }

                results.append(r)
            all_results.append(results)
        return all_results

    

class compute_config_position(compute_config):

    """
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_table]]:
        """
        compute_config_position::run
        """
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.POSITIONS
        tid = assignments.POSITIONS

        # for k, v in self.keys.items():
            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.10g} to {v:.10g} d={v-mm.chemical_system_get_value(csys, k):.10g}")
            # mm.chemical_system_set_value(csys, k, v)
            # mm.physical_system_set_value(psys, k, v)
        # if not self.keys:
        #     print(f"-> Obj Setting no values")
        # need to make this work for multi conformations
        all_results = []
        for keyset in self.keys:
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]] = []
            for eid, gde in gdb.entries.items():
                tbl = assignments.graph_db_table(topology.atom)
                for gid in gdb.graphs:
                    for rid in gde[tid][gid].rows:
                        pos0 = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
                        # print(f"MinPos Reuse={self.reuse}")
                        # print(f"MinPos pre torsions:")
                        # print(self.psys.models[2].labels)
                        # print(f"MinPos pre torsion values:")
                        # print(self.psys.models[2].values)
                        psys = mm.chemical_system_to_physical_system(csys, [pos0], ref=self.psys, reuse=self.reuse)
                        # print(f"MinPos Reuse={self.reuse}")
                        # print(f"MinPos torsions:")
                        # print(psys.models[2].labels)
                        # print(f"MinPos torsion values:")
                        # print(psys.models[2].values)
                        
                        reapply = set()
                        for k, v in keyset.items():
                            # mm.physical_system_get_value(psys, k, v)
                            # print(f"-> Obj Setting {k} from {mm.chemical_system_get_value(csys, k):.10g} to {v:.10g} d={v-mm.chemical_system_get_value(csys, k):.10g}")
                            # mm.chemical_system_set_value(csys, k, v)
                            mm.physical_system_set_value(psys, k, v)
                            reapply.add(k[0])
                        for m in reapply:
                            procs = csys.models[m].procedures
                            if len(procs) > 1:
                                # print(f"MinPos Reapply {m}")
                                for _ in range(1, len(psys.models[m].values)):
                                    psys.models[m].values.pop()
                                    psys.models[m].labels.pop()
                                for proc in procs[1:]:
                                    psys.models[m] = proc.assign(csys.models[m], psys.models[m], overrides={k[1:]: v for k, v in keyset.items() if k[0] == m})
                        # print(f"Initial xyz for EID {eid}:")
                        # print("\n".join(print_xyz(pos0)))
                        # pos = optimizers_scipy.optimize_positions_scipy(csys, psys, step_limit=self.step_limit, tol=self.tol)
                        pos = optimizers_openmm.optimize_positions_openmm(csys, psys, step_limit=self.step_limit, tol=self.tol)
                        # print(f"Optimized xyz for EID {eid}:")
                        # print("\n".join(print_xyz(pos)))
                        # pos = copy.deepcopy(pos0)
                        # struct = assignments.graph_assignment_to_graph_db_struct(pos, topology.atom)
                        gdg = assignments.graph_assignment_to_graph_db_graph(pos, topology.atom)
                        tbl.graphs[gid] = gdg
                r = {tbl_idx: tbl}
                results.append(r)
            all_results.append(results)
        # print("Done calculating", self)
        return all_results

class compute_config_penalty(compute_config):

    def __init__(self, targets):
        self.targets = targets

    def run(self, keys) -> dict:
        """
        compute_config_penalty::run
        """
        # csys = self.csys

        results = {}
        for k, ref in self.targets.items():
            val = keys.get(k)
            if val is not None:
                results[k] = round(val - ref, 12)
        return results


class objective_config:

    def __init__(self, addr, include=True, scale=1, coordinate_system="C"):
        self.addr = addr
        self.scale = scale
        self.step_limit = 200
        self.tol = 1e-5
        self.coordinate_system = coordinate_system
        self.include = include
        self.cache = {}
        self.batch_size = None

    def get_task(self, GDB, csys, keys=None, psys=None, reuse=None) -> compute_config:
        cc = compute_config(self.addr)

        cc.GDB = {i: GDB[i] for i in self.addr.eid}
        cc.csys = csys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.psys = psys
        cc.reuse = reuse
        return cc

    def compute_gradient_2pt(self, ref, X0, E, D: List[Dict[int, assignments.graph_db_table]], h):
        
        """
        D is the results for some objective idx
        so it is a list of tables (with their own tid)
                                   
        which now means that each is a table with structs which have indices
        this is the parameter deriv, which should be 2(X-X0)dX/dp
        and this is in particular dX/dpj
        """
        dxa = []
        dxb = []
        for etbls, dtbls  in zip(E, D):
            for tid, dtbl in dtbls.items():
                etbl = etbls[tid]
                for gid, dgdg in dtbl.graphs.items():
                    egdg = etbl.graphs[gid]
                    for rid, dgdr in dgdg.rows.items():
                        egdr = egdg.rows[rid]
                        for cid, dgdc in dgdr.columns.items():
                            egdc = egdr.columns[cid]
                            for ev, dv in zip(egdc.selections.values(), dgdc.selections.values()):
                                # vstr = "["+",".join([f"{vi:10.5f}" for vi in v]) + "]"
                                # v0str = "["+",".join([f"{vi:10.5f}" for vi in v0]) + "]"
                                # dxi = arrays.array_difference(dv, ev)
                                dxa.extend(ev)
                                dxb.extend(dv)
                                # dx += sum(x)
                                # print(f"{rid:3d} Position SSE: {sse:10.5f} {vstr} {vstr}")
        # print(f"Total Position SSE: {obj:10.5f} A^2")
        # dx = arrays.array_round(dx, 12)

        dxdp = arrays.array_difference(dxb, dxa)
        dxdp = arrays.array_round(dxdp, 12) 
        dxdp = arrays.array_scale(dxdp, 1/(2*h)) 

        d2xdp2 = arrays.array_add(
            arrays.array_scale(X0, -2),
            arrays.array_add(dxb, dxa)
        )
        d2xdp2 = arrays.array_round(d2xdp2, 12) 
        d2xdp2 = arrays.array_scale(d2xdp2, 1.0/(h*h))

        return dxdp, d2xdp2

    def compute_diff(self, GDB: assignments.graph_db, D: List[Dict[int, assignments.graph_db_table]], verbose=False):
        """
        D is the results for some objective idx
        so it is a list of tables (with their own tid)
                                   
        which now means that each is a table with structs which have indices
        """
        obj = []
        out = []
        for (eid, gde), tbls  in zip(GDB.entries.items(), D):
            rmse = []
            out.append(f"EID {eid} Position objective:")
            for tid, tbl in tbls.items():
                tbl0 = gde.tables[tid]
                for gid, gdg in tbl.graphs.items():
                    gdg0 = tbl0.graphs[gid]
                    for rid, gdr in gdg.rows.items():
                        gdr0 = gdg0.rows[rid]
                        for cid, gdc in gdr.columns.items():
                            gdc0 = gdr0.columns[cid]

                            # assume cartesian
                            x1l = [xyz for xyz in gdc.selections.values()]
                            x0l = [xyz for xyz in gdc0.selections.values()]

                            # print(f"{rid:3d} F= {vstr} F0= {v0str}")
                            for ic, x1, x0 in zip(gdc0.selections, x1l, x0l):
                                x = arrays.array_difference(x1, x0)
                                x = arrays.array_round(x, 12)
                                obj.extend(x)
                                xyz1 = " ".join([f"{xi:8.4f}" for xi in x1])
                                xyz0 = " ".join([f"{xi:8.4f}" for xi in x0])
                                dxyz = " ".join([f"{xi:8.4f}" for xi in  x])
                                rmse.append(arrays.array_inner_product(x, x))
                                out.append(f"    {eid:4d} {gid:2d} {rid:2d} {cid:2d} {ic[0]:3d} MM: {xyz1} QM: {xyz0} dP: {dxyz}")
            if len(rmse):
                rmse = (sum(rmse)/len(rmse))**.5
            else:
                rmse = 0
            out.append(f"Total Position SSE: {arrays.array_inner_product(obj, obj):10.5f} A^2 RMSE: {rmse:10.5f} A")
        return returns.success(arrays.array_round(obj, 12), out=out, err=[])

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
                                # vstr = "["+",".join([f"{vi:10.5f}" for vi in v]) + "]"
                                # v0str = "["+",".join([f"{vi:10.5f}" for vi in v0]) + "]"
                                x = arrays.array_difference(v, v0)
                                sse = arrays.array_inner_product(x, x)
                                obj += sse
                                # print(f"{rid:3d} Position SSE: {sse:10.5f} {vstr} {vstr}")
        # print(f"Total Position SSE: {obj:10.5f} A^2")
        return obj

class objective_config_gradient(objective_config):

    def get_task(self, GDB: assignments.graph_db, csys, keys=None, psys=None, reuse=None) -> compute_config:
        cc = compute_config_gradient(self.addr)
        cc.GDB = assignments.graph_db_get(GDB, self.addr)

        cc.csys = csys
        cc.psys = psys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.reuse = reuse
        return cc

    def compute_gradient_2pt(self, ref, X0, E0, E1, h):

        dxa = []
        dxb = []
        gde = ref.entries[self.addr.eid[0]]
        for etbls, dtbls  in zip(E0, E1):
            for tid, dtbl in dtbls.items():
                x0 = etbls[tid].values
                x1 = dtbl.values
                if self.coordinate_system == "D":
                    gid = list(gde.tables[assignments.POSITIONS].graphs)[0]
                    g = ref.graphs[gid]
                    smiles = ref.smiles[gid]
                    sel = gde.tables[assignments.POSITIONS].graphs[gid].rows[0].columns[0].selections
                    sel = {k: [v] for k, v in sel.items()}
                    pos = assignments.graph_assignment(smiles, sel, g)
                    x0 = array_numpy.dlcmatrix_project_gradients(pos, x0)
                    x1 = array_numpy.dlcmatrix_project_gradients(pos, x1)

                dxa.extend(x0)
                dxb.extend(x1)
                # x = arrays.array_difference(x1, x0)
                # dx.extend(x)
        dxdp = arrays.array_difference(dxb, dxa)
        dxdp = arrays.array_round(dxdp, 12) 
        dxdp = arrays.array_scale(dxdp, 1/(2*h)) 

        d2xdp2 = arrays.array_add(
            arrays.array_scale(X0, -2),
            arrays.array_add(dxb, dxa)
        )
        d2xdp2 = arrays.array_round(d2xdp2, 12) 
        d2xdp2 = arrays.array_scale(d2xdp2, 1.0/(h*h))

        return dxdp, d2xdp2

    def compute_diff(self, GDB: assignments.graph_db, D: List[Dict[int, assignments.graph_db_table]], verbose=False):
        """
        """

        out = []
        X = []
        for (eid, gde), tbls  in zip(GDB.entries.items(), D):
            rmse = []
            for tid, tbl in tbls.items():
                x0 = gde.tables[tid].values
                x1 = tbl.values

                if self.coordinate_system == "D":
                    gid = list(gde.tables[assignments.POSITIONS].graphs)[0]
                    g = GDB.graphs[gid]
                    smiles = GDB.smiles[gid]
                    sel = gde.tables[assignments.POSITIONS].graphs[gid].rows[0].columns[0].selections
                    sel = {k: [v] for k, v in sel.items()}
                    pos = assignments.graph_assignment(smiles, sel, g)

                    x1 = array_numpy.dlcmatrix_project_gradients(pos, x1)

                    if "D:x0" not in self.cache:
                        x0 = array_numpy.dlcmatrix_project_gradients(pos, x0)
                        self.cache['D:x0'] = x0
                    else:
                        x0 = self.cache['D:x0']


                x = arrays.array_difference(x1, x0)
                nlines = len(out)
                out.append(f" EID {eid} Gradient objective:")
                for i in range(len(x1)//3):
                    xx = x[3*i:3*i+3]
                    xyz1 = " ".join([f"{xi:12.4e}" for xi in x1[3*i:3*i+3]])
                    xyz0 = " ".join([f"{xi:12.4e}" for xi in x0[3*i:3*i+3]])
                    dxyz = " ".join([f"{xi:12.4e}" for xi in xx])
                    out.append(f"    {i+1:4d} MM: {xyz1} QM: {xyz0} dG: {dxyz}")
                    rmse.append(arrays.array_inner_product(xx, xx))
                if verbose:
                    for i in range(len(out) - nlines):
                        print(out[i])
                X.extend(x)
            if rmse:
                rmse = (sum(rmse)/len(rmse))**.5
            else:
                rmse = 0

            out.append(f" Total Gradient SSE: {arrays.array_inner_product(X,X):10.5f} (kJ/mol/A)^2 RMSE: {rmse:10.5f} kJ/mol/A")
        return returns.success(arrays.array_round(X, 12), out=out, err=[])

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



class objective_config_hessian(objective_config):


    def get_task(self, GDB: assignments.graph_db, csys, keys=None, psys=None, reuse=None) -> compute_config:
        cc = compute_config_hessian(self.addr)
        cc.GDB = assignments.graph_db_get(GDB, self.addr)

        cc.csys = csys
        cc.psys = psys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.reuse = reuse
        return cc

    def compute_gradient_2pt(self, ref, X0, E0, E1, h):

        dxa = []
        dxb = []
        gde = ref.entries[self.addr.eid[0]]

        for etbls, dtbls  in zip(E0, E1):
            # for tid, dtbl in dtbls.items():

            # x0 = etbls[tid].values
            # x1 = dtbl.values

            gid = list(gde.tables[assignments.POSITIONS].graphs)[0]
            g = ref.graphs[gid]
            DL = self.cache['IC:dl']
            (ics, B, B2) = self.cache['IC:b']

            hess0 = dtbls[assignments.HESSIANS].values
            grad0 = list(dtbls[assignments.GRADIENTS].values)
            x0 =  hessians.hessian_frequencies(g, hess0, grad0, DL, ics, B, B2)

            hess1 = etbls[assignments.HESSIANS].values
            grad1 = list(etbls[assignments.GRADIENTS].values)
            x1 =  hessians.hessian_frequencies(g, hess1, grad1, DL, ics, B, B2)

            # keep this for when I want to try internal hess
            # x0 = array_numpy.dlcmatrix_project_gradients(pos, x0)
            # x1 = array_numpy.dlcmatrix_project_gradients(pos, x1)

            dxa.extend(x0)
            dxb.extend(x1)
            # x = arrays.array_difference(x1, x0)
            # dx.extend(x)
        dxdp = arrays.array_difference(dxb, dxa)
        dxdp = arrays.array_round(dxdp, 12) 
        dxdp = arrays.array_scale(dxdp, 1/(2*h)) 

        d2xdp2 = arrays.array_add(
            arrays.array_scale(X0, -2),
            arrays.array_add(dxb, dxa)
        )
        d2xdp2 = arrays.array_round(d2xdp2, 12) 
        d2xdp2 = arrays.array_scale(d2xdp2, 1.0/(h*h))

        return dxdp, d2xdp2

    def compute_diff(self, GDB: assignments.graph_db, D: List[Dict[int, assignments.graph_db_table]], verbose=False):
        """
        """

        out = []
        X = []
        for gde, tbls  in zip(GDB.entries.values(), D):

            gid = list(gde.tables[assignments.POSITIONS].graphs)[0]
            g = GDB.graphs[gid]
            smiles = GDB.smiles[gid]

            hess_qm = gde.tables[assignments.HESSIANS].values
            hess_qm_sel = gde.tables[assignments.POSITIONS].graphs[gid].rows[0].columns[0].selections
            # shim until I get off the fence on the dimensions of sel
            hess_qm_sel = {k: [v] for k, v in hess_qm_sel.items()}

            hess_qm_pos = assignments.graph_assignment(smiles, hess_qm_sel, g)

            xyz = list([x[0] for x in hess_qm_sel.values()])
            # xyz = None
            sym = graphs.graph_symbols(g)
            mass = [[vibration.mass_table[sym[n]]]*3 for n in sym]
            sym = list(sym.values())

            if True or self.coordinate_system == "IC":
                if "IC:qm_freq" not in self.cache:

                    torsions = False
                    pairs = False
                    ics, B = hessians.bmatrix(hess_qm_pos, torsions=torsions, pairs=pairs, remove1_3=False) # pairs=False, torsions=False, outofplanes=False)
                    ics, B2 = hessians.b2matrix(hess_qm_pos, torsions=torsions, pairs=pairs,remove1_3=False) # pairs=False, torsions=False, outofplanes=False)

                    hess_qm_freq, hess_qm_modes, DL = vibration.hessian_modes(hess_qm, sym, xyz, mass, 0, remove=0, stdtr=True, return_DL=True)

                    # hess_qm_ic = project_ics(B, hess_qm)
                    # ic_qm_fcs = dict(zip(ics, np.diag(hess_qm_ic)))
                    # x0 = array_numpy.dlcmatrix_project_gradients(pos, x0)
                    self.cache['IC:qm_freq'] = hess_qm_freq
                    self.cache['IC:b'] = (ics, B, B2)
                    self.cache['IC:dl'] = DL

                x0 = list(self.cache['IC:qm_freq'])

                hess_mm = tbls[assignments.HESSIANS].values
                grad_mm = list(tbls[assignments.GRADIENTS].values)
                DL = self.cache['IC:dl']
                (ics, B, B2) = self.cache['IC:b']

                x1 =  list(hessians.hessian_frequencies(g, hess_mm, grad_mm, DL, ics, B, B2))


            x = arrays.array_difference(x1, x0)
            N = len(x)

            nlines = len(out)
            mae = [abs(xi) for xi in x]
            out.append("\nComputed Hessian Frequencies (cm-1):")
            for i, (f1, f0) in enumerate(zip(x1, x0), 1):
                out.append(f"     {i:4d} MM: {f1: 8.1f} QM: {f0:8.1f} Diff: {f1-f0:8.1f}")
            out.append("          ----------------------------------------")
            out.append(f"     SAE  MM: {sum(x1): 8.1f} QM: {sum(x0):8.1f} Diff: {sum(mae):8.1f}")
            out.append(f"     MAE  MM: {sum(x1)/N: 8.1f} QM: {sum(x0)/N:8.1f} Diff: {sum(mae)/N:8.1f}")
            out.append(f"Int. SAE  MM: {sum(x1[6:]): 8.1f} QM: {sum(x0):8.1f} Diff: {sum(mae):8.1f}")
            out.append(f"Int. MAE  MM: {sum(x1[6:])/(N-6): 8.1f} QM: {sum(x0[6:])/(N-6):8.1f} Diff: {sum(mae[6:])/(N-6):8.1f}")
            if verbose:
                for i in range(len(out) - nlines):
                    print(out[i])
            X.extend(x)

        return returns.success(arrays.array_round(X, 12), out=out, err=[])

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
        cc.GDB = assignments.graph_db_get(GDB, self.addr)
        cc.csys = csys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
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
        cc.GDB = assignments.graph_db_get(GDB, self.addr)
        cc.csys = csys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.psys = psys
        cc.reuse = reuse
        cc.step_limit = self.step_limit
        cc.tol = self.tol
        return cc

class objective_config_penalty(objective_config):

    def __init__(self, keys=None, scale=1, polynomial={2: 1.0}):
        """

        """
        if keys is None:
            keys = {}

        self.keys = keys
        self.reference = {}
        self.scale = scale
        assert all([type(n) is int for n in polynomial])
        self.polynomial = polynomial
        self.include = True
        
    def get_task(self) -> compute_config:

        # keys is the reference
        ref = self.reference

        # based on self.keys, set the target values from ref
        targets = {}
        for rk, rv in ref.items():
            for pk, pv in reversed(self.keys.items()):
                match = all([x is None or x == y for x, y in zip(pk, rk)])
                if match:
                    if pv is None:
                        targets[rk] = rv
                    else:
                        targets[rk] = pv
                    break

        cc = compute_config_penalty(targets)
        return cc

    def compute_gradient(self, result):
        
        """

        """
        DX = 0
        for dx in result.values():
            for n, m in self.polynomial.items():
                if n == 1:
                    if dx < 0:
                        DX -= m
                    else:
                        DX += m
                elif n % 2:
                    DX += m*n*dx**(n-1)
                else:
                    DX += m*n*abs(dx**(n-1))

        return round(DX, 12)

    def compute_hessian(self, result):
        
        """

        """
        DX = 0
        for dx in result.values():
            for n, m in self.polynomial.items():
                if n == 1:
                    pass
                elif n == 2:
                    DX += 2*m
                elif n % 2:
                    DX += m*n*(n-1)*dx**(n-2)
                else:
                    DX += m*n*(n-1)*abs(dx**(n-2))

        return round(DX, 12)

    def compute_diff(self, result, verbose=False):
        """
        """
        DX = 0
        for dx in result.values():
            for n, m in self.polynomial.items():
                if n % 2:
                    DX += m*abs(dx)**n
                else:
                    DX += m*dx**n

        return returns.success(round(DX, 12), out=[], err=[])

    def compute_total(self, result):
        """
        """
        obj = 0
        return obj
        


class objective_tier:
    def __init__(self):
        self.objectives: Dict[int, objective_config] = {}
        h = 1e-6
        self.h = h
        self.bounds = {
            "s": (0,None), # scale
            "c": (0,None), # cutoff
            "r": (0,None), # sigma
            "e": (h,None), # epsilon
            "k": (0,None), # 
            ("k", "b"): (0,None), # scale
            ("l", "b"): (0,None), # scale
            ("k", "a"): (0,None), # scale
            ("l", "a"): (0, 3.1415), # scale
            ("k", "t"): (-10,10), # scale
            ("k", "i"): (-10,10), # scale
            ("s", "q"): (0,None), # scale
            ("s", "n"): (0,None), # scale
            None: (None, None)
        }
        self.priors = {
            "s": (None, None), # scale
            "c": (None, None), # cutoff
            "r": (None, .1), # sigma
            "e": (None, 1), # epsilon
            "k": (None, 1e9), # 
            ("k", "b"): (None, 1e9), # scale
            ("l", "b"): (None, .1), # scale
            ("k", "a"): (None, 1e9), # scale
            ("l", "a"): (None, .2), # scale
            ("k", "t"): (None, 1.0), # scale
            ("k", "i"): (None, 1.0), # scale
            ("s", "q"): (None, 1.0), # scale
            ("s", "n"): (None, 1.0), # scale
            None: (None, None)
        }
        self.fit_models = None
        self.fit_symbols = None
        self.fit_names = None
        self.fit_names_exclude = None
        self.step_limit = None
        self.maxls = 20
        self.penalties: List[objective_config_penalty] = []
        self.accept: int = 0 # 0 keeps all (essentially disables)
        self.anneal = False

    def key_filter(self, x):
        if self.fit_models is None or self.fit_symbols is None:
            return True

        r = (x[0] in self.fit_models) and (x[1] in self.fit_symbols)
        if r and self.fit_names:
            r &= x[2] in self.fit_names
        if r and self.fit_names_exclude:
            r &= x[2] not in self.fit_names_exclude

        return r


def objective_tier_get_keys(ot, csys):
    keys = *filter(ot.key_filter, mm.chemical_system_iter_keys(csys)),
    return keys
    
def gdb_to_physical_systems(gdb, csys, ref=None, reuse=None):
    psysref = {}
    failures = 0
    if len(gdb.entries) > 100:
        print(f"{datetime.datetime.now()} Starting parameterization")
    for i, (eid, gde) in enumerate(gdb.entries.items(), 1):
        if len(gdb.entries) > 100:
            print(f"\r{datetime.datetime.now()} Parameterizing.. {i:6d}/{len(gdb.entries)}", end="")
        tid = assignments.POSITIONS
        gid = list(gde.tables[tid].graphs)[0]
        rid = 0
        pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
        try:
            psysref[eid] = mm.chemical_system_to_physical_system(csys, [pos], ref=ref, reuse=reuse)
        except Exception as e:
            print(f"\n\nWarning: could not parameterize {eid}, skipping.")
            print(f"SMILES: {pos.smiles}")
            print("Exception:")
            print(e)
            failures += 1
            raise e
    if len(gdb.entries) > 100:
        print()
    if failures:
        print(f"There were {failures} failures that will be skipped")
    return psysref

def objective_tier_run(
    ot, gdb: assignments.graph_db, csys: mm.chemical_system, keys,
    oid=None, psysref=None, reuse=None, wq=None, verbose=False
):

    # need this to avoid race conditions
    csys = copy.deepcopy(csys)
    # build the dataset and input ff
    ws = None
    if wq and configs.remote_compute_enable:
        ws = compute.workqueue_new_workspace(wq, shm={"csys": csys})
    elif configs.processors > 1:
        ws = compute.workspace_local('127.0.0.1', 0, shm={"csys": csys})

    

    # build the initial psys (for charges)
    # if verbose:
    #     print(datetime.datetime.now(), "Building physical systems")
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
    # keys = [k for k in keys if ot.key_filter(k)]
    # assigned = [(i, k, v) for psys in psysref.values() for i,m in enumerate(psys.models) for a in m.labels[0].values() for k, v in a.items()]
    # keys = [k for k in keys if tuple(k[:3]) in assigned or k[1] in "s"]
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
    
    # we need history to identify work for each iteration
    # otherwise we might accidentally use work from step n-1 that gets sent in
    # when we are already at step n
    history = []

    kv = dict(zip(keys, x0))
    
    for penalty in ot.penalties:
        penalty.reference.clear()
        penalty.reference.update(kv)


    # y0 = optimizers_scipy.fit_gdb(x0, *args) 

    bounds = []
    for k in keys:
        b = ot.bounds.get((k[1], k[2].lower()), False)

        if b is False:
            b = ot.bounds.get((k[1], k[2][0].lower()), False)

        if b is False:
            b = ot.bounds.get(k[1], False)

        if b is False:
            b = (None, None)

        bounds.append(b)

    priors = []
    for i, k in enumerate(keys):
        b = ot.priors.get((k[1], k[2].lower()), False)

        if b is False:
            b = ot.priors.get((k[1], k[2][0].lower()), False)

        if b is False:
            b = ot.priors.get(k[1], False)

        if b is False:
            b = (None, None)

        if b[0] is None:
            b = (x0[i], b[1])
        if b[1] is None:
            b = (b[0], 1.0)

        priors.append(b)
        if bounds[i][0] is not None:
            bounds[i] = ((bounds[i][0] - b[0])/b[1], bounds[i][1])
        if bounds[i][1] is not None:
            bounds[i] = (bounds[i][0], (bounds[i][1] - b[0])/b[1])

    x0 = [(x-p[0])/p[1] for x, p in zip(x0, priors)]

    args = (keys, csys, gdb, objectives, priors, ot.penalties, history, psysref, reuse, ws, verbose)

    ret = optimizers_scipy.optimize_forcefield_gdb_scipy(x0, args, bounds=bounds, step_limit=ot.step_limit, maxls=ot.maxls, anneal=ot.anneal)

    result, y0, y1, gx = ret.value

    ret.out.append(f">>> Initial Objective {y0:10.5g}")
    ret.out.append(f">>> Final Objective   {y1:10.5g}")
    ret.out.append(f">>> Percent change    {(y1-y0)/y0*100:10.5g}%")
    if verbose:
        for i in [-3, -2, -1]:
            print(ret.out[i])

    if ws:
        ws.close()
        ws = None


    kv = {k: v*p[1] + p[0] for k,v,p in zip(keys, result, priors)}

    return returns.success((kv, y0, y1, gx), out=ret.out, err=[])

def objective_run_distributed(obj, shm=None):
    if shm and "csys" in shm.__dict__:
        obj.csys = shm.csys
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
    if models is None:
        models = {}
    elif type(models) is not dict:
        models = {i: None for i in models}

    bounds = {}
    strat = forcefield_optimization_strategy(bounds)
    for m, cm in enumerate(csys.models):
        if m not in models:
            continue


        for h, hidx in enumerate(mm.chemical_model_iter_smarts_hierarchies(cm)):

            nodes = [n.name for n in hidx.index.nodes.values() if n.type == "parameter"]
            strat.reference_list.extend(nodes)

            assert type(models[m]) is not str
            if models[m] is not None:
                strat.target_list.extend(models[m])
            else:
                strat.target_list.extend(nodes)


            # splitter = configs.smarts_splitter_config(
            #     1, 2, 0, 0, 0, 0, False, True, 0, True, True, False, False
            # )
            splitter = configs.smarts_splitter_config(
                1, 2, 0, 0, 0, 0, True, True, 0, True, True, True, True
            )
            extender = configs.smarts_extender_config(
                0, 0, True
            )
            bounds[m] = configs.smarts_perception_config(splitter, extender)
            break

    strat.bounds.update(bounds)
    return strat


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

        # when modifying dihedrals, set the frequency limit
        self.modify_outofplane_frequency_limit = 12
        self.modify_torsion_frequency_limit = 12

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

        # Set the bond and angle lengths to whatever the inputs are based on
        # FF assignment
        self.enable_reset_bond_lengths = False
        self.enable_reset_angle_lengths = False

        # Set the bond and angle lengths to whatever the inputs are based on
        # FF assignment. Requires hessians.
        self.enable_reset_bond_stiffness = False
        self.enable_reset_angle_stiffness = False
        self.enable_reset_torsion_stiffness = False
        self.enable_reset_outofplane_stiffness = False
        self.hessian_projection_method = "native"
        
        # After splitting and modifying periodicities, reset to a more 
        self.enable_dihedral_periodicity_reset = True
        self.dihedral_periodicity_reset_max_n = 6
        self.dihedral_periodicity_reset_alpha = -.5
        self.dihedral_periodicity_reset_min_k = 1e-3
        self.dihedral_periodicity_reset_max_k = 5.0

        self.reference_list = []

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
        if self.steps is None:
            self.build_steps()
        return optimization.optimization_strategy_iteration_next(self, clusters)

    def build_steps(self):
        if self.steps is None:
            self.steps = []
        self.steps.extend(
            forcefield_optimization_strategy_build_macro_iterations(self)
        )

def forcefield_optimization_strategy_build_macro_iterations(strat: forcefield_optimization_strategy):

    macro_iters = []
    boundlist = [x.splitter for x in strat.bounds.values()]
    if not boundlist:
        return macro_iters
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
                bits = bounds.bit_search_limit

                search_cursor += 1
                if search_cursor < strat.cursor:
                    continue

                if strat.enable_merge:
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
                        if m not in [2,3]:
                            continue

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
                        s.models.append(m)
                        s.cluster = None
                        s.overlap = [0]
                        s.direct_enable = False
                        s.direct_limit = False
                        s.operation = strat.MODIFY
                        s.iterative_enable = False
                        s.modify_torsion_frequency_limit = strat.modify_torsion_frequency_limit
                        s.modify_outofplane_frequency_limit = strat.modify_outofplane_frequency_limit

                        splitter = configs.smarts_splitter_config(
                            mbounds.splitter.bit_search_min,
                            bits,
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
                            mbounds.splitter.bit_search_min,
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


    strat.cursor = 0

    return macro_iters

def fit(csys, gdb, objective, psystems, nodes, wq=None, verbose=False):

    assigned_nodes = sorted(set([
        (m, l) for psys in psystems.values()
            for m, pm in enumerate(psys.models)
                for proc in pm.labels
                    for glbl in proc.values()
                        for t, l in glbl.items()
    ]))

    fitkeys = [
        x for x in objective_tier_get_keys(objective, csys)
            if (x[0], x[2]) in assigned_nodes
    ]

    fitting_models = set((k[0] for k in fitkeys))
    reuse0=[k for k,_ in enumerate(csys.models) if k not in fitting_models]

    ret = objective_tier_run(
        objective,
        gdb,
        csys,
        fitkeys,
        psysref=psystems,
        reuse=reuse0,
        wq=wq,
        verbose=True
    )
    return ret


eid_t = int
def ff_optimize(
    csys0: mm.chemical_system,
    gdb: assignments.graph_db,
    psystems: Dict[eid_t, mm.physical_system],
    strategy: forcefield_optimization_strategy,
    chemical_objective,
    initial_objective: objective_tier,
    tiers: List[objective_tier],
    final_objective: objective_tier,
) -> mm.chemical_system:

    started = datetime.datetime.now()
    max_line = 0

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

    gcd = csys0.perception.gcd
    labeler = csys0.perception.labeler
    icd = codecs.intvec_codec(gcd.primitive_codecs, gcd.atom_primitives, gcd.bond_primitives)

    reset_config = {
        "bond_l": strategy.enable_reset_bond_lengths,
        "bond_k": strategy.enable_reset_bond_stiffness,
        "angle_l": strategy.enable_reset_angle_lengths,
        "angle_k": strategy.enable_reset_angle_stiffness,
        "torsion_k": strategy.enable_reset_torsion_stiffness,
        "outofplane_k": strategy.enable_reset_outofplane_stiffness,
        "dihedral_p": strategy.enable_dihedral_periodicity_reset,
        "dihedral_max_n": strategy.dihedral_periodicity_reset_max_n,
        "dihedral_alpha": strategy.dihedral_periodicity_reset_alpha,
        "dihedral_min_k": strategy.dihedral_periodicity_reset_min_k,
        "dihedral_max_k": strategy.dihedral_periodicity_reset_max_k
    }
    wq = compute.workqueue_local('0.0.0.0', configs.workqueue_port)

    # print(datetime.datetime.now(), "Resetting bonds and angles to reference values")
    ws = compute.workqueue_new_workspace(wq)
    psystems = reset(
        reset_config,
        csys0,
        gdb,
        psystems,
        verbose=True,
        ws=ws
    )
    ws.close()

    csys = copy.deepcopy(csys0)

    G0 = {i: icd.graph_encode(g) for i, g in gdb.graphs.items()}

    N = 1
    n_ics = 1

    repeat = set()
    visited = set()
    iteration = 0
    N_str = "{:" + str(len(str(n_ics))) + "d}"
    success = False

    keys = mm.chemical_system_iter_keys(csys)
    keys = [k for k in keys if initial_objective.key_filter(k)]

    assigned_nodes = sorted(set([
        (m, l) for psys in psystems.values()
            for m, pm in enumerate(psys.models)
                for proc in pm.labels
                    for glbl in proc.values()
                        for t, l in glbl.items()
    ]))

    print("Initial parameter assignments of dataset:")
    print_chemical_system(csys, show_parameters=[x[1] for x in assigned_nodes])

    assigned_nodes = [x for x in assigned_nodes if x[0] in strategy.bounds and x[1] not in "sc"] 
    assigned_nodes = set([x[1] for x in assigned_nodes])

    print("### BESMARTS chemical perception on the following assignments ###")
    print_chemical_system(csys, show_parameters=assigned_nodes)
    print("#################################################################")

    reuse0 = set(range(len(csys.models))).difference((k[0] for k in keys))
    print(f"{datetime.datetime.now()} Will be caching models {reuse0}")

    # if strategy.enable_dihedral_periodicity_reset:
    #     print("Resetting torsion periodicities")
    #     psystems = reset_project_torsions(
    #         csys,
    #         gdb,
    #         psystems,
    #         max_n=strategy.dihedral_periodicity_reset_max_n,
    #         alpha=strategy.dihedral_periodicity_reset_alpha,
    #         m=2,
    #         verbose=True
    #     )

    kv0 = mm.chemical_system_iter_keys(csys)
    print(f"{datetime.datetime.now()} Computing physical objective")

    fitwq = wq
    if configs.processors == 1:
        fitwq = None

    ret = fit(csys, gdb, initial_objective, psystems, assigned_nodes, wq=fitwq, verbose=True)
    kv, P00, P0, gp0 = ret.value

    print(f"{datetime.datetime.now()} Computing chemical objective")

    CP0 = P0
    CX0 = mm.chemical_system_smarts_complexity(csys)

    C00 = chemical_objective(csys, P0=len(G0), C0=CX0)
    C0 = C00
    C = C00
    print(f"{datetime.datetime.now()} C0={C0}")
    X0 = P0 + C0
    P = P0
    print(datetime.datetime.now(), f"Initial objective: X={P0+C0:13.6g} P={P0:13.6g} C={C0:13.6g}")
    for k, v in kv.items():
        v0 = kv0[k]
        mm.chemical_system_set_value(csys, k, v)
        print(f"{str(k):20s} | New: {v:12.6g} Ref {v0:12.6g} Diff {v-v0:12.6g}")

    for psys in psystems.values():
        reapply = set()
        for k, v in kv.items():
            # print(f"Setting pval to {k}={v}")
            # mm.chemical_system_set_value(csys, k, v)
            mm.physical_system_set_value(psys, k, v)
            reapply.add(k[0])
        for m in reapply:
            procs = csys.models[m].procedures
            if len(procs) > 1:

                for _ in range(1, len(psys.models[m].values)):
                    psys.models[m].values.pop()
                    psys.models[m].labels.pop()
                for proc in procs[1:]:
                    psys.models[m] = proc.assign(
                        csys.models[m],
                        psys.models[m],
                        overrides={
                            k[1:]: v for k, v in kv.items() if k[0] == m
                        }
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

        nodes = [
            x
            for x in strategy.tree_iterator(csys)
            if x.type == "parameter"
            and x.name in assigned_nodes
            # if strategy.cursor == -1
            # or strategy.cursor >= step_tracker.get((x.category, x.name), 0)
            # and x.type == "parameter"
        ]
        # remove any that are not in the models

        print(f"Targets for this macro step {strategy.cursor+1}:")
        for nidx, n in enumerate(nodes, 1):
            print(nidx, n.category, n.name)
        print(f"N Targets: {len(nodes)}")

        if len(nodes) == 0:
            print("Warning, no targets returned. Skipping parameter search.")
            break

        print(f"Step tracker for current macro step {strategy.cursor+1}")
        for n, v in step_tracker.items():
            print(n, v + 1)
        print()

        macro: optimization.optimization_iteration = strategy.macro_iteration(
            nodes
        )

        candidates = {}
        pq = []
        n_added = 0
        n_macro = len(strategy.steps)

        mdls = set([m for mstep in strategy.steps[strategy.cursor-1].steps for m in mstep.models])
        mds = " ".join([f"{i}:{csys.models[i].name}" for i in mdls])
        t = datetime.datetime.now()
        print("*******************")
        print(
            f" iteration={iteration:4d}"
            f" macro={strategy.cursor:3d}/{n_macro}"
            f" micro={len(macro.steps)}"
            f" X={X0:9.5g} P={P0:9.5g} C={C0:9.5g}"
            f" models={mds}"
        )
        print("*******************")
        print()

        # compute the current psystems

        print(f"{datetime.datetime.now()} Initial parameterization using reuse={reuse0}")

        reuse = reuse0
        psys = {}
        for eid, gde in gdb.entries.items():
            tid = assignments.POSITIONS
            gid = list(gde.tables[tid].graphs)[0]
            rid = 0
            pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
            psys[eid] = mm.chemical_system_to_physical_system(
                csys,
                [pos],
                ref=psystems[eid],
                reuse=reuse
            )
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
            cfg.depth_max = d
            cfg.depth_min = S0_depth

            t = datetime.datetime.now()

            print(
                f"{t} Collecting SMARTS for {S.name} and setting to depth={S0_depth}"
            )

            selected_ics = []
            selected_graphs = set()
            for eid, ps in psystems.items():
                gids = [*gdb.entries[eid].tables[assignments.POSITIONS].graphs]
                glbls = ps.models[int(S.category[0])].labels
                for gid, lbls in zip(gids, glbls):
                    for ic, term_lbls in lbls.items():
                        if S.name in term_lbls.values():
                            selected_ics.append((gid, ic))
                            selected_graphs.add(gid)
            # aa = cst.mappings[S.name]
            # selected_graphs = set((x[0] for x in aa))
            aa = selected_ics
            G = {k: v for k, v in G0.items() if k in selected_graphs}
            del selected_graphs

            assn_s = aa
            # assn_s = {i: assn[i] for i in cst.mappings[S.name]}

            iteration += 1

            t = datetime.datetime.now()
            print(
                f"\n =="
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
                f"\n"
            )

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

                # note that if the parameter list is flat, deleting won't be
                # possible
                for p_j, jidx in enumerate(hidx.index.below[S.index]):
                    J = hidx.index.nodes[jidx]
                    if J.type != "parameter":
                        continue

                    for overlap in step.overlap:
                        key = (overlap, macro.cursor, p_j)
                        cnd = (S, J, step, None, None, None, None)
                        candidates[key] = cnd

            elif step.operation == strategy.MODIFY:

                if S.type == "parameter" and S.category[0] in [2,3]:
                    s_per = csys.models[S.category[0]].topology_terms['n'].values[S.name] 
                    f_max = 6
                    if S.category[0] == 2:
                        f_max = step.modify_torsion_frequency_limit
                    elif S.category[0] == 3:
                        f_max = step.modify_outofplane_frequency_limit
                    print(f"Considering frequencies up to {f_max}")
                    mod_max = step.pcp.splitter.bit_search_limit
                    mod_min = step.pcp.splitter.bit_search_min


                    p_j = 0
                    s_add = tuple(sorted(set(range(1, f_max+1)).difference(s_per)))
                    s_per = tuple(s_per)

                    # key = (tuple(s_add), macro.cursor, p_j)
                    # cnd = (S, S, step, None, None, None, None)
                    # candidates[key] = cnd
                    # p_j = 1
                    # print(f"Adding modifactions min={mod_min} max={mod_max}")
                    # print(f"Adding modification {p_j} {S.name}: {s_add}")

                    for to_remove in range(0, min(mod_max, len(s_per))+1):
                        for rem_combo in map(list, itertools.combinations(s_per, to_remove)):
                            for to_add in range(0, min(mod_max, len(s_add))+1):
                                for add_combo in map(list, itertools.combinations(s_add, to_add)):
                                    modify = tuple(sorted([-x for x in rem_combo] + add_combo, key=lambda x: abs(x)))
                                    # print(f"candidate modification {p_j} {S.name}: {modify} {add_combo} {rem_combo}")
                                    if not set(add_combo).symmetric_difference(rem_combo):
                                        continue
                                    if not (set(s_per).symmetric_difference(rem_combo) or add_combo):
                                        continue
                                    if not modify:
                                        continue
                                    if len(modify) < mod_min or len(modify) > mod_max:
                                        continue
                                    key = (modify, macro.cursor, p_j)
                                    cnd = (S, S, step, None, None, None, None)
                                    candidates[key] = cnd
                                    p_j += 1
                                    print(f"Adding modification {p_j} {S.name}: {modify}")

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
                    graphs.subgraph_as_structure(
                        mm.chemical_system_get_node_hierarchy(csys, x[1]).subgraphs[x[1].index],
                        mm.chemical_system_get_node_hierarchy(csys, x[1]).topology
                    )
                    for x in candidates.values()
                ]
                # Sj_lst = [
                #     graphs.subgraph_as_structure(cst.hierarchy.subgraphs[x[1].index], topo)
                #     for x in candidates.values()
                # ]
            elif step.operation == strategy.MODIFY:
                Sj_lst = [
                    graphs.subgraph_as_structure(
                        mm.chemical_system_get_node_hierarchy(csys, x[0]).subgraphs[x[0].index],
                        mm.chemical_system_get_node_hierarchy(csys, x[0]).topology
                    )
                    for x in candidates.values()
                ]
                # Sj_lst = [
                #     graphs.subgraph_as_structure(cst.hierarchy.subgraphs[x[1].index], topo)
                #     for x in candidates.values()
                # ]
            Sj_sma = pool.map_async(gcd.smarts_encode, Sj_lst).get()
            del Sj_lst

        cnd_n = len(candidates)
        t = datetime.datetime.now()

        print_chemical_system(csys, show_parameters=assigned_nodes)

        visited.clear()
        repeat.clear()

        procs = (
            os.cpu_count() if configs.processors is None else configs.processors
        )
        
        print(f"Scoring and filtering {len(candidates)} candidates for operation={step.operation}")
        for t, tier in enumerate(tiers):
            print(f"Tier {t}: Scoring and filtering {len(candidates)} candidates for operation={step.operation}")
            if tier.accept == 0:
                print(f"Tier {t}: Accepting all so we skip")
                continue
            elif len(candidates) <= tier.accept:
                print(f"Tier {t}: Accepting all candidates so we skip")
                continue
            cnd_keys = {i: k for i, k in enumerate(candidates, 1)}

            fitkeys = objective_tier_get_keys(tier, csys)
            fitkeys = [k for k in fitkeys if k[1] in "skeler" and tier.key_filter(k)]

            fitting_models = set([x[0] for x in fitkeys])
            fitting_models.update(strategy.bounds)

            reuse=[k for k,_ in enumerate(csys.models) if k not in fitting_models]

            tier_psystems = psystems
            reuse = [x for x in range(len(csys.models))]
            reset_config_search = reset_config
            # reset_config_search = {
            #     "bond_l": False,
            #     "bond_k": False,
            #     "angle_l": False,
            #     "angle_k": False
            # }
            shm = compute.shm_local(0, data={
                "objective": tier,
                "csys": csys,
                "gdb": gdb,
                "reuse": reuse,
                "psysref": tier_psystems,
                "reset_config": reset_config_search,
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

            if configs.processors == 1 and not configs.remote_compute_enable:
                work = {}
                for k, v in iterable.items():
                    r = calc_tier_distributed(*v[0], **v[1], verbose=True, shm=shm)
                    work[k] = r
            elif len(tier.objectives) < 500 and len(candidates) > 1 and configs.remote_compute_enable:
                print(logs.timestamp(), f"Each worker will compute a full candidate N={len(iterable)}")
                ws = compute.workqueue_new_workspace(wq, address=addr, nproc=procs, shm=shm)
                # # this modifies the csys, relabels and computes objective
                chunksize = 1
                # need to loop through candidates for large fits rather than one worker per candidate
                work = compute.workspace_submit_and_flush(
                    ws,
                    calc_tier_distributed,
                    iterable,
                    chunksize,
                    math.log(len(iterable)//10, 10),
                    len(iterable),
                    verbose=True
                )
                ws.close()
                ws = None
            else:
                # the means we have candidates with lots of things to compute,
                # so do each one at a time
                print(logs.timestamp(), f"Dispatching candidate tasks= {len(iterable)} in serial")
                work = {}
                for i, unit in iterable.items():
                    print(logs.timestamp(), f"Running candidate task {i[0]}/{len(iterable)}")
                    args = unit[0]
                    kwds = unit[1]
                    kwds["wq"] = wq
                    kwds["verbose"] = True
                    work[i] = calc_tier_distributed(*args, **kwds, shm=shm)
            # now just sum over the jobs
            # return keep, X, obj, match_len
            work_new = {}
            for i, _ in enumerate(candidates, 1):
                if i not in work_new:
                    work_new[i] = [0, 0, 0, 0, {}, []]
                for ij, j in work:
                    if i == ij:
                        ret = work[(i,j)]
                        line = ret.value
                        work_new[i][0] |= int(line[0])
                        work_new[i][1] += line[1]
                        work_new[i][2]  = line[2]
                        work_new[i][3] += line[3]
                        # probably use average or something else
                        work_new[i][4].update(line[4])
                        work_new[i][5].extend(ret.out)
                    
            work_full = work
            work = work_new

            tier_scores = []
            max_line = 0
            for j, cnd_i in enumerate(sorted(work), 1):
                (keep, cP, c, match_len, kv, out) = work[cnd_i]
                # cnd_i, key, unit = unit
                (S, Sj, step, _, _, _, _) = candidates[cnd_keys[cnd_i]]

                cC = chemical_objective(csys, P0=len(G0), C0=CX0, c=c)
                cX = cP + cC
                if keep:
                    heapq.heappush(tier_scores, (cX, cnd_i))
            accept = tier.accept
            if accept:
                if type(accept) is float and accept > 0 and accept < 1:
                    accept = max(1, int(len(tier_scores)*accept))
                    print(f"Fraction acceptance {tier.accept*100}% N={accept}/{len(tier_scores)}")
                accepted_keys = [
                    x[1] for x in heapq.nsmallest(accept, tier_scores)
                ]
            else:
                accepted_keys = [
                    x[1] for x in heapq.nsmallest(len(tier_scores), tier_scores)
                ]

            cout_line = (
                f" Initial objectives: "
                f" X= {C0+P0:10.5f}"
                f" P= {P0:10.5f}"
                f" C= {C0:10.5f}"
            )
            print(cout_line)
            print(f"Accepted {len(accepted_keys)} candidates from tier summary")
            for j, cnd_i in enumerate(accepted_keys + list(set(work).difference(accepted_keys)), 1):
                (keep, cP, c, match_len, kv, out) = work[cnd_i]
                # cnd_i, key, unit = unit
                (S, Sj, step, _, _, _, _) = candidates[cnd_keys[cnd_i]]
                edits = cnd_keys[cnd_i][0]
                if edits:
                    edits = str(edits)
                else:
                    edits = f"{Sj_sma[cnd_i-1]}"
                cC = chemical_objective(csys, P0=len(G0), C0=CX0, c=c)
                dP = cP - P0
                dC = cC - C0
                cX = cP + cC
                dX = dP + dC
                K = "Y" if j <= len(accepted_keys) else "N"
                F = ">" if j <= len(accepted_keys) else " "
                cout_line = (
                    f"{F} Cnd. {cnd_i:4d}/{len(work)}"
                    f" N= {match_len:6d}"
                    f" dP= {dP:14.5f}"
                    f" dC= {dC:14.5f}"
                    # f" X0= {X0:10.5f}"
                    f" d(P+C)= {dX:14.5f}"
                    f" d%= {100*(cX-X0)/(X0):10.3f}%"
                    f" {S.name:6s}  " 
                    f" {edits}"
                )
                max_line = max(len(cout_line), max_line)
                # print(datetime.datetime.now())
                print(cout_line)
                if j == len(accepted_keys):
                    print("-"*max_line)
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


            shm = compute.shm_local(0, data={
                "objective": initial_objective,
                "csys": csys,
                "gdb": gdb,
                "reuse": reuse0,
                "psysref": psystems,
                "reset_config": reset_config
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
            print(datetime.datetime.now(), f"Generated {len(candidates)} x {len(initial_objective.objectives)//len(j)} = {len(iterable)} candidate evalulation tasks")

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
                procs = len(iterable)


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

            if step.operation != strategy.MERGE and (macroamt + macroampc + microamt + microampc == 0):
                # use X0 so later dX will be 0 and kept
                # if we do this for merges, every merge will be taken..
                work = {i: (1, X0, 0.0, 1) for i in cnd_keys}

            else:
                if configs.processors == 1 and not configs.remote_compute_enable:
                    work = {}
                    for k, v in iterable.items():
                        r = calc_tier_distributed(*v[0], **v[1], verbose=True, shm=shm)
                        work[k] = r
                elif len(initial_objective.objectives) < 500 and len(candidates) > 1 and configs.remote_compute_enable:
                    # this means each candidate has relatively few targets to compute, so we can let each worker handle one candidate
                    print(logs.timestamp(), f"Each worker will compute a full candidate N={len(iterable)}")
                    ws = compute.workqueue_new_workspace(wq, address=addr, nproc=procs, shm=shm)
                    # # this modifies the csys, relabels and computes objective
                    # so i should use objective_tier_run instead and loop through the iterable
                    chunksize = 1
                    work = compute.workspace_submit_and_flush(
                        ws,
                        calc_tier_distributed,
                        iterable,
                        chunksize,
                        0.0,
                        len(iterable),
                        verbose=True,
                    )
                    ws.close()
                    ws = None
                else:
                    # the means we have candidates with lots of things to compute,
                    # so do each one at a time
                    print(logs.timestamp(), f"Dispatching candidate tasks= {len(iterable)} in serial")
                    work = {}
                    for i, unit in iterable.items():
                        print(logs.timestamp(), f"Running candidate task {i[0]}/{len(iterable)}")
                        args = unit[0]
                        kwds = unit[1]
                        kwds["wq"] = wq
                        kwds["verbose"] = True
                        work[i] = calc_tier_distributed(*args, **kwds, shm=shm)

                # now just sum over the jobs
                # return keep, X, obj, match_len
                work_new = {}
                for i, _ in enumerate(candidates, 1):
                    if i not in work_new:
                        work_new[i] = [0, 0, 0, 0, {}, []]
                    for ij, j in work:
                        if i == ij:
                            ret = work[(i,j)]
                            line = ret.value
                            work_new[i][0] |= int(line[0])
                            work_new[i][1] += line[1]
                            work_new[i][2]  = line[2]
                            work_new[i][3] += line[3]
                            work_new[i][4].update(line[4])
                            work_new[i][5].extend(ret.out)

                work_full = work
                work = work_new

            print(f"The unfiltered results of the candidate scan N={len(work)} total={len(iterable)} oper={step.operation}:")


            best_reuse = None
            if reuse_cnd:
                # just append the best to work and let the loop figure it out
                best_reuse = sorted(reuse_cnd.items(), key=lambda y: (-y[1][0], y[1][1], y[1][2], y[1][3]))[0]
                work[best_reuse[0]] = best_reuse[1]

            cout_line = (
                f" Initial objectives: "
                f" X= {C0+P0:10.5f}"
                f" P= {P0:10.5f}"
                f" C= {C0:10.5f}"
            )
            print(cout_line)
            cnd_kv = {}
            for j, cnd_i in enumerate(sorted(work), 1):
                (keep, cP, c, match_len, kv, out) = work[cnd_i]
                # cnd_i, key, unit = unit
                (S, Sj, step, _, _, _, _) = candidates[cnd_keys[cnd_i]]

                dP = cP - P0
                # cC = C0 + dcC
                cC = chemical_objective(csys, P0=len(G0), C0=CX0, c=c)
                cX = cP + cC
                dP = cP - P0
                dC = cC - C0
                dX = dP + dC
                # dX = cX - X0
                keep = keep and dX <= 0.0

                if step.operation == strategy.SPLIT:
                    visited.add(S.name)
                elif step.operation == strategy.MERGE:
                    visited.add(Sj.name)
                elif step.operation == strategy.MODIFY:
                    visited.add(S.name)


                reused_line = ""
                if best_reuse is not None and cnd_i == best_reuse[0]:
                    reused_line="*"
                K = "Y" if keep else "N"
                edits = cnd_keys[cnd_i][0]
                if edits:
                    edits = str(edits)
                else:
                    edits = f"{Sj_sma[cnd_i-1]}"
                cout_line = (
                    f"Cnd. {cnd_i:4d}/{len(work)}"
                    f" N= {match_len:6d} K= {K}"
                    f" dP= {dP:14.5f}"
                    f" dC= {dC:14.5f}"
                    f" d(P+C)= {dX:14.5f}"
                    f" d%= {100*(cX-X0)/(X0):10.3f}%"
                    f" {S.name:6s} {reused_line}" 
                    f" {edits}"
                )
                max_line = max(len(cout_line), max_line)
                # print(datetime.datetime.now())
                print(cout_line, end=" " * (max_line - len(cout_line)) + '\n')
                sys.stdout.flush()

                if match_len == 0:
                    if step.operation in [strategy.SPLIT, strategy.MODIFY]:
                        keep = False
                        ignore.add(cnd_i)
                        continue



                # We prefer to add in this order
                cout_key = None

                # print sorted at the end but only for new
                # this is to speed things up
                cout_key = (-int(keep), cX, match_len, cnd_i, S.name)
                cout[cout_key] = cout_line

                # use these below to determine the best ones to keep
                heapq.heappush(cout_sorted_keys, cout_key)
                cnd_kv[cout_key] = kv

                if not keep:
                    ignore.add(cnd_i)
                    continue


                if cnd_i in kept:
                    ignore.add(cnd_i)
                    continue

            print("\r" + " " * max_line)


            # print sorted at the end
            print(f"Nanostep {n_nano}: The filtered results of the candidate scan N={len(cout)} total={len(iterable)} oper={step.operation}:")
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
                    case5 = ck[3] not in ignore
                    if case3 and case4 and case5:
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

            success = False
            added = False
            if len(cnd_keep) == 1 and len(nodes) == 1:
                print("Only one modification, keeping result and printing output:")
                (keep, cP, c, match_len, kv, out) = work[cnd_keep[0][3]]
                print(f"\n".join(out))

                C = chemical_objective(csys, P0=len(G0), C0=CX0)
                X = cnd_keep[0][1]
                P = X - C
                print(datetime.datetime.now(), f"Macro objective: P={P:13.6g} C={C:13.6g} DX={P+C-X0:13.6g}")

                kv = cnd_kv[cnd_keep[0]]
                for k, v in kv.items():
                    v0 = mm.chemical_system_get_value(csys, k)
                    mm.chemical_system_set_value(csys, k, v)
                    print(f"{str(k):20s} | New: {v:12.6g} Ref {v0:12.6g} Diff {v-v0:12.6g}")

                psysref = {
                    i: mm.chemical_system_to_physical_system(
                        csys,
                        psystems[i].models[0].positions,
                        ref=psystems[i],
                        reuse=reuse0
                    ) for i in psystems
                }

            else:


                # print("*just force testing this code path")

                psysref = {
                    i: mm.chemical_system_to_physical_system(
                        csys,
                        psystems[i].models[0].positions,
                        ref=psystems[i],
                        reuse=reuse0
                    ) for i in psystems
                }

                ws = None
                if wq:
                    ws = compute.workqueue_new_workspace(wq)
                psysref = reset(reset_config, csys, gdb, psysref, verbose=True, ws=ws)
                if ws:
                    ws.close()
                print("Multiple modifications, doing another fit with all accepted*")

                reuse = [x for x in range(len(csys.models))]

                fitkeys = objective_tier_get_keys(initial_objective, csys)
                fitting_models = set((k[0] for k in fitkeys))
                fitkeys = [k for k in fitkeys if (k[1] in "skeler" and k[2] in assigned_nodes) or k[0] in fitting_models]
                reuse=[k for k,_ in enumerate(csys0.models) if k not in fitting_models]
                
                kv0 = {k: mm.chemical_system_get_value(csys, k) for k in fitkeys}
                # for k, v in kv0.items():
                #     print(f"{str(k):20s} | v0 {v:12.6g}")
                print_chemical_system(csys, show_parameters=assigned_nodes.union([x.name for x in nodes.values()]))
                ret = objective_tier_run(
                    initial_objective,
                    gdb,
                    csys,
                    fitkeys,
                    psysref=psysref,
                    reuse=reuse,
                    wq=wq,
                    verbose=True
                )
                kv, _, P, gp = ret.value
                print("\n".join(ret.out))
                C = chemical_objective(csys, P0=len(G0), C0=CX0)
                X = P + C
                dX = X - X0
                print(datetime.datetime.now(), f"Macro objective: {P:13.6g} C={C:13.6g} DX={P+C-X0:13.6g}")
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

                    success = True
                    added = True
                    repeat.add(hent.name)
                    assigned_nodes.add(hent.name)
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

                    if hent.name in assigned_nodes:
                        assigned_nodes.remove(hent.name)
                    success = True
                    added = True

                elif step.operation == strategy.MODIFY:
                    if (S.category, S.name) in step_tracker:
                        step_tracker.pop((S.category, S.name))
                    repeat.add((S.category, S.name))
                    visited.add((S.category, S.name))
                    success = True
                    added = True

            print(datetime.datetime.now(), "Chemical system after nanostep:")
            # print the tree
            print_chemical_system(csys, show_parameters=assigned_nodes)

            X0 = X
            P0 = P
            C0 = C
            psystems = psysref


        if strategy.macro_accept_max_total > 0 and n_added > 0:
            strategy.repeat_step()

        print(f"There were {n_added} successful operations")

        print(f"{datetime.datetime.now()} Visited", visited)

        print(datetime.datetime.now(), "Saving chk.cs.p")
        pickle.dump([gdb, csys, strategy, psystems], open("chk.cst.p", "wb"))

        print(datetime.datetime.now(), "Macro step done.")

        print()
        print("#"*120)
        print()

    print(f"{datetime.datetime.now()} Strategy done.")

    if final_objective and (initial_objective or tiers):

        if strategy.enable_dihedral_periodicity_reset:
            print("Resetting torsion periodicities")
            psystems = reset_project_torsions(
                csys,
                gdb,
                psystems,
                max_n=strategy.dihedral_periodicity_reset_max_n,
                alpha=strategy.dihedral_periodicity_reset_alpha,
                m=2,
                verbose=True
            )

        print(f"{datetime.datetime.now()} Computing final fit")

        if configs.processors == 1:
            psystems = reset(reset_config, csys, gdb, psystems, verbose=True, ws=None)
            kv, _, P, gp = fit(csys, gdb, final_objective, psystems, assigned_nodes, wq=None, verbose=True).value
        else:
            ws = compute.workqueue_new_workspace(wq)
            psystems = reset(reset_config, csys, gdb, psystems, verbose=True, ws=ws)
            ws.close()
            kv, _, P, gp = fit(csys, gdb, final_objective, psystems, assigned_nodes, wq=wq, verbose=True).value

        C = chemical_objective(csys, P0=len(G0), C0=CX0)
        print(f"{datetime.datetime.now()} C0={C0}")
        X = P + C
        print(datetime.datetime.now(), f"Final objective: X={X:13.6g} P={P:13.6g} C={C:13.6g}")
        for k, v in kv.items():
            v0 = mm.chemical_system_get_value(csys, k)
            mm.chemical_system_set_value(csys, k, v)
            print(f"{str(k):20s} | New: {v:12.6g} Ref {v0:12.6g} Diff {v-v0:12.6g}")

        for psys in psystems.values():
            reapply = set()
            for k, v in kv.items():
                # mm.chemical_system_set_value(csys, k, v)
                mm.physical_system_set_value(psys, k, v)
                reapply.add(k[0])
            for m in reapply:
                procs = csys.models[m].procedures
                if len(procs) > 1:
                    for _ in range(1, len(psys.models[m].values)):
                        psys.models[m].values.pop()
                        psys.models[m].labels.pop()
                    for proc in procs[1:]:
                        psys.models[m] = proc.assign(csys.models[m], psys.models[m], overrides={k[1:]: v for k, v in kv.items() if k[0] == m})

    ended = datetime.datetime.now()

    print(f"Start time: {started}")
    print(f"End   time: {ended}")

    # return csys
    wq.close()
    wq = None

    print_chemical_system(csys)

    print(f"{datetime.datetime.now()} Saving final checkpoint to chk.cst.p")
    pickle.dump([gdb, csys, strategy, psystems], open("chk.cst.p", "wb"))

    print(f"{datetime.datetime.now()} Saving final (csys, (P0, P), (C0, C)) to csys.p")
    pickle.dump((csys, (P00, P), (C00, C)), open("csys.p", "wb"))

    return csys, (P00, P), (C00, C)

def calc_tier_distributed(S, Sj, operation, edits, oid, verbose=False, wq=None, shm=None):

    # copy once
    csys = copy.deepcopy(shm.csys)
    # csys = shm.csys

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
    kv, y0, P, gx, C = {}, 0, 0, [], 0
    # need to perform the operation and then add to keys
    # would also need to add the node to the FF
    if operation == optimization.optimization_strategy.SPLIT:
        node = mm.chemical_model_smarts_hierarchy_copy_node(cm, cid, pid, uid, S, None)
        hidx.subgraphs[node.index] = Sj
        sma = gcd.smarts_encode(Sj)
        hidx.smarts[node.index] = sma
        # dC = mm.graph_complexity(Sj, scale=.01, offset=-len(hidx.topology.primary)*.01)


    elif operation == optimization.optimization_strategy.MERGE:
        node = Sj
        sma = hidx.smarts[node.index]
        g = hidx.subgraphs[node.index]
        # dC = -mm.graph_complexity(Sj, scale=.01, offset=-len(hidx.topology.primary)*.01)
        mm.chemical_model_smarts_hierarchy_remove_node(cm, cid, pid, uid, Sj)

    elif operation == optimization.optimization_strategy.MODIFY:
        node = S
        # dC = -mm.graph_complexity(Sj, scale=.01, offset=-len(hidx.topology.primary)*.01)
        pers = cm.topology_terms['n'].values[S.name]
        for n in edits:
            if n < 0 and -n in pers:
                i = cm.topology_terms['n'].values[S.name].index(-n)
                for t in 'kpn':
                    cm.topology_terms[t].values[S.name].pop(i)
            elif n > 0 and n not in pers:
                cm.topology_terms['n'].values[S.name].append(n)
                cm.topology_terms['p'].values[S.name].append(0.0)
                cm.topology_terms['k'].values[S.name].append(0.0)
        if not cm.topology_terms['n'].values[S.name]:
            print(pers, edits)
            assert False

            
        # mm.chemical_model_smarts_hierarchy_remove_node(cm, cid, pid, uid, S)


    # since we only changed by Sj
    reuse = [x for x in range(len(csys.models)) if x != cid]

    # this does the param refresh after the modification
    # so it will reSMARTS cid
    psysref = {
        i: mm.chemical_system_to_physical_system(
            csys,
            psysref[i].models[0].positions,
            ref=psysref[i],
            reuse=reuse
        ) for i in psysref
    }

    # reset and potentially relabel psys
    ws = None
    if wq and configs.remote_compute_enable:
        ws = compute.workqueue_new_workspace(wq, shm={})
    elif configs.processors > 1 and len(objective.objectives) > 5:
        ws = compute.workspace_local('127.0.0.1', 0, shm={})
    psysref = reset(shm.reset_config, csys, gdb, psysref, verbose=True, ws=ws)
    if ws:
        ws.close()

    # now that we have all reSMARTS, just reuse everything
    # reuse = [x for x in range(len(csys.models))]

    # keys = objective_tier_get_keys(objective, csys)
    # fitting_models = set((k[0] for k in keys))
    # assigned = [(i, k, v) for psys in psysref.values() for i,m in enumerate(psys.models) for a in m.labels[0].values() for k, v in a.items()]
    # keys = [k for k in keys if tuple(k[:3]) in assigned or k[1] in "s"]
    # reuse=[m for m in csys.models if m not in fitting_models]
    # print("Fitting keys are:")
    # print(keys)
    # kv0 = {k: mm.chemical_system_get_value(csys, k) for k in keys}
    # for k, v in kv0.items():
    #     print(f"{str(k):20s} | v0 {v:12.6g}")
    assigned_nodes = sorted(set([
        (m, l) for psys in psysref.values()
            for m, pm in enumerate(psys.models)
                for proc in pm.labels
                    for glbl in proc.values()
                        for t, l in glbl.items()
    ]))

    fitkeys = [
        x for x in objective_tier_get_keys(objective, csys)
            if (x[0], x[2]) in assigned_nodes
    ]

    fitting_models = set((k[0] for k in fitkeys))
    reuse=[k for k,_ in enumerate(csys.models) if k not in fitting_models]

    if operation == optimization.optimization_strategy.SPLIT:

        match_len = 0
        old_match = 0
        for psys in psysref.values():
            pm = psys.models[cid]
            for ic, lbls in pm.labels[0].items():
                lbls = set(lbls.values())
                if node.name in lbls:
                    match_len += 1
                elif S.name in lbls:
                    old_match += 1

        # new_keys = mm.chemical_system_smarts_hierarchy_get_node_keys(cm, cid, pid, uid, node)
        if old_match == 0 or match_len == 0:   
            keep = False
    elif operation in [optimization.optimization_strategy.MERGE, optimization.optimization_strategy.MODIFY]: 

        match_len = 0
        for psys in psysref.values():
            pm = psys.models[cid]
            for ic, lbls in pm.labels[0].items():
                lbls = set(lbls.values())
                if S.name in lbls:
                    match_len += 1

    # print("Matches", match_len, "Old matches", old_match)
    # C = graphs.graph_bits(Sj) / len(Sj.nodes) /1000 + len(Sj.nodes)
    # configs.remote_compute_enable = False
    # wq = compute.workqueue_local('127.0.0.1', 0)
    out = []
    kv = {}
    if keep:
        ret = objective_tier_run(objective, gdb, csys, fitkeys, oid=oid, psysref=psysref, reuse=reuse, wq=wq, verbose=verbose)
        kv, y0, P, gx = ret.value
        out = ret.out

    # print(f"{S.name}->{sma:40s} OID={oid} {keep} {X} {C} {match_len}")

    c = mm.chemical_system_smarts_complexity(csys)

    # if wq:
    #     wq.close()
    return returns.success((keep, P, c, match_len, kv), out=out)

def print_chemical_system(csys, show_parameters=None):
    mm.chemical_system_print(csys, show_parameters=show_parameters)


def chemical_objective(csys, P0=1.0, C0=1.0, A=0.01, B=1.0, C=1.0, c=None):

    if c is None:
        c = mm.chemical_system_smarts_complexity(csys, B=B, C=C)
    CC = A * P0 * math.exp(A * (c - C0))
    return CC

def perform_operations(
        csys: mm.chemical_system,
        candidates,
        keys,
        Sj_sma,
    ):


    nodes = {}
    ignore = set()
    for cnd_i, key in keys.items():
        (S, Sj, step, _, _, _, _) = candidates[key]
        (edits, _, p_j) = key
        # param_name = "p."
        sma = ""
        added = False
        cid, pid, uid = S.category
        cm = csys.models[cid]


        if step.operation == optimization.optimization_strategy.SPLIT:
            hidx = mm.chemical_system_get_node_hierarchy(csys, S)
            if hidx is None:
                breakpoint()
                print(f"Invalid node for operation: {S.name}")
                continue
            topo = hidx.topology
            # param_name = "p" + str(group_number)
            node = mm.chemical_model_smarts_hierarchy_copy_node(cm, cid, pid, uid, S, None)
            node.type = "parameter"
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
            # S might have been deleted previously
            # and Sj is a graph in splitting
            if Sj.name in ignore:
                continue
            hidx = mm.chemical_system_get_node_hierarchy(csys, Sj)
            if hidx is None:
                breakpoint()
                print(f"Invalid node for operation: {Sj.name}")
                continue
            topo = hidx.topology
            if Sj.index in hidx.subgraphs:
                hidx.subgraphs.pop(Sj.index)
            if Sj.index in hidx.smarts:
                hidx.smarts.pop(Sj.index)

            mm.chemical_model_smarts_hierarchy_remove_node(cm, cid, pid, uid, Sj)
            ignore.add(Sj.name)
            nodes[cnd_i] = Sj

        elif step.operation == optimization.optimization_strategy.MODIFY:

            if (S.name, 1) in ignore:
                continue

            pers = cm.topology_terms['n'].values[S.name]
            for n in edits:
                if n < 0 and -n in pers:
                    i = cm.topology_terms['n'].values[S.name].index(-n)
                    for t in 'kpn':
                        cm.topology_terms[t].values[S.name].pop(i)
                elif n > 0 and n not in pers:
                    cm.topology_terms['n'].values[S.name].append(n)
                    cm.topology_terms['p'].values[S.name].append(0.0)
                    cm.topology_terms['k'].values[S.name].append(0.0)
            nodes[cnd_i] = S
            ignore.add((S.name, 1))

    return csys, nodes

def print_xyz(pos, comment="") -> List[str]:
    lines = []
    lines.append(str(len(pos.selections)))
    lines.append(comment)
    for ic, xyz in pos.selections.items():
        n = pos.graph.nodes[ic[0]]
        sym = primitives.element_tr[str(n.primitives['element'].on()[0])]
        try:
            x, y, z = xyz[0][:3]
        except TypeError:
            x, y, z = xyz[:3]
        lines.append(f"{sym:8s} {x:.6f} {y:.6f} {z:.6f}")
    return lines

def chemical_system_get_hessian(csys, psystems, m, fn, names=None) -> dict:
    """
    Look for hessians, do the projections
    Need to calculate the MM hessians
    set alpha to sum(QMic)/sum(MMic), set MMic = alpha * MMic 
    """
    kv = {(k[0], 'l', k[2], None): [] for k in mm.chemical_system_iter_keys(csys) if k[0] == m and k[1] == 'l'}
    for psys in psystems:
        pm: mm.physical_model = psys.models[m]
        pos = pm.positions[0]
        measure = fn(pos)
        for ic, ic_terms in pm.labels[0].items():
            lbl = ic_terms['l']
            if names and lbl not in names:
                continue
            x = measure.selections[ic][0]
            kv[(m, 'l', lbl, None)].extend(x)
    return kv

def chemical_system_cluster_data(csys, m, sag, objective, strategy=None):
    if strategy is None:
        splitter = configs.smarts_splitter_config(
            1, 2, 0, 0, 0, 0, True, True, 0, False, True, True, True
        )
        extender = configs.smarts_extender_config(
            0, 0, True
        )
        cfg = configs.smarts_perception_config(splitter, extender)
        optimization = cluster_optimization.optimization_strategy_default(cfg)
        optimization.overlaps=[100]
    else:
        optimization = strategy

    gcd = csys.perception.gcd
    labeler = csys.perception.labeler

    hidx = csys.models[m].procedures[0].smarts_hierarchies[0].copy()

    optimization.target_list = [n.name for n in hidx.index.nodes.values()]
    optimization.reference_list = [n.name for n in hidx.index.nodes.values()]
    optimization.merge_protect_list = [n.name for n in hidx.index.nodes.values()]
    initial_conditions = clusters.clustering_initial_conditions(gcd, sag, hidx=hidx, labeler=labeler, prefix=csys.models[m].symbol)

    if objective.is_discrete():
        cst = cluster_optimization.cluster_classifications(
            gcd, labeler, sag, objective, optimization=optimization, initial_conditions=initial_conditions
        )
    else:
        cst = cluster_optimization.cluster_means(
            gcd, labeler, sag, objective, optimization=optimization, initial_conditions=initial_conditions
        )

    cm = csys.models[m]
    proc = cm.procedures[0]

    newhidx = cst.hierarchy
    cst.hierarchy = None

    uid = 0
    refhidx = proc.smarts_hierarchies[uid]

    roots = trees.tree_index_roots(refhidx.index)
    assert len(roots) == 1
    root = roots[0]
    assert root.type == "hierarchy"

    for idx in [n.index for n in refhidx.index.nodes.values() if n.index != root.index]:
        if idx in refhidx.smarts:
            refhidx.smarts.pop(idx)
        name = refhidx.index.nodes[idx].name
        for (ui, namei) in list(proc.topology_parameters):
            if namei == name:
                proc.topology_parameters.pop((ui, namei))
        for sym, terms in cm.topology_terms.items():
            if name in terms.values:
                terms.values.pop(name)

        trees.tree_index_node_remove(refhidx.index, idx)

    nodes = []
    for newroot in trees.tree_index_roots(newhidx.index):
        if newroot.type == "parameter":
            nodes.append(newroot)
        for n in tree_iterators.tree_iter_breadth_first(newhidx.index, newroot):
            if n.type == "parameter":
                nodes.append(n)

    node_map = {}
    for cstnode in nodes:
        above = None
        if newhidx.index.above[cstnode.index] is None:
            above = root.index
        elif newhidx.index.above[cstnode.index] in node_map:
            above = node_map[newhidx.index.above[cstnode.index]]
        else:
            above = root.index
        # cstnode.name = cm.symbol + cstnode.name
        cstnode.type = "parameter"
        cstnode.category = (m, 0, 0)

        old_index = cstnode.index
        new_node = trees.tree_index_node_add(refhidx.index, above, cstnode)

        node_map[old_index] = cstnode.index

        pkey = (0, new_node.name)

        if m == 0:
            terms = {"k": new_node.name, "l": new_node.name}
            proc.topology_parameters[pkey] = terms
            kval = 10.0
            lval = 1.3
            cm.topology_terms["k"].values[new_node.name] = [kval]
            cm.topology_terms["l"].values[new_node.name] = [lval]
        elif m == 1:
            terms = {"k": new_node.name, "l": new_node.name}
            proc.topology_parameters[pkey] = terms
            kval = 50.0
            lval = 109*3.14/180
            cm.topology_terms["k"].values[new_node.name] = [kval]
            cm.topology_terms["l"].values[new_node.name] = [lval]
        elif m == 2:
            terms = {"n": new_node.name, "p": new_node.name, "k": new_node.name}
            proc.topology_parameters[pkey] = terms
            cm.topology_terms["n"].values[new_node.name] = [1,2,3]
            cm.topology_terms["p"].values[new_node.name] = [0,0,0]
            cm.topology_terms["k"].values[new_node.name] = [0,0,0]

        sma = newhidx.smarts[old_index]
        refhidx.smarts[new_node.index] = sma
        if old_index in newhidx.subgraphs:
            refhidx.subgraphs[new_node.index] = newhidx.subgraphs[old_index]

    refhidx.topology = cm.topology

    return csys

def chemical_system_cluster_force_constants(csys, gdb, sep, topo):
    sag = []
    
    measure = None
    m = 0
    if topo == topology.bond:
        measure = assignments.graph_assignment_geometry_bonds
        m = 0
    elif topo == topology.angle:
        measure = assignments.graph_assignment_geometry_angles
        m = 1
    assert measure

    for gid in gdb.graphs:
        measurements = {}
        g = gdb.graphs[gid]
        atoms = {atom: [] for atom in graphs.graph_atoms(g)}
        smiles = gdb.smiles[gid]
        for eid, gde in gdb.entries.items():
            gdt = gde.tables[assignments.POSITIONS]
            gdg = gdt.graphs.get(gid)

            if gdg is None:
                continue

            for rid, row in gdg.rows.items():
                for cid, col in row.columns.items():
                    for ic, xyz in col.selections.items():
                        atoms[ic].append(xyz)

        pos = assignments.graph_assignment(smiles, atoms, g)
        r = measure(pos)
        for ic, rv in r.selections.items():
            measurements[ic] = [y for x in rv for y in x]
        sag.append(assignments.smiles_assignment_float(smiles, measurements))
    sag = assignments.smiles_assignment_group(sag, topo)
    objective = cluster_objective.clustering_objective_mean_separation(sep, sep)

    gcd = csys.perception.gcd
    labeler = csys.perception.labeler

    splitter = configs.smarts_splitter_config(
        1, 2, 0, 0, 0, 0, True, True, 0, False, True, True, True
    )
    extender = configs.smarts_extender_config(
        0, 0, True
    )
    cfg = configs.smarts_perception_config(splitter, extender)
    optimization = cluster_optimization.optimization_strategy_default(cfg)
    
    cst = cluster_optimization.cluster_means(
        gcd, labeler, sag, objective, optimization=optimization, initial_conditions=None
    )

    cst.hierarchy

    cm = csys.models[m]
    proc = cm.procedures[0]

    newhidx = cst.hierarchy
    cst.hierarchy = None

    uid = 0
    refhidx = proc.smarts_hierarchies[uid]

    roots = trees.tree_index_roots(refhidx.index)
    assert len(roots) == 1
    root = roots[0]
    assert root.type == "hierarchy"
    

    for idx in [n.index for n in refhidx.index.nodes.values() if n.index != root.index]:
        if idx in refhidx.smarts:
            refhidx.smarts.pop(idx)
        name = refhidx.index.nodes[idx].name
        for (ui, namei) in list(proc.topology_parameters):
            if namei == name:
                proc.topology_parameters.pop((ui, namei))
        for sym, terms in cm.topology_terms.items():
            if name in terms.values:
                terms.values.pop(name)

        trees.tree_index_node_remove(refhidx.index, idx)


    nodes = []
    for newroot in trees.tree_index_roots(newhidx.index):
        nodes.append(newroot)
        for n in tree_iterators.tree_iter_breadth_first(newhidx.index, newroot):
            nodes.append(n)

    node_map = {}
    for cstnode in nodes:
        above = None
        if newhidx.index.above[cstnode.index] is None:
            above = root.index
        else:
            above = node_map[newhidx.index.above[cstnode.index]]
        cstnode.name = cm.symbol + cstnode.name
        cstnode.type = "parameter"
        cstnode.category = (m, 0, 0)

        old_index = cstnode.index
        new_node = trees.tree_index_node_add(refhidx.index, above, cstnode)

        node_map[old_index] = cstnode.index

        pkey = (0, new_node.name)
        terms = {"k": new_node.name, "l": new_node.name}
        proc.topology_parameters[pkey] = terms

        if m == 0:
            kval = 10.0
            lval = 1.3
        else:
            kval = 50.0
            lval = 109*3.14/180
        cm.topology_terms["k"].values[new_node.name] = [kval]
        cm.topology_terms["l"].values[new_node.name] = [lval]

        sma = newhidx.smarts[old_index]
        refhidx.smarts[new_node.index] = sma
        refhidx.subgraphs[new_node.index] = newhidx.subgraphs[old_index]

    refhidx.topology = cm.topology

    return csys


def chemical_system_cluster_geom(csys, gdb, sep, topo):
    sag = []
    
    measure = None
    m = 0
    if topo == topology.bond:
        measure = assignments.graph_assignment_geometry_bonds
        m = 0
    elif topo == topology.angle:
        measure = assignments.graph_assignment_geometry_angles
        m = 1
    assert measure

    for gid in gdb.graphs:
        measurements = {}
        g = gdb.graphs[gid]
        atoms = {atom: [] for atom in graphs.graph_atoms(g)}
        smiles = gdb.smiles[gid]
        for eid, gde in gdb.entries.items():
            gdt = gde.tables[assignments.POSITIONS]
            gdg = gdt.graphs.get(gid)

            if gdg is None:
                continue

            for rid, row in gdg.rows.items():
                for cid, col in row.columns.items():
                    for ic, xyz in col.selections.items():
                        atoms[ic].append(xyz)

        pos = assignments.graph_assignment(smiles, atoms, g)
        r = measure(pos)
        for ic, rv in r.selections.items():
            measurements[ic] = [y for x in rv for y in x]
        sag.append(assignments.smiles_assignment_float(smiles, measurements))
    sag = assignments.smiles_assignment_group(sag, topo)
    objective = cluster_objective.clustering_objective_mean_separation(sep, sep)
    

    
    gcd = csys.perception.gcd
    labeler = csys.perception.labeler

    splitter = configs.smarts_splitter_config(
        1, 2, 0, 0, 0, 0, True, True, 0, False, True, True, True
    )
    extender = configs.smarts_extender_config(
        0, 0, True
    )
    cfg = configs.smarts_perception_config(splitter, extender)
    optimization = cluster_optimization.optimization_strategy_default(cfg)
    
    cst = cluster_optimization.cluster_means(
        gcd, labeler, sag, objective, optimization=optimization, initial_conditions=None
    )

    cst.hierarchy

    cm = csys.models[m]
    proc = cm.procedures[0]


    newhidx = cst.hierarchy
    cst.hierarchy = None

    uid = 0
    refhidx = proc.smarts_hierarchies[uid]

    roots = trees.tree_index_roots(refhidx.index)
    assert len(roots) == 1
    root = roots[0]
    assert root.type == "hierarchy"
    

    for idx in [n.index for n in refhidx.index.nodes.values() if n.index != root.index]:
        if idx in refhidx.smarts:
            refhidx.smarts.pop(idx)
        name = refhidx.index.nodes[idx].name
        for (ui, namei) in list(proc.topology_parameters):
            if namei == name:
                proc.topology_parameters.pop((ui, namei))
        for sym, terms in cm.topology_terms.items():
            if name in terms.values:
                terms.values.pop(name)

        trees.tree_index_node_remove(refhidx.index, idx)


    nodes = []
    for newroot in trees.tree_index_roots(newhidx.index):
        nodes.append(newroot)
        for n in tree_iterators.tree_iter_breadth_first(newhidx.index, newroot):
            nodes.append(n)

    node_map = {}
    for cstnode in nodes:
        above = None
        if newhidx.index.above[cstnode.index] is None:
            above = root.index
        else:
            above = node_map[newhidx.index.above[cstnode.index]]
        cstnode.name = cm.symbol + cstnode.name
        cstnode.type = "parameter"
        cstnode.category = (m, 0, 0)

        old_index = cstnode.index
        new_node = trees.tree_index_node_add(refhidx.index, above, cstnode)

        node_map[old_index] = cstnode.index

        pkey = (0, new_node.name)
        terms = {"k": new_node.name, "l": new_node.name}
        proc.topology_parameters[pkey] = terms

        if m == 0:
            kval = 10.0
            lval = 1.3
        else:
            kval = 50.0
            lval = 109*3.14/180
        cm.topology_terms["k"].values[new_node.name] = [kval]
        cm.topology_terms["l"].values[new_node.name] = [lval]

        sma = newhidx.smarts[old_index]
        refhidx.smarts[new_node.index] = sma
        refhidx.subgraphs[new_node.index] = newhidx.subgraphs[old_index]

    refhidx.topology = cm.topology

    return csys

def chemical_system_cluster_angles(csys: mm.chemical_system, gdb, sep=0.01):
    return chemical_system_cluster_geom(csys, gdb, sep, topology.angle)

def chemical_system_cluster_bonds(csys: mm.chemical_system, gdb, sep=0.001):
    return chemical_system_cluster_geom(csys, gdb, sep, topology.bond)

def smiles_assignment_force_constants(gdb, alpha=-.25, guess_periodicity=True, max_n=3, min_dihedral_k=100, max_dihedral_k=100):

    sag_map = {
        "bond_k": [],
        "bond_l": [],
        "angle_k": [],
        "angle_l": [],
        "torsion_n": [],
        "torsion_p": [],
        "torsion_k": [],
    }

    eid_hess = [eid for eid, e in gdb.entries.items() if assignments.HESSIANS in e.tables]

    for i, eid in enumerate(eid_hess, 1):
        # psys = psystems[eid]
        hessian = gdb.entries[eid].tables[assignments.HESSIANS].values

        # labeled_ic = [x for model in psystems[eid].models for x in model.labels[0]]
        # B = {ic:x for ic, x in zip(icb, B) if ic in labeled_ic}
        # B = list(B.keys()), list(B.values())

        pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, assignments.POSITIONS)
        g = pos.graph

        icb, B = hessians.bmatrix(pos, torsions=True, outofplanes=False, pairs=False ,remove1_3=True, linear_torsions=145.0)

        # xyz = np.vstack([x[0] for x in pos.selections.values()], dtype=float)
        # xyz = xyz.round(12)
        # sym = graphs.graph_symbols(pos.graph)
        # mass = np.array([[vibration.mass_table[sym[n]]]*3 for n in sym])

        # sym = list(sym.values())

        # hess_qm_freq, hess_qm_modes = vibration.hessian_modes(hess_qm, sym, xyz, mass, 0, remove=0, stdtr=True, verbose=False)
        # omega = np.round(hess_qm_freq, 12)
        # omega_qm = omega
        hess_qm_ic = hessians.project_ics(B, hessian)
        hess_qm_ic = [hess_qm_ic[i][i] for i in range(len(hess_qm_ic))]

        ic_qm_fcs = dict(zip(icb, hess_qm_ic))

        sel = {}
        for ic in graphs.graph_bonds(pos.graph):
            if ic in ic_qm_fcs:
                sel[ic] = [ic_qm_fcs[ic]]

        ga = assignments.graph_assignment(pos.smiles, sel, pos.graph)
        sag_map["bond_k"].append(ga)
        bond_l = assignments.graph_assignment_geometry_bonds(pos)
        bond_l.selections = {ic: [x for y in v for x in y] for ic, v in bond_l.selections.items()}
        ga = assignments.graph_assignment(pos.smiles, bond_l.selections, pos.graph)
        sag_map["bond_l"].append(ga)

        sel = {}
        for ic in graphs.graph_angles(pos.graph):
            if ic in ic_qm_fcs:
                sel[ic] = [ic_qm_fcs[ic]]
        ga = assignments.graph_assignment(pos.smiles, sel, pos.graph)
        sag_map["angle_k"].append(ga)

        angle_l = assignments.graph_assignment_geometry_angles(pos)
        angle_l.selections = {ic: [x for y in v for x in y] for ic, v in angle_l.selections.items()}
        ga = assignments.graph_assignment(pos.smiles, angle_l.selections, pos.graph)
        sag_map["angle_l"].append(ga)

        dih_angles = assignments.smiles_assignment_geometry_torsions_nonlinear(pos)

        sel_n = {}
        sel_k = {}
        sel_p = {}
        for angle_ic, angle_val in dih_angles.selections.items():
            bond = geometry.bond(angle_ic[1:3])
            all_dihed = [v for ic, v in dih_angles.selections.items() if geometry.bond(ic[1:3]) == bond]
            angles = [x for y in all_dihed for z in y for x in z]
            vals = [ic_qm_fcs[angle_ic]]
            csys_n0 = [1, 2, 3]
            csys_p0 = [0, 0, 0]
            csys_k0 = [0.0, 0.0, 0.0]
            npk0 = csys_n0, csys_p0, csys_k0

            if guess_periodicity:
                csys_n = []
                csys_p = []
                for i in [*range(1, max_n+1)]:
                    if all(math.cos(i*t) < alpha for t in angles) or all(math.cos(i*t)  >  -alpha for t in angles):
                        csys_n.append(i)
                        csys_p.append(0)
                        break
                        # print(f"Consider {i}")
            else:
                csys_n = [n for n in csys_n0]
                csys_p = [n for n in csys_p0]

            # csys_n = [*range(1,121)]
            # csys_p = [*[0]*120]
            new_k_lst = []
            # new_k = sum(vals)/len(vals)
            deriv = [math.cos, math.sin, math.cos, math.sin]
            sign = [1, -1, -1, 1]
            hq = sum(vals) / len(vals)

            if not csys_n:
                print(f"Could not find any appropriate periodicities for ic={angle_ic} max_n={max_n} alpha={alpha}")
                csys_n = [n for n in range(1, max_n+1)]
                csys_p = [0 for n in range(1, max_n+1)]
                print("Angles:")
                print(angles)
                print("cosines:")
                for i in range(1, max_n+1):
                    print(i, [*(math.cos(i*t) for t in angles)])
                print(angles)

            changed = True
            while changed:
                A = []
                b = []
                changed = False

                if angles:
                    for i in range(len(csys_n)):
                        row = []
                        if i == 0:
                            b.append(hq)
                        else:
                            b.append(0)
                        for n, p in zip(csys_n, csys_p):
                            x = [sign[(i+2)%4] * n**(i+2)*deriv[(i+2) % 4](n*t - p) for t in angles]
                            x = sum(x) / len(x)
                            row.append(x)
                        A.append(row)
                    new_k_lst = np.linalg.solve(A, b)
                    # print("Calculated", new_k_lst)
                else:
                    new_k_lst = [0 for _ in csys_n]
                # at this point we have our new_n and new_p; project onto npk0

                if guess_periodicity:
                    # new_k_lst = project_torsions(npk0, (csys_n, csys_p), angles)
                    new_k = []
                    new_p = []
                    new_n = []
                    changed = False
                    max_k = max_dihedral_k
                    min_k = 1e-4
                    for n,p,k in zip(csys_n, csys_p, new_k_lst):
                        if abs(k) > min_k and abs(k) < max_k:
                            new_n.append(n)
                            new_p.append(p)
                            new_k.append(k)
                            # print("Added n=", n, p, k)

                    if not new_n:
                        print("Warning, all fitted values were out of range")
                        i = arrays.argmin([abs(x) for x in new_k_lst])
                        new_n.append(csys_n[i])
                        k = new_k_lst[i]
                        if k < -max_k:
                            k = -max_k
                        elif k > max_k:
                            k = max_k
                        if k < 0:
                            new_p.append(math.pi)
                            new_k.append(-k)
                        else:
                            new_p.append(0)
                            new_k.append(k)
                    if set(new_n).symmetric_difference(csys_n):
                        changed = True
                    # if changed:
                    csys_n = new_n
                    csys_p = new_p
                    new_k_lst = new_k
                else:
                    new_k = []
                    for k in new_k_lst:
                        if k < -max_k:
                            k = -max_k
                        elif k > max_k:
                            k = max_k
                        new_k.append(k)
                    new_k_lst = new_k

            sel_n[angle_ic] = new_n
            sel_p[angle_ic] = new_p
            sel_k[angle_ic] = new_k

        ga = assignments.graph_assignment(pos.smiles, sel_k, pos.graph)
        sag_map["torsion_k"].append(ga)

        ga = assignments.graph_assignment(pos.smiles, sel_n, pos.graph)
        sag_map["torsion_n"].append(ga)

        ga = assignments.graph_assignment(pos.smiles, sel_p, pos.graph)
        sag_map["torsion_p"].append(ga)

    return sag_map

def reset(reset_config, csys, gdb, psystems=None, verbose=False, ws=None):

    assert type(reset_config) is dict

    reset_bond_k = reset_config.get("bond_k", False)
    reset_bond_l = reset_config.get("bond_l", False)
    reset_angle_k = reset_config.get("angle_k", False)
    reset_angle_l = reset_config.get("angle_l", False)
    reset_torsion_k = reset_config.get("torsion_k", False)
    reset_outofplane_k = reset_config.get("outofplane_k", False)

    guess_periodicity = reset_config.get("dihedral_p", False)
    max_n = reset_config.get("dihedral_max_n", 3)
    alpha = reset_config.get("dihedral_alpha", -.25)
    max_k = reset_config.get("dihedral_max_k", 5.0)
    min_k = reset_config.get("dihedral_max_k", 1e-3)

    # print(datetime.datetime.now(), f"Resetting bonds and angles: bk={reset_bond_k} bl={reset_bond_l} ak={reset_angle_k} al={reset_angle_l} tk={reset_torsion_k} ok={reset_outofplane_k}")
    if psystems is None:
        psystems = gdb_to_physical_systems(gdb, csys)
    if not any([reset_bond_k, reset_bond_l, reset_angle_k, reset_angle_l, reset_torsion_k]):
        return psystems
    eid_hess = [eid for eid, e in gdb.entries.items() if assignments.HESSIANS in e.tables]
    psys_hess = [psystems[eid] for eid in eid_hess]

    reset = set()
    if reset_bond_l:
        mod = mm.chemical_system_reset_bond_lengths(csys, psys_hess)
        if verbose:
            print_chemical_system(csys, show_parameters=[x[2] for x in mod])
        reset.add(0)
    if reset_angle_l:
        mod = mm.chemical_system_reset_angles(csys, psys_hess)
        if verbose:
            print_chemical_system(csys, show_parameters=[x[2] for x in mod])
        reset.add(1)


    t = None
    hvals = {}
    psys_hic_all = {}
    if reset_bond_k or reset_angle_k or reset_torsion_k or reset_outofplane_k:
        psys_hic = {}
        if ws:
            iterable = {}
            for eid in eid_hess:
                icb, B = hessians.bmatrix(psystems[eid].models[0].positions[0], torsions=reset_torsion_k, outofplanes=reset_outofplane_k, pairs=False ,remove1_3=True)
                labeled_ic = [x for model in psystems[eid].models for x in model.labels[0]]
                B = {ic:x for ic, x in zip(icb, B) if ic in labeled_ic}
                B = list(B.keys()), list(B.values())

                iterable[eid] = [[csys, psystems[eid],gdb.entries[eid].tables[assignments.HESSIANS].values], {"B": B}]
            results =  compute.workspace_submit_and_flush(ws, hessians.hessian_project_onto_ics, iterable, verbose=True)
            psys_hic_all.update(results)
        else:
            for i, eid in enumerate(eid_hess, 1):
                psys = psystems[eid]
                hessian = gdb.entries[eid].tables[assignments.HESSIANS].values
                if verbose:
                    print(f"Projecting hessian for EID {eid} {i:8d}/{len(eid_hess)} {psys.models[0].positions[0].smiles}")
                icb, B = hessians.bmatrix(psystems[eid].models[0].positions[0], torsions=reset_torsion_k, outofplanes=reset_outofplane_k, pairs=False ,remove1_3=True)

                labeled_ic = [x for model in psystems[eid].models for x in model.labels[0]]
                B = {ic:x for ic, x in zip(icb, B) if ic in labeled_ic}
                B = list(B.keys()), list(B.values())

                hic = hessians.hessian_project_onto_ics(csys, psys, hessian, verbose=verbose, B=B)
                psys_hic_all[eid] = hic

    if reset_bond_k:
        psys_hics = []
        psyss = []
        for eid in eid_hess:
            psys = psystems[eid]
            hic = psys_hic_all[eid]
            psys_hics.append({k: [[hic.selections[k]]] for k in graphs.graph_bonds(hic.graph) if k in hic.selections})
            psyss.append(psys)

        hvals[0] = mm.chemical_system_groupby_names(csys, 0, psyss, psys_hics)
        psyss.clear()
        psys_hics.clear()

    if reset_angle_k:
        psys_hics = []
        psyss = []
        for eid in eid_hess:
            psys = psystems[eid]
            hic = psys_hic_all[eid]
            psys_hics.append({k: [[hic.selections[k]]] for k in graphs.graph_angles(hic.graph) if k in hic.selections})
            psyss.append(psys)

        hvals[1] = mm.chemical_system_groupby_names(csys, 1, psyss, psys_hics)
        psyss.clear()
        psys_hics.clear()

    if reset_torsion_k:
        psys_hics = []
        psyss = []
        for eid in eid_hess:
            psys = psystems[eid]
            hic = psys_hic_all[eid]
            psys_hics.append({k: [[hic.selections[k]]] for k in graphs.graph_torsions(hic.graph) if k in hic.selections})
            psyss.append(psys)

        hvals[2] = mm.chemical_system_groupby_names(csys, 2, psyss, psys_hics)
        psyss.clear()
        psys_hics.clear()

    if reset_outofplane_k:
        psys_hics = []
        psyss = []
        for eid in eid_hess:
            psys = psystems[eid]
            hic = psys_hic_all[eid]
            psys_hics.append({k: [[hic.selections[k]]] for k in graphs.graph_outofplanes(hic.graph) if k in hic.selections})
            psyss.append(psys)

        hvals[3] = mm.chemical_system_groupby_names(csys, 3, psyss, psys_hics)
        psyss.clear()
        psys_hics.clear()

    psys_hic_all.clear()

    for m, hic in hvals.items():
        for lbl, vals in hic.items():
            if not vals:
                continue
            try:
                csys_k = mm.chemical_system_get_value_list(csys, (m, 'k', lbl))
            except KeyError:
                continue
            if m == 1:
                # new_k = min(vals)
                new_k = sum(vals)/len(vals)
                if reset_outofplane_k or reset_torsion_k:
                    new_k *= 1.0
                new_k_lst = [new_k for _ in csys_k]
            elif m == 0:
                new_k = sum(vals)/len(vals)
                new_k_lst = [new_k for _ in csys_k]
            elif m == 2 or m == 3:

                oldtorsionmethod = False
                if oldtorsionmethod:
                    csys_n = mm.chemical_system_get_value_list(csys, (m, 'n', lbl))
                    csys_p = mm.chemical_system_get_value_list(csys, (m, 'p', lbl))
                    # go into the psys and get the angles
                    angles = []
                    for psys in psys_hess:
                        indices = [ic for ic, terms in psys.models[m].labels[0].items() if lbl in terms.values()]
                        if m == 2:
                            psys_angles = assignments.smiles_assignment_geometry_torsions(psys.models[m].positions[0], indices=indices)
                        elif m == 3:
                            psys_angles = assignments.smiles_assignment_geometry_outofplanes(psys.models[m].positions[0], indices=indices)
                        angles.extend([t for y in psys_angles.selections.values() for x in y for t in x])
                    new_k_lst = []
                    # new_k = sum(vals)/len(vals)
                    for n, p in zip(csys_n, csys_p):
                        k_vals = [hq/(-n**2*math.cos(n*t - p)) for hq, t in zip(vals, angles)]
                        new_k_lst.append(sum(k_vals)/len(k_vals))
                    # new_k = 1.3943109100240065
                else:
                    # these are the original we want to fit to
                    csys_n0 = mm.chemical_system_get_value_list(csys, (m, 'n', lbl))
                    csys_p0 = mm.chemical_system_get_value_list(csys, (m, 'p', lbl))
                    csys_k0 = mm.chemical_system_get_value_list(csys, (m, 'k', lbl))
                    npk0 = csys_n0, csys_p0, csys_k0

                    # go into the psys and get the angles
                    angles = []
                    idiv = []
                    for psys in psys_hess:
                        indices = [ic for ic, terms in psys.models[m].labels[0].items() if lbl in terms.values()]
                        inner_bonds = [geometry.bond(ic[1:3]) for ic in indices]
                        all_dihed = [geometry.bond(ic[1:3]) for ic in psys.models[m].labels[0] if geometry.bond(ic[1:3]) in inner_bonds]
                        idiv.extend([all_dihed.count(ic) for ic in set(all_dihed)])
                        if m == 2:
                            # psys_angles = assignments.smiles_assignment_geometry_torsions(psys.models[m].positions[0], indices=indices)
                            psys_angles = assignments.smiles_assignment_geometry_torsions_nonlinear(psys.models[m].positions[0], indices=indices)
                        elif m == 3:
                            psys_angles = assignments.smiles_assignment_geometry_outofplanes(psys.models[m].positions[0], indices=indices)
                        angles.extend([t for y in psys_angles.selections.values() for x in y for t in x])
                    print("idiv:", idiv)
                    idiv_expected = sum(idiv)/len(idiv)
                    idiv_expected = 1
                    print("Expected idiv (set to 1):", idiv_expected)
                    # print(angles)

                    print("Reference npk0", npk0)
                    if guess_periodicity:
                        csys_n = []
                        csys_p = []
                        for i in [*range(1, max_n+1)]:
                            if all(math.cos(i*t) < alpha for t in angles) or all(math.cos(i*t)  >  -alpha for t in angles):
                                csys_n.append(i)
                                csys_p.append(0)
                                break
                                # print(f"Consider {i}")
                    else:
                        csys_n = [n for n in csys_n0]
                        csys_p = [n for n in csys_p0]

                    # csys_n = [*range(1,121)]
                    # csys_p = [*[0]*120]
                    new_k_lst = []
                    # new_k = sum(vals)/len(vals)
                    deriv = [math.cos, math.sin, math.cos, math.sin]
                    sign = [1, -1, -1, 1]
                    hq = sum(vals) / len(vals)

                    if not csys_n:
                        print(f"Could not find any appropriate periodicities for label={lbl} max_n={max_n} alpha={alpha}")
                        csys_n = [n for n in range(1, max_n+1)]
                        csys_p = [0 for n in range(1, max_n+1)]
                        print("Angles:")
                        print(angles)
                        print("cosines:")
                        for i in range(1, max_n+1):
                            print(i, [*(math.cos(i*t) for t in angles)])
                        print(angles)

                    changed = True
                    while changed:
                        A = []
                        b = []
                        changed = False

                        if angles:
                            for i in range(len(csys_n)):
                                row = []
                                if i == 0:
                                    b.append(hq)
                                else:
                                    b.append(0)
                                for n, p in zip(csys_n, csys_p):
                                    x = [sign[(i+2)%4] * n**(i+2)*deriv[(i+2) % 4](n*t - p) for t in angles]
                                    x = idiv_expected * sum(x) / len(x)
                                    row.append(x)
                                A.append(row)
                            new_k_lst = np.linalg.solve(A, b)
                            print("Calculated", new_k_lst)
                        else:
                            new_k_lst = [0 for _ in csys_n]
                        # at this point we have our new_n and new_p; project onto npk0

                        if guess_periodicity:
                            # new_k_lst = project_torsions(npk0, (csys_n, csys_p), angles)
                            new_k = []
                            new_p = []
                            new_n = []
                            changed = False
                            # max_k = 5
                            # min_k = 1e-3
                            for n,p,k in zip(csys_n, csys_p, new_k_lst):
                                if abs(k) > min_k and abs(k) < max_k:
                                    new_n.append(n)
                                    new_p.append(p)
                                    new_k.append(k)
                                    print("Added n=", n, p, k)

                            if not new_n:
                                print("Warning, all fitted values were out of range")
                                i = arrays.argmin([abs(x) for x in new_k_lst])
                                new_n.append(csys_n[i])
                                new_p.append(csys_p[i])
                                k = new_k_lst[i]
                                if k < -max_k:
                                    k = -max_k
                                elif k > max_k:
                                    k = max_k
                                new_k.append(k)
                            if set(new_n).symmetric_difference(csys_n):
                                changed = True
                            # if changed:
                            csys_n = new_n
                            csys_p = new_p
                            new_k_lst = new_k
                        else:
                            new_k = []
                            for k in new_k_lst:
                                if k < -max_k:
                                    k = -max_k
                                elif k > max_k:
                                    k = max_k
                                new_k.append(k)
                            new_k_lst = new_k


                    if verbose:
                        if guess_periodicity:
                            print(f"Setting {lbl} n from {npk0[0]} to {new_n}")
                            print(f"Setting {lbl} p from {npk0[1]} to {new_p}")
                        print(f"Setting {lbl} k from {csys_k} to {new_k_lst}")
                    if guess_periodicity:
                        mm.chemical_system_set_value_list(csys, (m, 'n', lbl), new_n)
                        mm.chemical_system_set_value_list(csys, (m, 'p', lbl), new_p)
                    # mm.chemical_system_set_value_list(csys, (m, 'k', lbl), new_k_lst)
                    reset.add(m)

            if verbose:
                print(f"Setting {lbl} k from {csys_k} to {new_k_lst}")
            mm.chemical_system_set_value_list(csys, (m, 'k', lbl), new_k_lst)
            reset.add(m)

    hvals.clear()

    if reset:
        psys = {}
        reuse = list(set(range(len(csys.models))).difference(reset))
        for eid, gde in gdb.entries.items():
            tid = assignments.POSITIONS
            gid = list(gde.tables[tid].graphs)[0]
            rid = 0
            pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
            psys[eid] = mm.chemical_system_to_physical_system(csys, [pos], ref=psystems[eid], reuse=list(reuse))
        psystems = psys

    return psystems


def project_torsions(npk0, np1, angles):
    """
    Find the k values
    """
    deriv = [math.cos, math.sin, math.cos, math.sin]
    sign = [1, -1, -1, 1]
    A = []
    b = []
    for i in range(1, 1+len(np1[0])):
        row = [sum(
            (sign[i % 4] * n**i * deriv[i % 4](n*t - p) for t in angles)
        ) for n, p in zip(*np1)]
        A.append(row)
        b.append(sum([
            sign[i % 4] * n**i * k * deriv[i % 4](n*t - p)
            for t in angles for n, p, k in zip(*npk0)]))

    A = np.array(A)
    b = np.array(b)
    x = np.linalg.solve(A, b)
    return x

def reset_project_torsions(csys, gdb, psystems, max_n=6, alpha=-.25, m=2, verbose=False, ws=None):

    reset = set()
    hvals = {}
    psys_hic_all = {}
    psys_hics = []
    psyss = []
    eid_hess = [eid for eid, e in gdb.entries.items() if assignments.HESSIANS in e.tables]
    psys_hess = [psystems[eid] for eid in eid_hess]

    psys_hic = {}
    if ws:
        iterable = {eid: [[csys, psystems[eid],gdb.entries[eid].tables[assignments.HESSIANS].values], {"B": hessians.bmatrix(psystems[eid].models[0].positions[0], torsions=True, outofplanes=False, pairs=False ,remove1_3=True)}] for eid in eid_hess}
        results =  compute.workspace_submit_and_flush(ws, hessians.hessian_project_onto_ics, iterable, verbose=True)
        psys_hic_all.update(results)
    else:
        for i, eid in enumerate(eid_hess, 1):
            psys = psystems[eid]
            hessian = gdb.entries[eid].tables[assignments.HESSIANS].values
            if verbose:
                print(f"Projecting hessian for EID {eid} {i:8d}/{len(eid_hess)} {psys.models[0].positions[0].smiles}")
            hic = hessians.hessian_project_onto_ics(csys, psys, hessian, verbose=verbose, B=hessians.bmatrix(psystems[eid].models[0].positions[0], torsions=True, outofplanes=False, pairs=False ,remove1_3=True))
            psys_hic_all[eid] = hic

    for eid in eid_hess:
        psys = psystems[eid]
        hic = psys_hic_all[eid]
        if m == 2:
            psys_hics.append({k: [[hic.selections[k]]] for k in graphs.graph_torsions(hic.graph)})
        elif m == 3:
            psys_hics.append({k: [[hic.selections[k]]] for k in graphs.graph_outofplanes(hic.graph)})
        else:
            assert False
        psyss.append(psys)

    hvals[m] = mm.chemical_system_groupby_names(csys, m, psyss, psys_hics)
    psyss.clear()
    psys_hics.clear()
    for m, hic in hvals.items():
        for lbl, vals in hic.items():
            if not vals:
                continue
            csys_k = mm.chemical_system_get_value_list(csys, (m, 'k', lbl))

            # these are the original we want to fit to
            csys_n0 = mm.chemical_system_get_value_list(csys, (m, 'n', lbl))
            csys_p0 = mm.chemical_system_get_value_list(csys, (m, 'p', lbl))
            csys_k0 = mm.chemical_system_get_value_list(csys, (m, 'k', lbl))
            npk0 = csys_n0, csys_p0, csys_k0

            # go into the psys and get the angles
            angles = []
            for psys in psys_hess:
                indices = [ic for ic, terms in psys.models[m].labels[0].items() if lbl in terms.values()]
                if m == 2:
                    psys_angles = assignments.smiles_assignment_geometry_torsions(psys.models[m].positions[0], indices=indices)
                elif m == 3:
                    psys_angles = assignments.smiles_assignment_geometry_outofplanes(psys.models[m].positions[0], indices=indices)
                angles.extend([t for y in psys_angles.selections.values() for x in y for t in x])
            # print(angles)

            csys_n = []
            csys_p = []
            for i in range(1, max_n+1):
                if all(math.cos(i*t) < alpha for t in angles):
                    csys_n.append(i)
                    csys_p.append(0)
                    # print(f"Consider {i}")
            # csys_n = [*range(1,121)]
            # csys_p = [*[0]*120]
            new_k_lst = []
            # new_k = sum(vals)/len(vals)
            deriv = [math.cos, math.sin, math.cos, math.sin]
            sign = [1, -1, -1, 1]
            hq = sum(vals) / len(vals)

            if not csys_n:
                print(f"Could not find any appropriate periodicities for label={lbl} max_n={max_n} alpha={alpha}")
                continue
            changed = True
            while changed:
                A = []
                b = []

                for i in range(len(csys_n)):
                    row = []
                    if i == 0:
                        b.append(hq)
                    else:
                        b.append(0)
                    for n, p in zip(csys_n, csys_p):
                        x = [sign[(i+2)%4] * n**(i+2)*deriv[(i+2) % 4](n*t - p) for t in angles]
                        x = sum(x) / len(x) * 2
                        row.append(x)
                    A.append(row)
                new_k_lst = np.linalg.solve(A, b)

                new_k = []
                new_p = []
                new_n = []
                changed = False
                for n,p,k in zip(csys_n, csys_p, new_k_lst):
                    if abs(k) > 1e-4:
                        new_n.append(n)
                        new_p.append(p)
                        new_k.append(k)
                        print("Added n=", n, p, k)
                    else:
                        changed = True
                if changed:
                    csys_n = new_n
                    csys_p = new_p

            # print("Final npk:")
            # for n,p,k in zip(new_n, new_p, new_k):
            #     print(n, p, k)
            # new_k_lst = new_k

            # at this point we have our new_n and new_p; project onto npk0
            new_k_lst = project_torsions(npk0, (new_n, new_p), angles)

            if verbose:
                print(f"Setting {lbl} k from {csys_k} to {new_k_lst}")
                print(f"Setting {lbl} n from {npk0[0]} to {new_n}")
                print(f"Setting {lbl} p from {npk0[1]} to {new_p}")
            mm.chemical_system_set_value_list(csys, (m, 'n', lbl), new_n)
            mm.chemical_system_set_value_list(csys, (m, 'p', lbl), new_p)
            mm.chemical_system_set_value_list(csys, (m, 'k', lbl), new_k_lst)
            reset.add(m)

    hvals.clear()

    if reset:
        psys = {}
        reuse = list(set(range(len(csys.models))).difference(reset))
        for eid, gde in gdb.entries.items():
            tid = assignments.POSITIONS
            gid = list(gde.tables[tid].graphs)[0]
            rid = 0
            pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid, rid)
            psys[eid] = mm.chemical_system_to_physical_system(csys, [pos], ref=psystems[eid], reuse=list(reuse))
        psystems = psys
    return psystems

def parse_xyz(xyzdata):

    N = None
    lines = xyzdata.split('\n')

    if N is None:
        N = int(lines[0].split()[0])

    assert N == int(lines[0].split()[0])

    syms = []
    xyzs = []
    for chunk in arrays.batched(lines, N+2):
        sym = [None]*N
        xyz = [None]*N
        if chunk and chunk[0]:
            for i, line in enumerate(chunk, -2):
                if i >= 0:
                    s, x, y, z = line.split()[:4]
                    sym[i] = s
                    xyz[i] = [[*map(float, (x, y, z))]]
        if all(sym) and all(xyz):
            syms.append(sym)
            xyzs.append(xyz)
    return syms, xyzs

def xyz_to_graph_assignment(gcd, smi, xyzdata: List, indices=None) -> assignments.graph_assignment:
    """
    
    """
    g = graphs.subgraph_as_graph(gcd.smiles_decode(smi))
    sel = {}

    for xyzs in xyzdata:
        for ic, xyz in enumerate(xyzs, 1):
            if (ic,) not in sel:
                sel[ic,] = []
            sel[ic,].extend(xyz)

    return assignments.graph_assignment(smi, sel, g)

