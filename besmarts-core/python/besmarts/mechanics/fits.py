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
from typing import List, Dict, Callable
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
from besmarts.mechanics import smirnoff_models

import scipy.optimize
# import scipy.stats
eid_t = int

PRECISION = configs.precision

POSITIONS = assignments.POSITIONS
HESSIANS = assignments.HESSIANS
GRADIENTS = assignments.GRADIENTS

au2kcal = 627.51
kj2kcal = 1/4.184


class compute_config:
    """
    Configuration for computing some property of an entry in a graph_db. The
    addr is the address of the entry in the graph_db. Keys are from a chemical
    system, and are used to temporarility update the parameters. This is
    generally used to computer numerical gradients and hessians.
    """

    def __init__(self, addr):
        self.addr = addr
        self.keys = []
        self.enable_minimization = True
        self.fc_reset_config = None
        self.method_vib_modes = "min"
        self.analytic = False
        self.analytic_use_gradients = True
        self.dimer_start = 0

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_table]]:

        return []

    def prepare_systems(self, keyset=None):

        gdb = self.GDB
        csys = self.csys
        psys = None
        traj = []
        if keyset is None:
            keyset = {}

        for eid, gde in gdb.entries.items():

            system = assignments.graph_db_graph_to_graph_assignment(
                gdb,
                eid,
                POSITIONS,
            )
            traj.append(system)

            if psys is None:
                psys = mm.chemical_system_to_physical_system(
                    csys,
                    system,
                    ref=self.psys[eid],
                    reuse=self.reuse
                )
                reapply = set()
                for k, v in keyset.items():
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
                                    k[1:]: v for k, v in keyset.items()
                                    if k[0] == m
                                }
                            )

                if self.fc_reset_config:
                    if self.fc_reset_config.get("unique", False):
                        psys = reset_unique(
                            self.fc_reset_config,
                            gde,
                            psys
                        )
                    else:
                        psys = reset(
                            self.fc_reset_config,
                            csys,
                            assignments.graph_db_get_entries(gdb, eid),
                            psystems={eid: psys}
                        )[eid].value
        return traj, psys

class compute_config_energy_total(compute_config):
    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        """
        Compute the energy of the entries and store them in graph_db tables.
        """
        # print("Starting calculation", self)
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.ENERGY
        tid = assignments.POSITIONS

        all_results = []
        for keyset in self.keys:
            results = []
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]]
            psys = None
            traj = []
            traj, psys0 = self.prepare_systems(keyset)
            # for eid, gde in gdb.entries.items():

            #     system = assignments.graph_db_graph_to_graph_assignment(
            #         gdb,
            #         eid,
            #         tid,
            #     )
            #     traj.append(system)

            #     if psys is None:
            #         psys = mm.chemical_system_to_physical_system(
            #             csys,
            #             system,
            #             ref=self.psys[eid],
            #             reuse=self.reuse
            #         )
            #         reapply = set()
            #         for k, v in keyset.items():
            #             mm.physical_system_set_value(psys, k, v)
            #             reapply.add(k[0])
            #         for m in reapply:
            #             procs = csys.models[m].procedures
            #             if len(procs) > 1:
            #                 for _ in range(1, len(psys.models[m].values)):
            #                     psys.models[m].values.pop()
            #                     psys.models[m].labels.pop()
            #                 for proc in procs[1:]:
            #                     psys.models[m] = proc.assign(
            #                         csys.models[m],
            #                         psys.models[m],
            #                         overrides={
            #                             k[1:]: v for k, v in keyset.items()
            #                             if k[0] == m
            #                         }
            #                     )

            #         if self.fc_reset_config:
            #             if self.fc_reset_config.get("unique", False):
            #                 psys = reset_unique(
            #                     self.fc_reset_config,
            #                     gde,
            #                     psys
            #                 )
            #             else:
            #                 psys = reset(
            #                     self.fc_reset_config,
            #                     csys,
            #                     assignments.graph_db_get_entries(gdb, eid),
            #                     psystems={eid: psys}
            #                 )[eid].value

                # ene = objectives.physical_system_energy(psys, csys)
            ene_traj = optimizers_openmm.physical_system_energy_openmm(psys0, csys, pos_traj=traj)
            for ene in ene_traj:
                tbl = assignments.graph_db_table(topology.null)

                # ene = objectives.physical_model_internal_energy(psys, csys)
                # ene_dict = {}
                # inter_dict = {}
                # for m, enes in ene.items():
                    # dimer_interaction = {ic: v for ic, v in enes.items() if ic[0][1] <= n and ic[1][1] > n}
                    # inter_dict[m] = sum([x for y in dimer_interaction.values() for x in y])
                    # total = {ic: v for ic, v in enes.items()}
                    # ene_dict[m] = sum([x for y in total.values() for x in y])

                tbl.values.append(ene)
                    # print("Calculated energy is", ene)
                    # ene = objectives.physical_system_energy(psys, csys)

                r = {tbl_idx: tbl}
                # print("Result is\n", r)
                results.append(r)
            all_results.append(results)
        # print("Done calculating", self)
        return all_results

class compute_config_energy_sapt(compute_config):
    """
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        """
        Compute the energy of the entries and store them in graph_db tables.
        """
        # print("Starting calculation", self)
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.SAPT
        tid = assignments.POSITIONS
        n = self.dimer_start

        all_results = []
        for keyset in self.keys:
            results = []
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]]
            psys = None
            traj = []
            traj, psys0 = self.prepare_systems(keyset)
            for (eid, gde), system in zip(gdb.entries.items(), traj):
                psys = mm.chemical_system_to_physical_system(
                    csys,
                    system,
                    ref=psys0,
                    reuse=self.reuse
                )
                ene = objectives.physical_model_internal_energy(psys, csys)
                ene_dict = {}
                inter_dict = {}
                for m, enes in ene.items():
                    dimer_interaction = {ic: v for ic, v in enes.items() if ic[0][1] <= n and ic[1][1] > n}
                    inter_dict[m] = sum([x for y in dimer_interaction.values() for x in y])
                    total = {ic: v for ic, v in enes.items()}
                    ene_dict[m] = sum([x for y in total.values() for x in y])

                tbl = assignments.graph_db_table(topology.null)
                tbl.values.append([inter_dict, ene_dict])
                    # print("Calculated energy is", ene)
                    # ene = objectives.physical_system_energy(psys, csys)

                r = {tbl_idx: tbl}
                # print("Result is\n", r)
                results.append(r)

            # for eid, gde in gdb.entries.items():
            #     tbl = assignments.graph_db_table(topology.null)
            #     rids = assignments.graph_db_table_get_row_ids(gde[tid])

            #     system = assignments.graph_db_graph_to_graph_assignment(
            #         gdb,
            #         eid,
            #         tid,
            #     )

            #     psys = mm.chemical_system_to_physical_system(
            #         csys,
            #         system,
            #         ref=self.psys[eid],
            #         reuse=self.reuse
            #     )
            #     reapply = set()
            #     for k, v in keyset.items():
            #         mm.physical_system_set_value(psys, k, v)
            #         reapply.add(k[0])
            #     for m in reapply:
            #         procs = csys.models[m].procedures
            #         if len(procs) > 1:
            #             for _ in range(1, len(psys.models[m].values)):
            #                 psys.models[m].values.pop()
            #                 psys.models[m].labels.pop()
            #             for proc in procs[1:]:
            #                 psys.models[m] = proc.assign(
            #                     csys.models[m],
            #                     psys.models[m],
            #                     overrides={
            #                         k[1:]: v for k, v in keyset.items()
            #                         if k[0] == m
            #                     }
            #                 )

            #     if self.fc_reset_config:
            #         if self.fc_reset_config.get("unique", False):
            #             psys = reset_unique(
            #                 self.fc_reset_config,
            #                 gde,
            #                 psys
            #             )
            #         else:
            #             psys = reset(
            #                 self.fc_reset_config,
            #                 csys,
            #                 assignments.graph_db_get_entries(gdb, eid),
            #                 psystems={eid: psys}
            #             )[eid].value

            #     # ene = objectives.physical_system_energy(psys, csys)
            #     ene = objectives.physical_model_internal_energy(psys, csys)
            #     ene_dict = {}
            #     inter_dict = {}
            #     for m, enes in ene.items():
            #         dimer_interaction = {ic: v for ic, v in enes.items() if ic[0][1] <= n and ic[1][1] > n}
            #         inter_dict[m] = sum([x for y in dimer_interaction.values() for x in y])
            #         total = {ic: v for ic, v in enes.items()}
            #         ene_dict[m] = sum([x for y in total.values() for x in y])

            #     tbl.values.append([inter_dict, ene_dict])
            #         # print("Calculated energy is", ene)
            #         # ene = objectives.physical_system_energy(psys, csys)

            #     r = {tbl_idx: tbl}
            #     # print("Result is\n", r)
            #     results.append(r)
            all_results.append(results)
        # print("Done calculating", self)
        return all_results


class compute_config_gradient(compute_config):
    """
    Configuration for computing the MM gradient
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        """
        Compute the gradient of the entries and store them in graph_db tables.
        compute_config_gradient::run
        """
        csys = self.csys
        gdb = self.GDB
        tbl_idx = GRADIENTS
        tid = assignments.POSITIONS

        # assume all have same topo

        all_results = []
        for keyset in self.keys:
            results = []
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]]
            traj = []
            psys = None
            traj, psys = self.prepare_systems(keyset)
            # for eid, gde in gdb.entries.items():

            #     system = assignments.graph_db_graph_to_graph_assignment(
            #         gdb,
            #         eid,
            #         tid,
            #     )
            #     traj.append(system)

            #     if psys is None:
            #         psys = mm.chemical_system_to_physical_system(
            #             csys,
            #             system,
            #             ref=self.psys[eid],
            #             reuse=self.reuse
            #         )

            #         reapply = set()
            #         for k, v in keyset.items():
            #             mm.physical_system_set_value(psys, k, v)
            #             reapply.add(k[0])
            #         for m in reapply:
            #             procs = csys.models[m].procedures
            #             if len(procs) > 1:
            #                 for _ in range(1, len(psys.models[m].values)):
            #                     psys.models[m].values.pop()
            #                     psys.models[m].labels.pop()
            #                 for proc in procs[1:]:
            #                     psys.models[m] = proc.assign(
            #                         csys.models[m],
            #                         psys.models[m],
            #                         overrides={
            #                             k[1:]: v for k, v in keyset.items()
            #                             if k[0] == m
            #                         }
            #                     )
            #         if self.fc_reset_config:
            #             if self.fc_reset_config.get("unique", False):
            #                 psys = reset_unique(
            #                     self.fc_reset_config,
            #                     gde,
            #                     psys
            #                 )
            #             else:
            #                 psys = reset(
            #                     self.fc_reset_config,
            #                     csys,
            #                     assignments.graph_db_get_entries(gdb, [eid]),
            #                     psystems={eid: psys}
            #                 ).value[eid]

            traj_gx = optimizers_openmm.physical_system_gradient_openmm(
                psys,
                csys,
                pos_traj=traj
            )
            for gx0 in traj_gx:
                # gx0 = objectives.physical_system_gradient(psys, csys)
                # dgx = arrays.array_difference(gx0, gx0)
                # assert all((abs(i) < 1e-10 for i in dgx))
                # e0 = optimizers_openmm.physical_system_energy_openmm(
                #     psys,
                #     csys
                # )
                # e1 = objectives.physical_system_energy(psys, csys)
                # breakpoint()

                gx = arrays.array_round(gx0, PRECISION)
                # print(f"GRADIENT EID {eid}", gx)
                tbl = assignments.graph_db_table(topology.atom)
                tbl.values = []
                tbl.values.extend(gx)


                    # gq = objectives.physical_system_gradient_internal(psys, csys)
                    # tbl.values['q'] = gq

                    # qic = objectives.physical_system_internals(psys, csys)
                    # tbl.values['ic'] = qic

                    # jac = objectives.physical_system_bmatrix(psys, csys)
                    # tbl.values['b'] = jac

                r = {tbl_idx: tbl}
                # print("Result is\n", gx)
                results.append(r)
            all_results.append(results)
        # print("Done calculating", self)
        return all_results


class compute_config_hessian(compute_config):
    """
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_entry]]:
        """
        Compute the hessian of the entries and store them in graph_db tables.
        compute_config_hessian::run
        """
        # print("Starting calculation", self)
        csys = self.csys
        gdb = self.GDB
        tbl_idx = HESSIANS
        tid = assignments.POSITIONS

        all_results = []
        for ki, keyset in enumerate(self.keys, 1):
            results = []
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]]
            traj, psys0 = self.prepare_systems(keyset)

            for (eid, gde), system in zip(gdb.entries.items(), traj):
                tbl_hess = assignments.graph_db_table(topology.null)
                tbl_grad = assignments.graph_db_table(topology.null)
                tbl_pos = assignments.graph_db_table(topology.atom)
                # rids = assignments.graph_db_table_get_row_ids(gde[tid])
                psys = mm.chemical_system_to_physical_system(
                    csys,
                    system,
                    ref=psys0,
                    reuse=range(len(psys0.models))
                )
                # reapply = set()
                # for k, v in keyset.items():
                #     mm.physical_system_set_value(psys, k, v)
                #     reapply.add(k[0])
                # for m in reapply:
                #     procs = csys.models[m].procedures
                #     if len(procs) > 1:
                #         # this is for post procedures like mixing rules
                #         for _ in range(1, len(psys.models[m].values)):
                #             psys.models[m].values.pop()
                #             psys.models[m].labels.pop()
                #         for proc in procs[1:]:
                #             psys.models[m] = proc.assign(
                #                 csys.models[m],
                #                 psys.models[m],
                #                 overrides={
                #                     k[1:]: v for k, v in keyset.items()
                #                     if k[0] == m
                #                 }
                #             )

                # pos = psys.models[0].positions
                N3 = sum([len(x) for posi in system for x in posi.selections.values()])*3
                # if self.fc_reset_config:
                #     if self.fc_reset_config.get("unique", False):
                #         psys = reset_unique(
                #             self.fc_reset_config,
                #             gde,
                #             psys,
                #             verbose=True
                #         )
                #     else:
                #         psys = reset(
                #             self.fc_reset_config,
                #             csys,
                #             assignments.graph_db_get_entries(gdb, [eid]),
                #             psystems={eid: psys}
                #         ).value[eid]

                pos = psys.models[0].positions
                if self.enable_minimization or self.method_vib_modes == "min":
                    minpos = optimizers_openmm.optimize_positions_openmm(csys, psys, tol=1e-12)
                    pos = minpos
                    tbl_pos.values = minpos
                    psys.models[0].positions = minpos


                # print("PSYS:")
                # for m, pm in enumerate(psys.models):
                #     print(m)
                #     for ic, terms in pm.values[0].items():
                #         print(m, ic, terms)

                # hess_mm = objectives.physical_system_hessian(
                #     psys,
                #     csys,
                #     h=1e-4
                # )
                hess_mm = None
                try:
                    if self.analytic:
                        hess_mm = objectives.physical_system_hessian_analytic(psys, csys, use_gradients=self.analytic_use_gradients)
                    else:
                        hess_mm = optimizers_openmm.physical_system_hessian_openmm(
                            psys,
                            csys,
                            h=1e-4
                        )
                except ArithmeticError as e:
                    pass
                # xyz = list([x[0] for x in pos.selections.values()])

                # sym = graphs.graph_symbols(pos.graph)
                # mass = [[vibration.mass_table[sym[n]]]*3 for n in sym]
                # sym = list(sym.values())

                # for i, row in enumerate(hess_mm):
                #     hess_mm[i] = arrays.array_scale(row, 1/4.184)

                # freq, modes, DL = vibration.hessian_modes(
                #     hess_mm,
                #     sym,
                #     xyz,
                #     mass,
                #     0,
                #     remove=0,
                #     stdtr=True,
                #     return_DL=True,
                #     verbose=False
                # )
                # print(freq)

                # print(list(hessians.hessian_frequencies(
                #     pos.graph, hess_mm, None, DL, None, None, None)
                # ))
                s = 1/4.184
                if hess_mm is None:
                    s = 0.0
                    hess_mm = [*[[*[0.0]*N3]]*N3]

                hx = []
                # print("MM Hessian")
                for i, row in enumerate(hess_mm):
                    hx.append(arrays.array_scale(row, 1/4.184))
                    # hx.append(arrays.array_scale(row, 1))
                    # print(hx[-1])
                tbl_hess.values.extend(hx)

                # gx = optimizers_openmm.physical_system_gradient_openmm(
                #     psys,
                #     csys
                # )
                gx = list([0.0] * N3)
                tbl_grad.values.extend(gx)

                r = {
                    POSITIONS: tbl_pos,
                    HESSIANS: tbl_hess,
                    GRADIENTS: tbl_grad
                }
                if not self.enable_minimization:
                    r.pop(POSITIONS)

                results.append(r)
            all_results.append(results)
        return all_results


class compute_config_position(compute_config):
    """
    """

    def run(self) -> List[Dict[assignments.tid_t, assignments.graph_db_table]]:
        """
        Compute the positions that minimize the energy of the entries and store
        them in graph_db tables.
        compute_config_position::run
        """
        csys = self.csys
        gdb = self.GDB
        tbl_idx = assignments.POSITIONS
        tid = assignments.POSITIONS

        all_results = []
        for keyset in self.keys:
            results = []
            results: List[Dict[assignments.tid_t, assignments.graph_db_table]]
            traj, psys0 = self.prepare_systems(keyset)
            for (eid, gde), system in zip(gdb.entries.items(), traj):
                gdt = gdb.entries[eid].tables[tid]
                gids = gdt.graphs
                tbl = assignments.graph_db_table(topology.atom)
                # pos0 = assignments.graph_db_graph_to_graph_assignment(
                #     gdb,
                #     eid,
                #     tid,
                # )
                # psys = mm.chemical_system_to_physical_system(
                #     csys,
                #     pos0,
                #     ref=psys0,
                #     reuse=range(len(psys0.models))
                # )
                psys = mm.chemical_system_to_physical_system(
                    csys,
                    system,
                    ref=psys0,
                    reuse=range(len(psys0.models))
                )
                # psys = copy.deepcopy(psys)

                # reapply = set()
                # for k, v in keyset.items():
                #     mm.physical_system_set_value(psys, k, v)
                #     reapply.add(k[0])
                # for m in reapply:
                #     procs = csys.models[m].procedures
                #     if len(procs) > 1:
                #         for _ in range(1, len(psys.models[m].values)):
                #             psys.models[m].values.pop()
                #             psys.models[m].labels.pop()
                #         for proc in procs[1:]:
                #             psys.models[m] = proc.assign(
                #                 csys.models[m],
                #                 psys.models[m],
                #                 overrides={
                #                     k[1:]: v for k, v in keyset.items()
                #                     if k[0] == m
                #                 }
                #             )
                # if self.fc_reset_config:
                #     if self.fc_reset_config.get("unique", False):
                #         psys = reset_unique(
                #             self.fc_reset_config,
                #             gde,
                #             psys
                #         )
                #     else:
                #         psys = reset(
                #             self.fc_reset_config,
                #             csys,
                #             assignments.graph_db_get_entries(gdb, [eid]),
                #             psystems={eid: psys}
                #         ).value[eid]


                # print(f"Initial xyz for EID {eid}:")
                # print("\n".join(print_xyz(pos0)))

                opt = optimizers_openmm.optimize_positions_openmm
                # opt = optimizers_scipy.optimize_positions_scipy
                pos = opt(
                    csys,
                    psys,
                    step_limit=self.step_limit,
                    tol=self.tol
                )
                # print(f"Optimized xyz for EID {eid}:")
                # print("\n".join(["\n".join(print_xyz(posi)) for posi in pos]))
                for pi, (posi, gid) in enumerate(zip(pos, gids)):

                    gdg = assignments.graph_assignment_to_graph_db_graph(
                        posi,
                        topology.atom
                    )
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
        Compute the penalty function of set of parameters.
        compute_config_penalty::run
        """
        # csys = self.csys

        results = {}
        for k, ref in self.targets.items():
            val = keys.get(k)
            if val is not None:
                results[k] = round(val - ref, PRECISION)
        return results


class objective_config:
    """
    Configuration of an objective function for force field fitting. The
    objective has two parts. The first part is how the properties are
    calculated and is configured with a compute_config. The next is how the
    objective is computed from the properties and is configured here.
    """

    def __init__(self, addr, include=True, scale=1, coordinate_system="C"):
        self.addr = addr
        self.scale = scale
        self.step_limit = 200
        self.tol = 1e-5
        self.coordinate_system = coordinate_system
        self.include = include
        self.cache = {}
        self.batch_size = None
        self.enable_minimization = True
        self.fc_reset_config = None
        self.method_vib_modes = "min"
        # forward difference 1st order
        self.grad_mode = "f1"
        self.verbose = 0
        self.weights = None
        self.gdb = None



        self.fit_models = {}
        self.fit_symbols = {}

    def get_task(
        self,
        GDB,
        csys,
        keys=None,
        psys=None,
        reuse=None
    ) -> compute_config:
        """
        Build and return a compute_config that will compute the needed
        properties.
        """
        cc = compute_config(self.addr)

        if self.gdb is None:
            self.gdb = assignments.graph_db_get(GDB, self.addr)
        cc.GDB = self.gdb
        cc.csys = csys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.psys = psys
        cc.reuse = reuse
        cc.fc_reset_config = self.fc_reset_config
        return cc

    def compute_gradient(
        self,
        ref,
        evals,
        h
    ):
        """
        Compute the parameter gradients of this objective.
        """

        f = []
        for D in evals:
            dxa = []
            for dtbls in D:
                for tid, dtbl in dtbls.items():
                    for gid, dgdg in dtbl.graphs.items():
                        for rid, dgdr in dgdg.rows.items():
                            for cid, dgdc in dgdr.columns.items():
                                for dv in dgdc.selections.values():
                                    dxa.extend(dv)
            f.append(dxa)

        if self.grad_mode == "f1":
            dxdp, d2xdp2 = finite_difference_forward_1(f, h)
        elif self.grad_mode == "c2":
            dxdp, d2xdp2 = finite_difference_central_2(f, h)

        return dxdp, d2xdp2

    # def compute_gradient_2pt(
    #     self,
    #     ref,
    #     X0,
    #     E,
    #     D: List[Dict[int, assignments.graph_db_table]],
    #     h
    # ):
    #     """
    #     Compute the parameter gradients of this objective.
    #     """

    #     dxa = []
    #     dxb = []
    #     for etbls, dtbls in zip(E, D):
    #         for tid, dtbl in dtbls.items():
    #             etbl = etbls[tid]
    #             for gid, dgdg in dtbl.graphs.items():
    #                 egdg = etbl.graphs[gid]
    #                 for rid, dgdr in dgdg.rows.items():
    #                     egdr = egdg.rows[rid]
    #                     for cid, dgdc in dgdr.columns.items():
    #                         egdc = egdr.columns[cid]
    #                         for ev, dv in zip(
    #                             egdc.selections.values(),
    #                             dgdc.selections.values()
    #                         ):
    #                             dxa.extend(ev)
    #                             dxb.extend(dv)

    #     dxdp = arrays.array_difference(dxb, dxa)
    #     dxdp = arrays.array_round(dxdp, PRECISION)
    #     dxdp = arrays.array_scale(dxdp, 1/(2*h))

    #     d2xdp2 = arrays.array_add(
    #         arrays.array_scale(X0, -2),
    #         arrays.array_add(dxb, dxa)
    #     )
    #     d2xdp2 = arrays.array_round(d2xdp2, PRECISION)
    #     d2xdp2 = arrays.array_scale(d2xdp2, 1.0/(h*h))

    #     return dxdp, d2xdp2

    def compute_diff(
        self,
        GDB: assignments.graph_db,
        D: List[Dict[int, assignments.graph_db_table]],
        verbose=False
    ):
        """
        Compute the difference of the computed properties versus a reference
        """

        obj = []
        out = []
        for (eid, gde), tbls in zip(GDB.entries.items(), D):
            rmse = []
            if self.verbose > 1:
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

                            for ic, x1, x0 in zip(gdc0.selections, x1l, x0l):
                                x = arrays.array_difference(x1, x0)
                                x = arrays.array_round(x, PRECISION)
                                obj.extend(x)
                                rmse.append(arrays.array_inner_product(x, x))
                                if self.verbose > 1:
                                    xyz1 = " ".join([f"{xi:8.4f}" for xi in x1])
                                    xyz0 = " ".join([f"{xi:8.4f}" for xi in x0])
                                    dxyz = " ".join([f"{xi:8.4f}" for xi in  x])
                                    addr = f"    {eid:4d} {gid:2d} {rid:2d} "
                                    addr += f"{cid:2d} {ic[0]:3d}"
                                    out_str  = f" MM: {xyz1} QM: {xyz0} dP: {dxyz}"
                                    out.append(addr + out_str)
            N = len(rmse)
            if len(rmse):
                rmse = (sum(rmse)/len(rmse))**.5
            else:
                rmse = 0
            X2 = arrays.array_inner_product(obj, obj)
            if self.verbose > 0:
                out.append(" ".join([
                    f"Total Position SSE: {X2:10.5f} A^2",
                    f"RMSE: {rmse:10.5f} A N={N}"
                ]))
        return returns.success(arrays.array_round(obj, PRECISION), out=out, err=[])


class objective_config_gradient(objective_config):

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.enable_minimization = False
        self.force_mode = "default"

    def get_task(
        self,
        GDB: assignments.graph_db,
        csys,
        keys=None,
        psys=None,
        reuse=None
    ) -> compute_config_gradient:

        cc = compute_config_gradient(self.addr)
        cc.GDB = assignments.graph_db_get(GDB, self.addr)

        cc.csys = csys
        cc.psys = psys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.reuse = reuse
        cc.fc_reset_config = self.fc_reset_config
        return cc

    def compute_gradient(self, ref, evals, h):


        dxdp = None
        d2xdp2= None

        # for eid in self.addr.eid:
        #     f = []
        #     gde = ref.entries[eid]
        #     # parameter steps
        f = []
        for E0 in evals:
            dxa = []
            # one for each eid
            for gde, dtbls in zip(ref.entries.values(), E0):
                TID = assignments.GRADIENTS
                dtbl = dtbls[TID]
                # for tid, dtbl in dtbls.items():
                x0 = dtbl.values
                if self.force_mode == "magnitude":
                    x0 = [arrays.array_magnitude(x0[i:i+3]) for i in range(0, len(x0), 3)]
                # if self.coordinate_system == "D":
                #     gid = list(gde[POSITIONS].graphs)[0]
                #     g = ref.graphs[gid]
                #     smiles = ref.smiles[gid]
                #     sel = gde[POSITIONS][gid][0][0].selections
                #     sel = {k: [v] for k, v in sel.items()}
                #     pos = assignments.graph_assignment(smiles, sel, g)
                #     x0 = array_numpy.dlcmatrix_project_gradients(pos, x0)

                dxa.extend(x0)
            f.append(dxa)

        if self.grad_mode == "f1":
            dxdp, d2xdp2 = finite_difference_forward_1(f, h)
        elif self.grad_mode == "c2":
            dxdp, d2xdp2 = finite_difference_central_2(f, h)

        # if dxdp is None:
        #     dxdp = dxdpi
        #     d2xdp2 = d2xdp2i
        # else:
        #     dxdp = arrays.array_add(dxdp, dxdpi)
        #     d2xdp2 = arrays.array_add(d2xdp2, d2xdp2i)


        return dxdp, d2xdp2

    def compute_diff(
        self,
        GDB: assignments.graph_db,
        D: List[Dict[int, assignments.graph_db_table]],
        verbose=False
    ):
        """
        """

        verbose=True
        out_all = []
        X2 = 0
        N = 0
        XX = []
        all_rmse = []
        coord = "DLC" if self.coordinate_system == "D" else "CART"
        for (eid, gde), tbls in zip(GDB.entries.items(), D):
            X = []
            out = []
            rmse = []

            for tid, tbl in tbls.items():
                x0cart = x0 = gde.tables[tid].values
                x1cart = x1 = tbl.values

                if self.coordinate_system == "D":
                    gid = list(gde.tables[POSITIONS].graphs)[0]
                    g = GDB.graphs[gid]
                    smiles = GDB.smiles[gid]
                    sel = gde.tables[POSITIONS][gid][0][0].selections
                    sel = {k: [v] for k, v in sel.items()}
                    pos = assignments.graph_assignment(smiles, sel, g)

                    x1 = array_numpy.dlcmatrix_project_gradients(pos, x1cart)

                    if "D:x0" not in self.cache:
                        self.cache['D:x0cart'] = x0cart
                        # x0 = array_numpy.dlcmatrix_project_gradients(pos, x0)
                        self.cache['D:x0'] = (
                            array_numpy.dlcmatrix_project_gradients(pos, x0cart)
                        )
                    else:
                        x0 = self.cache['D:x0']
                        x0cart = self.cache['D:x0cart']
                else:
                    x0 = x0cart
                    x1 = x1cart
                    if self.force_mode == "magnitude":
                        x0 = [arrays.array_magnitude(x0[i:i+3]) for i in range(0, len(x0), 3)]
                        x1 = [arrays.array_magnitude(x1[i:i+3]) for i in range(0, len(x1), 3)]

                x = arrays.array_difference(x1, x0)
                xcart = arrays.array_difference(x1cart, x0cart)
                nlines = len(out)
                if self.verbose > 0:
                    mode = "magnitude" if self.force_mode == "magnitude" else ""
                    if self.verbose > 1:
                        out.append(f" EID {eid} TID {tid} Gradient {mode} objective:")
                    if mode == "magnitude":

                        for i in range(len(x0)):
                            xyz1 = f"{x1[i]:12.4e}"
                            xyz0 = f"{x0[i]:12.4e}"
                            dxyz = f"{x[i]:12.4e}"
                            if self.verbose > 1:
                                out.append(f"   {i+1:4d} MM: {xyz1} QM: {xyz0} dG: {dxyz}")
                            rmse.append((x1[i] - x0[i])**2)
                    else:
                        totalx0 = [0.0, 0.0, 0.0]
                        totalx1 = [0.0, 0.0, 0.0]
                        totalxx = [0.0, 0.0, 0.0]
                        for i in range(len(x1cart)//3):
                            xx = xcart[3*i:3*i+3]
                            totalx0 = arrays.array_add(totalx0, x0cart[3*i:3*i+3])
                            totalx1 = arrays.array_add(totalx1, x1cart[3*i:3*i+3])
                            totalxx = arrays.array_add(totalxx, xx)
                            xyz1 = " ".join([f"{xi:12.4e}" for xi in x1cart[3*i:3*i+3]])
                            xyz0 = " ".join([f"{xi:12.4e}" for xi in x0cart[3*i:3*i+3]])
                            dxyz = " ".join([f"{xi:12.4e}" for xi in xx])
                            if self.verbose > 1:
                                out.append(f"   {i+1:4d} MM: {xyz1} QM: {xyz0} dG: {dxyz}")
                            rmse.append(arrays.array_inner_product(xx, xx))
                    # if self.verbose > 1:
                    #     xyz1 = " ".join([f"{xi:12.4e}" for xi in totalx1])
                    #     xyz0 = " ".join([f"{xi:12.4e}" for xi in totalx0])
                    #     dxyz = " ".join([f"{xi:12.4e}" for xi in totalxx])
                    #     out.append(f"   Sum: MM: {xyz1} QM: {xyz0} dG: {dxyz}")

                    # show the IC grads
                    # i = len(out)
                    # out.append(f" EID {eid} IC gradients:")
                    # gq = tbl.values['q']
                    # gqic = tbl.values['ic']
                    # jac = tbl.values['b']
                    # for m, vals in gq.items():
                    #     for ic, v in vals.items():
                    #         q0 = " ".join([f"{xi:12.4e}" for xi in v])
                    #         icq = " ".join([f"{xi:12.4e}" for xi in gqic[m][ic]])
                    #         b = " ".join([f"{xi:12.4e}" for x in jac[m].selections[ic] for xi in x[0]])
                    #         out.append(f"{str(ic)} {q0} kJ/mol/q q: {icq} b: {b}")
                    if verbose:
                        for i in range(len(out)-nlines, len(out)):
                            print(out[i])
                X.extend(x)
                XX.extend(x)
            N = len(rmse)
            if rmse:
                all_rmse.extend(rmse)
                rmse = (sum(rmse)/len(rmse))**.5
            else:
                rmse = 0


            X = arrays.array_inner_product(X, X)
            X2 += X
            if self.verbose > 0:
                out.append(" ".join([
                    f"EID {eid} {coord} Gradient SSE: {X:10.5f} (kJ/mol/A)^2",
                    f"CART RMSE: {rmse:10.5f} kJ/mol/A N: {N}"
                ]))
            out_all.extend(out)

        out = out_all
        # X2 = arrays.array_inner_product(X, X)
        rmse = all_rmse
        N = len(rmse)
        if rmse:
            rmse = (sum(rmse)/len(rmse))**.5
        else:
            rmse = 0
        if self.verbose > 0:
            out.append(" ".join([
                f"Total {coord} Gradient SSE: {X2:10.5f} (kJ/mol/A)^2",
                f"CART RMSE: {rmse:10.5f} kJ/mol/A N: {N}"
            ]))
        return returns.success(arrays.array_round(XX, PRECISION), out=out, err=[])


class objective_config_hessian(objective_config):

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.enable_minimization = kwds.get("enable_minimization", False)
        self.method_vib_modes = kwds.get("method_vib_modes", "min")
        if self.method_vib_modes == "min":
            self.enable_minimization = True

        self.freq_range = (-5000, 5000)
        self.analytic = False
        self.analytic_use_gradients = True

    def get_task(
        self,
        GDB: assignments.graph_db,
        csys,
        keys=None,
        psys=None,
        reuse=None
    ) -> compute_config_hessian:
        cc = compute_config_hessian(self.addr)
        cc.GDB = assignments.graph_db_get(GDB, self.addr)

        cc.csys = csys
        cc.psys = psys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.reuse = reuse
        cc.enable_minimization = self.enable_minimization
        cc.method_vib_modes = self.method_vib_modes
        cc.fc_reset_config = self.fc_reset_config

        cc.analytic = self.analytic
        cc.analytic_use_gradients = self.analytic_use_gradients

        return cc

    def compute_gradient(self, ref, evals, h):

        f = []
        eid = self.addr.eid[0]
        gde = ref.entries[eid]

        pos = assignments.graph_db_graph_to_graph_assignment(
            ref,
            eid,
            POSITIONS,
        )


        for E in evals:
            dxa = []
            for etbls in E:
                # for tid, dtbl in dtbls.items():

                # x0 = etbls[tid].values
                # x1 = dtbl.values

                # breakpoint()
                # etbls[POSITIONS]

                # gid = list(gde.tables[POSITIONS].graphs)[0]
                # g = ref.graphs[gid]
                DL = self.cache['IC:dl']
                (ics, B, B2) = self.cache['IC:b']


                xyz = []
                for pi, posi in enumerate(pos):
                    for xyzi in posi.selections.values():
                        xyz.extend(xyzi)

                sym = [s for posi in pos for s in graphs.graph_symbols(posi.graph).values()]
                mass = np.array([[vibration.mass_table[s]]*3 for s in sym])

                hess_mm = etbls[HESSIANS].values

                if self.method_vib_modes in ["min", "mm", "mm+"]:

                    remove = 0
                    if self.method_vib_modes == "mm":
                        remove = 0

                    x1, mm_modes, mm_dl = vibration.hessian_modes(
                        hess_mm,
                        sym,
                        xyz,
                        mass,
                        0,
                        remove=remove,
                        stdtr=True,
                        return_DL=True
                    )
                    if x1 is None:
                        dxa.extend([0 for _ in range(len(hess_mm))])

                    else:
                        cost = np.abs(np.dot(np.array(DL).T, mm_dl))
                        # s_row, s_col = scipy.optimize.linear_sum_assignment(1 - cost)
                        # x1 = [x1[i] for i in s_row]
                        if self.method_vib_modes != "mm+":
                            x1 = list(sorted(x1))
                        else:
                            DL = self.cache['IC:dl']
                            cost = np.abs(np.dot(np.array(DL).T, mm_dl))
                            costheta = np.degrees(np.arccos(np.dot(np.array(DL).T, mm_dl)))
                            costheta[costheta > 90] -= 180
                            # s_row, s_col = scipy.optimize.linear_sum_assignment(cost)
                            # x1 = [x1[i] for i in s_col]
                            x1 = list(x1)
                            # angle = np.arccos(
                            #     np.abs(cost[s_row, s_col]).sum()/len(cost)
                            # )
                            angle = np.arccos(np.diag(cost))
                            cost = 1 - np.abs(np.dot(np.array(DL).T, mm_dl))
                            costheta = np.degrees(np.arccos(np.dot(np.array(DL).T, mm_dl)))
                            costheta[costheta > 90] -= 180
                            s_row, s_col = scipy.optimize.linear_sum_assignment(cost)
                            x1 = [x1[i] for i in s_col]
                            # x1 = list(x1)

                elif self.method_vib_modes == "qm":

                    hess_mm = etbls[HESSIANS].values
                    grad_mm = list(etbls[GRADIENTS].values)

                    # x1mat = list(hessians.hessian_transform(
                    #     g, hess_mm, grad_mm, DL, ics, B, B2)
                    # )
                    # x1 = [x1mat[i][i] for i in range(len(x1mat))]
                    # x1offdiag = [
                    #     x1mat[i][j] for i in range(len(x1mat)) for j in range(i)
                    # ]
                    # x1.extend(x1offdiag)
                    x1 = list(hessians.hessian_frequencies(
                        mass, hess_mm, grad_mm, DL, ics, B, B2)
                    )
                    # x1 = arrays.array_scale(x1, 2)
                elif self.method_vib_modes == "qm+":
                    grad_mm = list(etbls[GRADIENTS].values)
                    # Need Q' which is Q and the extra R modes
                    # DL2 = vibration.gram_schmidt(DL.T).T
                    hess_qm = gde.tables[HESSIANS].values
                    hess_mm = etbls[HESSIANS].values

                    x1 = list(hessians.hessian_frequencies(
                        mass, hess_mm, grad_mm, DL, ics, B, B2)
                    )
                    dH = [arrays.array_difference(hess_mm[i], hess_qm[i]) for i in range(len(hess_qm))]
                    # dHw, _, _= vibration.hessian_modes(
                    #     dH,
                    #     sym,
                    #     xyz,
                    #     mass,
                    #     0,
                    #     remove=0,
                    #     stdtr=True,
                    #     return_DL=True,
                    #     verbose=False
                    # )

                    dHw = list(hessians.hessian_frequencies(
                        mass, dH, grad_mm, DL, ics, B, B2)
                    )
                    # x1.extend(x1)
                    x1.extend(dHw)
                    # x0.extend([0]*len(dHw))
                    # x0.extend([0]*(len(x0) - len(x1)))

                    # keep this for when I want to try internal hess
                    # x0 = array_numpy.dlcmatrix_project_gradients(pos, x0)
                    # x1 = array_numpy.dlcmatrix_project_gradients(pos, x1)

                dxa.extend(x1)
            f.append(dxa)

        if self.grad_mode == "f1":
            dxdp, d2xdp2 = finite_difference_forward_1(f, h)
        elif self.grad_mode == "c2":
            dxdp, d2xdp2 = finite_difference_central_2(f, h)

        return dxdp, d2xdp2

    # def compute_gradient_2pt(self, ref, X0, E0, E1, h):

    #     dxa = []
    #     dxb = []
    #     gde = ref.entries[self.addr.eid[0]]

    #     for etbls, dtbls in zip(E0, E1):
    #         # for tid, dtbl in dtbls.items():

    #         # x0 = etbls[tid].values
    #         # x1 = dtbl.values

    #         gid = list(gde.tables[POSITIONS].graphs)[0]
    #         g = ref.graphs[gid]
    #         DL = self.cache['IC:dl']
    #         (ics, B, B2) = self.cache['IC:b']

    #         hess0 = dtbls[HESSIANS].values
    #         grad0 = list(dtbls[GRADIENTS].values)
    #         x0 = hessians.hessian_frequencies(g, hess0, grad0, DL, ics, B, B2)

    #         hess_mm = etbls[HESSIANS].values

    #         if self.enable_minimization or self.method_vib_modes == "min":
    #             pos1 = etbls[POSITIONS].values
    #             xyz = list([x[0] for x in pos1.selections.values()])

    #         if False and self.method_vib_modes in ["min", "mm"]:
    #             sym = graphs.graph_symbols(g)
    #             mass = [[vibration.mass_table[sym[n]]]*3 for n in sym]
    #             sym = list(sym.values())

    #             remove = 0
    #             if self.method_vib_modes == "mm":
    #                 remove = 0

    #             x1, mm_modes, mm_dl = vibration.hessian_modes(
    #                 hess_mm,
    #                 sym,
    #                 xyz,
    #                 mass,
    #                 0,
    #                 remove=remove,
    #                 stdtr=True,
    #                 return_DL=True
    #             )
    #             cost = np.abs(np.dot(np.array(DL).T, mm_dl))
    #             # s_row, s_col = scipy.optimize.linear_sum_assignment(1 - cost)
    #             # x1 = [x1[i] for i in s_row]
    #             x1 = list(x1)

    #         elif self.method_vib_modes == "qm":
    #             grad_mm = list(etbls[GRADIENTS].values)
    #             x1 = list(hessians.hessian_frequencies(
    #                 g, hess_mm, grad_mm, DL, ics, B, B2)
    #             )
    #         # elif self.method_vib_modes == "qm+":
    #         #     grad_mm = list(etbls[GRADIENTS].values)
    #         #     # Need Q' which is Q and the extra R modes
    #         #     DL2 = vibration.gram_schmidt(DL2.T).T

    #         #     x1 = list(hessians.hessian_frequencies(
    #         #         g, hess_mm, grad_mm, DL2, ics, B, B2)
    #         #     )
    #         #     x0.extend([0]*(len(x0) - len(x1)))
                

    #         # keep this for when I want to try internal hess
    #         # x0 = array_numpy.dlcmatrix_project_gradients(pos, x0)
    #         # x1 = array_numpy.dlcmatrix_project_gradients(pos, x1)

    #         dxa.extend(x0)
    #         dxb.extend(x1)

    #     dxdp = arrays.array_difference(dxb, dxa)
    #     dxdp = arrays.array_round(dxdp, PRECISION)
    #     dxdp = arrays.array_scale(dxdp, 1/(2*h))

    #     d2xdp2 = arrays.array_add(
    #         arrays.array_scale(X0, -2),
    #         arrays.array_add(dxb, dxa)
    #     )
    #     d2xdp2 = arrays.array_round(d2xdp2, PRECISION)
    #     d2xdp2 = arrays.array_scale(d2xdp2, 1.0/(h*h))

    #     return dxdp, d2xdp2

    def compute_diff(
        self,
        GDB: assignments.graph_db,
        D: List[Dict[int, assignments.graph_db_table]],
        verbose=False
    ):
        """
        """

        out = []
        X = []
        for (eid, gde), tbls in zip(GDB.entries.items(), D):

            gid = list(gde.tables[assignments.POSITIONS].graphs)[0]
            g = GDB.graphs[gid]
            # smiles = GDB.smiles[gid]

            hess_qm = gde.tables[HESSIANS].values
            # hess_qm_sel = gde.tables[POSITIONS][gid][0][0].selections
            # # shim until I get off the fence on the dimensions of sel
            # hess_qm_sel = {k: [v] for k, v in hess_qm_sel.items()}

            # hess_qm_pos = assignments.graph_assignment(smiles, hess_qm_sel, g)

            hess_qm_pos = assignments.graph_db_graph_to_graph_assignment(
                GDB,
                eid,
                POSITIONS,
            )

            # xyz = list([x[0] for x in hess_qm_sel.values()])
            xyz = []
            for pi, posi in enumerate(hess_qm_pos):
                for xyzi in posi.selections.values():
                    xyz.extend(xyzi)

            sym = [s for posi in hess_qm_pos for s in graphs.graph_symbols(posi.graph).values()]
            mass = np.array([[vibration.mass_table[s]]*3 for s in sym])

            if True or self.coordinate_system == "IC":
                if (
                    ("IC:qm_freq" not in self.cache)
                ):

                    torsions = True
                    pairs = False
                    ics, B = assignments.bmatrix(
                        hess_qm_pos,
                        torsions=torsions,
                        pairs=pairs,
                        remove1_3=True
                    )
                    ics, B2 = assignments.b2matrix(
                        hess_qm_pos,
                        torsions=torsions,
                        pairs=pairs,
                        remove1_3=True
                    )

                    hess_qm_freq, hess_qm_modes, DL = vibration.hessian_modes(
                        hess_qm,
                        sym,
                        xyz,
                        mass,
                        0,
                        remove=0,
                        stdtr=True,
                        return_DL=True,
                        verbose=False
                    )
                    keep_modes = [
                        i
                        for i, w in enumerate(hess_qm_freq)
                        if w > self.freq_range[0] and w < self.freq_range[1]
                    ]
                    DL = DL[:, keep_modes]
                    hess_qm_freq = [hess_qm_freq[i] for i in keep_modes]


                    self.cache['IC:qm_freq'] = hess_qm_freq
                    self.cache['IC:b'] = (ics, B, B2)
                    self.cache['IC:dl'] = DL

                x0 = list(self.cache['IC:qm_freq'])
                hess_mm = tbls[HESSIANS].values
                grad_mm = list(tbls[GRADIENTS].values)

                # grad1 = list(tbls[GRADIENTS].values)

                if self.verbose > 2:
                    out.append(
                        f"EID {self.addr.eid[0]} "
                        f"FC reset mode is {self.fc_reset_config}"
                    )

                if self.method_vib_modes in ["min", "mm", "mm+"]:
                    # sym = graphs.graph_symbols(g)
                    # mass = [[vibration.mass_table[sym[n]]]*3 for n in sym]
                    # sym = list(sym.values())

                    remove = 0
                    if self.method_vib_modes == "mm":
                        remove = 0

                    x1, mm_modes, mm_dl = vibration.hessian_modes(
                        hess_mm,
                        sym,
                        xyz,
                        mass,
                        1,
                        remove=remove,
                        stdtr=True,
                        return_DL=True,
                        verbose=False
                    )
                    # if self.method_vib_modes in ["min", "mm"]:
                    #     DL = mm_dl
                    if self.method_vib_modes != "mm+":
                        x1 = list(sorted(x1))
                    else:
                        DL = self.cache['IC:dl']
                        cost = np.abs(np.dot(np.array(DL).T, mm_dl))
                        costheta = np.degrees(np.arccos(np.dot(np.array(DL).T, mm_dl)))
                        costheta[costheta > 90] -= 180
                        # s_row, s_col = scipy.optimize.linear_sum_assignment(cost)
                        # x1 = [x1[i] for i in s_col]
                        x1 = list(x1)
                        # angle = np.arccos(
                        #     np.abs(cost[s_row, s_col]).sum()/len(cost)
                        # )
                        angle = np.arccos(np.diag(cost))
                        if self.verbose > 2:
                            out.append(
                                f"1Cost for EID {self.addr.eid[0]} "
                                f"is total {np.sum(angle)} "
                                f"mean {np.mean(angle)} "
                                f"degrees {np.degrees(np.mean(angle))}"
                            )
                            for i in range(len(cost)):
                                out.append(
                                    f"Unmapped Vib. assignment {i+1} -> {i+1} "
                                    f"QMfreq {x0[i]} "
                                    f"cost {cost[i][i]} "
                                    f"angle {costheta[i][i]} "
                                )
                            # print("\n".join(out))

                        cost = 1 - np.abs(np.dot(np.array(DL).T, mm_dl))
                        costheta = np.degrees(np.arccos(np.dot(np.array(DL).T, mm_dl)))
                        costheta[costheta > 90] -= 180
                        s_row, s_col = scipy.optimize.linear_sum_assignment(cost)
                        x1 = [x1[i] for i in s_col]
                        # x1 = list(x1)
                        angle = np.arccos(
                            np.abs(cost[s_row, s_col]).sum()/len(cost)
                        )
                        if self.verbose > 2:
                            out.append(
                                f"2Cost for EID {self.addr.eid[0]} "
                                f"is total {np.sum(angle)} "
                                f"mean {np.mean(angle)} "
                                f"degrees {np.degrees(np.mean(angle))}"
                            )
                            for i, j in zip(s_row, s_col):
                                out.append(
                                    f"Mapped Vib. assignment {i+1} -> {j+1} "
                                    f"QMfreq {x0[i]} "
                                    f"cost {cost[i][j]} "
                                    f"angle {costheta[i][j]} "
                                )
                            # print("\n".join(out))

                        # x0 = [x0[i] for i in s_row]
                        # d = np.abs(np.array(x1) - np.array(x0))


                elif self.method_vib_modes == "qm":
                    DL = self.cache['IC:dl']
                    (ics, B, B2) = self.cache['IC:b']

                    x1mat = list(hessians.hessian_transform(
                        mass, hess_mm, grad_mm, DL, ics, B, B2)
                    )
                    x1 = [x1mat[i][i] for i in range(len(x1mat))]
                    x1offdiag = [
                        x1mat[i][j] for i in range(len(x1mat)) for j in range(i)
                    ]
                    # x1 = list(hessians.hessian_frequencies(
                    #     g, hess_mm, grad_mm, DL, ics, B, B2)
                    # )
                    # x1 = arrays.array_scale(x1, 2)
                elif self.method_vib_modes == "qm+":
                    DL = self.cache['IC:dl']
                    (ics, B, B2) = self.cache['IC:b']

                    dH = [arrays.array_difference(hess_mm[i], hess_qm[i]) for i in range(len(hess_qm))]
                    # dHw, _, _= vibration.hessian_modes(
                    #     dH,
                    #     sym,
                    #     xyz,
                    #     mass,
                    #     0,
                    #     remove=0,
                    #     stdtr=True,
                    #     return_DL=True,
                    #     verbose=False
                    # )

                    x1mat = list(hessians.hessian_transform(
                        g, hess_mm, grad_mm, DL, ics, B, B2)
                    )
                    x1 = [x1mat[i][i] for i in range(len(x1mat))]
                    # dHw = list(hessians.hessian_frequencies(
                    #     g, dH, grad_mm, DL, ics, B, B2)
                    # )
                    dmat = list(hessians.hessian_transform(
                        g, dH, grad_mm, DL, ics, B, B2)
                    )
                    dHw = [dmat[i][i] for i in range(len(dmat))]

                    x1.extend(dHw)
                    # qw = list(hessians.hessian_frequencies(
                    #     g, hess_qm, grad_mm, DL, ics, B, B2)
                    # )
                    # x1.extend(x1)
                    x0.extend([0]*len(dHw))
                    # x0.extend(x0)

                    # x1mat = list(hessians.hessian_transform(
                    #     g, hess_mm, grad_mm, DL, ics, B, B2)
                    # )
                    # x1 = [x1mat[i][i] for i in range(len(x1mat))]
                    # x1offdiag = [
                    #     x1mat[i][j] for i in range(len(x1mat)) for j in range(i)
                    # ]
                    x1offdiag = [0]
                    # x1 = list(hessians.hessian_frequencies(
                    #     g, hess_mm, grad_mm, DL, ics, B, B2)
                    # )
                    # x1 = arrays.array_scale(x1, 2)
                    # x0.extend([0]*(len(x0) - len(x1)))

            x = arrays.array_difference(x1, x0)
            N = len(x)

            mae = [abs(xi) for xi in x]
            t = "Numerical"
            if self.analytic:
                t = "Analytic"
            if self.verbose > 0:
                out.append(f"\nEID {eid} Computed {t} Hessian Frequencies (cm-1):")
            if self.verbose > 1:
                for i, (f1, f0) in enumerate(zip(x1, x0), 1):
                    out.append(f"      {i:4d} MM: {f1: 8.1f} QM: {f0:8.1f} Diff: {f1-f0:8.1f}")
                out.append("           ----------------------------------------")
            if self.verbose > 0:
                xx0 = [x*x for x in x0]
                xx1 = [x*x for x in x1]
                sse = [x*x for x in mae]
                out.append(f"      SAE  MM: {sum(x1): 8.1f} QM: {sum(x0):8.1f} Diff: {sum(mae):8.1f}")
                out.append(f"      MAE  MM: {sum(x1)/N: 8.1f} QM: {sum(x0)/N:8.1f} Diff: {sum(mae)/N:8.1f}")
                out.append(f"      SSE  MM: {sum(xx1): 8.1e} QM: {sum(xx0):8.1e} Diff: {sum(sse):8.1e}")
                out.append(f"     RMSE  MM: {(sum(xx1)/N)**.5: 8.1f} QM: {(sum(xx0)/N)**.5:8.1f} Diff: {(sum(sse)/N)**.5:8.1f}")
                M = 6
                if N == 6:
                    M = 5
                out.append(f"Int.  SAE  MM: {sum(x1[M:]): 8.1f} QM: {sum(x0):8.1f} Diff: {sum(mae):8.1f}")
                out.append(f"Int.  MAE  MM: {sum(x1[M:])/(N-M): 8.1f} QM: {sum(x0[M:])/(N-M):8.1f} Diff: {sum(mae[M:])/(N-M):8.1f}")
                out.append(f"Int.  SSE  MM: {sum(xx1[M:])/(N-M): 8.1e} QM: {sum(xx0[M:])/(N-M):8.1e} Diff: {sum(sse[M:]):8.1e} N: {N-M}")
                out.append(f"Int. RMSE  MM: {(sum(xx1[M:])/(N-M))**.5: 8.1f} QM: {(sum(xx0[M:])/(N-M))**.5:8.1f} Diff: {(sum(sse[M:])/(N-M))**.5:8.1f}")
            # if verbose:
            #     for i in range(len(out) - nlines):
            #         print(out[i])
            X.extend(x)

            if self.method_vib_modes == "qm":
                if verbose > 0:
                    out.append(f"Off-diag MAE: {sum(map(abs, x1offdiag)): 8.1f}")
                # X.extend(x1offdiag)

        return returns.success(arrays.array_round(X, PRECISION), out=out, err=[])


class objective_config_energy_sapt(objective_config):

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.dimer_start = 0

    def get_task(
        self,
        GDB: assignments.graph_db,
        csys,
        keys=None,
        psys=None,
        reuse=None
    ) -> compute_config:
        cc = compute_config_energy_sapt(self.addr)
        if self.gdb is None:
            self.gdb = assignments.graph_db_get(GDB, self.addr)
        cc.GDB = self.gdb

        cc.csys = csys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.psys = psys
        cc.reuse = reuse
        cc.dimer_start = self.dimer_start
        cc.fc_reset_config = self.fc_reset_config
        return cc

    def compute_gradient(
        self,
        ref,
        evals,
        h
    ):
        """
        Compute the parameter gradients of this objective.
        """

        dxa = []
        f = []
        m = 5
        for E0 in evals:
            dxa = []
            for dtbls in E0:
                tid = assignments.SAPT
                dtbl = dtbls[tid]

                x0 = dtbl.values
                el = dtbl.values[0][0][4]
                lj = dtbl.values[0][0][5]
                # dimer_interaction = {ic: v for ic, v in LJ.items() if ic[0][1] <= n and ic[1][1] > n}
                # ene = sum([x for y in dimer_interaction.values() for x in y])
                # dxa.extend([lj+el])
                dxa.extend([lj])
            f.append(dxa)

        if self.grad_mode == "f1":
            dxdp, d2xdp2 = finite_difference_forward_1(f, h)
        elif self.grad_mode == "c2":
            dxdp, d2xdp2 = finite_difference_central_2(f, h)

        return dxdp, d2xdp2


    def compute_diff(
        self,
        GDB: assignments.graph_db,
        D: List[Dict[int, assignments.graph_db_table]],
        verbose=False
    ):
        """
        Compute the difference of the computed properties versus a reference
        """

        obj = []
        out = []
        n = self.dimer_start
        m = 5
        for (eid, gde), tbls in zip(GDB.entries.items(), D):
            rmse = []
            if self.verbose > 1:
                out.append(f"EID {eid} Energy objective:")
            tid = assignments.SAPT
            tbl = tbls[tid]
            tbl0 = gde.tables[tid]

            ene_dict = tbl.values[0][0]
            tot_dict = tbl.values[0][1]

            # for m, enes in :
            #     LJ = tbl.values[0][m]
            #     dimer_interaction = {ic: v for ic, v in enes.items() if ic[0][1] <= n and ic[1][1] > n}
            #     ene_dict[m] = sum([x for y in dimer_interaction.values() for x in y])

            el_qm = tbl0.values[0][0]
            lj_qm = tbl0.values[0][1]
            el_mm = ene_dict[4]
            lj_mm = ene_dict[5]
            obj.append(el_mm - el_qm)
            obj.append(lj_mm - lj_qm)
            # for tid, tbl in tbls.items():
            #     tbl0 = gde.tables[tid]

                # for gid, gdg in tbl.graphs.items():
                #     gdg0 = tbl0.graphs[gid]
                #     for rid, gdr in gdg.rows.items():
                #         gdr0 = gdg0.rows[rid]
                #         for cid, gdc in gdr.columns.items():
                #             gdc0 = gdr0.columns[cid]

                #             # assume cartesian
                #             x1l = [xyz for xyz in gdc.selections.values()]
                #             x0l = [xyz for xyz in gdc0.selections.values()]

                #             for ic, x1, x0 in zip(gdc0.selections, x1l, x0l):
                #                 x = arrays.array_difference(x1, x0)
                #                 x = arrays.array_round(x, PRECISION)
                #                 obj.extend(x)
                #                 rmse.append(arrays.array_inner_product(x, x))
                #                 if self.verbose > 1:
                #                     xyz1 = " ".join([f"{xi:8.4f}" for xi in x1])
                #                     xyz0 = " ".join([f"{xi:8.4f}" for xi in x0])
                #                     dxyz = " ".join([f"{xi:8.4f}" for xi in  x])
                #                     addr = f"    {eid:4d} {gid:2d} {rid:2d} "
                #                     addr += f"{cid:2d} {ic[0]:3d}"
                #                     out_str  = f" MM: {xyz1} QM: {xyz0} dP: {dxyz}"
                #                     out.append(addr + out_str)
            # if len(rmse):
            #     rmse = (sum(rmse)/len(rmse))**.5
            # else:
            #     rmse = 0
            # X2 = arrays.array_inner_product(obj, obj)
            el = obj[0]**2
            lj = obj[1]**2
            tot = lj_mm + el_mm - (lj_qm + el_qm)
            tot2 = tot**2
            if self.verbose > 0:
                out.extend([" ".join([
                    f"Model {m} Energy: {e:10.5f} (kJ/mol)"
                ]) for m, e in tot_dict.items()])
                out.append(" ".join([
                    f"Dimer Electrostatics  QM: {el_qm:10.5f} MM: {el_mm:10.5f} SSE: {el:10.5f} (kJ/mol)^2 N: 1",
                    f"RMSE: {el**.5:10.5f} kJ/mol "
                    f"Delta: {obj[0]:10.5f} kJ/mol"
                ]))
                out.append(" ".join([
                    f"Dimer Exchange+Induction+Dispersion+dHF QM: {lj_qm:10.5f} MM: {lj_mm:10.5f} SSE: {lj:10.5f} (kJ/mol)^2 N: 1",
                    f"RMSE: {lj**.5:10.5f} kJ/mol "
                    f"Delta: {obj[1]:10.5f} kJ/mol"
                ]))
                out.append(" ".join([
                    f"Dimer Total QM: {lj_qm+el_qm:10.5f} MM: {lj_mm+el_mm:10.5f} SSE: {tot2:10.5f} (kJ/mol)^2 N: 1",
                    f"RMSE: {tot2**.5:10.5f} kJ/mol "
                    f"Delta: {tot:10.5f} kJ/mol"
                ]))

        return returns.success(arrays.array_round([obj[1]], PRECISION), out=out, err=[])
        # return returns.success(arrays.array_round([tot], PRECISION), out=out, err=[])

class objective_config_energy_total(objective_config):

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.ene_mode = "single"

    def get_task(
        self,
        GDB: assignments.graph_db,
        csys,
        keys=None,
        psys=None,
        reuse=None
    ) -> compute_config:
        cc = compute_config_energy_total(self.addr)
        if self.gdb is None:
            self.gdb = assignments.graph_db_get(GDB, self.addr)
        cc.GDB = self.gdb
        cc.csys = csys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.psys = psys
        cc.reuse = reuse
        cc.fc_reset_config = self.fc_reset_config
        return cc

    def compute_gradient(
        self,
        ref,
        evals,
        h
    ):
        """
        Compute the parameter gradients of this objective.
        """

        TID = assignments.ENERGY

        f = []
        for E0 in evals:
            dxa = []
            # one for each eid
            for gde, dtbls in zip(ref.entries.values(), E0):
                dtbl = dtbls[TID]
                # for tid, dtbl in dtbls.items():
                x0 = dtbl.values
                if self.ene_mode == "single":
                    dxa.extend(x0)
                elif self.ene_mode == "pairs":
                    dxa.extend(x0*len(x0))
            f.append(dxa)

        if self.grad_mode == "f1":
            dxdp, d2xdp2 = finite_difference_forward_1(f, h)
        elif self.grad_mode == "c2":
            dxdp, d2xdp2 = finite_difference_central_2(f, h)
            
        return dxdp, d2xdp2


    def compute_diff(
        self,
        GDB: assignments.graph_db,
        D: List[Dict[int, assignments.graph_db_table]],
        verbose=False
    ):
        """
        Compute the difference of the computed properties versus a reference
        """

        obj = []
        out = []
        tid = assignments.ENERGY

        min_idx = None
        min_ref_ene = None
        min_ene = None

        X1 = []
        X0 = []
        for (eid, gde), tbls in zip(GDB.entries.items(), D):

            tbl0 = gde.tables[tid]
            x0 = tbl0.values

            tbl = tbls[tid]
            x1 = tbl.values

            min_idx = arrays.argmin(x0)
            if min_ref_ene is None:
                min_ref_ene = min(x0)
            if min_ene is None or x0[min_idx] < min_ref_ene:
                min_ene = x1[min_idx]

            X0.extend(x0)
            X1.extend(x1)

        X1 = arrays.array_translate(X1, -min_ene)

        mode = self.ene_mode
        weight = [1 for x in range(len(X1))]
        if self.weights is not None:
            weight = self.weights
        if self.verbose > 0:
            out.extend(
                    f"EID {m:4d} {mode} W={w:10.5f} MM Energy: {x1:10.5f} (kJ/mol) QM Energy {x0:10.5f} (kJ/mol) Delta {x1-x0:10.5f} (kJ/mol)"
             for m, (x1, x0, w) in enumerate(zip(X1, X0, weight)))
        if self.ene_mode == "single":
            obj = arrays.array_difference(X1, X0)
            if self.verbose > 0:
                X2 = arrays.array_inner_product(obj, obj)
                out.append(" ".join([
                    f"Energy SSE: {X2:10.5f} (kJ/mol)^2 N: {len(x1)}",
                    f"RMSE: {(X2/len(x1))**.5:10.5f} kJ/mol "
                ]))
        elif self.ene_mode == "pairs":
            for xi, x0 in enumerate(X0):
                obji = arrays.array_translate(X1, -x0)
                obj.extend(obji)
                if self.verbose > 0:
                    X2 = arrays.array_inner_product(obji, obji)
                    out.append(" ".join([
                        f"Energy from QM {xi:4d} SSE: {X2:10.5f} (kJ/mol)^2 N: {len(x1)}",
                        f"RMSE: {(X2/len(x1))**.5:10.5f} kJ/mol "
                    ]))
            if self.verbose > 0:
                X2 = arrays.array_inner_product(obj, obj)
                out.append(" ".join([
                    f"Total Pairwise Energy SSE: {X2:10.5f} (kJ/mol)^2 N: {len(obj)}",
                    f"RMSE: {(X2/len(obj))**.5:10.5f} kJ/mol "
                ]))


        # for (eid, gde), tbls in zip(GDB.entries.items(), D):
        #     rmse = []
        #     if self.verbose > 1:
        #         out.append(f"EID {eid} Energy objective:")
        #     tbl0 = gde.tables[tid]

        #     x0 = tbl0.values

        #     tbl = tbls[tid]
        #     x1 = tbl.values
        #     x1 = arrays.array_translate(x1, -min_ene)
            
        #     obj.extend(arrays.array_difference(x1, x0))
            # ene_dict = tbl.values[0][0]
            # tot_dict = tbl.values[0][1]

            # for m, enes in :
            #     LJ = tbl.values[0][m]
            #     dimer_interaction = {ic: v for ic, v in enes.items() if ic[0][1] <= n and ic[1][1] > n}
            #     ene_dict[m] = sum([x for y in dimer_interaction.values() for x in y])

            # el_qm = tbl0.values[0][0]
            # lj_qm = tbl0.values[0][1]
            # el_mm = ene_dict[4]
            # lj_mm = ene_dict[5]
            # obj.append(el_mm - el_qm)
            # obj.append(lj_mm - lj_qm)
            # for tid, tbl in tbls.items():
            #     tbl0 = gde.tables[tid]

                # for gid, gdg in tbl.graphs.items():
                #     gdg0 = tbl0.graphs[gid]
                #     for rid, gdr in gdg.rows.items():
                #         gdr0 = gdg0.rows[rid]
                #         for cid, gdc in gdr.columns.items():
                #             gdc0 = gdr0.columns[cid]

                #             # assume cartesian
                #             x1l = [xyz for xyz in gdc.selections.values()]
                #             x0l = [xyz for xyz in gdc0.selections.values()]

                #             for ic, x1, x0 in zip(gdc0.selections, x1l, x0l):
                #                 x = arrays.array_difference(x1, x0)
                #                 x = arrays.array_round(x, PRECISION)
                #                 obj.extend(x)
                #                 rmse.append(arrays.array_inner_product(x, x))
                #                 if self.verbose > 1:
                #                     xyz1 = " ".join([f"{xi:8.4f}" for xi in x1])
                #                     xyz0 = " ".join([f"{xi:8.4f}" for xi in x0])
                #                     dxyz = " ".join([f"{xi:8.4f}" for xi in  x])
                #                     addr = f"    {eid:4d} {gid:2d} {rid:2d} "
                #                     addr += f"{cid:2d} {ic[0]:3d}"
                #                     out_str  = f" MM: {xyz1} QM: {xyz0} dP: {dxyz}"
                #                     out.append(addr + out_str)
            # if len(rmse):
            #     rmse = (sum(rmse)/len(rmse))**.5
            # else:
            #     rmse = 0
            # X2 = arrays.array_inner_product(obj, obj)
            # el = obj[0]**2
            # lj = obj[1]**2
            # tot = lj_mm + el_mm - (lj_qm + el_qm)
            # tot2 = tot**2
            # if self.verbose > 0:
            #     out.extend([" ".join([
            #         f"Model {m} Energy: {e:10.5f} (kJ/mol)"
            #     ]) for m, e in tot_dict.items()])
            #     out.append(" ".join([
            #         f"Dimer Electrostatics  QM: {el_qm:10.5f} MM: {el_mm:10.5f} SSE: {el:10.5f} (kJ/mol)^2 N: 1",
            #         f"RMSE: {el**.5:10.5f} kJ/mol "
            #         f"Delta: {obj[0]:10.5f} kJ/mol"
            #     ]))
            #     out.append(" ".join([
            #         f"Dimer Exchange+Induction+Dispersion+dHF QM: {lj_qm:10.5f} MM: {lj_mm:10.5f} SSE: {lj:10.5f} (kJ/mol)^2 N: 1",
            #         f"RMSE: {lj**.5:10.5f} kJ/mol "
            #         f"Delta: {obj[1]:10.5f} kJ/mol"
            #     ]))
            #     out.append(" ".join([
            #         f"Dimer Total QM: {lj_qm+el_qm:10.5f} MM: {lj_mm+el_mm:10.5f} SSE: {tot2:10.5f} (kJ/mol)^2 N: 1",
            #         f"RMSE: {tot2**.5:10.5f} kJ/mol "
            #         f"Delta: {tot:10.5f} kJ/mol"
            #     ]))

        return returns.success(arrays.array_round(obj, PRECISION), out=out, err=[])
        # return returns.success(arrays.array_round([tot], PRECISION), out=out, err=[])


class objective_config_position(objective_config):

    def get_task(
        self,
        GDB,
        csys,
        keys=None,
        psys=None,
        reuse=None
    ) -> compute_config:
        cc = compute_config_position(self.addr)
        if self.gdb is None:
            self.gdb = assignments.graph_db_get(GDB, self.addr)
        cc.GDB = self.gdb
        cc.csys = csys
        if keys:
            cc.keys.extend(keys)
        else:
            cc.keys.append({})
        cc.psys = psys
        cc.reuse = reuse
        cc.step_limit = self.step_limit
        cc.tol = self.tol
        cc.fc_reset_config = self.fc_reset_config
        return cc

class objective_config_penalty(objective_config):

    def __init__(self, keys=None, scale=1, polynomial={2: 1.0}):
        if keys is None:
            keys = {}

        self.keys = keys
        self.reference = {}
        self.scale = scale
        assert all([type(n) is int for n in polynomial])
        self.polynomial = polynomial
        self.include = True
        self.reference_from_reset = False

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

        return round(DX, PRECISION)

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

        return round(DX, PRECISION)

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

        return returns.success(round(DX, PRECISION), out=[], err=[])


class objective_tier:
    """
    A configuration for how to compute a total objective over a list of
    objectives. A tier is meant to be calculated for 1 or more force fields,
    which is what happens during parameter searching. Tiers have a fitting
    aspect that is configured here. For example, an objective might only
    perform 1 or 2 fitting steps for each force field, and the top N force
    fields would be accepted (and passed to the next tier.
    """

    def __init__(self):
        self.objectives: Dict[int, objective_config] = {}
        h = 1e-6
        self.h = h
        self.bounds = {
            "s": (0, None),
            "c": (0, None),
            "r": (0, None),
            "e": (h, None),
            "k": (0, None),
            ("k", "b"): (0, None),
            ("l", "b"): (0, None),
            ("k", "a"): (0, None),
            ("l", "a"): (0, 3.1415),
            ("k", "t"): (-10, 10),
            ("k", "i"): (-10, 10),
            ("s", "q"): (0, None),
            ("s", "n"): (0, None),
            None: (None, None)
        }
        self.priors = {
            "s": (None, None),
            "c": (None, None),
            "r": (None, .1),
            "e": (None, 1),
            "k": (None, 1e9),
            ("k", "b"): (None, 1e9),
            ("l", "b"): (None, .1),
            ("k", "a"): (None, 1e9),
            ("l", "a"): (None, .2),
            ("k", "t"): (None, 1.0),
            ("k", "i"): (None, 1.0),
            ("s", "q"): (None, 1.0),
            ("s", "n"): (None, 1.0),
            None: (None, None)
        }
        self.fit_models = None
        self.fit_symbols = None
        self.fit_names = None
        self.fit_names_exclude = None
        self.method = "trust-ncg"
        self.method_config = {}
        self.step_limit = None
        self.maxls = 20
        self.penalties: List[objective_config_penalty] = []
        self.accept: int = 0 # 0 keeps all (essentially disables)
        self.anneal = False
        self.minstep = 10**(-PRECISION)
        self.ftol = 1e-5
        self.gtol = 1e-5

    def key_filter(self, x):
        if self.fit_models is None or self.fit_symbols is None:
            return True

        if type(self.fit_symbols) is dict:
            r = self.fit_symbols.get(x[0]) == x[1]
        else:
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
    """
    Build a physical system for each entry in the graph_db
    """

    psysref = {}
    failures = 0
    if len(gdb.entries) > 100:
        print(f"{datetime.datetime.now()} Starting parameterization")
    for i, (eid, gde) in enumerate(gdb.entries.items(), 1):
        if len(gdb.entries) > 100:
            print(f"\r{datetime.datetime.now()} Parameterizing.. {i:6d}/{len(gdb.entries)}", end="")
        tid = assignments.POSITIONS
        # gid = list(gde.tables[tid].graphs)
        # rid = 0
        pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, tid, gid=None, rid=None)
        try:
            if type(ref) is dict:
                ref0 = ref[eid]
            else:
                ref0 = ref
            psysref[eid] = mm.chemical_system_to_physical_system(csys, pos, ref=ref0, reuse=reuse)
        except Exception as e:
            print(f"\n\nWarning: could not parameterize {eid}, skipping.")
            print(f"SMILES:")
            for pi in pos:
                print(f"{pi.smiles}")
            print("Exception:")
            print(e)
            failures += 1
    if len(gdb.entries) > 100:
        print()
    if failures:
        print(f"There were {failures} failures that will be skipped")
    return psysref


def objective_tier_run(
    ot,
    gdb: assignments.graph_db,
    csys: mm.chemical_system,
    keys,
    oid=None,
    psysref=None,
    reuse=None,
    wq=None,
    verbose=False
):
    """
    Run the given objective and return the resulting parameters, and the
    initial and final objective and gradients.

    Parameters
    ----------
    ot: The objective_tier object
    gdb: The complete graph_db
    csys: The initial chemical system
    keys: The keys in the chemical_system that are going to be fit
    oid: The objective IDs to use in the computation
    psysref: The reference physical systems to use when applied parameters
        should be reused (charges etc)
    reuse: Models that be reused during parameterization
    wq: The wq to use for distributing computing
    verbose: Whether to make a lot of noise
    """

    # need this to avoid race conditions
    csys = copy.deepcopy(csys)
    # build the dataset and input ff
    ws = None
    addr = False
    shm = {"csys": csys}
    if wq and configs.remote_compute_enable:
        shm["compute_verbosity"] = 2
        ws = compute.workqueue_new_workspace(wq, shm=shm)
    elif configs.processors > 1:
        addr = ('127.0.0.1', 0)
        ws = compute.workspace_local(*addr, shm=shm)

    # build the initial psys (for charges)
    # if verbose:
    #     print(datetime.datetime.now(), "Building physical systems")
    if psysref is None:
        psysref = gdb_to_physical_systems(gdb, csys)

    if oid:
        ot.objectives = {i: ot.objectives[i] for i in oid}

    x0 = [mm.chemical_system_get_value(csys, k) for k in keys]

    # we need history to identify work for each iteration
    # otherwise we might accidentally use work from step n-1 that gets sent in
    # when we are already at step n
    history = []

    kv = dict(zip(keys, x0))
    for penalty in ot.penalties:
        penalty.reference.clear()
        penalty.reference.update(kv)

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

    if not bounds:
        bounds = None

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

    minstep = ot.minstep

    args = (
        keys,
        csys,
        gdb,
        ot,
        priors,
        ot.penalties,
        history,
        psysref,
        reuse,
        ws,
        verbose,
        minstep
    )

    ret = optimizers_scipy.optimize_forcefield_gdb_scipy(
        x0,
        args,
        bounds=bounds,
        step_limit=ot.step_limit,
        maxls=ot.maxls,
        ftol=ot.ftol,
        gtol=ot.gtol,
        anneal=ot.anneal
    )

    result, y0, y1, gx = ret.value

    ret.out.append(f">>> Initial Objective {y0:10.5g}")
    ret.out.append(f">>> Final Objective   {y1:10.5g}")
    if y0 > 0.0:
        ret.out.append(f">>> Percent change    {(y1-y0)/y0*100:10.5g}%")

    if verbose:
        for i in [-3, -2, -1]:
            print(ret.out[i])

    if ws:
        if wq:
            compute.workqueue_remove_workspace(wq, ws)
        ws.close()
        ws = None

    kv = {k: v*p[1] + p[0] for k, v, p in zip(keys, result, priors)}

    return returns.success((kv, y0, y1, gx), out=ret.out, err=[])


def objective_run_distributed(obj, shm=None):
    if shm and "csys" in shm.__dict__:
        obj.csys = shm.csys
    return obj.run()


def forcefield_optimization_strategy_default(csys, models=None):
    """
    Configure a force field optimiation strategy with default settings. The
    strategy determines how to step the optimization forward, choosing which
    hyperparameter to try next in the search, such as bit depth and which
    model to focus on.
    """

    if models is None:
        models = {}
    elif type(models) is not dict:
        models = {i: None for i in models}

    bounds = {}
    strat = forcefield_optimization_strategy(bounds)

    strat.primitives = list(csys.perception.gcd.atom_primitives + csys.perception.gcd.bond_primitives)

    for m, cm in enumerate(csys.models):
        if m not in models:
            continue

        hiers = mm.chemical_model_iter_smarts_hierarchies(cm)

        for h, hidx in enumerate(hiers):

            nodes = [
                n.name for n in hidx.index.nodes.values()
                if n.type == "parameter"
            ]

            strat.reference_list.extend(nodes)

            assert type(models[m]) is not str
            if models[m] is not None:
                strat.target_list.extend(models[m])
            else:
                strat.target_list.extend(nodes)

            splitter = configs.smarts_splitter_config(
                1, 2, 0, 0, 0, 0, True, True, 0, True, True, True, True
            )
            extender = configs.smarts_extender_config(
                0, 0, True
            )
            bounds[m] = configs.smarts_perception_config(splitter, extender)

            # Currently, all SMARTS hierarchies are first so we cheat
            break

    strat.bounds.update(bounds)
    return strat


class forcefield_optimization_strategy(optimization.optimization_strategy):
    """
    Configuration for fitting a force field using BESMARTS (this package).
    """

    def __init__(self, bounds):

        self.bounds = bounds

        self.objective_accept_total = [0]

        # Only consider the top N clusters of each objective state
        self.objective_accept_clusters = [0]

        # Update objective on each evaluation. Some objectives change if new
        # clusters are added. This option determines whether accepting causes a
        # refresh
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
        self.modify_outofplane_frequency_limit = 3
        self.modify_torsion_frequency_limit = 3

        self.steps: List[optimization.optimization_iteration] = []
        self.tree_iterator: Callable = mm.chemical_system_iter_smarts_hierarchies_nodes

        self.step_tracker = {}

        # Number of operations to accept per macro step
        # Relabeling is done here
        # self.accept_max = 1 will give best performance, but will cost the
        # most self.accept_max = 0 is no max
        self.macro_accept_max_total: int = 1

        # Number of operations to accept per micro step
        # We do not relabel here but instead just keep this many
        # self.accept_max = 1 will give best performance, but will cost the
        # most self.accept_max = 0 is no max
        self.micro_accept_max_total: int = 1

        # Number of operations to accept per step per cluster
        # self.accept_max = 1 will give best performance, but will cost the
        # most self.accept_max = 0 is no max
        self.macro_accept_max_per_cluster: int = 1

        # Number of operations to accept per step per cluster
        # self.accept_max = 1 will give best performance, but will cost the
        # most self.accept_max = 0 is no max
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

        self.bond_k_skip = []
        self.bond_l_skip = []

        self.angle_k_skip = []
        self.angle_l_skip = []

        self.torsion_k_skip = []
        self.torsion_n_skip = []

        self.outofplane_k_skip = []
        self.outofplane_n_skip = []

        # The reference list to use for merge protect and target splitting.
        # Because parameters can be added later, we only try to protect
        # the original (reference) set.
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
        return optimization.optimization_strategy_iteration_next(
            self,
            clusters
        )

    def build_steps(self):
        if self.steps is None:
            self.steps = []
        self.steps.extend(
            forcefield_optimization_strategy_build_macro_iterations(self)
        )


def forcefield_optimization_strategy_build_macro_iterations(
    strat: forcefield_optimization_strategy
):
    """
    Build the macro iterations of a strategy, which define each step in the
    optimization. Each macro step will be a list of micro steps, where each
    micro step represents e.g. a parameter in the hierarchy.
    """

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
                bits += 1

                search_cursor += 1
                if search_cursor < strat.cursor:
                    continue

                steps = []
                if strat.enable_split:
                    for m, mbounds in strat.bounds.items():

                        mbounds: configs.smarts_perception_config

                        if mbounds.splitter.bit_search_limit + 1 < bits:
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
                            mbounds.splitter.primitives
                        )
                        extender = configs.smarts_extender_config(
                            branches, branches, mbounds.extender.include_hydrogen
                        )
                        config = configs.smarts_perception_config(
                            splitter, extender
                        )
                        s.pcp = config

                        steps.append(s)

                if strat.enable_merge:
                    for m, mbounds in strat.bounds.items():

                        mbounds: configs.smarts_perception_config

                        # go through the model bounds and add if current
                        # settings are a subset of bounds
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
                            mbounds.splitter.primitives
                        )
                        extender = configs.smarts_extender_config(0, 0, True)
                        config = configs.smarts_perception_config(
                            splitter, extender
                        )
                        s.pcp = config
                        steps.append(s)

                if strat.enable_modify:
                    for m, mbounds in strat.bounds.items():

                        mbounds: configs.smarts_perception_config
                        if m not in [2, 3]:
                            continue

                        if mbounds.splitter.bit_search_limit + 1 < bits:
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
                            mbounds.splitter.primitives
                        )
                        extender = configs.smarts_extender_config(0, 0, True)
                        config = configs.smarts_perception_config(
                            splitter, extender
                        )
                        s.pcp = config
                        steps.append(s)

                macro_step = optimization.optimization_iteration(steps)
                macro_iters.append(macro_step)


    strat.cursor = 0

    return macro_iters


def fit(csys, gdb, objective, psystems, nodes, wq=None, verbose=False):
    """
    Helper function to easily run a force field fit

    Parameters
    ----------
    csys: The reference chemical system to fit
    gdb: The graph_db
    objective: The objective tier to run
    psystems: The parameterized systems
    nodes: unused
    wq: The workqueue_local to use for distributing computing
    """

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
    reuse0 = [k for k, _ in enumerate(csys.models) if k not in fitting_models]

    print("REUSE IS", reuse0)
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
    """
    The toplevel function to run a full BESMARTS force field fit.

    Parameters
    ----------
    csys0: The initial chemical system to fit
    gdb: The graph_db with all data
    psystems: The initial parameterized systems
    strategy: The fitting strategy
    chemical_objective: The function that will compute the chemical objective
    initial_objective: The objective_tier that will perform the initial fit
    tiers: The tiers that will score each parameter candidate
    final_objective: The objective_tier that will score the remaining
    candidates passed by the tiers
    """

    started = datetime.datetime.now()
    max_line = 0

    if not strategy.steps:
        print("Optimization strategy is building steps...")
        strategy.build_steps()

    print(" ".join([
        f"{datetime.datetime.now()}",
        "The optimization strategy has the following iterations:"
    ]))

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
    icd = codecs.intvec_codec(
        gcd.primitive_codecs,
        gcd.atom_primitives,
        gcd.bond_primitives
    )

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
        "dihedral_max_k": strategy.dihedral_periodicity_reset_max_k,
        "bond_l_skip": strategy.bond_l_skip,
        "bond_k_skip": strategy.bond_k_skip,
        "angle_l_skip": strategy.angle_l_skip,
        "angle_k_skip": strategy.angle_k_skip,
        "torsion_n_skip": strategy.torsion_n_skip,
        "torsion_k_skip": strategy.torsion_k_skip,
        "outofplane_n_skip": strategy.outofplane_n_skip,
        "outofplane_k_skip": strategy.outofplane_k_skip,

    }
    wq = compute.workqueue_local('0.0.0.0', configs.workqueue_port)

    ws = compute.workqueue_new_workspace(wq)
    ret = reset(
        reset_config,
        csys0,
        gdb,
        psystems,
        verbose=True,
        ws=ws
    )
    psystems = ret.value
    print("\n".join(ret.out))
    compute.workqueue_remove_workspace(wq, ws)
    ws.close()
    ws = None

    csys = copy.deepcopy(csys0)

    G0 = {i: icd.graph_encode(g) for i, g in gdb.graphs.items()}

    n_ics = 1

    repeat = set()
    visited = set()
    iteration = 0
    N_str = "{:" + str(len(str(n_ics))) + "d}"
    success = False

    keys = mm.chemical_system_iter_keys(csys)
    keys = [k for k in keys if initial_objective.key_filter(k)]

    assigned_nodes = sorted(set([
        (m, l)
        for psys in psystems.values()
        for m, pm in enumerate(psys.models)
        for proc in pm.labels
        for glbl in proc.values()
        for t, l in glbl.items()
    ]))

    print("Initial parameter assignments of dataset:")
    mm.chemical_system_print(csys, show_parameters=[x[1] for x in assigned_nodes])

    assigned_nodes = [
        x for x in assigned_nodes
        if x[0] in strategy.bounds and x[1] not in "sc"
    ]
    assigned_nodes = set([x[1] for x in assigned_nodes])

    print("### BESMARTS chemical perception on the following assignments ###")
    mm.chemical_system_print(csys, show_parameters=assigned_nodes)
    print("#################################################################")

    reuse0 = set(range(len(csys.models))).difference((k[0] for k in keys))
    print(f"{datetime.datetime.now()} Will be caching models {reuse0}")

    kv0 = mm.chemical_system_iter_keys(csys)
    print(f"{datetime.datetime.now()} Computing physical objective")

    fitwq = wq
    if configs.processors == 1:
        fitwq = None

    restarting = 0
    if strategy.cursor > 0 or strategy.step_tracker:
        restarting = initial_objective.step_limit
        initial_objective.step_limit = 0

    ret = fit(
        csys,
        gdb,
        initial_objective,
        psystems,
        assigned_nodes,
        wq=fitwq,
        verbose=True
    )
    kv, P00, P0, gp0 = ret.value

    if strategy.cursor > 0 or strategy.step_tracker:
        initial_objective.step_limit = restarting

    print(f"{datetime.datetime.now()} Computing chemical objective")

    CX0 = mm.chemical_system_smarts_complexity(csys)

    C00 = chemical_objective(csys, P0=math.log(len(psystems)+1), c=CX0)
    C0 = C00
    C = C00
    print(f"{datetime.datetime.now()} C0={C0}")
    X0 = P0 + C0
    P = P0
    print(" ".join([
        f"{datetime.datetime.now()}",
        f"Initial objective: X={P0+C0:13.6g} P={P0:13.6g} C={C0:13.6g}"
    ]))
    for k, v in kv.items():
        v0 = kv0[k]
        mm.chemical_system_set_value(csys, k, v)
        print(f"{str(k):20s} | New: {v:12.6g} Ref {v0:12.6g} Diff {v-v0:12.6g}")

    smirnoff_models.smirnoff_write_version_0p3(csys, "csys.initial.offxml")

    for psys in psystems.values():
        reapply = set()

        for k, v in kv.items():
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

    union_cache = {}
    step_tracker = strategy.step_tracker

    while strategy.steps:
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
            and any(strategy.cursor >= y for y in step_tracker.get((x.category, x.name), {
                    strategy.SPLIT: -1,
                    strategy.MERGE: -1,
                    strategy.MODIFY: -1
                }.values()))
            and ((strategy.reference_list and x.name not in strategy.reference_list) or
                (strategy.reference_list and strategy.target_list and x.name in strategy.target_list))
            # and x.type == "parameter"
        ]

        for n in nodes:
            tkey = n.category, n.name
            if tkey not in step_tracker:
                step_tracker[tkey] = {
                    strategy.SPLIT: 0,
                    strategy.MERGE: 0,
                    strategy.MODIFY: 0
                }
        # remove any that are not in the models

        print(f"Targets for this macro step {strategy.cursor+1}:")
        for nidx, n in enumerate(nodes, 1):
            print(nidx, n.category, n.name)
        print(f"N Targets: {len(nodes)}")

        if len(nodes) == 0:
            print("Warning, no targets returned. Skipping parameter search.")
            break

        print(f"Step tracker for current macro step {strategy.cursor+1}")
        for n, x in step_tracker.items():
            print(n, {k: v+1 for k, v in x.items()})
        print()

        macro: optimization.optimization_iteration = strategy.macro_iteration(
            nodes
        )

        candidates = {}

        n_added = 0
        n_macro = len(strategy.steps)

        mdls = set([
            m
            for mstep in strategy.steps[strategy.cursor-1].steps
            for m in mstep.models
        ])
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

        print(" ".join([
            f"{datetime.datetime.now()}",
            f"Initial parameterization using reuse={reuse0}"
        ]))

        reuse = reuse0
        psys = {}
        for eid, gde in gdb.entries.items():
            tid = assignments.POSITIONS
            pos = assignments.graph_db_graph_to_graph_assignment(
                gdb,
                eid,
                tid,
            )
            psys[eid] = mm.chemical_system_to_physical_system(
                csys,
                pos,
                ref=psystems[eid],
                reuse=reuse
            )
        psystems = psys

        print(f"{datetime.datetime.now()} Saving checkpoint to chk.cst.p")
        pickle.dump([gdb, csys, strategy, psystems], open("chk.cst.p", "wb"))


        #######################################################################
        # go through the strategy and generate all candidates
        candidates, iters = generate_candidates(
            csys,
            psystems,
            gdb,
            G0,
            macro,
            strategy,
            union_cache,
            wq
        )

        if iters == 0:
            print(
                f"{datetime.datetime.now()} Warning, this macro step had no micro steps"
            )
            continue

        step = macro.steps[-1]

        iteration += iters

        print(f"{datetime.datetime.now()} Scanning done.")
        # *********************************************************************

        # print(datetime.datetime.now())
        # print(f"\n\nGenerating SMARTS on {len(candidates)}")
        Sj_sma = []

        with multiprocessing.pool.Pool(configs.processors) as pool:
            Sj_lst = []
            for (_,_,_,oper), x in candidates.items():
                
                if oper == strategy.SPLIT:
                    # Sj_lst = [candidates[x[1]][1] for x in pq]
                    # Sj_lst = [graphs.subgraph_as_structure(x[1], topo) for x in candidates.values()]
                    Sj_lst.append(x[1])
                elif oper == strategy.MERGE:
                    Sj_lst.append(
                        graphs.subgraph_as_structure(
                            mm.chemical_system_get_node_hierarchy(csys, x[1]).subgraphs[x[1].index],
                            mm.chemical_system_get_node_hierarchy(csys, x[1]).topology
                        )
                    )
                    
                    # Sj_lst = [
                    #     graphs.subgraph_as_structure(cst.hierarchy.subgraphs[x[1].index], topo)
                    #     for x in candidates.values()
                    # ]
                elif oper == strategy.MODIFY:
                    Sj_lst.append(
                        graphs.subgraph_as_structure(
                            mm.chemical_system_get_node_hierarchy(csys, x[0]).subgraphs[x[0].index],
                            mm.chemical_system_get_node_hierarchy(csys, x[0]).topology
                        )
                    )
                    # Sj_lst = [
                    #     graphs.subgraph_as_structure(cst.hierarchy.subgraphs[x[1].index], topo)
                    #     for x in candidates.values()
                    # ]
            Sj_sma = pool.map_async(gcd.smarts_encode, Sj_lst).get()
            del Sj_lst

        cnd_n = len(candidates)
        t = datetime.datetime.now()

        mm.chemical_system_print(csys, show_parameters=assigned_nodes)

        visited.clear()
        repeat.clear()

        procs = (
            os.cpu_count() if configs.processors is None else configs.processors
        )

        print(f"Scoring and filtering {len(candidates)} candidates for operation={step.operation}")

        candidates, Sj_sma = process_tiers(
            tiers,
            candidates,
            gdb,
            csys,
            psystems,
            chemical_objective,
            P00,
            CX0,
            X0,
            C0,
            P0,
            Sj_sma,
            strategy,
            step,
            reset_config,
            wq
        )

        print(f"Scanning {len(candidates)} candidates for operation={step.operation}")
        repeat, visited, success, n_added, csys, psystems, X0, P0, C0 = insert_candidates(
            candidates,
            Sj_sma,
            strategy,
            step,
            step_tracker,
            csys,
            psystems,
            gdb,
            assigned_nodes,
            reuse0,
            reset_config,
            chemical_objective,
            initial_objective,
            tiers,
            C0,
            P0,
            X0,
            CX0,
            G0,
            union_cache,
            visited,
            wq
        )

        if n_added > 0:

            strategy.repeat_step()
            print((
                f"{datetime.datetime.now()} "
                "Saving restart checkpoint to "
                f"restart_{iteration}_{macro.cursor}.p"
            ))
            pickle.dump(
                [gdb, csys, strategy, psystems],
                open(f"restart_{iteration}_{macro.cursor}.p", "wb")
            )
            smirnoff_models.smirnoff_write_version_0p3(
                csys,
                f"restart_{iteration}_{macro.cursor}.offxml"
            )

        print(f"There were {n_added} successful operations")

        print(f"{datetime.datetime.now()} Visited", visited)

        print(datetime.datetime.now(), "Saving chk.cs.p")
        pickle.dump([gdb, csys, strategy, psystems], open("chk.cst.p", "wb"))


        print(datetime.datetime.now(), "Macro step done.")

        print()
        print("#"*120)
        print()

    print(f"{datetime.datetime.now()} Strategy done.")
    CX0 = mm.chemical_system_smarts_complexity(csys)
    C = chemical_objective(csys, P0=math.log(len(psystems)+1), c=CX0)

    if final_objective and (initial_objective or tiers):

        print(f"{datetime.datetime.now()} Computing final fit")

        if ws:
            compute.workqueue_remove_workspace(wq, ws)
            ws.close()
            ws = None
        ws = None
        if configs.processors > 1:
            ws = compute.workqueue_new_workspace(wq)
        ret = reset(reset_config, csys, gdb, psystems, verbose=True, ws=ws)
        psystems = ret.value
        print("\n".join(ret.out))
        kv, _, P, gp = fit(csys, gdb, final_objective, psystems, assigned_nodes, wq=wq, verbose=True).value
        if ws:
            compute.workqueue_remove_workspace(wq, ws)
            ws.close()
            ws = None

        # C = chemical_objective(csys, P0=len(G0), C0=CX0)
        # print(f"{datetime.datetime.now()} C0={C0}")
        # X = P + C
        print(datetime.datetime.now(), f"Final objective: P={P:13.6g}")
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

    mm.chemical_system_print(csys)

    print(f"{datetime.datetime.now()} Saving final checkpoint to chk.cst.p")
    pickle.dump([gdb, csys, strategy, psystems], open("chk.cst.p", "wb"))

    print(f"{datetime.datetime.now()} Saving final (csys, (P0, P), (C0, C)) to csys.p")
    pickle.dump((csys, (P00, P), (C00, C)), open("csys.p", "wb"))

    return csys, (P00, P), (C00, C)


def calc_tier_distributed(S, Sj, operation, edits, oid, verbose=False, wq=None, shm=None):

    # copy once
    csys = copy.deepcopy(shm.csys)
    # csys = shm.csys
    print(f"{verbose=}")

    hidx = mm.chemical_system_get_node_hierarchy(csys, S)
    cm = mm.chemical_system_get_node_model(csys, S)
    cid = S.category[0]
    pid = S.category[1]
    uid = S.category[2]

    # node_ref = trees.tree_node_copy(S)

    gcd = csys.perception.gcd
    objective = shm.objective

    reuse = shm.reuse
    gdb = assignments.graph_db_get_entries(
        shm.gdb,
        [*set((e for o in oid for e in objective.objectives[o].addr.eid))]
    )

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
        node = mm.chemical_model_smarts_hierarchy_copy_node(cm, pid, uid, S, None)
        hidx.subgraphs[node.index] = Sj
        sma = gcd.smarts_encode(Sj)
        hidx.smarts[node.index] = sma

    elif operation == optimization.optimization_strategy.MERGE:
        node = Sj
        sma = hidx.smarts[node.index]
        g = hidx.subgraphs[node.index]
        mm.chemical_model_smarts_hierarchy_remove_node(cm, cid, pid, uid, Sj)

    elif operation == optimization.optimization_strategy.MODIFY:
        node = S
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
        ws = compute.workqueue_new_workspace(wq, address=None, shm={})
    elif configs.processors > 1 and len(objective.objectives) > 5:
        addr = ('127.0.0.1', 0)
        ws = compute.workspace_local(*addr, shm={})

    out = []
    ret = reset(shm.reset_config, csys, gdb, psysref, verbose=verbose, ws=ws)
    out.extend(ret.out)

    # keep a list of keys in case we reset but don't fit
    kvreset = mm.chemical_system_iter_keys(csys)

    out.append("RESETTING")
    for k, v in shm.reset_config.items():
        out.append(f"{k} : {v}")
    psysref = ret.value

    if verbose:
        print("\n".join(ret.out))
    if ws:
        if wq:
            compute.workqueue_remove_workspace(wq, ws)
        ws.close()
        ws = None

    assigned_nodes = sorted(set([
        (m, l) for psys in psysref.values()
        for m, pm in enumerate(psys.models)
        for proc in pm.labels
        for glbl in proc.values()
        for t, l in glbl.items()
    ]))

    fitkeys = [
        x
        for x in objective_tier_get_keys(objective, csys)
        if (x[0], x[2]) in assigned_nodes
    ]

    fitting_models = set((k[0] for k in fitkeys))
    reuse = [k for k, _ in enumerate(csys.models) if k not in fitting_models]

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

    elif operation in [
        optimization.optimization_strategy.MERGE,
        optimization.optimization_strategy.MODIFY
    ]:

        match_len = 0
        for psys in psysref.values():
            pm = psys.models[cid]
            for ic, lbls in pm.labels[0].items():
                lbls = set(lbls.values())
                if S.name in lbls:
                    match_len += 1

    # print("Matches", match_len, "Old matches", old_match)
    # C = graphs.graph_bits(Sj) / len(Sj.nodes) /1000 + len(Sj.nodes)
    kv = {}
    if keep:
        ret = objective_tier_run(
            objective,
            gdb,
            csys,
            fitkeys,
            oid=oid,
            psysref=psysref,
            reuse=reuse,
            wq=wq,
            verbose=verbose
        )
        kv, y0, P, gx = ret.value
        kvreset.update(kv)
        out.extend(ret.out)

    # print(f"{S.name}->{sma:40s} OID={oid} {keep} {X} {C} {match_len}")

    c = mm.chemical_system_smarts_complexity(csys)
    kv = kvreset

    return returns.success((keep, P, c, match_len, kv), out=out)


def print_chemical_system(csys, show_parameters=None):
    mm.chemical_system_print(csys, show_parameters=show_parameters)


def chemical_objective(csys, P0=1.0, C0=1.0, A=1.0, B=1.0, C=1.0, c=None):

    if c is None:
        c = mm.chemical_system_smarts_complexity(csys, B=B, C=C)
    CC = A * P0 * c**2
    # print(f"CHEM OBJECTIVE {CC=} {A=} {c=}")
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
        (edits, _, p_j, oper) = key
        # param_name = "p."
        cid, pid, uid = S.category
        cm = csys.models[cid]

        if oper == optimization.optimization_strategy.SPLIT:
            hidx = mm.chemical_system_get_node_hierarchy(csys, S)
            if hidx is None:
                breakpoint()
                print(f"Invalid node for operation: {S.name}")
                continue

            # param_name = "p" + str(group_number)
            node = mm.chemical_model_smarts_hierarchy_copy_node(cm, pid, uid, S, None)
            node.type = "parameter"
            # print(datetime.datetime.now(), '*** 2')
            assert Sj.select
            Sj = graphs.subgraph_relabel_nodes(
                Sj, {n: i for i, n in enumerate(Sj.select, 1)}
            )
            print(f" ::: Splitting {S.name} :::")
            # Sj = graphs.subgraph_to_structure(Sj, topo)

            hidx.subgraphs[node.index] = graphs.subgraph_copy(Sj)
            hidx.smarts[node.index] = Sj_sma[cnd_i-1]

            nodes[cnd_i] = node

        elif oper == optimization.optimization_strategy.MERGE:

            # S might have been deleted previously
            # and Sj is a graph in splitting
            if Sj.name in ignore:
                continue

            hidx = mm.chemical_system_get_node_hierarchy(csys, Sj)
            if hidx is None:
                breakpoint()
                print(f"Invalid node for operation: {Sj.name}")
                continue

            print(f" ::: Removing {Sj.name} from parent {S.name} :::")
            mm.chemical_model_smarts_hierarchy_remove_node(cm, cid, pid, uid, Sj)
            ignore.add(Sj.name)
            nodes[cnd_i] = Sj

        elif oper == optimization.optimization_strategy.MODIFY:

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
    return assignments.graph_assignment_to_format_xyz(pos, comment=comment)


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

    lbls = [n.name for n in hidx.index.nodes.values()]
    optimization.target_list = lbls
    optimization.reference_list = lbls
    optimization.merge_protect_list = lbls
    initial_conditions = clusters.clustering_initial_conditions(
        gcd,
        sag,
        hidx=hidx,
        labeler=labeler,
        prefix=csys.models[m].symbol
    )

    cluster_fn = cluster_optimization.cluster_means
    if objective.is_discrete():
        cluster_optimization.cluster_classifications

    cst = cluster_fn(
        gcd,
        labeler,
        sag,
        objective,
        optimization=optimization,
        initial_conditions=initial_conditions
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

    for idx in [
        n.index
        for n in refhidx.index.nodes.values()
        if n.index != root.index
    ]:
    #     if idx in refhidx.smarts:
    #         refhidx.smarts.pop(idx)
    #     name = refhidx.index.nodes[idx].name
    #     for (ui, namei) in list(proc.topology_parameters):
    #         if namei == name:
    #             proc.topology_parameters.pop((ui, namei))
    #     for sym, terms in cm.topology_terms.items():
    #         if name in terms.values:
    #             terms.values.pop(name)

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

        if new_node.name not in lbls:
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

def smiles_assignment_force_constants(
    gdb,
    alpha=-.25,
    guess_periodicity=True,
    max_n=3,
    min_dihedral_k=1e-3,
    max_dihedral_k=100
):

    sag_map = {
        "bond_k": [],
        "bond_l": [],
        "angle_k": [],
        "angle_l": [],
        "torsion_n": [],
        "torsion_p": [],
        "torsion_k": [],
        "pair_l": [],
    }

    eid_hess = [eid for eid, e in gdb.entries.items() if HESSIANS in e.tables]

    for i, eid in enumerate(eid_hess, 1):
        # psys = psystems[eid]
        hessian = gdb.entries[eid].tables[HESSIANS].values

        # labeled_ic = [x for model in psystems[eid].models for x in model.labels[0]]
        # B = {ic:x for ic, x in zip(icb, B) if ic in labeled_ic}
        # B = list(B.keys()), list(B.values())

        pos = assignments.graph_db_graph_to_graph_assignment(gdb, eid, assignments.POSITIONS)
        # g = pos.graph

        icb, B = assignments.bmatrix(pos, torsions=True, outofplanes=False, pairs=False ,remove1_3=True, linear_torsions=145.0)

        # xyz = np.vstack([x[0] for x in pos.selections.values()], dtype=float)
        # xyz = xyz.round(PRECISION)
        # sym = graphs.graph_symbols(pos.graph)
        # mass = np.array([[vibration.mass_table[sym[n]]]*3 for n in sym])

        # sym = list(sym.values())

        # hess_qm_freq, hess_qm_modes = vibration.hessian_modes(hess_qm, sym, xyz, mass, 0, remove=0, stdtr=True, verbose=False)
        # omega = np.round(hess_qm_freq, PRECISION)
        # omega_qm = omega
        ic_qm_fcs = {}
        if B:

            hess_qm_ic = hessians.project_ics(B, hessian)
            hess_qm_ic = [hess_qm_ic[i][i] for i in range(len(hess_qm_ic))]

            ic_qm_fcs = dict(zip(icb, hess_qm_ic))

        sel = {}
        bidx = assignments.graph_assignment_matrix_bond_indices(pos)
        for ic in bidx:
            if ic in ic_qm_fcs:
                sel[ic] = [ic_qm_fcs[ic]]

        # ga = assignments.graph_assignment(pos.smiles, sel, pos.graph)
        sag_map["bond_k"].append([pos, sel])

        bond_l = assignments.graph_assignment_geometry_bond_matrix(pos, bidx)
        bond_l.selections = {
            ic: [x for y in v for x in y]
            for ic, v in bond_l.selections.items()
        }
        # ga = assignments.graph_assignment(pos.smiles, bond_l.selections, pos.graph)
        sag_map["bond_l"].append([pos, bond_l])

        pair_l = assignments.graph_assignment_geometry_pair_matrix(pos)
        pair_l.selections = {
            ic: [x for y in v for x in y]
            for ic, v in pair_l.selections.items()
        }
        # ga = assignments.graph_assignment(pos.smiles, pair_l.selections, pos.graph)
        sag_map["pair_l"].append([pos, sel])

        sel = {}
        aidx = assignments.graph_assignment_matrix_angle_indices(pos)
        for ic in aidx:
            if ic in ic_qm_fcs:
                sel[ic] = [ic_qm_fcs[ic]]
        # ga = assignments.graph_assignment(pos.smiles, sel, pos.graph)
        sag_map["angle_k"].append([pos, sel])

        angle_l = assignments.graph_assignment_geometry_angle_matrix(pos, aidx)
        angle_l.selections = {
            ic: [x for y in v for x in y]
            for ic, v in angle_l.selections.items()
        }
        # ga = assignments.graph_assignment(pos.smiles, angle_l.selections, pos.graph)
        sag_map["angle_l"].append([pos, angle_l])

        tidx = assignments.graph_assignment_matrix_torsion_indices(pos)
        dih_angles = assignments.smiles_assignment_geometry_torsion_matrix_nonlinear(pos)

        sel_n = {}
        sel_k = {}
        sel_p = {}
        for angle_ic, angle_val in dih_angles.selections.items():
            bond = geometry.bond(angle_ic[1:3])
            all_dihed = [v for ic, v in dih_angles.selections.items() if ic[1:3] == bond]
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
                    if (
                        all(math.cos(i*t) < alpha for t in angles)
                        or all(math.cos(i*t) > -alpha for t in angles)
                    ):
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
                if angles:
                    csys_n = [n for n in range(1, max_n+1)]
                    csys_p = [0 for n in range(1, max_n+1)]
                else:
                    csys_n = [n for n in csys_n0]
                    csys_p = [n for n in csys_p0]
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
            print(f"Calculated final", new_n, new_k, new_p)

        for ic in tidx:
            if ic not in sel_k:
                sel_k[ic] = []
                sel_n[ic] = []
                sel_p[ic] = []
        # ga = assignments.graph_assignment(pos.smiles, sel_k, pos.graph)
        sag_map["torsion_k"].append([pos, sel_k])

        # ga = assignments.graph_assignment(pos.smiles, sel_n, pos.graph)
        sag_map["torsion_n"].append([pos, sel_n])

        # ga = assignments.graph_assignment(pos.smiles, sel_p, pos.graph)
        sag_map["torsion_p"].append([pos, sel_p])

    return sag_map


def reset_project_dihedrals(
    csys,
    psystems,
    psys_hess,
    vals,
    psys_grad,
    eid_grad,
    gvals,
    m,
    lbl,
    guess_periodicity=True,
    max_n=3,
    max_k=10,
    alpha=-.5,
    min_k=1e-3,
    verbose=False
):

    out = []

    csys_n0 = mm.chemical_system_get_value_list(csys, (m, 'n', lbl))
    csys_p0 = mm.chemical_system_get_value_list(csys, (m, 'p', lbl))
    csys_k0 = mm.chemical_system_get_value_list(csys, (m, 'k', lbl))
    npk0 = csys_n0, csys_p0, csys_k0

    # go into the psys and get the angles
    grad_angles = []
    grad_gq = []
    angles = []
    idiv = []
    for psys in psys_hess:
        pos = psys.models[m].positions
        indices = [
            ic for ic, terms in psys.models[m].labels[0].items()
            if lbl in terms.values()
        ]
        inner_bonds = [geometry.bond(ic[1:3]) for ic in indices]
        idiv.extend([(
            (max(1, len(graphs.graph_connection(posi.graph, a[1]))-1)) *
            (max(1, len(graphs.graph_connection(posi.graph, b[1]))-1)))
            for posi in pos
            for a, b in inner_bonds
        ])
        if m == 2:
            psys_angles = assignments.smiles_assignment_geometry_torsion_matrix_nonlinear(pos, indices=indices)
        elif m == 3:
            psys_angles = assignments.smiles_assignment_geometry_outofplane_matrix(pos, indices=indices)
        psys_angles = [
            t for y in psys_angles.selections.values()
            for x in y for t in x
        ]
        angles.extend(psys_angles)
    if gvals[m][lbl]:
        grad_gq = gvals[m][lbl]
    for i, eid in enumerate(eid_grad, 1):
        # if eid not in psys_grad:
        #     continue
        psys = psystems[eid]
        pos = psys.models[m].positions
        indices = [
            ic for ic, terms in psys.models[m].labels[0].items()
            if lbl in terms.values()
        ]
        # inner_bonds = [geometry.bond(ic[1:3]) for ic in indices]
        # idiv.extend([(
        #     (max(1, len(graphs.graph_connection(pos.graph, a))-1)) *
        #     (max(1, len(graphs.graph_connection(pos.graph, b))-1)))
        #     for a, b in inner_bonds
        # ])
        if m == 2:
            psys_angles_sel = assignments.smiles_assignment_geometry_torsion_matrix_nonlinear(pos, indices=indices)
        elif m == 3:
            psys_angles_sel = assignments.smiles_assignment_geometry_outofplane_matrix(pos, indices=indices)
        psys_angles = [
            t for y in psys_angles_sel.selections.values()
            for x in y for t in x
        ]
        grad_angles.extend(psys_angles)
    #     icb, B = assignments.bmatrix(pos, **bmat_config)
    #     gx = arrays.array_scale(gdb.entries[eid].tables[GRADIENTS].values, 1/4.184)
    #     gq = dict(zip(icb, hessians.project_gradient(B, gx)))

    #     # grad_angle and grad_gq are (should be) in order
    #     for ic in indices:
    #         v = gq.get(ic)
    #         if v is not None:
    #             # grad_gq.append(v * au2kcal)
    #             grad_gq.append(v * kj2kcal)
    #             # ang = psys_angles_sel.selections.get(ic)
    #             # if ang is None:
    #             #     continue
    #             # ang = ang[0][0]
    #             # indicator = "  "
    #             # if eid in eid_hess:
    #             #     indicator = "G0"
    #             # else:
    #             #     indicator = "G "
    #             # line = (
    #             #     f"{lbl:6s} {eid:5d} "
    #             #     f"{'-'.join(map(str, ic)):12s} "
    #             #     f"{ang:12.6f} {indicator} {v:12.6f}"
    #             #     "\n"
    #             # )
    #             # data_file.write(line)
    #             # data_file.flush()

    idiv_expected = (sum(idiv)/len(idiv))/4
    # SCALE
    # idiv_expected = 1.0
    line = f"idiv unique: {set(idiv)}"
    out.append(line)

    # idiv_expected = 1
    line = f"Expected idiv: {idiv_expected}"
    out.append(line)
    # print("Expected idiv:", idiv_expected)
    # print(angles)


    # line = f"Expected idiv: {idiv_expected}"
    line = f"Reference npk0 {npk0}"
    out.append(line)

    if guess_periodicity:
        csys_n = []
        csys_p = []
        for i in [*range(1, max_n+1)]:
            if (
                all(math.cos(i*t) < alpha for t in angles) or
                all(math.cos(i*t) > -alpha for t in angles)
            ):
                csys_n.append(i)
                csys_p.append(0)
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
    hq = sum(vals) # / len(vals)

    gq = None
    if grad_gq:
        gq = sum(grad_gq) # /len(grad_gq)
        line = f"Expected gq: {gq:.6g} kcal/mol/rad"
        out.append(line)
    if gq:
        gq_err = arrays.array_translate(grad_gq, -gq)
        gq_sse = arrays.array_inner_product(gq_err, gq_err)
        gq_err = None
        line = f"Expected gq var: {gq_sse/len(grad_gq):.6g} (kcal/mol/rad)^2"
        out.append(line)
        line = f"Expected gq minmax: {min(grad_gq):.6g}, {max(grad_gq):.6g}"
        out.append(line)

    if not csys_n:
        line = (
            "Could not find any appropriate periodicities for "
            f"label={lbl} max_n={max_n} alpha={alpha}"
        )
        out.append(line)
        csys_n = [n for n in range(1, max_n+1)]
        csys_p = [0 for n in range(1, max_n+1)]
        line = "First 10 Angles:"
        out.append(line)
        line = " ".join(map("{:6g}".format, angles[:10]))
        out.append(line)
        line = "Cosines:"
        out.append(line)
        for i in range(1, max_n+1):
            line = f"{i} " + " ".join(
                map(lambda t: f"{math.cos(i*t):6g}", angles[:10])
            )
            out.append(line)
            # print(i, [*(math.cos(i*t) for t in angles[:10])])

    changed = True
    while changed:
        A = []
        b = []
        changed = False

        if angles:
            # for hqi, angi in zip(vals, angles):
            #     out.append(f"Angle: {angi:12.4f} h_dih {hqi:12.6f}")
            if grad_gq and grad_angles:
                # for gqi, angi in zip(grad_gq, grad_angles):
                #     out.append(f"Angle: {angi:12.4f} g_dih {gqi:12.6f}")
                b.append(gq)
                row = []
                # angs = [scipy.stats.circmean(grad_angles, high=np.pi, low=-np.pi)]
                angs = grad_angles
                for n, p in zip(csys_n, csys_p):
                    i = -1
                    x = [sign[(i+2)%4] * n**(i+2)*deriv[(i+2) % 4](n*t - p) for t in angs]
                    x = sum(x) * idiv_expected
                    row.append(x)

                out.append("Using gradients in fit; first derivative row is:")
                out.append(" ".join(map("{:12.6f}".format, row)))
                A.append(row)
            for i in range(len(csys_n)):
                if i == 0:
                    b.append(hq)
                else:
                    b.append(0)
                row = []
                for n, p in zip(csys_n, csys_p):
                    # angs = [scipy.stats.circmean(angles, high=np.pi, low=-np.pi)]
                    angs = angles
                    x = [sign[(i+2)%4] * n**(i+2)*deriv[(i+2) % 4](n*t - p) for t in angs]
                    x = sum(x) * idiv_expected
                    row.append(x)
                A.append(row)
            out.append("A is:")
            for r in A:
                out.append(" ".join(map("{:12.5f}".format, r)))
            out.append("B is:")
            out.append(" ".join(map("{:12.5f}".format, b)))

            # new_k_lst = np.linalg.solve(A, b)
            new_k_lst, residual, rank, s = np.linalg.lstsq(A, b, rcond=None)
            line = "Calculated " + " ".join(map(lambda x: "n={:d} k={:12.5f}".format(*x),zip(csys_n, new_k_lst)))
            out.append(line)
            line = f"Error {residual}"
            out.append(line)
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
            for n, p, k in zip(csys_n, csys_p, new_k_lst):
                if abs(k) > min_k and abs(k) < max_k:
                    new_n.append(n)
                    new_p.append(p)
                    new_k.append(k)
                    # line = f"Added n= {n:d} {p:.2f} {k:12.5f}"
                    # out.append(line)
                # else:
                #     new_k = []
                #     new_p = []
                #     new_n = []

            if set(new_n).symmetric_difference(csys_n):
                line = "Warning, some fitted values were out of range"
                out.append(line)
                new_k = []
                new_p = []
                new_n = []

                s = arrays.argsort([abs(x) for x in new_k_lst])
                if len(s) > 1:
                    line = f"Removing n= {csys_n[s[-1]]:d}"
                    out.append(line)
                    s = s[:-1]
                new_n.extend([csys_n[i] for i in s])
                new_p.extend([csys_p[i] for i in s])
                for i in s:
                    k = new_k_lst[i]
                    if k < -max_k:
                        k = -max_k
                    elif k > max_k:
                        k = max_k
                    new_k.append(k)
                new_k_lst = new_k
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

    if guess_periodicity:

        mm.chemical_system_set_value_list(csys, (m, 'n', lbl), new_n)
        line = (
            f"Setting {lbl} N_g= {len(grad_gq):5d} n "
            f"from {npk0[0]} to {new_n}"
        )
        out.append(line)

        mm.chemical_system_set_value_list(csys, (m, 'p', lbl), new_p)
        line = (
            f"Setting {lbl} N_g= {len(grad_gq):5d} p "
            f"from {npk0[1]} to {new_p}"
        )
        out.append(line)
    line = (
        f"Setting {lbl} N= {len(vals):5d} k "
        f"from {csys_k0} to {new_k_lst}"
    )
    out.append(line)

    return returns.success(new_k_lst, out=out)


def reset(
    reset_config,
    csys,
    gdb,
    psystems=None,
    use_min_gradients=True,
    verbose=False,
    ws=None,
    shm=None
):

    assert type(reset_config) is dict
    out = []
    bond_k_skip = reset_config.get("bond_k_skip", [])
    bond_l_skip = reset_config.get("bond_l_skip", [])
    angle_k_skip = reset_config.get("angle_k_skip", [])
    angle_l_skip = reset_config.get("angle_l_skip", [])
    torsion_k_skip = reset_config.get("torsion_k_skip", [])
    torsion_n_skip = reset_config.get("torsion_n_skip", [])
    outofplane_k_skip = reset_config.get("outofplane_k_skip", [])
    outofplane_n_skip = reset_config.get("outofplane_n_skip", [])

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
    min_k = reset_config.get("dihedral_min_k", 1e-3)

    if psystems is None:
        psystems = gdb_to_physical_systems(gdb, csys)
    if not any([
        reset_bond_k,
        reset_bond_l,
        reset_angle_k,
        reset_angle_l,
        reset_torsion_k,
        reset_outofplane_k
    ]):
        out.append("No reset terms enabled. Skipping.")
        return returns.success(psystems)

    if reset_config.get("unique", False):
        return returns.success({
            eid: reset_unique(
                reset_config,
                gdb.entries[eid],
                ps, verbose=verbose
            )
            for eid, ps in psystems.items()
        })

    if ws:
        ws.reset()

    eid_hess = [eid for eid, e in gdb.entries.items() if HESSIANS in e.tables]
    eid_grad = []
    if use_min_gradients:
        eid_grad.extend([
            eid for eid, e in gdb.entries.items() if GRADIENTS in e.tables
        ])
    out.append(f"There are {len(eid_grad)} gradients")
    psys_hess = [psystems[eid] for eid in eid_hess]
    if use_min_gradients:
        psys_grad = [psystems[eid] for eid in eid_grad]
        out.append(f"Setting to use gradients in fits N={len(psys_grad)}")
    else:
        psys_grad = [psystems[eid] for eid in eid_grad if eid not in eid_hess]
        out.append(f"Setting to use not gradients in fits N={len(psys_grad)}")

    reset = set()
    if reset_bond_l:
        mod = mm.chemical_system_reset_bond_lengths(csys, psys_hess, skip=bond_l_skip)
        # if verbose:
        #     mm.chemical_system_print(csys, show_parameters=[x[2] for x in mod])
        reset.add(0)
    if reset_angle_l:
        mod = mm.chemical_system_reset_angles(csys, psys_hess, skip=angle_l_skip)
        if verbose:
            mm.chemical_system_print(csys, show_parameters=[x[2] for x in mod])
        reset.add(1)
    if reset_torsion_k:
        hierarchy = csys.models[2].procedures[0].smarts_hierarchies[0]
        names = [x.name for x in hierarchy.index.nodes.values()]
        if verbose:
            mm.chemical_system_print(csys, show_parameters=names)
        reset.add(2)
    if reset_outofplane_k:
        hierarchy = csys.models[3].procedures[0].smarts_hierarchies[0]
        names = [x.name for x in hierarchy.index.nodes.values()]
        if verbose:
            mm.chemical_system_print(csys, show_parameters=names)
        reset.add(3)

    bmat_config = dict(
        torsions=True,
        outofplanes=True,
        pairs=False,
        remove1_3=True
    )
    hvals = {}
    psys_hic_all = {}

    gvals = {}
    psys_gic_all = {}
    if reset_bond_k or reset_angle_k or reset_torsion_k or reset_outofplane_k:
        if ws:
            iterable = {}
            bmats = {}
            for eid in eid_hess:
                pos = psystems[eid].models[0].positions
                if eid in bmats:
                    B = bmats[eid]
                else:
                    icb, B = assignments.bmatrix(pos, **bmat_config)
                    labeled_ic = [
                        x
                        for model in psystems[eid].models
                        for x in model.labels[0]
                    ]
                    B = {ic:x for ic, x in zip(icb, B) if ic in labeled_ic}
                    B = list(B.keys()), list(B.values())
                    bmats[eid] = B

                psys = psystems[eid]
                # models = list(psys.models)
                # psys.models = [None, None, None, None, psys.models[4], psys.models[5]]
                # hess_mm = optimizers_openmm.physical_system_hessian_openmm(
                #     psys,
                #     csys,
                #     h=1e-4
                # )
                # psys.models = models
                hessian = gdb.entries[eid].tables[HESSIANS].values
                # hessian = [arrays.array_difference(a, arrays.array_scale(b, 1/4.184)) for a, b in zip(hessian, hess_mm)]
                iterable[(eid, "hess")] = [
                    [
                        csys,
                        psystems[eid],
                        hessian
                    ],
                    {"B": B}
                ]
            results = compute.workspace_submit_and_flush(
                ws,
                hessians.hessian_project_onto_ics,
                iterable,
                verbose=verbose,
                chunksize=100
            )
            psys_hic_all.update({
                k[0]: v for k, v in results.items() if k[1] == "hess"
            })
            iterable.clear()
            results.clear()
            for eid in eid_grad:
                pos = psystems[eid].models[0].positions
                if eid in bmats:
                    B = bmats[eid]
                else:
                    icb, B = assignments.bmatrix(pos, **bmat_config)
                    labeled_ic = [
                        x
                        for model in psystems[eid].models
                        for x in model.labels[0]
                    ]
                    B = {ic: x for ic, x in zip(icb, B) if ic in labeled_ic}
                    B = list(B.keys()), list(B.values())
                    bmats[eid] = B

                gx = arrays.array_scale(
                    gdb.entries[eid].tables[GRADIENTS].values,
                    kj2kcal
                )
                iterable[(eid, "grad")] = [[B[1], gx], {}]
            results = compute.workspace_submit_and_flush(
                ws,
                hessians.project_gradient,
                iterable,
                verbose=verbose,
                chunksize=100
            )
            iterable.clear()
            for eid in eid_grad:
                pos = psystems[eid].models[0].positions
                icb = bmats[eid][0]
                sel = dict(zip(icb, results[(eid, "grad")]))
                psys_gic_all[eid] = sel
            results.clear()
            bmats.clear()
        else:
            for i, eid in enumerate(eid_hess, 1):
                pos = psystems[eid].models[0].positions
                psys = psystems[eid]
                # models = list(psys.models)
                # psys.models = [None, None, None, None, psys.models[4], psys.models[5]]
                # hess_mm = optimizers_openmm.physical_system_hessian_openmm(
                #     psys,
                #     csys,
                #     h=1e-4
                # )
                # psys.models = models
                hessian = gdb.entries[eid].tables[HESSIANS].values
                # hessian = [arrays.array_difference(a, arrays.array_scale(b, 1/4.184)) for a, b in zip(hessian, hess_mm)]
                line = (
                    f"Projecting hessian for EID {eid} "
                    f"{i:8d}/{len(eid_hess)} "
                    # f"{psys.models[0].positions[0].smiles}"
                )
                out.append(line)
                icb, B = assignments.bmatrix(pos, **bmat_config)

                labeled_ic = [x for model in psystems[eid].models for x in model.labels[0]]
                B = {ic:x for ic, x in zip(icb, B) if ic in labeled_ic}
                B = list(B.keys()), list(B.values())

                hic = hessians.hessian_project_onto_ics(csys, psys, hessian, verbose=verbose, B=B)

                psys_hic_all[eid] = hic
            for i, eid in enumerate(eid_grad, 1):
                pos = psystems[eid].models[0].positions
                psys = psystems[eid]
                gx = arrays.array_scale(gdb.entries[eid].tables[GRADIENTS].values, 1/4.184)
                line = (
                    f"Projecting gradient for EID {eid} "
                    f"{i:8d}/{len(eid_hess)} "
                    # f"{psys.models[0].positions[0].smiles}"
                )
                out.append(line)
                icb, B = assignments.bmatrix(pos, **bmat_config)

                sel = {}
                labeled_ic = [x for model in psystems[eid].models for x in model.labels[0]]
                B = [x for ic, x in zip(icb, B) if ic in labeled_ic]
                if B:
                    # B = list(B.keys()), list(B.values())
                    gq = hessians.project_gradient(B, gx)
                    sel = dict(zip(icb, gq))
                    psys_gic_all[eid] = sel

    if reset_bond_k:
        psys_hics = []
        psyss = []
        for eid in eid_hess:
            bidx = assignments.graph_assignment_matrix_bond_indices(psys.models[0].positions)
            if not bidx:
                continue
            if eid not in psystems or eid not in psys_hic_all:
                continue
            psys = psystems[eid]
            hic = psys_hic_all[eid]
            bidx = assignments.graph_assignment_matrix_bond_indices(psys.models[0].positions)
            psys_hics.append({
                k: [[hic[k]]]
                for k in bidx
                if k in hic
            })
            psyss.append(psys)

        hvals[0] = mm.chemical_system_groupby_names(csys, 0, psyss, psys_hics)
        psyss.clear()
        psys_hics.clear()
        reset.add(0)

    if reset_angle_k:
        psys_hics = []
        psyss = []
        for eid in eid_hess:
            aidx = assignments.graph_assignment_matrix_angle_indices(psys.models[0].positions)
            if not aidx:
                continue
            if eid not in psystems or eid not in psys_hic_all:
                continue
            psys = psystems[eid]
            hic = psys_hic_all[eid]
            psys_hics.append({
                k: [[hic[k]]]
                for k in aidx
                if k in hic
            })
            psyss.append(psys)

        hvals[1] = mm.chemical_system_groupby_names(csys, 1, psyss, psys_hics)
        psyss.clear()
        psys_hics.clear()
        reset.add(1)

    if reset_torsion_k:
        psys_hics = []
        psyss = []
        for eid in eid_hess:
            tidx = assignments.graph_assignment_matrix_torsion_indices(psys.models[0].positions)
            if not tidx:
                continue
            if eid not in psystems or eid not in psys_hic_all:
                continue

            hic = psys_hic_all[eid]
            psys = psystems[eid]
            psys_hics.append({
                k: [[hic[k]]]
                for k in tidx
                if k in hic
            })
            psyss.append(psys)
        hvals[2] = mm.chemical_system_groupby_names(csys, 2, psyss, psys_hics)
        reset.add(2)

        psyss.clear()
        psys_hics.clear()

        for eid in eid_grad:
            indices = assignments.smiles_assignment_geometry_torsion_matrix_nonlinear(pos)
            if not indices:
                continue
            if eid not in psystems or eid not in psys_gic_all:
                continue
            psys = psystems[eid]
            pos = psys.models[0].positions
            hic = psys_gic_all[eid]
            psys_hics.append({
                k: [[hic[k]]]
                for k in indices.selections
                if k in hic
            })
            psyss.append(psys)
        gvals[2] = mm.chemical_system_groupby_names(csys, 2, psyss, psys_hics)

    if reset_outofplane_k:
        psys_hics = []
        psyss = []
        for eid in eid_hess:
            oidx = assignments.graph_assignment_matrix_outofplane_indices(psys.models[0].positions)
            if not oidx:
                continue
            if eid not in psystems or eid not in psys_hic_all:
                continue
            psys = psystems[eid]
            hic = psys_hic_all[eid]
            psys_hics.append({
                k: [[hic[k]]]
                for k in oidx
                if k in hic
            })
            psyss.append(psys)

        hvals[3] = mm.chemical_system_groupby_names(csys, 3, psyss, psys_hics)
        psyss.clear()
        psys_hics.clear()

        for eid in eid_grad:
            oidx = assignments.graph_assignment_matrix_outofplane_indices(psys.models[0].positions)
            if not oidx:
                continue
            if eid not in psystems or eid not in psys_gic_all:
                continue
            psys = psystems[eid]
            hic = psys_gic_all[eid]
            psys_hics.append({
                k: [[hic[k]]]
                for k in oidx
                if k in hic
            })
            psyss.append(psys)
        gvals[3] = mm.chemical_system_groupby_names(csys, 3, psyss, psys_hics)
        reset.add(3)

    psys_hic_all.clear()
    psys_gic_all.clear()

    # data_file = open("ic_grads.dat", 'a')


    for m, hic in hvals.items():
        for lbl, vals in hic.items():
            line = f"Resetting {lbl} N= {len(vals):5d}"
            out.append(line)
            if not vals:
                continue
            if (
                (m == 0 and lbl in bond_k_skip) or
                (m == 1 and lbl in angle_k_skip) or
                (m == 2 and lbl in torsion_k_skip) or
                (m == 3 and lbl in outofplane_k_skip)
            ):
                continue
            try:
                csys_k = mm.chemical_system_get_value_list(csys, (m, 'k', lbl))
            except KeyError:
                continue
            if m == 1:
                # new_k = min(vals)
                new_k = sum(vals)/len(vals)
                angles = []
                idiv = []
                for psys in psys_hess:
                    pos = psys.models[m].positions
                    indices = [
                        ic[1] for ic, terms in psys.models[m].labels[0].items()
                        if lbl in terms.values()
                    ]
                    nbr = [
                        len(graphs.graph_connection(pos[a[0]].graph, a[1]))
                        for a in indices
                    ]
                    idiv.extend(nbr)
                idiv_expected = sum(idiv)/len(idiv)
                # idiv_expected = 1.0
                # SCALE
                # if reset_outofplane_k or reset_torsion_k:
                #     new_k *= idiv_expected * 1.3 / 3
                if reset_outofplane_k or reset_torsion_k:
                    new_k *= idiv_expected/ 2
                new_k_lst = [new_k for _ in csys_k]
                line = (
                    f"Setting {lbl} N= {len(vals):5d} k "
                    f"from {csys_k} to {new_k_lst}"
                )
                out.append(line)
            elif m == 0:
                new_k = sum(vals)/len(vals)
                new_k_lst = [new_k for _ in csys_k]
                line = (
                    f"Setting {lbl} N= {len(vals):5d} k "
                    f"from {csys_k} to {new_k_lst}"
                )
                out.append(line)
            elif (
                (reset_torsion_k and m == 2) or
                (reset_outofplane_k and m == 3)
            ):

                p = guess_periodicity
                if m == 2 and torsion_n_skip:
                    p = False
                if m == 3 and outofplane_n_skip:
                    p = False
                ret = reset_project_dihedrals(
                    csys,
                    psystems,
                    psys_hess,
                    vals,
                    psys_grad,
                    eid_grad,
                    gvals,
                    m,
                    lbl,
                    guess_periodicity=p,
                    max_n=max_n,
                    max_k=max_k,
                    alpha=alpha,
                    min_k=min_k,
                    verbose=verbose
                )
                new_k_lst = ret.value
                out.extend(ret.out)
                reset.add(m)

            # if verbose:
            mm.chemical_system_set_value_list(csys, (m, 'k', lbl), new_k_lst)
            reset.add(m)

    # data_file.close()
    hvals.clear()

    if reset:
        line = "Resetting systems with new values"
        out.append(line)
        psys = {}
        reuse = list(set(range(len(csys.models))).difference(reset))
        for eid, gde in gdb.entries.items():
            tid = POSITIONS
            gid = list(gde.tables[tid].graphs)
            rid = None
            pos = assignments.graph_db_graph_to_graph_assignment(
                gdb,
                eid,
                tid,
                gid,
                rid
            )
            psys[eid] = mm.chemical_system_to_physical_system(
                csys,
                pos,
                ref=psystems[eid],
                reuse=list(reuse)
            )
        psystems = psys

    return returns.success(psystems, out=out)


def reset_unique_dihedral(
    psystem,
    hic,
    m,
    guess_periodicity,
    alpha,
    max_n,
    min_k,
    max_k,
    gradients=None,
    verbose=False
):

    pos = psystem.models[m].positions
    for ic in psystem.models[m].values:
        idiv_expected = 1
        if m == 2:
            angles = assignments.smiles_assignment_geometry_torsion_matrix_nonlinear(psystem.models[m].positions, indices=[ic])
            nba = len([x for x in pos.graph.edges if ic[1] in x])
            nbb = len([x for x in pos.graph.edges if ic[2] in x])
            # SCALE UNIQUE
            idiv_expected = ((nba-1) * (nbb-1)) / 4
            # idiv_expected = 100
            # print("IDIV IS", idiv_expected)
        elif m == 3:
            angles = assignments.smiles_assignment_geometry_outofplane_matrix(psystem.models[m].positions, indices=[ic])
        angles = tuple((z for y in angles.selections.values() for x in y for z in x))

        csys_n0 = psystem.models[m].values[0][ic]['n']
        csys_p0 = psystem.models[m].values[0][ic]['p']
        csys_k0 = psystem.models[m].values[0][ic]['k']

        npk0 = csys_n0, csys_p0, csys_k0

        if verbose:
            print("Reference npk0", npk0)
        if guess_periodicity:
            csys_n = []
            csys_p = []
            for i in [*range(1, max_n+1)]:
                if all(math.cos(i*t) < alpha for t in angles) or all(math.cos(i*t)  >  -alpha for t in angles):
                    csys_n.append(i)
                    csys_p.append(0)
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
        hq = hic.selections.get(ic, 0)

        out = []
        grad_gq = None
        if gradients:
            grad_gq = gradients.selections.get(ic)
        gq = None
        if grad_gq:
            gq = sum(grad_gq) # /len(grad_gq)
            line = f"Expected gq: {gq:.6g} kcal/mol/rad"
            out.append(line)
        if gq:
            gq_err = arrays.array_translate(grad_gq, -gq)
            gq_sse = arrays.array_inner_product(gq_err, gq_err)
            gq_err = None
            line = f"Expected gq var: {gq_sse/len(grad_gq):.6g} (kcal/mol/rad)^2"
            out.append(line)
            line = f"Expected gq minmax: {min(grad_gq):.6g}, {max(grad_gq):.6g}"
            out.append(line)
        if out:
            print("\n".join(out))
            out.clear()

        if not csys_n:
            if angles:
                csys_n = [n for n in range(1, max_n+1)]
                csys_p = [0 for n in range(1, max_n+1)]
            else:
                csys_n = [n for n in csys_n0]
                csys_p = [n for n in csys_p0]

            if verbose:
                print(f"Could not find any appropriate periodicities for ic={ic} max_n={max_n} alpha={alpha}")
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
                if gq:
                    b.append(gq)
                    row = []
                    # angs = [scipy.stats.circmean(grad_angles, high=np.pi, low=-np.pi)]
                    for n, p in zip(csys_n, csys_p):
                        i = -1
                        x = [sign[(i+2)%4] * n**(i+2)*deriv[(i+2) % 4](n*t - p) for t in angles]
                        x = sum(x) * idiv_expected
                        row.append(x)

                    out.append("Using gradients in fit; first derivative row is:")
                    out.append(" ".join(map("{:12.6f}".format, row)))
                    if out:
                        print("\n".join(out))
                        out.clear()
                    A.append(row)
                for i in range(len(csys_n)):
                    row = []
                    if i == 0:
                        b.append(hq)
                    else:
                        b.append(0)
                    for n, p in zip(csys_n, csys_p):
                        x = [sign[(i+2)%4] * n**(i+2)*deriv[(i+2) % 4](n*t - p) for t in angles]
                        # x = sum(x) / len(x) / idiv_expected
                        x = sum(x) * idiv_expected
                        # x = sum(x) / len(x)
                        row.append(x)
                    A.append(row)
                new_k_lst, residual, rank, s = np.linalg.lstsq(A, b, rcond=None)
                # new_k_lst = np.linalg.solve(A, b)
                if verbose:
                    print("A is:")
                    for r in A:
                        print(" ".join(map("{:12.5f}".format, r)))
                    print("B is:")
                    print(" ".join(map("{:12.5f}".format, b)))
                    print("Calculated npk", csys_n, csys_p, new_k_lst)
                    line = f"Error {residual}"
                    print(line)
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
                        if verbose:
                            print("Added n=", n, p, k)
                    else:
                        new_k = []
                        new_p = []
                        new_n = []
                        break

                if not new_n:
                    if verbose:
                        print("Warning, all fitted values were out of range")
                    i = arrays.argmin([abs(x) for x in new_k_lst])
                    s = arrays.argsort([abs(x) for x in new_k_lst])
                    s = [s[0]]
                    new_n.extend([csys_n[i] for i in s])
                    new_p.extend([csys_p[i] for i in s])
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
            print(f"Set m={m} ic={ic} to nkp= {csys_n} {new_k_lst} {csys_p}")
        psystem.models[m].values[0][ic]['n'].clear()
        psystem.models[m].values[0][ic]['n'].extend(csys_n)
        psystem.models[m].values[0][ic]['k'].clear()
        # new_k_lst = [x*angles[0] for x in new_k_lst]
        psystem.models[m].values[0][ic]['k'].extend(new_k_lst)
        psystem.models[m].values[0][ic]['p'].clear()
        psystem.models[m].values[0][ic]['p'].extend(csys_p)
    return psystem



def reset_unique(reset_config, gde, psystem, verbose=False):

    verbose=True
    assert type(reset_config) is dict

    reset_bond_k = reset_config.get("bond_k", True)
    reset_bond_l = reset_config.get("bond_l", True)
    reset_angle_k = reset_config.get("angle_k", True)
    reset_angle_l = reset_config.get("angle_l", True)
    reset_torsion_k = reset_config.get("torsion_k", True)
    reset_outofplane_k = reset_config.get("outofplane_k", True)

    guess_periodicity = reset_config.get("dihedral_p", False)
    max_n = reset_config.get("dihedral_max_n", 3)
    alpha = reset_config.get("dihedral_alpha", -.25)
    max_k = reset_config.get("dihedral_max_k", 10.0)
    min_k = reset_config.get("dihedral_min_k", 1e-3)

    if not any([
        reset_bond_k,
        reset_bond_l,
        reset_angle_k,
        reset_angle_l,
        reset_torsion_k
    ]):
        return psystem

    assert HESSIANS in gde.tables
    # eid_hess = [eid for eid, e in gdb.entries.items() if HESSIANS in e.tables]
    # psys_hess = [psystems[eid] for eid in eid_hess]

    pos = psystem.models[0].positions
    # reset = set()
    if reset_bond_l:
        lengths = assignments.graph_assignment_geometry_bond_matrix(pos)
        for ic, r in lengths.selections.items():
            r = r[0][0]
            l0 = psystem.models[0].values[0][ic]["l"][0]
            psystem.models[0].values[0][ic]["l"][0] = r
            if verbose:
                print(f"Set m=0 ic={ic} from l={l0} to l={r}")

    if reset_angle_l:
        lengths = assignments.graph_assignment_geometry_angles(pos)
        for ic, r in lengths.selections.items():
            r = r[0][0]
            l0 = psystem.models[1].values[0][ic]["l"][0]
            psystem.models[1].values[0][ic]["l"][0] = r
            if verbose:
                print(f"Set m=1 ic={ic} from l={l0} to l={r}")

    if reset_bond_k or reset_angle_k or reset_torsion_k or reset_outofplane_k:

        hessian = gde.tables[HESSIANS].values
        icb, B = assignments.bmatrix(pos, torsions=reset_torsion_k, outofplanes=reset_outofplane_k, pairs=False ,remove1_3=True)
        if B:
            gic = None
            if GRADIENTS in gde.tables:
                gx = arrays.array_scale(gde.tables[GRADIENTS].values, 1/4.184)
                gq = hessians.project_gradient(B, gx)
                sel = dict(zip(icb, [[x] for x in gq]))
                gic = assignments.graph_assignment(
                    pos.smiles,
                    sel,
                    pos.graph
                )

            labeled_ic = [x for model in psystem.models for x in model.labels[0]]
            B = {ic: x for ic, x in zip(icb, B) if ic in labeled_ic}
            B = list(B.keys()), list(B.values())


            hic = hessians.hessian_project_onto_ics(None, psystem, hessian, verbose=verbose, B=B)
            if reset_bond_k:
                for ic in graphs.graph_bonds(pos.graph):
                    if ic not in icb:
                        continue
                    k0 = psystem.models[0].values[0][ic]["k"][0]
                    k1 = hic.selections[ic]
                    psystem.models[0].values[0][ic]["k"][0] = k1
                    if verbose:
                        print(f"Set m=0 ic={ic} from k={k0} to k={k1}")
            if reset_angle_k:
                for ic in graphs.graph_angles(pos.graph):
                    if ic not in icb:
                        continue
                    k0 = psystem.models[1].values[0][ic]["k"][0]
                    k1 = hic.selections[ic]
                    nba = len([x for x in pos.graph.edges if ic[1] in x])
                    # psystem.models[1].values[0][ic]["k"][0] = k1 * 1.3 * nba / 3
                    # psystem.models[1].values[0][ic]["k"][0] = k1 * nba * 2/3
                    # psystem.models[1].values[0][ic]["k"][0] = k1 * nba  / 2
                    # SCALE UNIQUE
                    psystem.models[1].values[0][ic]["k"][0] = k1 / 2
                    if verbose:
                        print(f"Set m=1 ic={ic} from k={k0} to k={k1}")

            if reset_torsion_k:
                psystem = reset_unique_dihedral(
                    psystem, hic, 2, guess_periodicity, alpha, max_n, min_k, max_k, gradients=gic, verbose=verbose
                )

            if reset_outofplane_k:
                psystem = reset_unique_dihedral(
                    psystem, hic, 3, guess_periodicity, alpha, max_n, min_k, max_k, gradients=gic, verbose=verbose
                )

    return psystem

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
    eid_hess = [eid for eid, e in gdb.entries.items() if HESSIANS in e.tables]
    psys_hess = [psystems[eid] for eid in eid_hess]

    psys_hic = {}
    if ws:
        iterable = {
            eid: [[
                    csys,
                    psystems[eid],
                    gdb.entries[eid].tables[HESSIANS].values
                ],
                {"B": assignments.bmatrix(
                    psystems[eid].models[0].positions,
                    torsions=True,
                    outofplanes=False,
                    pairs=False,
                    remove1_3=True
                )}
            ] for eid in eid_hess
        }

        results = compute.workspace_submit_and_flush(
            ws,
            hessians.hessian_project_onto_ics,
            iterable,
            verbose=True
        )
        psys_hic_all.update(results)

    else:

        for i, eid in enumerate(eid_hess, 1):
            psys = psystems[eid]
            hessian = gdb.entries[eid].tables[HESSIANS].values
            if verbose:
                print(f"Projecting hessian for EID {eid} {i:8d}/{len(eid_hess)}")

            hic = hessians.hessian_project_onto_ics(
                csys,
                psys,
                hessian,
                verbose=verbose,
                B=assignments.bmatrix(
                    psystems[eid].models[0].positions,
                    torsions=True,
                    outofplanes=False,
                    pairs=False,
                    remove1_3=True
                )
            )
            psys_hic_all[eid] = hic

    for eid in eid_hess:
        psys = psystems[eid]
        hic = psys_hic_all[eid]
        if m == 2:
            psys_hics.append({
                k: [[hic.selections[k]]]
                for k in graphs.graph_torsions(hic.graph) if k in hic.selections
            })
        elif m == 3:
            psys_hics.append({
                k: [[hic.selections[k]]]
                for k in graphs.graph_outofplanes(hic.graph) if k in hic.selections
            })
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
                    psys_angles = assignments.smiles_assignment_geometry_torsions(psys.models[m].positions, indices=indices)
                elif m == 3:
                    psys_angles = assignments.smiles_assignment_geometry_outofplanes(psys.models[m].positions, indices=indices)
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
                    i2 = i + 2
                    ix = (i2) % 4
                    for n, p in zip(csys_n, csys_p):
                        x = [
                            sign[ix] * n**i2*deriv[ix](n*t - p)
                            for t in angles
                        ]
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
            tid = POSITIONS
            gid = list(gde.tables[tid].graphs)[0]
            rid = 0
            pos = assignments.graph_db_graph_to_graph_assignment(
                gdb,
                eid,
                tid,
                gid,
                rid
            )
            psys[eid] = mm.chemical_system_to_physical_system(
                csys,
                [pos],
                ref=psystems[eid],
                reuse=list(reuse)
            )
        psystems = psys

    return psystems


def generate_candidates(
    csys,
    psystems,
    gdb,
    G0,
    macro,
    strategy,
    union_cache,
    wq
):

    step = None
    gcd = csys.perception.gcd
    step_tracker = strategy.step_tracker
    icd = codecs.intvec_codec(
        gcd.primitive_codecs,
        gcd.atom_primitives,
        gcd.bond_primitives
    )
    iteration = 0
    # G0 = {i: icd.graph_encode(g) for i, g in gdb.graphs.items()}
    n_macro = len(strategy.steps)

    n_ics = 1
    N_str = "{:" + str(len(str(n_ics))) + "d}"

    candidates = {}

    while not optimization.optimization_iteration_is_done(macro):
        # t = datetime.datetime.now()
        # print(f"{t} Initializing new loop on macro {strategy.cursor}")

        step: optimization.optimization_step = (
            optimization.optimization_iteration_next(macro)
        )
        # step_tracker[tkey] = strategy.cursor

        n_micro = len(macro.steps)
        config: configs.smarts_perception_config = step.pcp
        S = step.cluster
        tkey = (S.category, S.name)
        hidx = mm.chemical_system_get_node_hierarchy(csys, S)
        topo = hidx.topology

        oper = step.operation

        if S.type == 'parameter' and S.index not in hidx.subgraphs:
            hidx.subgraphs[S.index] = gcd.smarts_decode(hidx.smarts[S.index])
        if type(hidx.subgraphs[S.index]) is str:
            # could not parse smarts (e.g. recursive) so we skip
            step_tracker[tkey][oper] = strategy.cursor
            continue
        S0 = graphs.subgraph_to_structure(hidx.subgraphs[S.index], topo)

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
        aa = []
        for eid, ps in psystems.items():
            gids = [*gdb.entries[eid].tables[assignments.POSITIONS].graphs]
            glbls = ps.models[int(S.category[0])].labels
            for gid, lbls in dict(zip(gids, glbls)).items():
                if gid not in selected_graphs:
                    for ic, term_lbls in lbls.items():
                        if type(ic[0]) is int:
                            ic = ic,
                        if S.name in term_lbls.values():
                            #TODO: this prunes all intermolecular pairs.
                            #We cannot generate [*:1].[*:2] splits on multigraph
                            #systems
                            ic = tuple((ici[1] for ici in ic))
                            selected_ics.append((gid, ic))
                            selected_graphs.add(gid)

        G = {k: v for k, v in G0.items() if k in selected_graphs}
        del selected_graphs

        aa = selected_ics
        assn_s = aa

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

        oper = step.operation

        if tkey in step_tracker:
            if step.index < step_tracker[tkey].get(oper, -1):
                continue

        if oper == strategy.SPLIT:
            new_candidates = {}
            new_candidates_direct = {}
            direct_success = False

            print(f"Attempting to split {S.name}:")
            s0split = graphs.structure_copy(S0)
            if config.splitter.primitives:
                graphs.graph_set_primitives_atom(s0split, config.splitter.primitives)
                graphs.graph_set_primitives_bond(s0split, config.splitter.primitives)
            print("S0:", gcd.smarts_encode(S0), "split_space:", gcd.smarts_encode(s0split))

            if not aa:
                print("No matches.")
                step_tracker[tkey][oper] = strategy.cursor
                continue

            print(f"Matched N={len(aa)}")
            seen = set()
            extend_config = config.extender.copy()
            extend_config.depth_max = config.splitter.branch_depth_limit
            extend_config.depth_min = config.splitter.branch_depth_min

            # For each node, I present just.. the chemical objective
            # until I can make a case for IC objective accounting

            for seen_i, i in enumerate(aa, 1):
                g = graphs.graph_to_structure(
                    icd.graph_decode(G[i[0]]),
                    i[1],
                    topo
                )
                graphs.structure_extend(extend_config, [g])
                g = graphs.structure_remove_unselected(g)
                gs = graphs.structure_copy(g)
                if config.splitter.primitives:
                    graphs.graph_set_primitives_atom(gs, config.splitter.primitives)
                    graphs.graph_set_primitives_bond(gs, config.splitter.primitives)

                if seen_i < 100:
                    print(
                        f"{seen_i:06d} {str(i):24s}",
                        # objective.report([x]),
                        gcd.smarts_encode(gs), "<",
                        gcd.smarts_encode(g),
                    )
                    seen.add(gs)
            print()
            if len(seen) < 2 and len(assn_s) < 100:
                print(f"Skipping {S.name} since all graphs are the same")
                step_tracker[tkey][oper] = strategy.cursor
                continue

            if len(seen) < 2 and len(assn_s) < 100:
                print(f"Skipping {S.name} since all graphs are the same")
                step_tracker[tkey][oper] = strategy.cursor
                continue

            if graphs.structure_max_depth(S0) > config.splitter.branch_depth_min:
                print("This parameter exceeds current depth. Skipping")
                step_tracker[tkey][oper] = strategy.cursor
                continue

            # this is where I generate all candidates
            if step.direct_enable and (
                    config.splitter.bit_search_limit < len(assn_s)
                ):
                assn_i = []
                if len(assn_s) < step.direct_limit:
                    # or form matches based on unique smarts
                    a = []
                    extend_config = config.extender.copy()
                    extend_config.depth_max = config.splitter.branch_depth_limit
                    extend_config.depth_min = config.splitter.branch_depth_min
                    for seen_i, (i, x) in enumerate(assn_s.items(), 1):
                        g = graphs.graph_to_structure(
                            icd.graph_decode(G[i[0]]),
                            i[1],
                            topo
                        )
                        graphs.structure_extend(extend_config, [g])
                        g = graphs.structure_remove_unselected(g)
                        a.append(g)

                    assn_i = list(range(len(a)))

                    pcp = step.pcp.copy()
                    pcp.extender = cfg
                    print("Direct splitting....")

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
                        matched_assn = tuple((assn_i[i] for i in matches))
                        unmatch_assn = tuple((assn_i[i] for i in unmatches))

                        new_candidates_direct[(step.overlap[0], None, p_j, oper)] = (
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

                (Q, new_candidates) = union_cache.get((S.index, S.name, strategy.cursor), (None, None))
                # here I need to get the graphs from the gdb
                # where aa refers to the indices
                if Q is None:
                    if len(aa) < 100:
                        Q = mapper.union_list_parallel(
                            G, aa, topo,
                            reference=S0,
                            max_depth=graphs.structure_max_depth(S0),
                            icd=icd
                        )
                    else:
                        Q = mapper.union_list_distributed(
                            G, aa, topo, wq,
                            reference=S0,
                            max_depth=graphs.structure_max_depth(S0),
                            icd=icd
                        )
                    union_cache[(S.index, S.name, strategy.cursor)] = (Q, new_candidates)
                else:
                    print(
                        f"{datetime.datetime.now()} Candidates N={len(new_candidates)} retreived from cache for node {S.index}:{S.name}"
                    )

                t = datetime.datetime.now()
                print(f"{t} Union is {gcd.smarts_encode(Q)}")

                if new_candidates is None:
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

                    print(
                        f"{datetime.datetime.now()} Collecting new candidates"
                    )
                    new_candidates = clusters.clustering_collect_split_candidates_serial(
                        S, ret, step, oper
                    )
                    for k in new_candidates:
                        v = list(new_candidates[k])
                        v[1] = graphs.subgraph_as_structure(
                            new_candidates[k][1],
                            topo
                        )
                        new_candidates[k] = tuple(v)
                    union_cache[(S.index, S.name, strategy.cursor)] = (Q, new_candidates)

            p_j_max = -1
            if candidates:
                p_j_max = max(x[2] for x in candidates) + 1
            for k, v in new_candidates.items():
                k = (k[0], k[1], k[2]+p_j_max, k[3])
                candidates[k] = v
            new_candidates = None

            p_j_max = -1
            if candidates:
                p_j_max = max(x[2] for x in candidates) + 1
            for k, v in new_candidates_direct.items():
                k = (k[0], k[1], k[2]+p_j_max, k[3])
                candidates[k] = v
            new_candidates_direct = None

        elif oper == strategy.MERGE:

            for p_j, jidx in enumerate(hidx.index.below[S.index]):
                J = hidx.index.nodes[jidx]
                if J.type != "parameter":
                    continue

                for overlap in step.overlap:
                    key = (overlap, macro.cursor, p_j, oper)
                    cnd = (S, J, step, None, None, None, None)
                    candidates[key] = cnd

        elif oper == strategy.MODIFY:

            if S.type == "parameter" and S.category[0] in [2, 3]:
                p_j = 0
                cm = csys.models[S.category[0]]
                s_per = cm.topology_terms['n'].values[S.name]
                f_max = 6
                if S.category[0] == 2:
                    f_max = step.modify_torsion_frequency_limit
                elif S.category[0] == 3:
                    f_max = step.modify_outofplane_frequency_limit
                print(f"Considering frequencies up to {f_max}")
                mod_max = step.pcp.splitter.bit_search_limit
                mod_min = step.pcp.splitter.bit_search_min

                s_add = tuple(sorted(set(range(1, f_max+1)).difference(s_per)))
                s_per = tuple(s_per)

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
                                key = (modify, macro.cursor, p_j, oper)
                                cnd = (S, S, step, None, None, None, None)
                                candidates[key] = cnd
                                p_j += 1
                                print(f"Adding modification {p_j} {S.name}: {modify}")
    return candidates, iteration


def process_tiers(
    tiers,
    candidates,
    gdb,
    csys,
    psystems,
    chemical_objective,
    P00,
    CX0,
    X0,
    C0,
    P0,
    Sj_sma,
    strategy,
    step,
    reset_config,
    wq
):

    n_ics = 1
    procs = (
        os.cpu_count() if configs.processors is None else configs.processors
    )
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
            (i, j): ((S, Sj, step.operation, edits, j), {"verbose": False})
            for i, (
                (edits, _, p_j, oper),
                (S, Sj, step, _, _, _, _),
            ) in enumerate(candidates.items(), 1)
            if mm.chemical_system_get_node_hierarchy(csys, S) is not None
        }
        print(
            datetime.datetime.now(),
            f"Generated {len(candidates)}",
            f"x {len(tiers[0].objectives)//len(j)}",
            f"= {len(iterable)} candidate evalulation tasks"
        )

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
                r = calc_tier_distributed(*v[0], **v[1], shm=shm)
                work[k] = r
        # elif len(tier.objectives) < 500 and len(candidates) > 1 and configs.remote_compute_enable:
        elif configs.processors > 1 or configs.remote_compute_enable:
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
                0, #math.log(len(iterable)//10, 10),
                len(iterable),
                verbose=True
            )
            compute.workqueue_remove_workspace(wq, ws)
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
                    ret = work[(i, j)]
                    line = ret.value
                    work_new[i][0] |= int(line[0])
                    work_new[i][1] += line[1]
                    work_new[i][2] = line[2]
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

            cC = chemical_objective(csys, P0=math.log(len(psystems)+1), c=c)
            cX = cP + (cC - C0)
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
            oper = cnd_keys[cnd_i][3]
            edits = cnd_keys[cnd_i][0]
            if edits:
                edits = str(edits)
            else:
                edits = f"{Sj_sma[cnd_i-1]}"
            cC = chemical_objective(csys, P0=math.log(len(psystems)+1), c=c)
            dP = cP - P0
            dC = cC - C0
            cX = cP + cC
            dX = dP + dC
            K = "Y" if j <= len(accepted_keys) else "N"
            F = ">" if j <= len(accepted_keys) else " "
            cout_line = (
                f"{F} Cnd. {cnd_i:4d}/{len(work)}"
                f" N= {match_len:6d} OP={oper:+d}"
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
    return candidates, Sj_sma


def process_accepted_candidates(
    candidates,
    cnd_keys,
    work,
    csys,
    strategy,
    Sj_sma,
    chemical_objective,
    C0,
    P0,
    X0,
    CX0,
    cscale,
    visited,
    ignore,
    kept
):

    cout_line = (
        f" Initial objectives: "
        f" X= {C0+P0:10.5f}"
        f" P= {P0:10.5f}"
        f" C= {C0:10.5f}"
    )
    print(cout_line)
    cout = {}
    cout_sorted_keys = []
    max_line = 0
    cnd_kv = {}
    for j, cnd_i in enumerate(sorted(work), 1):
        (keep, cP, c, match_len, kv, out) = work[cnd_i]
        # cnd_i, key, unit = unit
        (S, Sj, step, _, _, _, _) = candidates[cnd_keys[cnd_i]]

        dP = cP - P0
        # cC = C0 + dcC
        cC = chemical_objective(csys, P0=math.log(cscale+1), c=c)
        cX = cP + cC
        dP = cP - P0
        dC = cC - C0
        dX = dP + dC
        # dX = cX - X0
        keep = keep and dX <= 0.0

        oper = cnd_keys[cnd_i][3]
        if step.operation == strategy.SPLIT:
            visited.add((S.name, oper))
        elif step.operation == strategy.MERGE:
            visited.add((Sj.name, oper))
        elif step.operation == strategy.MODIFY:
            visited.add((S.name, oper))


        reused_line = ""
        K = "Y" if keep else "N"
        edits = cnd_keys[cnd_i][0]
        if edits:
            edits = str(edits)
        else:
            edits = f"{Sj_sma[cnd_i-1]}"
        cout_line = (
            f"Cnd. {cnd_i:4d}/{len(work)}"
            f" N= {match_len:6d} K= {K} OP={oper:+d}"
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

        # if this was valid but not accepted, we allow it to be
        # reconsidered if we repeat the step
        # if keep:
        #     step_tracker[(S.category, S.name)] = max(0, strategy.cursor - 1)

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
    return cout, cout_sorted_keys, cnd_kv


def select_best_accepted(
    cout,
    cout_sorted_keys,
    cnd_keys,
    strategy,
    X0,
    ignore,
    kept,
    macro_count,
):
    ck_i = 1

    n_added = sum(macro_count.values())
    cnd_keep = []
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

    return cnd_keep, micro_count


def insert_candidates(
    candidates,
    Sj_sma,
    strategy,
    step,
    step_tracker,
    csys,
    psystems,
    gdb,
    assigned_nodes,
    reuse0,
    reset_config,
    chemical_objective,
    initial_objective,
    tiers,
    C0,
    P0,
    X0,
    CX0,
    G0,
    union_cache,
    visited,
    wq
):

    visited = set()
    repeat = set()
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
    n_ics = 1
    procs = (
        os.cpu_count() if configs.processors is None else configs.processors
    )
    micro_added = 0
    success = False
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
            (i, j): ((S, Sj, step.operation, edits, j), dict(verbose=False))
            for i, (
                (edits, _, p_j, oper),
                (S, Sj, step, _, _, _, _),
            ) in enumerate(candidates.items(), 1)
            if mm.chemical_system_get_node_hierarchy(csys, S) is not None
            # ) in enumerate(candidates.items(), 1) for j in tiers[0].objectives
        }
        print(
            datetime.datetime.now(),
            f"Generated {len(candidates)} x "
            f"{len(initial_objective.objectives)//len(j)} = "
            f"{len(iterable)} candidate evalulation tasks"
        )

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

        if step.operation != strategy.MERGE and (macroamt + macroampc + microamt + microampc == 0):
            # use X0 so later dX will be 0 and kept
            # if we do this for merges, every merge will be taken..
            work = {i: (1, X0, 0.0, 1) for i in cnd_keys}

        else:
            if configs.processors == 1 and not configs.remote_compute_enable:
                work = {}
                print(logs.timestamp(), f"Dispatching candidate tasks= {len(iterable)} in serial")
                for i, (k, v) in enumerate(iterable.items()):
                    v[1]['verbose'] = True
                    print(logs.timestamp(), f"Running candidate task {i}/{len(iterable)}")
                    r = calc_tier_distributed(*v[0], **v[1], shm=shm)
                    work[k] = r
            elif configs.remote_compute_enable:
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
                compute.workqueue_remove_workspace(wq, ws)
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

        cout, cout_sorted_keys, cnd_kv = process_accepted_candidates(
            candidates,
            cnd_keys,
            work,
            csys,
            strategy,
            Sj_sma,
            chemical_objective,
            C0,
            P0,
            X0,
            CX0,
            len(psystems),
            visited,
            ignore,
            kept,
        )


        # print sorted at the end
        print(f"Nanostep {n_nano}: The filtered results of the candidate scan N={len(cout)} total={len(iterable)} oper={step.operation}:")
        if len(cout) == 0:
            continue

        cnd_keep, micro_count = select_best_accepted(
            cout,
            cout_sorted_keys,
            cnd_keys,
            strategy,
            X0,
            ignore,
            kept,
            macro_count,
        )
        micro_added = sum(micro_count.values())
        n_added += micro_added
        keys = {x[3]: cnd_keys[x[3]] for x in cnd_keep}

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
        print(f"Micro total: {sum(micro_count.values())}")

        print("Operations per parameter for this macro:")
        print(macro_count)
        print(f"Macro total: {sum(macro_count.values())}")

        if len(nodes) == 0:
            success = False
            added = False
            csys = csys_ref
            continue

        # print("Performing possibly requested reset")

        # ws = compute.workqueue_new_workspace(wq)
        # ret = reset(
        #     reset_config,
        #     csys,
        #     gdb,
        #     psystems,
        #     verbose=True,
        #     ws=ws
        # )
        # psystems = ret.value
        # print("\n".join(ret.out))
        # compute.workqueue_remove_workspace(wq, ws)
        # ws.close()
        # ws = None


        print("The resulting hierarchy is")
        mm.chemical_system_print(csys)

        success = False
        added = False
        if len(keys) == 1 and len(nodes) == 1:
            print("Only one modification, keeping result and printing output:")
            (keep, cP, c, match_len, kv, out) = work[cnd_keep[0][3]]
            print(f"\n".join(out))

            C = chemical_objective(csys, P0=math.log(len(psystems)+1), c=c)
            X = cnd_keep[0][1]
            P = X - C
            print(datetime.datetime.now(), f"Macro objective: P={P:13.6g} C={C:13.6g} DX={P+C-X0:13.6g}")

            kv = cnd_kv[cnd_keep[0]]
            kvgroups = {}
            for k, v in kv.items():
                if len(k) == 4:
                    if k[:3] not in kvgroups:
                        kvgroups[k[:3]] = {}
                    kvgroups[k[:3]][k[3]] = v

            # for k, v in list(kvgroups.items()):
            #     kvgroups[k] = [v[i] for i in sorted(v)]

            for k, vlst in kvgroups.items():
                vlst0 = mm.chemical_system_get_value_list(csys, k)
                for i, v in vlst.items():
                    v0 = 0
                    if i < len(vlst0):
                        v0 = vlst0[i]
                    print(f"{str(k):20s} | New: {v:12.6g} Ref {v0:12.6g} Diff {v-v0:12.6g}")
                mm.chemical_system_set_value_list(csys, k, [vlst[i] for i in sorted(vlst)])

            for k, v in kv.items():
                if len(k) == 3:
                    v0 = mm.chemical_system_get_value(csys, k, missing=0)
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
            S = list(nodes.values())[0]
            step_tracker[(S.category, S.name)] = {
                strategy.SPLIT: 0,
                strategy.MERGE: 0,
                strategy.MODIFY: 0
            }

        else:
            psysref = {
                i: mm.chemical_system_to_physical_system(
                    csys,
                    psystems[i].models[0].positions,
                    ref=psystems[i],
                    reuse=reuse0
                ) for i in psystems
            }

            if ws:
                compute.workqueue_remove_workspace(wq, ws)
                ws.close()
                ws = None
            ws = None
            if wq:
                ws = compute.workqueue_new_workspace(wq)
            psysref = reset(reset_config, csys, gdb, psysref, verbose=True, ws=ws).value
            if ws:
                compute.workqueue_remove_workspace(wq, ws)
                ws.close()
                ws = None
            print("Multiple modifications, doing another fit with all accepted*")

            reuse = [x for x in range(len(csys.models))]

            fitkeys = objective_tier_get_keys(initial_objective, csys)
            fitting_models = set((k[0] for k in fitkeys))
            fitkeys = [k for k in fitkeys if (k[1] in "skeler" and k[2] in assigned_nodes) or k[0] in fitting_models]
            reuse=[k for k,_ in enumerate(csys.models) if k not in fitting_models]
            
            kv0 = {k: mm.chemical_system_get_value(csys, k) for k in fitkeys}
            # for k, v in kv0.items():
            #     print(f"{str(k):20s} | v0 {v:12.6g}")
            mm.chemical_system_print(csys, show_parameters=assigned_nodes.union([x.name for x in nodes.values()]))
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
            C = chemical_objective(csys, P0=math.log(len(psystems)+1))
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

            hidx = mm.chemical_system_get_node_hierarchy(csys, S)
            if hidx is None:
                if S.name in repeat:
                    repeat.remove(S.name)
                if S.name in assigned_nodes:
                    assigned_nodes.remove(S.name)
                if (S.category, S.name) in step_tracker:
                    step_tracker[(S.category, S.name)] = dict.fromkeys(
                        step_tracker[(S.category, S.name)], 0
                    )
                if S.name in visited:
                    visited.remove(S.name)
                print("Warning, could not find node {S.name} in the hierarchy. Skipping")
                continue

            repeat.add(S.name)
            sma = Sj_sma[cnd_i-1]
            kept.add(cnd_i)
            # cnd_i = best[S.name][0] 
            # hent = Sj
            visited.add((hent.name, key[3]))
            edits = step.overlap[0]
            for union_idx in [(S.index, S.name), (hent.index, hent.name)]:
                for k in list(union_cache):
                    if union_idx == tuple(k[:2]):
                        print(f"Popping {k}")
                        union_cache.pop(k)

            if step.operation == strategy.SPLIT:

                success = True
                added = True
                repeat.add(hent.name)
                assigned_nodes.add(hent.name)
                step_tracker[(S.category, S.name)] = {
                    strategy.SPLIT: 0,
                    strategy.MERGE: 0,
                    strategy.MODIFY: 0
                }

                step_tracker[(hent.category, hent.name)] = {
                    strategy.SPLIT: 0,
                    strategy.MERGE: 0,
                    strategy.MODIFY: 0
                }

            elif step.operation == strategy.MERGE:

                if (hent.category, hent.name) in step_tracker:
                    step_tracker[(hent.category, hent.name)] = {
                        strategy.SPLIT: 0,
                        strategy.MERGE: 0,
                        strategy.MODIFY: 0
                    }
                else:
                    print("WARNING", hent.name, "missing from the tracker")

                step_tracker[(S.category, S.name)] = {
                    strategy.SPLIT: 0,
                    strategy.MERGE: 0,
                    strategy.MODIFY: 0
                }

                visited.add((S.category, S.name, key[3]))
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

                step_tracker[(S.category, S.name)] = {
                    strategy.SPLIT: 0,
                    strategy.MERGE: 0,
                    strategy.MODIFY: 0
                }
                repeat.add((S.category, S.name))
                visited.add((S.category, S.name, key[3]))
                success = True
                added = True

        print(datetime.datetime.now(), "Chemical system after nanostep:")
        # print the tree
        mm.chemical_system_print(csys, show_parameters=assigned_nodes)

        X0 = X
        P0 = P
        C0 = C
        psystems = psysref
    return repeat, visited, success, n_added, csys, psystems, X0, P0, C0


def finite_difference_forward_1(evals, h):

    a, b = evals
    dxdp = arrays.array_difference(b, a)
    dxdp = arrays.array_round(dxdp, PRECISION)
    dxdp = arrays.array_scale(dxdp, 1/h)
    d2xdp2 = list([0.0]*len(dxdp))
    return dxdp, d2xdp2


def finite_difference_central_2(evals, h):

    a, x, b = evals
    dxdp = arrays.array_difference(b, a)
    dxdp = arrays.array_round(dxdp, PRECISION)
    dxdp = arrays.array_scale(dxdp, 1/(2*h))
    d2xdp2 = arrays.array_add(
        arrays.array_scale(x, -2),
        arrays.array_add(b, a)
    )
    d2xdp2 = arrays.array_round(d2xdp2, PRECISION)
    d2xdp2 = arrays.array_scale(d2xdp2, 1.0/(h*h))

    return dxdp, d2xdp2
