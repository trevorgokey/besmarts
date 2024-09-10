"""
examples/11_aicp.py

Ab-Initio (valence) chemical perception (AICP)

This is a full valence fit of parameters starting from a base, empty set of
parameters. The initial parameters are first rapidly expanded by clustering
bonds, angles, and torsions based on Hessian projection, followed by a BESMARTS
parameter search which will refine the parameters to fit the objective. Here
we fit on the geometry, gradient, and frequencies of a single molecule. The
molecule is chemically diverse and has 25 atoms, so it presents a typical
challenge in force field design.

Estimated time to complete is 1 hour.
"""

import os
import pickle
import numpy as np
import base64
import tempfile

from besmarts.core import configs
from besmarts.core import perception
from besmarts.core import assignments
from besmarts.cluster import cluster_objective
from besmarts.cluster import cluster_optimization
from besmarts.mechanics import fits
from besmarts.mechanics import smirnoff_models

from besmarts.codecs import codec_rdkit
from besmarts.assign import hierarchy_assign_rdkit

configs.processors = 1
configs.remote_compute_enable = False
# To use additional workers, from a terminal run (on the separate work machine)
# python -m besmarts.worker <IP> 55555 np
# where ip is your publically visible IP where this is being run and np is the
# number of cores to assign to the worker. Recommended value for this example
# is 4-10 with a maximum number of 30 workers (due to tier acceptance).
configs.workqueue_port = 55555

# Show a little more output during the longer distributed compute
# to know it's doing something :)
# Currently will show an update every 10% progress or 60 seconds, whichever
# comes first.
configs.compute_verbosity = 1

# Known issue: if the command is abnormally terminated, some processes will
# not be killed and will not free the port. This will prevent the script from
# running again as it will try to bind the occupied port. To get around this
# (on Linux), run the following in a terminal
# pkill -f 11_aicp.py

# Boring constants yet to find a home. Need because gradient and hessian
# are in atomic units. Currently, we expect gradients to be in kJ/A, but
# Hessians in kcal/A/A. This will change in the future.
au2ang = 0.529177249
hess2kcal = 627.51 / (0.529177249**2)
au2kj = 627.51*4.184


def mystrategy():

    # The chemical perception strategy. Use a SMARTS search depth of up to 2
    # bits, and do not look at neighboring atoms.

    splitter = configs.smarts_splitter_config(
        1, 2, 0, 0, 0, 0, True, True, 0, False, True, True, True
    )

    # Disable general patterns like [!#6][!#8]. In general these will lead
    # to a hierarchy with fewer parameters as it greatly widens the search
    # space with almost no extra cost, but are harder to read and reason
    # with by SMARTS non-experts. We prefer specific splits here, like
    # [#6X4]-[#7]
    splitter.split_specific = True
    splitter.split_general = False
    extender = configs.smarts_extender_config(
        0, 0, True
    )
    cfg = configs.smarts_perception_config(splitter, extender)
    strategy = cluster_optimization.optimization_strategy_default(cfg)
    strategy.overlaps = [0]

    # Accept an unlimited number of parameters per iteration
    strategy.macro_accept_max_total = 0
    strategy.micro_accept_max_total = 0

    # But only accept at most one split per iteration
    strategy.macro_accept_max_per_cluster = 1
    strategy.micro_accept_max_per_cluster = 1

    return strategy


def pkl(self, fnm):
    with open(fnm, 'wb') as f:
        pickle.dump(self, f)


def pkl_load(fnm):
    with open(fnm, 'rb') as f:
        return pickle.load(f)


def load_xml(xml):

    gcd = codec_rdkit.graph_codec_rdkit()
    labeler = hierarchy_assign_rdkit.smarts_hierarchy_assignment_rdkit()
    pcp = perception.perception_model(gcd, labeler)

    fd, ff_fname = tempfile.mkstemp(suffix=".offxml")
    with open(fd, 'w') as f:
        f.write(xml)
    csys = smirnoff_models.smirnoff_load(ff_fname, pcp)
    os.remove(ff_fname)
    return csys


def configure_tiers(objs, penalties, fit_models, fit_symbols):

    # The initial fit to determine the reference
    initial = fits.objective_tier()
    final = fits.objective_tier()

    initial.objectives = objs

    # Don't perform any stochastic optimization
    initial.anneal = 0

    # line search in LBFGS
    initial.maxls = 20
    initial.fit_models = fit_models
    initial.fit_symbols = fit_symbols

    # Don't allow bonds and angles to have negative k values. Mostly due to the
    # fact that OpenMM refuses to take negative k values. Regardless, if any
    # ever get near 0, then the chemical model and/or objective config is
    # quite poor.
    initial.bounds[('k', 'b')] = (0, None)
    initial.bounds[('k', 'a')] = (0, None)
    initial.bounds[('k', 't')] = (None, None)
    initial.bounds[('k', 'i')] = (None, None)

    # The fits use parameter values in z-score space; these are the scale
    # factors, or previously know widths of the distribution. The parameters
    # in z-score space are (p-p0)/prior, where p0 is the expected value of
    # the parameter, but usually the initial value at the start of the
    # optimization.
    initial.priors = {
       "k": (None, 20),
       ("k", "b"): (None, 50),
       ("l", "b"): (None, .01),
       ("k", "a"): (None, 50),
       ("l", "a"): (None, .01),
       ("k", "t"): (None, 5.0),
       ("k", "i"): (None, 5.0),
       ("s", "s"): (None, .1),
       None: (None, None)
    }

    final.objectives = objs

    final.fit_models = fit_models
    final.fit_symbols = fit_symbols

    # There will be 10 (serial) iterations of full fits, each with a randomly
    # kicked set of parameters.
    final.anneal = 10
    final.maxls = 20

    final.bounds = initial.bounds
    final.priors = initial.priors

    tier = fits.objective_tier()
    tier.objectives = final.objectives

    # Only perform two optimization steps before scoring
    tier.step_limit = 2

    # Pass the best 5 candidates
    tier.accept = 5
    tier.maxls = 10
    tier.priors = final.priors
    tier.bounds = initial.bounds
    tier.fit_models = fit_models
    tier.fit_symbols = fit_symbols

    initial.penalties.extend(penalties)
    tier.penalties.extend(penalties)
    final.penalties.extend(penalties)

    tiers = [tier]

    return initial, tiers, final


def configure_chem_perception_strat(csys, split_models):

    strat = fits.forcefield_optimization_strategy_default(
        csys, models=split_models)

    strat.enable_split = 1
    strat.enable_merge = 1

    # Don't scan periodicity modifications since we estimate them each fit
    strat.enable_modify = 0

    # When scoring, reset the torsions
    strat.enable_dihedral_periodicity_reset = 1
    strat.dihedral_periodicity_reset_alpha = -.5
    strat.dihedral_periodicity_reset_max_n = 3
    strat.dihedral_periodicity_reset_max_k = 5
    strat.dihedral_periodicity_reset_min_k = 1e-3

    # For each iteration, accept only 1 split
    strat.macro_accept_max_total = 1

    # Reset torsions but keep bonds and angles during parameter search
    strat.enable_reset_bond_lengths = False
    strat.enable_reset_angle_lengths = False
    strat.enable_reset_bond_stiffness = False
    strat.enable_reset_angle_stiffness = False
    strat.enable_reset_torsion_stiffness = True
    strat.enable_reset_outofplane_stiffness = True

    # The default is to merge then split, we want split then merge
    strat.build_steps()
    strat.steps = list(reversed(strat.steps))

    return strat


def configure_objectives(opt_mols):

    # Configure the geometry, gradient, and frequency objectives. The scales
    # were determined by examining the first step, and scaling such
    # that they are relatively similar in scale. In this example, we put
    # slightly more emphasis on frequencies.
    objs = {}

    # Note that only one minimization step is performed. This is the slowest
    # part of each iteration so increasing this will cause a marked reduction
    # in speed.

    for eid in opt_mols:
        o = fits.objective_config_hessian(
                assignments.graph_db_address(
                    eid=[eid],
                ),
                scale=1e-4,
                include=True
        )
        o.batch_size = 1
        o.method_vib_modes = "qm"
        o.verbose = 2
        objs[len(objs)] = o

    for eid in opt_mols:
        o = fits.objective_config_position(
                assignments.graph_db_address(
                    eid=[eid],
                ),
                scale=1
        )
        o.step_limit = 1
        o.batch_size = 1
        o.tol = 1e-4
        o.verbose = 2
        objs[len(objs)] = o

    for eid in opt_mols:
        o = fits.objective_config_gradient(
            assignments.graph_db_address(
                eid=[eid],
            ),
            scale=1e-7,
            include=True
        )
        o.verbose = 2
        objs[len(objs)] = o

    return objs


def configure_penalties():

    # Examples of penalties. Since the parameters are reset each time,
    # we forego them.
    penalties = []

    # restrain equilibrium lengths of bonds and angles to the starting
    # values at the beginning of a fit
    # penalties.append(fits.objective_config_penalty(
    #     keys={
    #         (0, 'l', None, None): None,
    #         (1, 'l', None, None): None,
    #     },
    #     polynomial={1: 1.0},
    #     scale=1000.0,
    # ))
    # penalties.append(fits.objective_config_penalty(
    #     keys={
    #         (0, 'k', None, None): None,
    #         (1, 'k', None, None): None,
    #     },
    #     polynomial={2: 1.0},
    #     scale=0.01,
    # ))

    # restrain torsion k to 0
    # penalties.append(fits.objective_config_penalty(
    #     keys={
    #         (2, 'k', None, None): 0.0, # constant
    #         (3, 'k', None, None): 0.0
    #     },
    #     polynomial={2: 1.0},
    #     scale=0.1,
    # ))

    return penalties


def configure_molecule(csys, gdb: assignments.graph_db):

    # Load the molecule data and add to the graph_db object
    gcd = csys.perception.gcd

    smi = "[C:2]([S:6]1)([H:1])=[C:3]([C:16]([H:17])([H:18])[C:19](=[O:20])[O:21][C:22]([H:23])([H:24])[H:25])[N:4]=[C:5]1[N:7]([H:8])[S:9](=[O:10])(=[O:11])[C:12]([H:13])([H:14])[H:15]" 
    xyzdata, grad, hess = load_data()

    pos = assignments.xyz_to_graph_assignment(gcd, smi, xyzdata)

    grad *= au2kj/au2ang
    sel = {(i,): [x] for i, x in enumerate(grad, 1)}

    gx = assignments.graph_assignment(smi, sel, pos.graph)

    hess *= hess2kcal
    hx = hess.tolist()

    energy = 0

    eid, gid = assignments.graph_db_add_single_molecule_state(
        gdb,
        pos,
        gradients=gx,
        hessian=hx,
        energy=energy
    )
    return gdb


def main():

    prefix = "./11_aicp/"
    if not os.path.exists(prefix):
        os.mkdir(prefix)

    # build the dataset and input ff
    csys = load_xml(ai_offxml)

    # load the molecule
    gdb = assignments.graph_db()
    gdb = configure_molecule(csys, gdb)

    # Create the objectives for molecule. Here there will be 3 objectives:
    # geometry, gradients, and frequencies
    opt_mols = list(gdb.entries)
    objs = configure_objectives(opt_mols)

    # Search for parameters on bonds (0), angles (1), and torsions (2)
    split_models = [0, 1, 2]

    # Score the candidates on fit_models 0, 1, 2, i.e. these models are part
    # of the optimizations
    fit_models = [0, 1, 2]

    # These terms are fit (length and k). The term l is similar for bonds and
    # angles and represents to equilibrium value, and k is similar for all
    # bonds, angles, and torsions. The result is we fit to all bonds and
    # angles and all k values of torsions
    fit_symbols = ["l", "k"]

    # The penalties control the amount of objective increase as the parameters
    # deviate away from some initial value. In this example, there are no
    # penalties and so the parameters are free take on any value that minimizes
    # the molecular objective. Rather than use a penalty, we reset all
    # parameter values back to what the Hessian projection should be. This is a
    # softer way to control parameter drift, and if we are going to be stuck in
    # a local minima on the objective function surface, we would choose the
    # minima that is "close" to what the Hessian predicts.
    penalties = configure_penalties()

    # The initial tier will be run before any parameter search to create a
    # reference score. The final will be run when no further modifications of
    # the parameters were found.
    # The difference between initial and final is that the final fit includes
    # 10 steps of stochastic optimization using basinhopping. This is
    # controlled by the final.anneal member (set to 10).
    # The single tier will run two steps of parameter fitting, and then pass
    # the 30 best to the final scorer. Here, the initial tier will take the
    # final candidates and fit them before incorporating the best and repeating
    initial, tiers, final = configure_tiers(
        objs,
        penalties,
        fit_models,
        fit_symbols
    )

    # The function that penalizes complex SMARTS hierarchies. This objective
    # increases as more parameters are added
    co = fits.chemical_objective

    # Estimate the force constants for all bonds, angles, and torsions.
    # Furthermore, predict the periodicities up to n=6 and cap any k val to 
    # 100. Periodicities are preferred where all cos(q) < alpha or
    # cos(q) > -alpha, indicating all angles are in troughs of the
    # periodicitiy. We choose 100 to highlight distinct chemistries, but in 
    # the actual predicted values we cap to 5 kcal.
    sag_map = fits.smiles_assignment_force_constants(
        gdb,
        alpha=-.5,
        max_n=6,
        max_dihedral_k=100
    )

    # Cluster bond lengths where a split is accepted if the mean difference
    # is greater than 0.025 A.
    if os.path.exists(prefix + "csys.bond_l.p"):
        csys = pkl_load(prefix + "csys.bond_l.p")
    else:
        sag = assignments.smiles_assignment_group_bonds(sag_map["bond_l"])
        sep = 0.025
        obj = cluster_objective.clustering_objective_mean_separation(sep, sep)
        csys = fits.chemical_system_cluster_data(
            csys,
            0,
            sag,
            obj,
            strategy=mystrategy()
        )
        smirnoff_models.smirnoff_write_version_0p3(
            csys,
            prefix + "out.bond_l.offxml"
        )
        pkl(csys, prefix + "csys.bond_l.p")

    # Cluster bond force constants where a split is accepted if the mean
    # difference is greater than 10 kcal/A/A.
    if os.path.exists(prefix + "csys.bond_k.p"):
        csys = pkl_load(prefix + "csys.bond_k.p")
    else:
        sag = assignments.smiles_assignment_group_bonds(sag_map["bond_k"])
        sep = 10.0
        obj = cluster_objective.clustering_objective_mean_separation(sep, sep)
        csys = fits.chemical_system_cluster_data(
            csys,
            0,
            sag,
            obj,
            strategy=mystrategy())
        smirnoff_models.smirnoff_write_version_0p3(
            csys,
            prefix + "out.bond_k.offxml"
        )
        pkl(csys, prefix + "csys.bond_k.p")

    # Cluster angles where a split is accepted if the mean difference is
    # greater than 0.1 radians (5.7 degrees)
    if os.path.exists(prefix + "csys.angle_l.p"):
        csys = pkl_load(prefix + "csys.angle_l.p")
    else:
        sag = assignments.smiles_assignment_group_angles(sag_map["angle_l"])
        sep = 0.1
        obj = cluster_objective.clustering_objective_mean_separation(sep, sep)
        csys = fits.chemical_system_cluster_data(
            csys,
            1,
            sag,
            obj,
            strategy=mystrategy()
        )
        smirnoff_models.smirnoff_write_version_0p3(
            csys,
            prefix + "out.angle_l.offxml"
        )
        pkl(csys, prefix + "csys.angle_l.p")

    # Cluster angle force constants where a split is accepted if the mean
    # difference is greater than 5.0 kcal/rad/rad.
    if os.path.exists(prefix + "csys.angle_k.p"):
        csys = pkl_load(prefix + "csys.angle_k.p")
    else:
        sag = assignments.smiles_assignment_group_angles(sag_map["angle_k"])
        sep = 5.0
        obj = cluster_objective.clustering_objective_mean_separation(sep, sep)
        csys = fits.chemical_system_cluster_data(
            csys,
            1,
            sag,
            obj,
            strategy=mystrategy()
        )
        smirnoff_models.smirnoff_write_version_0p3(
            csys, prefix + "out.angle_k.offxml"
        )
        pkl(csys, prefix + "csys.angle_k.p")

    # Split on predicted cosine series. Each torsion will have a set of
    # periodicities assigned. We then build a descrete label for the given
    # set, and try to cluster torsions that have the same label, or the same
    # set of predicted cosine terms.
    if os.path.exists(prefix + "csys.torsion_p.p"):
        csys = pkl_load(prefix + "csys.torsion_p.p")
    else:
        assn = []
        for sag_p, sag_n in zip(sag_map["torsion_p"], sag_map["torsion_n"]):
            sel = {}
            for ic, plist in sag_p.selections.items():
                nlist = sag_n.selections[ic]
                pstr = []
                for ni in range(1, 7):
                    if ni in nlist:
                        idx = nlist.index(ni)
                        p = plist[idx]
                        if p == 0.0:
                            pstr.append("0")
                        elif abs((p - 3.14159)) < .0001:
                            pstr.append("1")
                        else:
                            pstr.append("2")
                    else:
                        pstr.append("x")
                sel[ic] = tuple(["".join(pstr)])
            ga = assignments.graph_assignment(sag_p.smiles, sel, sag_p.graph)
            assn.append(ga)

        sag = assignments.smiles_assignment_group_torsions(assn)
        obj = cluster_objective.clustering_objective_classification()

        strategy = mystrategy()
        strategy.overlaps = [0]
        csys = fits.chemical_system_cluster_data(
            csys,
            2,
            sag,
            obj,
            strategy=strategy
        )
        smirnoff_models.smirnoff_write_version_0p3(
            csys,
            prefix + "out.torsion_p.offxml"
        )
        fits.print_chemical_system(csys)
        pkl(csys, prefix + "csys.torsion_p.p")

    # Now that torsions are grouped as well as possible based on periodicities,
    # split further based on predicted k values. In this case, we translate
    # phases to zero, flipping signs of k when appropriate (when phase is pi).
    # We then split if the maximum difference in k values for any cosine term
    # is greater than 0.5.
    if os.path.exists(prefix + "csys.torsion_k.p"):
        csys = pkl_load(prefix + "csys.torsion_k.p")
    else:
        assn = []
        for sag_p, sag_n, sag_k in zip(
                sag_map["torsion_p"],
                sag_map["torsion_n"],
                sag_map["torsion_k"]
            ):
            sel = {}
            for ic, plist in sag_p.selections.items():
                nlist = sag_n.selections[ic]
                klist = sag_k.selections[ic]
                k_vals = []
                for ni in range(1, 7):
                    if ni in nlist:
                        idx = nlist.index(ni)
                        k = klist[idx]
                        if abs((p - 3.14159)) < .0001:
                            k_vals.append(-k)
                        else:
                            k_vals.append(k)
                    else:
                        k_vals.append(0.0)
                sel[ic] = tuple(k_vals)
            ga = assignments.graph_assignment(sag_p.smiles, sel, sag_p.graph)
            assn.append(ga)

        sag = assignments.smiles_assignment_group_torsions(assn)
        sep = .5
        obj = cluster_objective.clustering_objective_mean_separation(sep, sep)
        strategy = mystrategy()
        strategy.overlaps = [0]
        strategy.bounds.splitter.split_general = False
        csys = fits.chemical_system_cluster_data(
            csys,
            2,
            sag,
            obj,
            strategy=strategy
        )
        smirnoff_models.smirnoff_write_version_0p3(
            csys,
            prefix + "out.torsion_k.offxml"
        )
        fits.print_chemical_system(csys)
        pkl(csys, prefix + "csys.torsion_k.p")

    # parameterize all molecules
    psys = fits.gdb_to_physical_systems(gdb, csys)

    # Initial chemical perception is done. Perform a full reset of values
    # on each parameter. In this example, we split torsions using 6 cosine
    # terms, but for the actual force field we limit cosine periodicity to 3.
    reset_config = {
        "bond_l": True,
        "bond_k": True,
        "angle_l": True,
        "angle_k": True,
        "torsion_k": True,
        "outofplane_k": False,
        "dihedral_p": True,
        "dihedral_max_n": 3,
        "dihedral_alpha": -0.5,
        "dihedral_min_k": 1e-3,
        "dihedral_max_k": 5,
    }
    ret = fits.reset(reset_config, csys, gdb, psystems=psys, verbose=True)
    psys = ret.value
    print("\n".join(ret.out))

    print("After resetting and initializing, the parameters are:")
    fits.print_chemical_system(csys)
    smirnoff_models.smirnoff_write_version_0p3(csys, prefix + "out.offxml")

    # Configure the strategy for the chemical perception
    strat = configure_chem_perception_strat(csys, split_models)

    # Do the full parameter search using the objectives defined earlier. This
    # part is the most expensive and would benefit from additional workers.
    newcsys, (P0, P), (C0, C) = fits.ff_optimize(
        csys,
        gdb,
        psys,
        strat,
        co,
        initial,
        tiers,
        final
    )
    print("Initial objectives:")
    X0 = P0 + C0
    X = P + C
    print(f"Total= {X0:15.8g} Physical {P0:15.8g} Chemical {C0:15.8g}")
    print("Final objectives:")
    print(f"Total= {X:15.8g} Physical {P:15.8g} Chemical {C:15.8g}")
    print("Differences:")
    print(
        f"Total= {100*(X-X0)/X0:14.2f}%",
        f"Physical {100*(P-P0)/P0:14.2f}%",
        f"Chemical {100*(C-C0)/C0:14.2f}%"
    )

    smirnoff_models.smirnoff_write_version_0p3(
        newcsys,
        prefix + "final.offxml"
    )
    print("wrote final.offxml")


def load_data():

    global xyz_file
    global grad_file
    global hess_bytes

    hx_b = base64.decodebytes(hess_bytes)

    fd, hx_fname = tempfile.mkstemp(suffix=".npy")
    with open(fd, 'wb') as f:
        f.write(hx_b)

    cx = xyz_file

    fd, gx_fname = tempfile.mkstemp(suffix=".txt")
    with open(fd, 'w') as f:
        f.write(grad_file)

    gx = np.loadtxt(gx_fname)
    hx = np.load(hx_fname)

    os.remove(hx_fname)
    os.remove(gx_fname)

    return cx, gx, hx

# Below this is just data: the force field, and the xyz, grad, and hess of
# the molecule to fit (used by load_data)


ai_offxml = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Bonds
        version="0.4"
        potential="harmonic"
        fractional_bondorder_method="AM1-Wiberg"
        fractional_bondorder_interpolation="linear"
    >
        <Bond
            id="b1"
            smirks="[*:1]~[*:2]"
            length="1.5 * angstrom"
            k="600.0 * angstrom**-2 * mole**-1 * kilocalorie"
        ></Bond>
        <Bond
            id="b2"
            smirks="[*:1]-[#1:2]"
            length="1.0 * angstrom" k="800.0 *
            angstrom**-2 * mole**-1 * kilocalorie"
        ></Bond>
    </Bonds>
    <Angles version="0.3" potential="harmonic">
        <Angle
            id="a1"
            smirks="[*:1]~[X2:2]~[*:3]"
            angle="180.0 * degree"
            k="0.0 * mole**-1 * radian**-2 * kilocalorie"
        ></Angle>
        <Angle
            id="a2"
            smirks="[*:1]~[X3:2]~[*:3]"
            angle="120.0 * degree"
            k="120.0 * mole**-1 * radian**-2 * kilocalorie"
        ></Angle>
        <Angle
            id="a3"
            smirks="[*:1]~[X4:2]~[*:3]"
            angle="109.5 * degree"
            k="120.0 * mole**-1 * radian**-2 * kilocalorie"
        ></Angle>
    </Angles>
    <ProperTorsions
        version="0.4"
        potential="k*(1+cos(periodicity*theta-phase))"
        default_idivf="auto"
        fractional_bondorder_method="AM1-Wiberg"
        fractional_bondorder_interpolation="linear"
    >
        <Proper
            id="t1"
            smirks="[*:1]~[X4:2]~[X4:3]~[*:4]"
            periodicity1="1"
            phase1="0.0 * degree"
            k1="0.0 * mole**-1 * kilocalorie"
            idivf1="1.0"
        ></Proper>
        <Proper
            id="t2"
            smirks="[*:1]~[X4:2]~[X3:3]~[*:4]"
            periodicity1="1"
            phase1="0.0 * degree"
            k1="0.0 * mole**-1 * kilocalorie"
            idivf1="1.0"
        ></Proper>
        <Proper
            id="t3"
            smirks="[*:1]~[X4:2]~[X2:3]~[*:4]"
            periodicity1="1"
            phase1="0.0 * degree"
            k1="0.0 * mole**-1 * kilocalorie"
            idivf1="1.0"
        ></Proper>
        <Proper
            id="t4"
            smirks="[*:1]~[X3:2]~[X3:3]~[*:4]"
            periodicity1="1"
            phase1="0.0 * degree"
            k1="0.0 * mole**-1 * kilocalorie"
            idivf1="1.0"
        ></Proper>
        <Proper
            id="t5"
            smirks="[*:1]~[X3:2]~[X2:3]~[*:4]"
            periodicity1="1"
            phase1="0.0 * degree"
            k1="0.0 * mole**-1 * kilocalorie"
            idivf1="1.0"
        ></Proper>
        <Proper
            id="t6"
            smirks="[*:1]~[X2:2]~[X2:3]~[*:4]"
            periodicity1="1"
            phase1="0.0 * degree"
            k1="0.0 * mole**-1 * kilocalorie"
            idivf1="1.0"
        ></Proper>
    </ProperTorsions>
    <ImproperTorsions
        version="0.3"
        potential="k*(1+cos(periodicity*theta-phase))"
        default_idivf="auto"
    ></ImproperTorsions>
    <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1.0" cutoff="9.0 * angstrom" switch_width="1.0 * angstrom" method="cutoff">
        <Atom smirks="[#1:1]" epsilon="0.0157 * mole**-1 * kilocalorie" id="n1" rmin_half="0.6 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X4]" epsilon="0.01577948280971 * mole**-1 * kilocalorie" id="n2" rmin_half="1.48419980825 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X4]-[#7,#8,#9,#16,#17,#35]" epsilon="0.01640924602775 * mole**-1 * kilocalorie" id="n3" rmin_half="1.449786411317 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]" epsilon="0.0157 * mole**-1 * kilocalorie" id="n4" rmin_half="1.287 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])(-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]" epsilon="0.0157 * mole**-1 * kilocalorie" id="n5" rmin_half="1.187 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X4]~[*+1,*+2]" epsilon="0.0157 * mole**-1 * kilocalorie" id="n6" rmin_half="1.1 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X3]" epsilon="0.01561134320353 * mole**-1 * kilocalorie" id="n7" rmin_half="1.443812569645 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X3]~[#7,#8,#9,#16,#17,#35]" epsilon="0.01310699839698 * mole**-1 * kilocalorie" id="n8" rmin_half="1.377051329051 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X3](~[#7,#8,#9,#16,#17,#35])~[#7,#8,#9,#16,#17,#35]" epsilon="0.01479744504464 * mole**-1 * kilocalorie" id="n9" rmin_half="1.370482808197 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#6X2]" epsilon="0.015 * mole**-1 * kilocalorie" id="n10" rmin_half="1.459 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#7]" epsilon="0.01409081474669 * mole**-1 * kilocalorie" id="n11" rmin_half="0.6192778454102 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#8]" epsilon="1.232599966667e-05 * mole**-1 * kilocalorie" id="n12" rmin_half="0.2999999999997 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#16]" epsilon="0.0157 * mole**-1 * kilocalorie" id="n13" rmin_half="0.6 * angstrom"></Atom>
        <Atom smirks="[#6:1]" epsilon="0.0868793154488 * mole**-1 * kilocalorie" id="n14" rmin_half="1.953447017081 * angstrom"></Atom>
        <Atom smirks="[#6X2:1]" epsilon="0.21 * mole**-1 * kilocalorie" id="n15" rmin_half="1.908 * angstrom"></Atom>
        <Atom smirks="[#6X4:1]" epsilon="0.1088406109251 * mole**-1 * kilocalorie" id="n16" rmin_half="1.896698071741 * angstrom"></Atom>
        <Atom smirks="[#8:1]" epsilon="0.2102061007896 * mole**-1 * kilocalorie" id="n17" rmin_half="1.706036917087 * angstrom"></Atom>
        <Atom smirks="[#8X2H0+0:1]" epsilon="0.1684651402602 * mole**-1 * kilocalorie" id="n18" rmin_half="1.697783613804 * angstrom"></Atom>
        <Atom smirks="[#8X2H1+0:1]" epsilon="0.2094735324129 * mole**-1 * kilocalorie" id="n19" rmin_half="1.682099169199 * angstrom"></Atom>
        <Atom smirks="[#7:1]" epsilon="0.1676915150424 * mole**-1 * kilocalorie" id="n20" rmin_half="1.799798315098 * angstrom"></Atom>
        <Atom smirks="[#16:1]" epsilon="0.25 * mole**-1 * kilocalorie" id="n21" rmin_half="2.0 * angstrom"></Atom>
        <Atom smirks="[#15:1]" epsilon="0.2 * mole**-1 * kilocalorie" id="n22" rmin_half="2.1 * angstrom"></Atom>
        <Atom smirks="[#9:1]" epsilon="0.061 * mole**-1 * kilocalorie" id="n23" rmin_half="1.75 * angstrom"></Atom>
        <Atom smirks="[#17:1]" epsilon="0.2656001046527 * mole**-1 * kilocalorie" id="n24" rmin_half="1.85628721824 * angstrom"></Atom>
        <Atom smirks="[#35:1]" epsilon="0.3218986365974 * mole**-1 * kilocalorie" id="n25" rmin_half="1.969806594135 * angstrom"></Atom>
        <Atom smirks="[#53:1]" epsilon="0.4 * mole**-1 * kilocalorie" id="n26" rmin_half="2.35 * angstrom"></Atom>
        <Atom smirks="[#3+1:1]" epsilon="0.0279896 * mole**-1 * kilocalorie" id="n27" rmin_half="1.025 * angstrom"></Atom>
        <Atom smirks="[#11+1:1]" epsilon="0.0874393 * mole**-1 * kilocalorie" id="n28" rmin_half="1.369 * angstrom"></Atom>
        <Atom smirks="[#19+1:1]" epsilon="0.1936829 * mole**-1 * kilocalorie" id="n29" rmin_half="1.705 * angstrom"></Atom>
        <Atom smirks="[#37+1:1]" epsilon="0.3278219 * mole**-1 * kilocalorie" id="n30" rmin_half="1.813 * angstrom"></Atom>
        <Atom smirks="[#55+1:1]" epsilon="0.4065394 * mole**-1 * kilocalorie" id="n31" rmin_half="1.976 * angstrom"></Atom>
        <Atom smirks="[#9X0-1:1]" epsilon="0.003364 * mole**-1 * kilocalorie" id="n32" rmin_half="2.303 * angstrom"></Atom>
        <Atom smirks="[#17X0-1:1]" epsilon="0.035591 * mole**-1 * kilocalorie" id="n33" rmin_half="2.513 * angstrom"></Atom>
        <Atom smirks="[#35X0-1:1]" epsilon="0.0586554 * mole**-1 * kilocalorie" id="n34" rmin_half="2.608 * angstrom"></Atom>
        <Atom smirks="[#53X0-1:1]" epsilon="0.0536816 * mole**-1 * kilocalorie" id="n35" rmin_half="2.86 * angstrom"></Atom>
        <Atom smirks="[#1]-[#8X2H2+0:1]-[#1]" epsilon="0.1521 * mole**-1 * kilocalorie" id="n-tip3p-O" sigma="3.1507 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#8X2H2+0]-[#1]" epsilon="0 * mole**-1 * kilocalorie" id="n-tip3p-H" sigma="1 * angstrom"></Atom>
    </vdW>
    <Electrostatics
        version="0.3"
        scale12="0.0"
        scale13="0.0"
        scale14="0.8333333333"
        scale15="1.0"
        cutoff="9.0 * angstrom"
        switch_width="0.0 * angstrom"
        method="PME"
    ></Electrostatics>
    <LibraryCharges version="0.3">
    </LibraryCharges>
    <ToolkitAM1BCC version="0.3"></ToolkitAM1BCC>
</SMIRNOFF>
"""

xyz_file = """25
-1479.5743077535349
H   4.386657679753 2.296618590076 -5.383224800225
C   4.729216447807 1.358560040748 -4.967362254153
C   4.029588688776 0.447907408265 -4.224512040660
N   4.734845554933 -0.685034097617 -3.866750065397
C   5.961935922526 -0.632234373951 -4.298094321983
S   6.369091824126 0.820007250917 -5.197666646292
N   6.890556075016 -1.638770307391 -4.021042640243
H   7.802034379634 -1.560094928194 -4.468145060261
S   7.122081137003 -2.104909864238 -2.361176267778
O   8.380049769967 -2.853115266701 -2.398753079954
O   6.950479534920 -0.952216676150 -1.479011246763
C   5.734910370415 -3.220017390725 -2.124168148822
H   5.799832197553 -3.568485218546 -1.091581069911
H   4.822878265265 -2.649514601226 -2.301120300296
H   5.841650929028 -4.045777301286 -2.827136505585
C   2.612552551867 0.606749494723 -3.740778752355
H   2.113297736036 1.409805387624 -4.283654522125
H   2.070045618044 -0.331286809497 -3.882060462234
C   2.622499947537 0.981035144203 -2.263410653530
O   2.598724673502 2.121323447016 -1.844837975147
O   2.703142141389 -0.114290361616 -1.479461296283
C   2.868890453131 0.147529775256 -0.068059793352
H   2.901627072451 -0.833013290244 0.406127231341
H   2.030097178299 0.734382963010 0.312405184624
H   3.802267718594 0.687846390666 0.103502459240"""

grad_file = """
8.867688338395541170e-07 -3.754424730134653387e-08 1.667365844800210335e-06
1.150239820518479499e-06 -1.308969444553003455e-08 1.470935122667991254e-06
2.495585186899733354e-07 -4.529787465364695093e-07 1.524777688498269411e-06
1.620139719383552335e-07 -3.753862358881572200e-07 1.369112753440244676e-06
6.140865947659897295e-07 -2.143645095951506735e-07 2.134195407769816366e-06
9.080181438010084469e-07 -7.061938836320797580e-08 2.385068961739786504e-06
3.506062477264806562e-07 -2.909669887114470568e-07 2.221312183989384518e-06
5.616629453866262868e-07 -1.665920658582389132e-07 2.694406361033254175e-06
-4.504623955714145978e-07 -5.867181296395325776e-07 2.317091295946874935e-06
-4.084128442001511150e-07 -4.661122166230841340e-07 2.767559324507319184e-06
-7.989899498902250201e-07 -6.284681072744272706e-07 2.432283707586198329e-06
-4.557806281724359736e-07 -5.481652876838851579e-07 1.483587498316745260e-06
-9.813373544210194288e-07 -7.163182901180826286e-07 1.447570345227577231e-06
-4.111560588414340904e-07 -5.490364127491969420e-07 1.166405269928410702e-06
-1.725702603907184642e-07 -4.608892141464131925e-07 1.423565696622218848e-06
1.898413269678855472e-07 -3.329544413262756518e-07 5.805128472328318929e-07
4.604483124045393241e-07 -2.006463122519860573e-07 4.414001722856804007e-07
2.955562622423599760e-07 -4.039544570922352527e-07 2.037502129408671342e-07
-5.051658064095460379e-07 -7.291781050373772645e-07 7.044626729575350667e-07
-6.545498960071447110e-07 -4.222640180806492768e-07 8.397893402255977684e-07
-7.677721197696927603e-07 -7.475459855236699547e-07 6.669102585248835174e-07
-1.613363045728001283e-06 -9.382742307804247002e-07 6.337329958980199990e-07
-1.681759256652367988e-06 -9.558845554081618970e-07 5.249236686440767949e-07
-1.646404029111806790e-06 -9.800791430232650174e-07 3.797825569016081114e-07
-1.577956875999597474e-06 -9.337138379907856185e-07 1.146424160262189858e-06"""

hess_bytes = \
    b'k05VTVBZAQB2AHsnZGVzY3InOiAnPGY4JywgJ2ZvcnRyYW5fb3JkZXInOiBGYWxz'\
    b'ZSwgJ3NoYXBl\nJzogKDc1LCA3NSksIH0gICAgICAgICAgICAgICAgICAgICAgICA'\
    b'gICAgICAgICAgICAgICAgICAg\nICAgICAgICAgICAgIAq4pO0U39azP1vBsY8MG7'\
    b'i/5powdcHWoD+Q4efKfAK0v05wt8qyWrc/rYvs\nbx0ZoL8s1l0C90CFP+eobHUBH'\
    b'Yq/wpz80uT9bD9idA1I4t1VP1EWPaa8i2i/PWxpIi4IbT9ugm6e\nwaJ8v5zav4qi'\
    b'IUu/qerbGFTLYD/v0qYVTCx1vzvlkPKCaJI/5/Z5ZZslgr9Dq4Rz6KQ7P3M2+sR+'\
    b'\nV2I/omAN+IgTU7/TCuNSD5MvPyCzjURbezy/zh3iRDyPML9DRuBQNd0iP3osx1s'\
    b'cQyw/XjGcl+PQ\nUz+EbvrpspIAv+5MBTXVlPK+AqTaslRmJL/5Pj2+se4FP2NLUl'\
    b'Zhewy/oPC1h1D/G79GQZfLWMkS\nvwkIOJFDVwS/QkToMehOIL/utyu8oz4VP9Geq'\
    b'w+CzOw+k6kughcJ274XbrcTIlEBP0vIeJANAww/\nCKHlqImfA7+SE9dwlNjivkpZ'\
    b'G7H8WPK+0S8sLDaT1L4Ao7wo6wpIP8+iwwflOEu/kBBU0OpcUb8E\n1uqmGSU4vzi'\
    b'3tLR6YiI/K258BvyjDT+ANcQaDU0hv5gX+WvrNh8/Txc4yY71ED+B2FBBpR0lP8N'\
    b'S\n7GD2Qio/aZtMt6BuGz9fje5vqm0avxbwaSun1xi/dgPuDdi0F7/7bK2dp64Bv9'\
    b'Zu8hYJffE+Wxbd\nPLHSIr9VEynaDbfcPtqKIvEasvo+j142UJTC0T7WGBFHbEbNv'\
    b'hcVobtQJte+FIswJEr+ub5KaZpy\nvACsPkXnuC7smsm+Tef87tDXs75kDlIzLVyu'\
    b'vmj0iYRgsNs+EDARPqr1475bwbGPDBu4v8Y397is\nZNI/q5in7DAfv7+sN8f1Rg6'\
    b'3Pyubb6zv5dG/OQmGvdkLvj9bL7o9qiaBP/tD9nM1dYC/cwYdlRH6\nWD/7YSKUoU'\
    b'xCvxXKnQwvjyS/xkSu3OoAcj8MRR6g+CxUv0MLhOygM06/lbJdPJGZQz+ocZp4sJ'\
    b'pe\nvxjleWTt4Sw/SOg46ldPUb9Vx5ONs4Y3v44gSnPir1s/Bxvmlp0YNT8/TEuT2'\
    b'p7uPg5eySaGKR6/\nQIUL1kNuHL95L6ixW9rgvv2nRe7gJSe/4kJ+prUfOD/wii8Y'\
    b'+WcZP+Fcv9opjxS/bRDsgENR8L5I\n+g2LRbvHviWZnsYuEwk/JVfio/Nd/b4hWLb'\
    b'5SLwIv+pZbpNYZPM+7HH24SwI4b504idEADsDP5fT\nP6fOqNG+fJS6vDU8pr6sFv'\
    b'O8YMUAP/1baGsvqvE+uNG4qFfeAL9VSejuyDbhPkNVcuQQZ+6+wonH\n8VBL1z6gL'\
    b'Ni7M08UP2TSmNak8U6/siCK53MSXr8IdotdwKY+v+FyAspd9Do/mIxFCW91FD8z4'\
    b'Sc2\nZDMmv17gqow6RhE/5ZYQDNMmID84H+PKvrw1P+LJ1UppXzk/6Ln4tTYFQj+x'\
    b'n4nCBPMNv19E0178\npDG/no/sR76KM7+bRhPAXGEZv95YVk2HJhs/grmu/Cm5Ob/'\
    b'6afhB5MSePtaYz4+U1Pk+MAICUJJq\nET9N7+o+MN+hPuDPiwPHh90+NkFbjAqm/b'\
    b'4L45LRSbJvvmq9lzeRBsS+3BahGgpI2b6khiyLGT/Z\nvv7idCS0MeO+1yTuRIDf9'\
    b'L7mmjB1wdagP6uYp+wwH7+/T/PLSfwasz+B5NAAIQGev53Nocw0Or4/\n+sXQHhmK'\
    b's7+WrbX5vzODv3poPHMBR3o/W+ebBjDQUb+c33DiiPVTPxpeTw6Z53c/O1R0vt+J'\
    b'bj9n\n9pEll6J0P32PGAOZcFq/OiC4NwzeV79pD+9zcHNWP28CyFhMy3q/PEIwse6'\
    b'2Nz9YGuFETQxWv91H\nmiguTls/YOO7P/87Xz/xK1XlV+ogvxViVyjmySi/epmicw'\
    b'6YML/O6BJKO7o8P503eeKNkjC/8HyL\n4pBdMz9JdoZV+/wrPxn+wMYlQSq/5yOfH'\
    b'uoXIz81toQsnLoSv6W81WmXFgU/DGsD/d1s8T6jniDW\n5N/qvqzKzu1IxPU+MubD'\
    b'yrK6tT4eDX/kcpgBv9dmsWZXFbW+5Tqq383W5b6T0A353Qn7vtQJMKvA\nFwq/1Ff'\
    b'2caSSEj8vBmxdFPPkvmLr9NxlrPG+vxb8SaqxvD7RwGzhBzJTv8SHRQ3q5mC/Q3c'\
    b'e2/3I\nZb8Kd0IFcFg3PwCmgNDkOys/yRnhk2o8Hj8jOhVb9EMzv154HLc3SxW/Vs'\
    b'4NtT1cPj8Dm3Q1mgBC\nP+WerSNLzFY/Ke9qRNR+Kj/7dy22l64zv8qhVPPPLkm/v'\
    b'UJGd/X9EL9cutBMXocCPxvWtDiyxja/\noC3cpYxB+b614n6X0H/9vuKD3CRetwM/'\
    b'COYGJLlkA79c47jqhcjSvqsc5dQCvPC+QZbx8Pomyb4k\nxRA5I5fUPovS8FMXBNE'\
    b'+7PPx82Q0w76IqgKLOczhvgNE9XzpzPA+vjpmofRW776Q4efKfAK0v6w3\nx/VGDr'\
    b'c/geTQACEBnr8rcIMT717bP/mBMbgd7qe/hXYf1P/Et7/32AbYLXTGvziHtOwhtr'\
    b'e/Voh9\nnjnHtz+NlWIkSaCQPwEEiaEzOo6/Dvs4xn2wgT9o8R1RyhSLv06eVkVry'\
    b'Zs/wzFew3Y+db+krRpi\nxajDvyuJMvjvZKk/rxam79Mhgj+ard+8ECRuv/eYYvyn'\
    b'4Xk/7sOQokduSr88wjlQpCkmv+zXejmA\n9lY/x/csCJUVM78GHbGAC4Bjv1zA5j2'\
    b'HM2y/9GdJKpcCcj+lBBvTUmtTv1XlLCxPR10/XUEmBPrb\nUr9T866W/mtIP+LC4/'\
    b'XpHFO/T0WykDZQY79mIfAiBP0Hv/8sgVZOxBo/gBoOsLT92D6Z7A0kZo4L\nP3C8G'\
    b'VjnO/e+OgyLVd5udj4weUS+ijsFvwwc5AeVuDE/imShHKwHIL/0T+17hgcMP6LTE'\
    b'qNvDwG/\nhF1SC3P9176ld0+a+HOPv1Kk8dSkOI+/VKYEpKLJkT+YBCZ81ZhWP4w0'\
    b'R86Fmk6/uUf12jaUMT+M\nZsWElgRhv6Qsybj8CPs+TgaGcQTrUD87HVykhPRZv8r'\
    b'LEx0RaWS/dr62nTyEZb/eBcN1bcw7P25d\n/Hgjcko/FrwA1RW4Kj8rSl51I446P9'\
    b'c6D3xCZ9I+YQoeE8n/Tj8T17pEp2n7vrEXNwFUIRi/W9qz\nPsvOA7/jDJllCnHsv'\
    b'tYhtEO1hfE+oBlPftqn/T4612t5xB3sPvby78G20YI+lXmx3Coc+j6olyMm\n+Rj6'\
    b'voxOcPzMj/I+4kQz11v9FT9OcLfKslq3Pyubb6zv5dG/m82hzDQ6vj/4gTG4He6n'\
    b'v5Jw7+y3\nk+Q/sFodX4dj07+FMSm4xAi9vxEm/tF9nMu/mkHznvmEwD/O6nLGrhG'\
    b'hPywuLn6Mn6+/bG+BZPCV\nnD9oj/noutdiPxrlJmSaNKK/upEsgiX4mD8fABhrT+'\
    b'ycP1BQZORBgaq/Vkso6sGeZD+wn37bLuNi\nP/B6H4b3gHK/ZvRha/IMdT+9RpS9u'\
    b'u1Bv/V8ZEvRHD2/MF9DqMKoRD86yDX0/Tlev4jT7V2l+zS/\nn97d69c/cr/ncE0A'\
    b'kLhRP5L7UCNb/Di/Fv4KvSXlQj+l2J93Asw2P4StuHUnV0M/XnpWUZPLTj+D\nXTc'\
    b'aORwUPxxg7Z9YXve+TpP+rmrySD+3WFQd3UTDvuoX5kDSVSS/nnLZppoTED/su+L'\
    b'X8B8oP0Ix\nzDDj2TE/7kWjC0pON79jfIAyO4kWP6BsLCL2ZgY/2XkdhEWa+b6aG/'\
    b'B89Jw2P1eW3Yx333k/p/78\nb3oDNL9qFx8uQuM0P56WcJ4XpSU/kOArkbyOMb+Ma'\
    b'yU1Qg5eP57x51meKRI/wilb62PyRr+JX19g\n5TtOP6VTR7skr1O/HaLgD7O3Rr+7'\
    b'gS+ksM4qP7U+KSB8j0I/Fxy1ckEGFj+QagAVI+sov6hi3uUr\nJ0I/KtMNOkqOPb9'\
    b'fThfNv1cSv9mYAWppVR2/U+fX4U6+JT9Hkb0daFL7PvKkMrvWEQc/RfjpCfj0\nD7'\
    b'/ebk6AEPn9vlX9qF9h3ew+9xS/8+vi675F6HIfQzvgvvV+iC/azAq/msaOoLLWAL'\
    b'+ti+xvHRmg\nvzkJhr3ZC74/+cXQHhmKs7+Fdh/U/8S3v7BaHV+HY9O/IQSRktPl0'\
    b'T815018/Lu6Pxt+uL6CVb8/\nL/P3p3a5xL8JskW3dQyUv6h6xGpkZqY/NBiwvpPP'\
    b'd78YxzwKfm97P5Y+TxNDT5A/zzfV9w8tYL+Y\nkgkxAbyYP+upxoybpHO/deHeDFp'\
    b'loL9U4DfPFZpEP8AZn29SaWM/muo8ZdiPYr+QUjpBpjcePw04\n03vpiUW/t2UCgs'\
    b'VRU79MpmWlLdlRP+0bAQjI8VA/sm7PC6JhVj+8zW9ZWSEevz5nWZiGejO/BYFv\nf'\
    b'5gA2T5pTAt6YY8ov/7STqlxci0/PC756a5kOj/3N8hrKIATP8L2vFIPrwm/WRewT'\
    b'+ExOL9Q7qmK\nDUL5vnijRt2Fywo/LbHUXA32Dr9gyc2nFcQNv3Nu/e4yYCK/+ORR'\
    b'1MxmLT+GJ1XreC4Rv0kr4ybE\nlwW/gKytGSwJ/D7z0OaFqTR/P6Dt6IW65HY/QaY'\
    b'JIlhJUr/CnqxgD3ojv/rFRCaAHBO/c9Iw1nu6\n+r5ynoxAcKQqv4T8gBOqfzy/oc'\
    b'kwfSjkLz+RloSWAOUuv79838MJ3T8/HREaShN+Tz8JMvvkoOwW\nP/aNLv7v7yI/Y'\
    b'TjO8OlULj/ywVDvHsPmvtcRVvwcxRE/mDVsRXFzGb+MzGBR1ZPjPsOMPazOBAw/\n'\
    b'INuFjKsdDT89Vim5WnTAPnEIDa+C47a+iIZpyr0n0r6KtteDlVjwPlknueGevPC+'\
    b'QFt8ngzv6L7a\n5xs6Iw3evhI8+zOSzOO+8Hg0jrTx1j4r1l0C90CFP1svuj2qJoE'\
    b'/lq21+b8zg7/32AbYLXTGv4Ux\nKbjECL2/NOdNfPy7uj/pPTZlZQLlP/1RciIXuJ'\
    b'K/4J1T4zWHxr8cGag9NxDAv663oYzxiLE/vOUJ\nCZrvkb8MGXaUEJy1v2dUXrruv'\
    b'aE/kuT+LCXCiT/NqfXHytqov5qYPiK8V4Q/wwoo7TLyiz9t2MRo\ndCtYP7A4n24Q'\
    b'5Wy/cNzcyHQccz/hNuak1jb0vopE0Zu8KlA/U5vjKKHnD78pYOS1z7Blv0pCfQB3'\
    b'\nh2I/2cn7fRA/ar8sUVBrMJNNP6N3z4ETmz+/Tmqmurb7MT+4tKLBY3JEPxt9Z6N'\
    b'DZSK/fyEX/U2Y\nUz/qwbbt3+YWv1PpyJqyUES/jWuaJwbmRj/g6oVePY06P6N8Ou'\
    b'0mRS+/Ggig/urIID/ccwtBcwM6\nPw8Am36bQkU/5bvUghR0Tb/QQvqIvaQBP06vK'\
    b'1H5QAE/EcgrvoV7Bb9fQpPA9G/JvxxxgWTY8ns/\nae2itvIxpD+kABIJnLGHv8Jz'\
    b'Uaasuze/4luHZ7WqdD91cVoQfPqTv/7V5fMfNGI/fadfbSbwgD/h\nw0uSID9ov65'\
    b'7oTKGsXM/Qer6U7txnj/MX1G9fXM7v+2HGfyp5He/S2RmppfTe78nDf76Ma5Uv3b'\
    b'1\nFb79fW0/565WNP8ffL9wsM9eBXI5P96XTEZx0DU/kRYrSbfJXD/FGV2BpMH+Po'\
    b'jf1KHE0R0/GnLu\nfX/eOb+2lt53cIQZv5IVuyJ95yK/z+WOMpIOMb+oJ/nDKhzAv'\
    b'sO6NU9lzS+/1s8/6JSdPb/oqGx1\nAR2Kv/tD9nM1dYC/e2g8cwFHej83h7TsIba3'\
    b'vxEm/tF9nMu/G364voJVvz//UXIiF7iSv/0jHkE2\nvuA/2K27avKax7/A5Sf4I6y'\
    b'yP379HEvCFsq/4VUAeD8hqD/Qzv3dTzuxP8ECFpbiqnK/H8NngcQE\nlL+5r4dsS0'\
    b'afv/yRNr5BFoS/IWsTViJrlT8QODV4461pv9y6BBOhKHC/+iAxB+bjaz8VUsrxse'\
    b'c9\nv8xZPIY5I1I/Ub5gECT8Lj/QJ11yb35xP+rX0CAq8j0/pfE0u//PTL+YX8xub'\
    b'uU5P1wFBJwHjVi/\nyULTQ08LUT+ZqUAincFYv2CBNgxMwVM/gaJ2eg5ZVT+/0xpl'\
    b'JWIbP77ntyPNJgm/cd4Bv4cI4L7z\nr9f8wJUvv7xg04+Gufw+TJmWgzTPEb+gJ3Z'\
    b'fIK0jP+bWeRKUVke/MRKg/S/GPz8Da6HI3YASv0/9\nr8T7cLe+jceMPtVFBT9tb6'\
    b'eYuS+LPzYd2A+K8bO/+I+SSbF0eT9mMw/HBf6XP8hgiYJxJmO/7b67\nki/Dhb+p/'\
    b'jR3xv6cv7z0N7Bvc3U/+OrPFFONhz+27X9uBm95P6eQE+XKJVC/u/A4JK4Yfr/Vm'\
    b'BMi\nuKxYv44qybOubDG/EUZYMvyzQb8QUz/ThU9QPwtXJ5YhcVi/QT9vVy7iYD9s'\
    b'sDM66Vk2P9xFyTee\nqxy/bHudtOzH4D5GI3xRTaX3vozNr4aF/bw+fiaC0QNoIz/'\
    b'iA1VnoK/2PsLOLIlwrRW//oWqLEGX\nFL/tOwzTHNLkPkCAoclILRA/wJulxVrdKz'\
    b'/DnPzS5P1sP3kGHZUR+lg/WeebBjDQUb9ViH2eOce3\nP5pB8575hMA/L/P3p3a5x'\
    b'L/gnVPjNYfGv9itu2rymse/QLf6o0WY1T80CKnb3BWMvwFdtXX/Jak/\nXkEqRjUa'\
    b'tr+G5jl2BvuJvysRXtqddWa/WSNvG1oKhT8sLLRD9ZChP1Cwib/ww2c/ijRfIocq'\
    b'f7+N\naDCcmd99P4eOfnJaT1o/EHMaqcWPaT/nnyd5P0ZMPzvp3OG09VU/xKx3ItD'\
    b'YWT+0rpGu3qJEv8pN\nPvMNLVK/Txys8y9LST8Rl19n64xCv3vyXZee60Q/rkwjUn'\
    b'dwN78ehPxooRZBP6SxAL6Lykm/IAfW\nUkJOUL//TCkq4e4Pv6kt3J91hig/yspo5'\
    b'AO6Mb/dXm5MyT0Qv+9Ls1WXkxM/n79OGC0j8b5IIRpJ\neFMyv4nZT71WYyc/z5Ds'\
    b'2j1z7z54RDeTZPgGP8HMM1nlgN4+A5Uk5NEz8L4w+To2bDiqP/XUa4od\nUXU/1Y+'\
    b'n0w9qtb/g5FyfOO6Qv6b9HeUsLDM/xvY9v0dpeD9fn677drxZv/ye6XhbwWA/Uuz'\
    b'IQCcr\nST9IFF+bj6ujPzLLRMfAvHm/UVjOhnVpkb8tZ7C8Dv9sv0A80nAzZj4/dU'\
    b'CJGsQ88r7RtouS8o5c\nv0S5IU9YhUA/9BaUtpyrFb/ceAkXMo5cv8kVqu+V8Q0/s'\
    b'wCBEIJdR7/vL6nAo1gTv5kIgItyssa+\nE6a8JH8wCz9EB+vNkzoGP3T5wN7QMfI+'\
    b'VXILQKAaOT/RTr8L20IUP+EGGkyC/R8/hxkM5AsFMr9g\ndA1I4t1VP/BhIpShTEK'\
    b'/nN9w4oj1Uz+LlWIkSaCQP8/qcsauEaE/CbJFt3UMlL8bGag9NxDAv8Dl\nJ/gjrL'\
    b'I/Mwip29wVjL+cfLFFRtXjP+VwGQM+4bi/E00iVjZUwL+zo4Juv9bYv54oo4UrwJ'\
    b'C/yYst\nPNOsvD8A1Mbpmienv2m5f4Iu5aa/YgSBQZsAqD9XwpnrY+mqv1IHZgr8m'\
    b'5M/F1hGPMH3Zj+nQSNN\nuVNhv4nqdX8vTVM/l34e+RhBZb8/Q8gLQDBWPxg0L9y5'\
    b'2HK/Nujlgphdjz+hMLEADIFNP1w0z3VX\nAD4/z8NmBhVeWb+P25v/rEU0P1tdRA/'\
    b'oXzK/xF/lb6OkYr9KqRxpvQcvv3gd5HY/vFQ/QvPhaIC0\nYL9a2J1/5vQpPz7uKq'\
    b's3fkI/JgyiGYalIb8K1mQ0Cf9Wvz1rxzM0Tz4//xnfOeb3Lb9CxV5l4LUA\nvymND'\
    b'HiJ8xK/sHZ1l2h5BL8Fw1FZPV6Xv9qMBtp7AKA/HGqTrgpvgr9mRqHeqnV4v406s'\
    b'ZF1VFU/\n6eSa+j2gaT+9Y4HUWiJwP5Il7bmYOEk/UZLYYqxaUb/z1hEEGftKPxdi'\
    b'G5aU3VO/o4ouzDnodr/J\n3X9RDQ8yP+r1UEcXTVw/uVzUgJwyYz/I+8D8U3pJv2m'\
    b'e3lDq0zQ/b4OSLjqzZT/nUy/fhHA/vzTe\nrGytTRw/JGIK55qTWb9+GccH+Rj8vo'\
    b'oe9cGKbCK/UjmwtNycID+53X8/qF4BP4HzmSPMzCA/wWQg\nnJZLMz9hMjpBw2sQP'\
    b'+f7ms2MviU/Jcd0YDyyJT9QFj2mvItov/nJnQwvjyS/HV5PDpnndz8DBImh\nMzqO'\
    b'vysuLn6Mn6+/qHrEamRmpj+ut6GM8YixP379HEvCFsq/AV21df8lqT/kcBkDPuG4'\
    b'v/usGCpk\nFtg/YE8WfPUwv78FuCj5DHuCv5gUHyY9Qbi/Zdy433Jdlz8T9mB9NKZ'\
    b'2P7yxeOTwCIy/NKz8HEyF\njT/claQM0NeiP3PvFJuyE4A/mOGayXm/kL/CiTA1wc'\
    b'ZwP4pm3e7Pbmq/EXoI9iM6az+7IdsChzJq\nv1/3fKXspFu/pAw/Zr/Da78uscBhV'\
    b'3gaPwVSJ5y1AFA/OorYL7iaGj8l9rFhNONDPy4D4MQG/Ek/\nDh0EUxWUPz9aURtC'\
    b'8FtJv2h1juJ6vFS/8jqK5mMNZz9J96nY+rBQPygOIKtGajq/HOAyorSxKT/5\nUmE'\
    b'yx+5aP84KluJYpFu/bBPb8yStKb80FI05XjsOP2AD7E5wCy8/QRnj+vB/ML+J1aF'\
    b'JlQqEP+gG\n48a3ZXS/Z0e4edXIY7+m0gNED8MVv5f6RHSg/Uo/s/llEFvqFb+rZ8'\
    b'CT2vU5v3y9qkPGxCc/PAZj\nWzDfJD/Da4bIms02P6233tGm2uc+KA/0i8OcZz9oN'\
    b'e2ojcAsvwKQEJXONzW/ksYLuHt9RL+Vbh4f\nnU4mP/nZB2pvSEQ/neDrB48BS7+e'\
    b'nMqd5tYQv9YvZKMQbRi/Rzm3RXhlMj+uJlbQkE//Plerud4K\niQg/x9kxMY4PFr/'\
    b'EIAFimt55Pl+hNh6FZfu+Mx5XUowiCb/fTkLAHVXvPg84AD+U3QO/UtW2LC8z\nEL'\
    b'89bGkiLghtP8NErtzqAHI/O1R0vt+Jbj8O+zjGfbCBP2xvgWTwlZw/NBiwvpPPd7'\
    b'+95QkJmu+R\nv+BVAHg/Iag/XUEqRjUatr8TTSJWNlTAv2BPFnz1ML+/NmPu3P1Ty'\
    b'T/y/3T+wjq+P+tJSvr2g5c/\n/hmfkpzxu7/YwB2uermRP5F92bjxUKE/S2WZ00gG'\
    b'cr/jYuoTeE91v16NVpuynoS/maz6KHQ6fT/i\nECmCV4JHv5FP7nU0c1I/5OQ9sYs'\
    b'fJj/aMozSr1RlP0QUtwRzFBi/XSsTmUccgL9lmDwdM1tFv7NO\n8w4jQAq/3vzvZr'\
    b'xMUz/Cepwsj1hnvx19rEne/GA/IOGcWTKGUD/fAnG0z8FFP5B2DnUOhFg/IRgt\nZ'\
    b'eF2IL9/lAKr1TlEvxtmGQFXhyo/CS08Hdt7Fj/4AjDIPlhTv9DjafQAkCm/PzSe7'\
    b'AP6Rb/M8sLH\nD0FCvxRb8VYAyD8/rnLO9VaeI7+1gx28q8pmPzHRcC+HlIK/7KAX'\
    b'D3EibT9fhRdwt1pdP89vjWu7\nQjE/1k2kThcJGz9LukboTj8sP3oy+lu3FQ2/v8f'\
    b'3UpK9Kz825yfxZJhMv5Wl73vBY2A/fL619OmN\nQb9EHpmoPyocP+Jg4XrrK1C/ii'\
    b'0ONCg7J7+2Ao3voDFgP4CciTf/BVi/wbGNMIVPSr9w6Ic0vp4R\nP3bxayJeaRS/6'\
    b'+KV9cj0UD+RJ+unPp4HvwdoUCf2wwc/UQHJP7sE9T6BbAjr/uwNP/Gh8b7SJBO/\n'\
    b'DnADt6/8JL/eDeKcLJkPv2MZx3O5rQ8/T1bcUikmFz9ugm6ewaJ8vwhFHqD4LFS/'\
    b'ZvaRJZeidD9o\n8R1RyhSLv2WP+ei612I/GMc8Cn5vez8LGXaUEJy1v9DO/d1PO7E'\
    b'/huY5dgb7ib+zo4Juv9bYvwW4\nKPkMe4K/8/90/sI6vj+bH01cNxPnP7OTnSWCRc'\
    b'S/4gm1m6Vgub8dmlputsukvwmYhus0dJu/bS0x\nh5d2hT9O28r65+bFv+IEnC/gl'\
    b'bY/TTOx3y/ZiT9LvWlVio2Kv1I2TQFdRpo/HfefUarGgr99MW6o\nDxRkv2KVdGlF'\
    b'q4k/pLhVNXEOqb84sfwPOcZkPwcQRN5BxlW/3ug7Yzgjcz8zwuonqLBov9P702e/'\
    b'\nu18/ecOwRqltfz/U0gfZ+LxfP7k39Bg2+jq/bWzDkdk1Zj/zzxHJ5vthv+dRvnD'\
    b'Ffjq/FZ3GI5iB\nJD/HAMNXu7Q9v1fYJwZySWG/zSSWb2F/Xz8CQMORtCQiv05m3B'\
    b'33aDc/GXrCEJgpJj8u5t9mz/RV\nPwfTNiXcBFs/XequhXd3U7/UvUdySdweP/2AQ'\
    b'wkp2Ci/cgIIxFnzPj8SbNferZ1QvyBMKZUb4yy/\n8jo4VrzD7b6G1crz0XBhv1yG'\
    b'NzVUi2I/F7naW2T+TD/d8MuNDXA8P3EevEewF1C/Fa4sefVxR7+t\n4UYUuthRP4Q'\
    b'+KRXYuUy/ALhp9SzMUb98Ryi1o3IyP0dR/c/VPCq/U2gNSJF0XD8nud7pgK8IP5I'\
    b'w\n2g4BJSM/w9xzvorQEL+Egs5cTPXqPsvNz7qMfRy/JRuNlX0HNL/+XzJgW5ISv3'\
    b'FkmQR1Tx2/vIFJ\n3wjy976Z2r+KoiFLvzoLhOygM06/d48YA5lwWr9OnlZFa8mbP'\
    b'xrlJmSaNKK/lT5PE0NPkD9nVF66\n7r2hP8MCFpbiqnK/LBFe2p11Zr+fKKOFK8CQ'\
    b'v5kUHyY9Qbi/60lK+vaDlz+yk50lgkXEv3V1mW/j\nLNs/8XtMZnr9wb+0FC90VJ1'\
    b'nvyPtOMqRk7i/iDv6dG9doD9Xw5dvZ6i+P2BNxbr/7MW/kGSstxyr\nqz8wflB4wN'\
    b'yCv+xe2LcHGUo/0dDfHfLre7/HzcilLvd8P0p7fo1tUZW/oY451/n5oj/bNafjbK'\
    b'xZ\nvx801NhkrXM/6pRBhEtmbL+McZthEoVEP4/D5GjnJGW/PIAOmI3xgr+0iyxrO'\
    b'eRSv9bhoePttUY/\nIfZ9YbPkWr/sxfpfOeI5P08Y4fX5Njs/JK/boF6XKb9diKAw'\
    b'fhNUv8OgQpJkgC8/prttzYiKRD/P\nJ0noh204PzBL/otNuT+/UCODjPu3GD/41sX'\
    b'GqNtcP3EPaRd98XS/zh09+LEtbz8fWCU3d2pYP8mp\nKiMOBxw/zdDiF9TBQr85/1'\
    b'37PWNQv6R8DI8XkDQ/KPypWKGHJz9tIahQhvFIv9uoFsedhUc/w/ef\nfVpRY78vh'\
    b'3Z5Bm2XvlxwuJm0diE/dpW7RjPPUD/9BKTYoRUwP+2s5MEkXlC/GljJK8MdYD+2Y'\
    b'zzk\nOIUcv7DM7hkTPyM/gbMVxjkyTr96+zfvWpwcv5gNniJ5eyC/+pEy5uueJz+z'\
    b'kf2gmpQXP4Hw8l2P\nyPs+fCg8SvReJT+8hOZRvVv3vmvsyyfTGik/oX6tSuBRID+'\
    b'p6tsYVMtgP5myXTyRmUM/OyC4Nwze\nV7/BMV7Ddj51v7iRLIIl+Jg/0TfV9w8tYL'\
    b'+T5P4sJcKJPx/DZ4HEBJS/WiNvG1oKhT/Iiy0806y8\nP2XcuN9yXZc//Rmfkpzxu'\
    b'7/iCbWbpWC5v/F7TGZ6/cG/X1vBt1xB0D9NnGPC/+VaP2MuxxTs/6Y/\n3r1kw6h0'\
    b'r7+tDvQntCCCvyb1Py2qYas/r8LibxZZtb/nYnt8KFJ9P/AH3k4fC5K//sGvvEmu'\
    b'cz/e\nE7Aa4V+Wv1smN8PmIps/+SybehFBhb96f9KpzLxmvw/EYQXg5WI/gR3Wj3E'\
    b'Bb78MHJjFty9MP25c\nXhWEcRq/ZyCvfSBLZz8qvzmFXqhEPzso7VgZaUq/mtFclV'\
    b'TAQr8SwxTb1xQ1PzIVW3RDnie/dPdo\n1J8GFD95bZEa/AhTP3bEFq+lnVI/yt68R'\
    b'ASVSr+18wlu2xACvxs/duOKTB2/yyEVFZ+GGj+w2hvT\nmBdUvz/eEUZXZGE/L41l'\
    b'tnNqSz9NwZ2tpsQ5vy7R8KQ+UA0/VrZBcw/0ND8KL4d/W7dHP33XPDS6\nFDc/+Or'\
    b'S/N/YJr8GvgpD7VpQPx76O5HbKUK/k3TIweuvVD91E+IaVg4qvx2LzACu5xM/szH'\
    b'1a0Du\nNb99bmlAg1FBv4S5NGdUz0c/L7l96JRNTL8LaSW9wzcDP3r4FIt/PhA/pM'\
    b'o7K3pb+j4fmRy1CuAG\nP5qE39oyGvk+l8/PyJ31F79GmtlTTSkRvzwoK2dhtNU+e'\
    b'ETz+ttW5L5QmUZQ8Yf6PqHtSEfjjRm/\nHBJOKtLbGb/w0qYVTCx1v5Nxmniwml6/'\
    b'aQ/vc3BzVj+krRpixajDvx0AGGtP7Jw/lpIJMQG8mD/N\nqfXHytqov7yvh2xLRp+'\
    b'/LCy0Q/WQoT//08bpmienvxf2YH00pnY/2cAdrnq5kT8gmlputsukv6sU\nL3RUnW'\
    b'e/Q5xjwv/lWj/WZhCqFE7SP0pBcQNFqIG/9ZCOENM8sr9t2oGDW8l/P796pQfa0I'\
    b'A/739Z\ndHnThr+wh+mhf61DPz8ZbMEyaim/s4tTMID2Rr9ybtK703ZRPzauUJnH4'\
    b'HA/Q84jzyKCVD9hyAQ9\nuatDv1eXVPQGojq/rRF5VH9OP79UnizDVJQfv73FwoK2'\
    b'0Sg/6BHUIMLDRz9sRMeEyjkGvwNZm4jL\nIBm/sTXPUQ7cRL+3cnwxUWMyPxA/NrJ'\
    b'Mlhs/IBCR+MJzBL9uUc+M4HnvvlhzAHHjADM/W2BdR0Gu\nJr/nXuY+GikLv1+soI'\
    b'6uCAW/KCTDArKZDj9XwtHYEVI5P0tZ8n18Cma/pAxfRLudbz/fWQf4+XcK\nv/oiP'\
    b'bf1LQo/q61R1G2JLr8wnbR4zyBMv4Z/SVJS2UI/TLQeY3/mQj975Qwkc7NcvwJXb'\
    b'IF2VGC/\nBeWzIzNpZ7+yuWAvn+hCP1+fDoC8g10/lsg2ARS5VD+9rTfJ2vAgvyS6'\
    b'/+ac7TA/HIMUrX02Vj9B\nXZxv69YEP0vJDtlN9xG/lcDi0CXKR78JhGUNx7D3Pky'\
    b'vTUsPRuu+d2ittkUYFj8UV3bEjGPhPoeL\npW8mBgg/ut8wItlwHT/yX/haq7LUPs'\
    b'u/tzQDwNi+7M7bFBGAJD885ZDygmiSP+DmeWTt4Sw/bALI\nWEzLer8riTL472SpP'\
    b'1NQZORBgaq/6qnGjJukc7+YmD4ivFeEP/iRNr5BFoS/O7CJv/DDZz9suX+C\nLuWm'\
    b'v7qxeOTwCIy/kH3ZuPFQoT8LmIbrNHSbvyHtOMqRk7i/Yi7HFOz/pj9KQXEDRaiB'\
    b'v3LIy4gq\nmMg/BGYXYITktL/6AKh8cJ5nP0oFFxUydpi/ghbP32h9ij+PrkcByBZ'\
    b'NP4tjFHU69hI/AoXL4zHY\nUT8+agjHegBlvwpgzTI27ka/vD62KAx7er9mT4VSo5'\
    b'oyPxIJHnP1njc/DBnYelElOT8eTrUjyDo+\nP0vbRxXEHU0/01Sb3EBKDT99HoBub'\
    b'GcdP6f5A9zO2Cq/wcy8mkkUSD+tOHuVZZYYP/+2QBMU9Sy/\n3l2LOPK2Iz8Cq8f1'\
    b'pfAYPymSRJeQjUY/RXX/7f/2Pr+WVVLsAH8dPxAx2odsqgY/c0eFYZ3z/75Z\nEeu'\
    b'kgXY5v/eRlTOhJWI/z35iLBjjaD+K0oZ52ywzv+PDqhuvAN0+4Q7euaGgBD/ikOD'\
    b'Cth1WP1Vu\nYDIT+sK+QuNtteqyO7+NjotSuDMpP8Rjae2A3Vm/6eiPA+TaVz/USn'\
    b'zoUPeYvqkAnRxcejY/wfuf\nBruyQL+D69Dp/y1Dv2CiAiwOQUw/8u/vQu/EP785a'\
    b'HIUhTT2vhRL+Ug8uPk+SwvoiN1tJb84dUcE\nEx39PikP3dTFG9G+3EBnxAB1DL/i'\
    b'MKzE6nAGvzPt5A5VbgA/xIKZORRDCj+iKV+xgRoDP/fqpgjd\nXwK/YvVns0ww+L7'\
    b'n9nllmyWCvznoOOpXT1G/V0Iwse62Nz+vFqbv0yGCP0JLKOrBnmQ/d+HeDFpl\noL'\
    b'/BCijtMvKLPyJrE1Yia5U/lTRfIocqf79iBIFBmwCoPzSs/BxMhY0/SWWZ00gGcr'\
    b'9vLTGHl3aF\nP4Y7+nRvXaA/3r1kw6h0r7/1kI4Q0zyyvwNmF2CE5LS/n1Tbv7ELu'\
    b'j+Zc+j88uZgvxU6kAhsp4I/\nxOiCohDHaL8utxCDB/Apv/MyIILVdkc/sPpREsWx'\
    b'RT9E+3+tdV5AP+qHIXgRzUo/maJG7bZvbT/n\nHuVUwCxGv1uEkjq8fjQ/XSvQ1Oy'\
    b'LRr+hX7zE0HopPxs26GcoXVm/IhXjH2sBRb9vlMriwNTsvmT7\n88BTeAI/88z9hF'\
    b'DaQr+XRArk0UYJvz3wT53/BB0/+6M8k5SeD78N1oKRNZj6vuBF75JOQyq/RbFz\nD'\
    b'+IXNT+nAkSA98H9vhncu7KJ8iC/NkWSwvsbAT9fbIooqsZkPyHUgBRHbXM/DHc2T'\
    b'SZhcD9PdyXM\ngd8nP8+MxkJnLTS/IUX4iW92Gr+VdSftkKNNP8g6Rj8P/z+/gWU/'\
    b'7xMeSr8Czxf8/11Wv5UhP/5A\nqEI/TeujQ74JC78nGVVgoCfzPkWvZteGci6/i5Y'\
    b'Jz4qKHj+q+5SqiGovPzqKQ2DErUS/ivd/Gmww\nQD/0DFyJsGUcP4t5kGpjVvO+mK'\
    b'PtKe9rET/DDEWaMc3yvnjCaseXvP6+w5NOJuZ8Cz9mAeFcqnLM\nvsgB0dW328E+P'\
    b'97suq6F/75nM7o68amavk6I+sjJbfw+wqxAHK8xxL48q4Rz6KQ7P0rHk42zhje/\n'\
    b'WRrhRE0MVr+Vrd+8ECRuv62fftsu42I/W+A3zxWaRD+G2MRodCtYPxU4NXjjrWm/'\
    b'j2gwnJnffT9U\nwpnrY+mqv9uVpAzQ16I/5GLqE3hPdb9O28r65+bFv1fDl29nqL4'\
    b'/rg70J7Qggr9Y2oGDW8l/P/UA\nqHxwnmc/nnPo/PLmYL/4QFLqWs/kP2Knp5aVws'\
    b'G/UeTr35KwxL9hvVmEmQvXv3k2OLvTgqC/YLkd\n0cgIwz/macwFHfq0vxBvSndOk'\
    b'WO/dGsSxxTMoj8NOp/830yBv8bEcZYSC5E/FbUc8kfJoL/b7njm\nhIKIPw7uwgm/'\
    b'0Xk/UTbu9F5ZZj9F49PnOjuAPwT+4ZEi84C/idkToXQSkj/kZBYYwBBtv8lGVinM'\
    b'\ns2u/1/Lc/01GOz9KAc3XRZZiPyYxVLMB8WE/750USpT1U78Phqx1hQ41vxS8NmU'\
    b'd5Bg/N8lFgaNQ\nI78sF3FtbvYqv1yG3c+pOls/X6O2MEuSSb+uzXoo2VZAv7pWGd'\
    b'k6XyM/sga4FvIRLD+VVG/WIcEv\nP3XCbl10diu/4BETU2E48D5U1mr93LUvPzOEZ'\
    b'rOAJCq/97mseQDjXL9yOLKZ6ucLv0aodaWEQjY/\nJvE9Y+ZzSj/OsY4r9JHzvqIH'\
    b'1L6nXy+/4v89gZ45Uj/uw169S7Ilv7JcJwsg9gK/UG27FhHSR798\nP0BCRyDyvkC'\
    b'mAzm5HhO/ByHhfO5WED+U/yVXLqzwPggWanRw7BA/uILE+nkCHj86GNH47wD0PoD'\
    b'E\nB7QDlxo/iC5FgUGqED90NvrEfldiP5EgSnPir1s/6EeaKC5OWz/3mGL8p+F5P/'\
    b'J6H4b3gHK/vxmf\nb1JpYz+zOJ9uEOVsv+K6BBOhKHC/gY5+clpPWj9SB2YK/JuTP'\
    b'3XvFJuyE4A/X41Wm7KehL/iBJwv\n4JW2P2FNxbr/7MW/JfU/Laphqz+8eqUH2tCA'\
    b'P0kFFxUydpi/EzqQCGyngj9ip6eWlcLBv9W4KJFK\nPNE/fAZV9XRxwL8O48ht7Hm'\
    b'Av0C3urH4F52/SkcuuiUaoT86tzBWVS6PP/W096vyj6u/nT1wn15B\niz9jU2PoDH'\
    b'SJP1uMK1ZpHnK/p0zTNpeqlj/Y1Lq70fkHP8tGFolHNIY/We9YdDr8iL/vNSab5w'\
    b'FX\nPynmK3/id0M/lh7YwSFrij8i125aTbdov90gwAfC1Ee/lU5e5rdnND9WG5Qoc'\
    b'vJJv2s3tcTv91a/\nIKwlhHxdRT8RG7wrfi9QP/zKERU3EVw//zmlFv18ID8DJnkC'\
    b'eBZDv5S8TRwIzhM/eOgLzNlwP78y\nc5xzTO8gv1Zvuar1CP6+CN4aupI/BL81zbg'\
    b'Xdas2P0zZLDfrpyQ/pmvjuFBSGb8WrJq1xtkQP45n\nAea3Q0S/s4THYauJaj+b30'\
    b'Xby7YgP6jnZSRg5xG/BAi0ZhT+VL+1bongOAsCvx50NJSTjVI/DR7c\n1GDoYr/NS'\
    b'9NA2WMRPxhfxTnjFSa/7SeoaZaBUj+kMeZy4HEUP4KhrOOUFCI/Q5+hOqSCK789A'\
    b'alX\nrKL9vu55SsNULwm/W7edy9kXKL/fkn1i84DkvqIuNOxjDCC/SdCPI5KNG7+i'\
    b'YA34iBNTvxsb5pad\nGDU/XeO7P/87Xz/pw5CiR25Kv2X0YWvyDHU/m+o8ZdiPYr9'\
    b'u3NzIdBxzP/cgMQfm42s/D3MaqcWP\naT8dWEY8wfdmP5jhmsl5v5C/maz6KHQ6fT'\
    b'9MM7HfL9mJP5JkrLccq6s/rsLibxZZtb/wf1l0edOG\nv4AWz99ofYo/xOiCohDHa'\
    b'L9S5OvfkrDEv3sGVfV0ccC/HLg0JTdH1T/qUukeM33EP9jtW2vjCas/\nqJCuaVYf'\
    b'vb9wbEICDY91v2mhpFop5HE/uWM+opGevb+dlSIoRXFgP0fEFIz942M/RwNYyRM1'\
    b'Kz8Y\nhKM6yZklv0gz0dMCAnY/gG4W/MCqkr8uOTh2dMx5v+zvvQsNJVE/0iZTAW7'\
    b'ygb/zT3HeWBBoP0Dj\ncKSSqVo/2RuUJ3MFRT+1NcY5TOczP1Oxxjd9LT4/TG+oNj'\
    b'pZSb/BS+2nh8xIP6PZTGBZoFQ/4mCp\nowMgQr80Ue9zeSBBP/S1T/WAdDG/5FNV9'\
    b'sBlVj/QuLoU5CY0P6aYYAEJwOa+uLnEkm5pHr9iin0c\nIvwhv3C3CDjwViG/fVPk'\
    b'L1Lt+T6sUJqyd6Y/v878904LglE/YUxZfy6vZ799X5fz8NsPPwGHJ9e1\nOiS/XR7'\
    b'U2/vCUT8PnI9kack2P4+IbmVcJ1K/yfHJGhUUYD+IUKj+KtMwv99tMmeyZQA/af7'\
    b'J9LG6\nT79KByOcrQkkv55tbd0SLCC/iKoVFeyWKD+Yn+AsJ38gP0ZgxYYbU/s+/5'\
    b'leQNNWJz9YHKQjXaD0\nvk04o1LkBDI/KhMyYV99Iz/UCuNSD5MvP21MS5Panu4+8'\
    b'itV5VfqIL81wjlQpCkmv75GlL267UG/\nmFI6QaY3Hj8VN+ak1jb0vhRSyvGx5z2/'\
    b'558neT9GTD+nQSNNuVNhv8KJMDXBxnA/4RApgleCR79K\nvWlVio2Kvy9+UHjA3IK'\
    b'/52J7fChSfT+xh+mhf61DP5KuRwHIFk0/NrcQgwfwKb9hvVmEmQvXvwrj\nyG3seY'\
    b'C/6lLpHjN9xD8qPx/q6YzXP6G5mXfOFIw/wHIC0AaCxL9QzuOd0K2DP6MTSVIC1F'\
    b'S/9UNz\nS936Zr+InU4iVTtZvwy/6SRDX0A/Lh0FeIGXNj+tufNfPPpUv0iLjEzEz'\
    b'Fa/mOZbuHQGa79x1+Jx\nprc+v1v8Dk4Yy1I/qB6nxsIFTr/gR/i+GoU4P4XGl4ML'\
    b'zjU/qsAlEm6I/b4zSjK2Uxsov/Z5LlEU\naRK/yzYXN9BpDr8RrWuyDkUeP9ua7Sr'\
    b'5BQ8/yVxrNoiQ5L6TudW9fU0Yv9E7kL+39DQ/iR4bmgVi\nI7+txy8H2YkXv0//KK'\
    b'/Zn/2+tyDHxgrdBz/a7fxVcrUGPz0dnzBP1rY+LfheAg69jD4lRs9ro0nM\nvuXmi'\
    b'gRe1w4/VnLOccy4/j6lUR3YqGQIP6bg4CAotBi/udvXB+hkF79NtlMIicLsPqiLV'\
    b'va7g/G+\na7VqLICJEb8AYZlukUm6PihIzTWY67++men7JQQTBD9fCXOz6HHpvgoO'\
    b'Rn3/6ss+YrOGl55Dzz5T\n1Lh8FnPtPkzgvbgAxuG+frCtQBXexj71rYJH3Bzpvr/'\
    b'ROjL1YuQ+H1aBPKUL5z4ds41EW3s8vx9e\nySaGKR6/EGJXKObJKL/t13o5gPZWP/'\
    b'l8ZEvRHD2/DDjTe+mJRb+LRNGbvCpQP8xZPIY5I1I/O+nc\n4bT1VT+H6nV/L01TP'\
    b'4pm3e7Pbmq/kU/udTRzUj9QNk0BXUaaP+Ve2LcHGUo/8QfeTh8Lkr+DGWzB\nMmop'\
    b'v1pjFHU69hI/8TIggtV2Rz96Nji704Kgvz+3urH4F52/2e1ba+MJqz+huZl3zhSM'\
    b'P0IdvAt3\n7qA/evpZayfoob/Ne5fiDQKKv9hFL+nmqWA/7P3w/yxXab997DwjudZ'\
    b'fP6mx6vbsnV+/vdf1d2Qk\nZr/1p7Auv0P8vgyXxh53Bku/s+UG/6yCYL+QTvey2b'\
    b'lYP0aOlMTHCk+/1J6hq69PdD/1R2kXfZFV\nvzLGSXwOBVG/UPU1FickLD/GP3KPZ'\
    b'809Py8AuGd+UiE/KyVsUeBZKL/vS1tQ43QOPwFbUA6jBxY/\na9RxnhYBAL+PPQoq'\
    b'1E0SP0JMynzJqDe/nE5TCx0jOT8uOn3u5j8YP2XwfYS6Zfu+FJd/qPmREb/g\nSBt'\
    b'qDWsWvyNBKAJExt8+Y2eebMiWBj+TC2dQw1oUv5rj3OR76RC/s5j8FNci977S+UP'\
    b'xluvHvok0\nYRcZ6hI/cAVErFYBx74o/Lp82cbtvohsXew3Mvo+oeu1O5iC+75lAC'\
    b'op5XsBP9Bpu5BStOW+A3tx\nutE+8D5xBa2yryDqPt2ldzbuprk+wBDithGi5L7Z9'\
    b'YgQQOTrvktWS2nkJuM+J+Jh+n1c3r500u6J\nHPHDPmEeD0AlE+a+GUjPX2qjzb7P'\
    b'HeJEPI8wv0SFC9ZDbhy/e5micw6YML/F9ywIlRUzvy9fQ6jC\nqEQ/t2UCgsVRU79'\
    b'6m+MooecPv0u+YBAk/C4/wqx3ItDYWT+Wfh75GEFlvxB6CPYjOms/4uQ9sYsf\nJj'\
    b'8c959RqsaCv87Q3x3y63u//cGvvEmucz+wi1MwgPZGvwCFy+Mx2FE/q/pREsWxRT'\
    b'9fuR3RyAjD\nP0tHLrolGqE/qJCuaVYfvb+/cgLQBoLEv3r6WWsn6KG/Qs1ztnYew'\
    b'D/sFO5tnGugP8RreEY4Lnk/\nGmrlxN2XkL9fs/GgFY1mv5XNiPOi32K/wRuS8ktT'\
    b'cT82ZbYs0pdnvzkq9qq6OWO/S49mumlRaL+5\nc00LEw1kvyc5Y7ecHU8/hPinX5m'\
    b'Rab9LKwRAzzBPP4lfHUl65UQ/Z+VYA6CnHj9rRTtVcHg0v6Zf\nY3JWlTa//ppmfT'\
    b'AE6j7Nuqx9KaUMP1tN6bNC6yY/oE8hzkRWDT9CeL7Objz/vs3sSb4PTDQ/s3qy\nI'\
    b'wUtA78o0aSoj4oWv7CX7hXjgv8+D0mjAQza6D5auvKGSKvzPitAL1/hRPu+SK4Yi'\
    b'+6n474ga5jf\nrcT5Pljw3o9EDfg+ibSWadLAM7/N3Wjw6m4Dvz8cnF9F8QQ/uBwL'\
    b'EbFOIj/HrmVWRffXPmVU7dBU\nahe/Gfqo9FitKD9ESpaSec9zvhlyciK/3/Q+nDg'\
    b'/4TKHGb/uqa+bdOeyvsggXH+C5e2+rduIRJtE\n7z6460j9eHe4vrezJFr+uNE+QF'\
    b'KueoF67T6CDcf01lTDPnpijcNAPNg+gDZfQFLPzj5eRuBQNd0i\nP+ItqLFb2uC+w'\
    b'+gSSju6PD8KHbGAC4Bjvz7INfT9OV6/S6ZlpS3ZUT8vYOS1z7Blv8snXXJvfnE/\n'\
    b'mq6Rrt6iRL9iQ8gLQDBWP78h2wKHMmq/3zKM0q9UZT96MW6oDxRkv8vNyKUu93w/'\
    b'3hOwGuFflr+C\nbtK703ZRPz9qCMd6AGW/R/t/rXVeQD/jacwFHfq0vzm3MFZVLo8'\
    b'/cWxCAg2Pdb9RzuOd0K2DP8x7\nl+INAoq/7BTubZxroD//QgUHZennP/3ag9BIeM'\
    b'm/jvqXUL/+sL8UeGBcQSzcv5AMwI+FG8o/mOIr\nyS8bjD8Q9FXboG23v6pKLETIv'\
    b'6Q/rAzPs/DVoT8yMQxlWV6/vxHrL0ozhaW/AORd8DIddT9d9Ka1\nRUJlP4IVWHfq'\
    b'Cmw/TmSOHMQLV7+t+rg2MpiVv6B/i5tYw5G/a+GrqYG0cj/a7WeVSvhiP6xBW7tj'\
    b'\nvXE/e/MOIoB+Jz/rFTl5iio3P3TB0oWIaRM/JaqluUEoJT+HhslDFuwZP5vJ/kz'\
    b'ZKiG/lJl3tWZ9\n577Q1rN5i8Ivv3P0+dK1KPk+bCIc4jV9F79S+bm92184v4+1Gv'\
    b'8P7WQ/87B/VhSHS7+UbkLDUcsD\nP3mbKRDNCle/i+Rzfx8AHr8Q4jy4lZVFP9HUi'\
    b'vkfjFq/JIFcLjcQPT+BX9wTdWG3vglTSXXRgz0/\ng0H2KIEVQD9uastVBmQtv04O'\
    b'tjnDgQK/4tleqW3IMD9r26HxNwEiP9zr790yziS/3AdOUAJo9j7b\n3iYN6Crmvtb'\
    b'jOupiXx0/33C2Ejb/DT9qLMdbHEMsP1CoRe7gJSe/fzd54o2SML9dwOY9hzNsv9H'\
    b'T\n7V2l+zS/7BsBCMjxUD9HQn0Ad4diP+LX0CAq8j0/zk0+8w0tUr8VNC/cudhyv1'\
    b'b3fKXspFu/jhS3\nBHMUGL9ilXRpRauJP0x7fo1tUZW/WyY3w+Yimz88rlCZx+BwP'\
    b'x5gzTI27ka/4YcheBHNSj8Qb0p3\nTpFjv/O096vyj6u/b6GkWinkcT+6E0lSAtRU'\
    b'v9JFL+nmqWA/xGt4RjgueT/92oPQSHjJv88A0tSM\nH+k/v3MgQCStyT9XR2F/yEX'\
    b'KP9WfWxYaEcq/eA2OkptDfr+wkezSrSCkP1Bkdwl/UNm/DCOrcEZu\nzr9LMtdGkg'\
    b'+kv7DGramncLq/rE1HHcivej93mdIJhmSAv0n70gNXMHu/i7oO1bnYFT/YRKk0lx'\
    b'aK\nPwnOoMK+zIw/ngbRPqgBcr+f79qMv/iSv/THPJEXwo6/VCROXr7Qcz/6CV7pg'\
    b'nI0P4yslI0QfDY/\njkCtHf7QMT/5ffqEIb4ZvxOXUsw2FgU/4MVNUj5zED8A2EdT'\
    b'gJMwv9VFHtf/tg6/TkfkGda19L5h\ntNCgxTlFP24UCkfK80Y/P9i4ugCrV7+ch6n'\
    b'dg+FCvwBSxdwCHzC/mTg3mRP3QT8JzD99f7lHv7iU\nND+pm1C/1jGyhg6OUT/cbC'\
    b'UrgtFKP2ZNpepGYk4/fMpnslj8OL+liCPCeBsRP5aHP9EeXhC/XVUA\nunxyKz8kH'\
    b'zV2+3cyv0q/qDiwTRm/Hnogza9J2b7Y3Xcq4IwLP5dNPdfIY0G/3XnNQWKfLr9aM'\
    b'ZyX\n49BTP0hDfqa1Hzg/+HyL4pBdMz/zZ0kqlwJyP6De3evXP3K/t27PC6JhVj/Q'\
    b'yft9ED9qv7XxNLv/\nz0y/SBys8y9LST8z6OWCmF2PP6cMP2a/w2u/XCsTmUccgL+'\
    b'luFU1cQ6pv6GOOdf5+aI/9iybehFB\nhb9WziPPIoJUP70+tigMe3q/n6JG7bZvbT'\
    b'93axLHFMyiP5w9cJ9eQYs/uGM+opGevb8HRHNL3fpm\nv/H98P8sV2m/G2rlxN2Xk'\
    b'L+O+pdQv/6wv79zIEAkrck/o8JqATpb4T/tBGjexAiRP1iV6ImIp4S/\nZfje44TL'\
    b'tL+PBV+ascSgP0LQ/zDhfs2/MaoAdvfVz7/TFIXNeMiCP9+BuQADPoE/rUE3sE/B'\
    b'sL/e\nmoaWhe2YP4moUUyuxpM/hvH12pO2cr+2uxMTZPdvv1CqEPUFf3q/JKy96/y'\
    b'AJD/wHrrTYjiSv5E4\n4t9uMIu/3x8HA+IWYD9ZR90t/rIwP0l5lXB7uEO/wbeymu'\
    b'qETL9xt2JudgAnP0Hb4nvuqR6/+8wC\nI7z2ML8PtjHqI70pP+e3ZsRk/SA/g8qWS'\
    b'rKhCb+kHwX8JyROP7/Wi/WzemG/AhcPGl8Ufj9Jgejt\nVnYwv7RIHFEeoSs/pDiP'\
    b'dl5BZr+nNX4ajzZRv6NPVFIDSWQ/4junwTbcc794lzP60/xOP4YDcdf6\nriA/okI'\
    b'z1v1UYz9AarKs3Tk7P1nxBDyumjQ/rJKqIGlOPL9EI1diQ7o6vwXTkq/8Fha/s1W'\
    b'UWfrG\nP7/K0lvYBFcHP8lArEcmnEu/P4mw+is0Pr+2bvrpspIAvzqLLxj5Zxk/OX'\
    b'aGVfv8Kz+aBBvTUmtT\nv/hwTQCQuFE/uM1vWVkhHr/UUFBrMJNNP1BfzG5u5Tk/E'\
    b'JdfZ+uMQr/CMLEADIFNP32xwGFXeBo/\nT5g8HTNbRb9BsfwPOcZkP881p+NsrFm/'\
    b'en/Sqcy8Zr91yAQ9uatDv11PhVKjmjI/4h7lVMAsRr8O\nOp/830yBv2NTY+gMdIk'\
    b'/oJUiKEVxYD+LnU4iVTtZv3vsPCO51l8/YLPxoBWNZr8TeGBcQSzcv1dH\nYX/IRc'\
    b'o/6wRo3sQIkT9YK17pCkjdP5dN0W6tQ86/tzwOACVnbL8vMfD24IaBPwMfdbRriD'\
    b'g/QBwQ\nmO8Rfr+IKi5wY9yMvxoPx8XQeZY/NR3VYAT7ab9+vFpND+5Cv/1tKAgki'\
    b'zM/8t8mXJAITD/71Faz\nTNlwv6Y0b/MHZHy/xFIO0PYaZD9AkjwrkF43P9bvCur6'\
    b'Fkc/hLHi5fT1QL8Z4Immr2v6PnIF4i3d\nxRI/tKn0xzyNGj9iqS4JSL7qvtVoyy/'\
    b'skhI/oOHqrXXb+T6kKYamrXkSPxiBdCVqTxO/CPBsOqEa\nDD8KkZTUnmokP+MOxS'\
    b'6BpVG/8ldQgDrDL7/zMvP7t64Lv/+10Fxhdkc/HhrPBGjJMD8WkmfTuyAy\nv56Vm'\
    b'mlvhEA/RmXQ4ov+Jz+z0poVIwUTPzdr51aQOSa/nXl/g0EEOr9YgkIEWGIXP5V7f'\
    b'1gfMAG/\nqD9jCkLBGL+BTC761foXv3jRWn4hpxQ/SQ4GHv2S8j6XDa3dDcwBP3N6'\
    b'/gDuyw+/S3rACVIrAb+i\nSgU11ZTyvnFdv9opjxS/Df7AxiVBKr835SwsT0ddP5T'\
    b'7UCNb/Di/OmdZmIZ6M7/dds+BE5s/v1gF\nBJwHjVi/ffJdl57rRD9GNM91VwA+Pw'\
    b'ZSJ5y1AFA/qk7zDiNACr8tEETeQcZVvyQ01NhkrXM/C8Rh\nBeDlYj8/l1T0BqI6v'\
    b'8kIHnP1njc/YoSSOrx+ND/IxHGWEguRP2CMK1ZpHnK/RsQUjP3jYz/uvukk\nQ19A'\
    b'P6ax6vbsnV+/mc2I86LfYr+QDMCPhRvKP9afWxYaEcq/XZXoiYinhL+XTdFurUPO'\
    b'v6i5cz6H\n5sg/hMYokM7sez9b8yH7qK2cP722Kc0j7lA/oXeCAxIaMT8acvRB8FG'\
    b'LvwZUUSQw1IM/3sCPa6zL\nUL9zjWbkSChWP8iOYznDckI/CmwyaFUASb8aDL0Dqe'\
    b'RlvzkF4nEykmG/V9kFO7fTVT/2XJrm5eNg\nP76ocsMqOWE/EobZZDbFRD/G2Jtp0'\
    b'Mwqv17NShYjPCG/NMKIuZ/AKL9li+XWG+L9PnGwocULZwu/\nE38O6e4g8r6oHP9f'\
    b'VqcTP06jMwr42Bg/J4pH/bG85L43/zkHfhgjv4psyPTLsgE/xi4u42fGSD+A\n7HR'\
    b'yuqkhP1f+ULO5Dya/OlafIEvCM79m6q53qpouP/RJ5BGndDI/1gpdAZajQL/A60O'\
    b'xDy8uv3Wn\nJKho0yi/5DrrvO21Mj+vYUmeS14Bv17Ps0xePBI/s3bKgrY+DL+FrE'\
    b'CMZ5AdP/mpQI024ei+Nn+t\noQrn+L4WB3IX+IgNv+ySr6P6XR8/6kv20+8FDT/8o'\
    b'9qyVGYkvysT7IBDUfC+oiOfHuoXIz9UQSYE\n+ttSv1n+Cr0l5UI/m4lvf5gA2T4O'\
    b'a6a6tvsxP9JC00NPC1E/30wjUndwN7+/w2YGFV5Zv+KJ2C+4\nmho/5fzvZrxMUz/'\
    b'g6DtjOCNzP9mUQYRLZmy/eh3Wj3EBb7+1EXlUf04/v6kY2HpRJTk/cCvQ1OyL\nRr'\
    b'8VtRzyR8mgv6VM0zaXqpY/+QFYyRM1Kz+PHQV4gZc2P7rX9XdkJGa/whuS8ktTcT'\
    b'+a4ivJLxuM\nP3wNjpKbQ36/Zfje44TLtL+3PA4AJWdsv4LGKJDO7Hs/5JEOPgaSr'\
    b'j+kfw9+qSySP/JJfpCDvJC/\nwU4ybDZnhT/g22MsKwJnP9yCqv+CE2m/pNfbXozw'\
    b'gT8l81hNUJxjv5o8v7TCDVm/N7GEuPe1RD8V\nSiO8zh5QP9plXxOXL10/D4Myoge'\
    b'PSD8iecYNX9NbPwhYvsmMhVk/q+gyZSJEHL+g/13/95kNvwWs\njVA8myA/xzv9LM'\
    b'9yJz9ie3BN5G38vp7uwdJ80Nk+rzFZEQwTDD+DWmVr09gMvwXhedVh/e2+nA8i\nW'\
    b'f5b4D4WfYy6NhQsv/Vx0uJi5UA/j/u80XS6V7/UBJMR4PMMP+ajnoDMUBy/ifxHd'\
    b'C12QT+/j+S2\nPpAuP/gkzneRFEG/HG2aXvOhTz/8hQEtNZopv3ti9NnMJuq+V1wc'\
    b'lvuKPb+76QQhlcEUv44YrRo8\n9QW/bQkR6bfmFT8UfjbGp9sXP2AxHTFaheQ+/nw'\
    b'pJGWiFj8zXIWdCoX8vmYpuJioLSQ/Uhayh/ZW\nGD+qPj2+se4FP9P6DYtFu8e+Hr'\
    b'aELJy6Er9i866W/mtIP9XYn3cCzDY/iUwLemGPKL/8tKLBY3JE\nP5ypQCKdwVi/E'\
    b'YT8aKEWQT/K25v/rEU0PyT2sWE040M/v3qcLI9YZ78xwuonqLBov51xm2EShUQ/\n'\
    b'DxyYxbcvTD/hnizDVJQfvwxOtSPIOj4/nl+8xNB6KT/X7njmhIKIP+XUurvR+Qc/'\
    b'9oOjOsmZJb+g\nufNfPPpUv1GosC6/Q/y+N2W2LNKXZ78S9FXboG23v6+R7NKtIKQ'\
    b'/jgVfmrHEoD8yMfD24IaBP1rz\nIfuorZw/pX8Pfqkskj9Obdo8/b6wPwwo2+2Qw6'\
    b'i/ooVKq1kgor/q/Wl89N2CP32QthD4a5G/EMp0\nJy+nhr/t4fFB+pA1v9si9OYHX'\
    b'xO/+kq5q781NL9sSQIC4lloP7sl8aXFOVk/alVy7LOuML/QrqgG\nF5duv4jid+ue'\
    b'o2+/R2RXflkGVj9PxBJG5VT6vlui1y4Lcuu+6OnPsgXXG79l4qikjwD1vt7mAsYB'\
    b'\nsuE+xy8Wj/N+9757rcoVN+ETPyhPSYmV5Pg+C+2gFrR8/T6IgoGemwE4PyR8mYF'\
    b'NylS/2V2Am4Ar\nUT9bmxW4l5UVvyOKB4rTT0I/Xu1NbBpvJr9x9sTiYOw9vz0tpj'\
    b'5TOkw/G5BWxkrKQ79P8HNS2xgx\nP8G5W4Um+BC/00W3iltpBb/E7E6dwrQlP/EDL'\
    b'5h36BE/LszvXNxGIb8q8Qv3pqojv7HGVFuQtP8+\nuKlVkvhHFL+2Z7vDKKjuPjoU'\
    b'55ozLDO/e4B8JTUlJb9mS1JWYXsMvzyYnsYuEwk/C73VaZcWBT/n\nwuP16RxTv4+'\
    b'tuHUnV0M/PdNOqXFyLT/IfWejQ2Uiv1SBNgxMwVM/sbEAvovKSb9jXUQP6F8yvzQ'\
    b'D\n4MQG/Ek/HX2sSd78YD/I+9Nnv7tfP5PD5GjnJGW/+1teFYRxGr9ZxcKCttEoP0'\
    b'7bRxXEHU0/Ezbo\nZyhdWb8S7sIJv9F5P8pGFolHNIY/SzPR0wICdj9Li4xMxMxWv'\
    b'w2Xxh53Bku/Oir2qro5Y7+pSixE\nyL+kP1Bkdwl/UNm/QdD/MOF+zb/yHnW0a4g4'\
    b'P9e2Kc0j7lA/8kl+kIO8kL8NKNvtkMOov3A2SOzF\nbdk/dZAi5nDB0D9o/Dr0M6t'\
    b'/P2Cig56EjoO/XC0st3AblL+Ipxt9oSpWP9VDvOXOWiY/Tq4JwYfn\nJb/u8rQrSf'\
    b'tBv9pjs/+mbFm/ctmCNlblIr9HD1rV2GB5v8sEV4KV82a/91jmiwaKYD90Z34W4r'\
    b'Yi\nv2mILlLHh+C+uJW7dr2k0L7on8QJtSHnPsJpt3rIV/k+EhouDh24Az/Y4ILNq'\
    b'fMSPzeqxZTXtwq/\nycGhRlD2+z7ODIPBeYlBv575YhFjLTG/AjhCBo2eTb9wJF5m'\
    b'qqw5P/XSnWRfrDM/tvhphmy/Nj9n\nefA5BIs8PzTi/iewsBY/euIuT02/PT9HATo'\
    b'5SkxEv91S+/JvdkS/ZjLY+rTgM783J3gONd/1vq9Z\nlK95wua+HUwynlFVAb/4Jc'\
    b'bXT7woP6gvRQCsjBs/0NKFAi5IFz/DLbItyX0Kv58dTqt0mjk/l351\nXbgmMT+a8'\
    b'LWHUP8bv8FW4qPzXf2+nmoD/d1s8T5FRbKQNlBjv0h6VlGTy04/RC756a5kOj93I'\
    b'Rf9\nTZhTP3+idnoOWVU/FwfWUkJOUL+3X+Vvo6RivyYdBFMVlD8/KuGcWTKGUD91'\
    b'w7BGqW1/PzqADpiN\n8YK/ZyCvfSBLZz/qEdQgwsNHPx1Vm9xASg0/OBXjH2sBRb9'\
    b'tNu70XllmP1rvWHQ6/Ii/gW4W/MCq\nkr+05lu4dAZrv7LlBv+sgmC/UY9mumlRaL'\
    b'+sDM+z8NWhPwsjq3BGbs6/MqoAdvfVz785HBCY7xF+\nv593ggMSGjE/t04ybDZnh'\
    b'T+ihUqrWSCiv3WQIuZwwdA/CshCxHgF0D+l2Ehb/a53P/88ErbNPFy/\njWvMASE9'\
    b'cT+JvZ9IVHtdv68RMzudDme/defC92cDTz9d2mXH9yhCP1hhOwB0XU4/m/ZNF5eq'\
    b'Rb9R\n1EkzxiQ3v57cxM4XiyC/mKB1nBnuXz97yyFjW2UUv9ggquwY/xw/cQl2a5t'\
    b'kMz/xP2rzbUIMv6z7\nN9Id9go/8O79QIOmHD8qEwZZ6lUSv0Q54u1phRO/vVD64K'\
    b'CP/T6jXj9rP3c3v7VM77LOzUk/sl2D\nMe1aZ78ZtK8N+ZMiP+hn+XxDeg2/hlt+z'\
    b'aLjUD+45CLmRps9PwiWv9aTEU2/Tqfle1ZEYD94jMQW\n5ZFBv3bZsSiJLS2/1H0I'\
    b'EhonU7+Ee3L1kLIuv6eUyKNjZie/f+CeHP9rJT8Ly5uo99wrP7rHGpAf\n3RA/rPb'\
    b'0XfcuMD9PB9VwKMQBv108iPLrTEE/r6VQ3/nSKz9bQZfLWMkSvx1YtvlIvAi/r58'\
    b'g1uTf\n6r4oIfAiBP0Hv+ldNxo5HBQ/TDjIayiAEz8dw7bt3+YWv33TGmUlYhs/YE'\
    b'0pKuHuD79QqRxpvQcv\nv1ZRG0LwW0m/zwJxtM/BRT/i0gfZ+LxfP7CLLGs55FK/I'\
    b'r85hV6oRD/iQ8eEyjkGv1segG5sZx0/\n2JHK4sDU7L5K49PnOjuAP+k1JpvnAVc/'\
    b'MDk4dnTMeb9r1+Jxprc+v49O97LZuVg/u3NNCxMNZL8x\nMQxlWV6/v0gy10aSD6S'\
    b'/1xSFzXjIgj+CKi5wY9yMvyBy9EHwUYu/2ttjLCsCZz/o/Wl89N2CP2f8\nOvQzq3'\
    b'8/q9hIW/2udz/4/04gU7TdP3pT9bFkj72/kuOCtccknT92YeAn0lKov5TyHrQ/WX'\
    b'k/sniD\nAdA7kb/p8t0fwM/Pv5jT8mTIi8A/Nc0fkIV4pL/HXSoII3Cpvz6p0DXaO'\
    b'pg/L7pIhxd9kj9DMj5c\nSA79Pv8TQeO90Nw+OYJ8+DUa6T6Vq3ePWNPFPm1uYgHT'\
    b'+8A+OlK8uzm7ur46op2zbHDxvodjtbcd\nMvG+piqX2lXzxj7P6lt1tbv1voqHL3+'\
    b'h3xk/zBsF6Wn2I7/1pcZWq7DtvnpyCtZoLfq+XO7vjOB5\n+T5yaXXPZjUdP387uR'\
    b'hRGiG/lrh215m1/T5XuGfNDELIPih6PsHnuAC/epbHgVmMID8SilOUaU7B\nvoqc0'\
    b'R6yHs8+R4BQV7zn5r4KPTxctbKHvrUOMDK0Yqi+r7daPp0T975qc11JEI+nvkDzz'\
    b'g8BGcg+\nUJRNsDx54b6ZBziRQ1cEv45ZbpNYZPM+w8vO7UjE9T4kLYFWTsQaP+Vb'\
    b'7Z9YXve+t/a8Ug+vCb9m\n6ciaslBEv8DltyPNJgm/si3cn3WGKD9yHeR2P7xUP21'\
    b'1juJ6vFS/jXYOdQ6EWD+3N/QYNvo6v8fh\noePttUY/RijtWBlpSr/9WJuIyyAZv3'\
    b'35A9zO2Cq/0vvzwFN4Aj8E/uGRIvOAvyPmK3/id0M/9u+9\nCw0lUT9c/A5OGMtSP'\
    b'0aOlMTHCk+/JDljt5wdTz8T6y9KM4Wlv7HGramncLq/3YG5AAM+gT8bD8fF\n0HmW'\
    b'PwVUUSQw1IM/4YKq/4ITab99kLYQ+GuRv2iig56EjoO/+TwSts08XL96U/WxZI+9'\
    b'v8/GbyDT\npOA/gmeT/VH3nD/dfSKl8nZ7P0oKbitMfrO/k3UbYyW3tT+S8wlRITr'\
    b'AP7SUtgTLV8C/SrIYbVVR\nlz85ZdHI6sCaP9wOtxvlJsu/y63tE5w9wr8TJbW/w6'\
    b'H+vt7fklq/NQu/Qhtmw5Rf7L4DOkTKwj3w\nPmcg+2YVzAM/irMgvLDS9b7YK3E7Z'\
    b'HEDP+rSl1XUfgK/uouPSwDmAT/ULbQ/860Gv3yqMqH0cjW/\n4G1qlpfiLD9jziCM'\
    b'+LkFP2rZSXc+yjA/WYnb+R000T6D0heDrxbvPsqc+70ZgzA/bSmEftSn8L6E\njUz'\
    b'lsRogv1YtQSTNlxC/ZA0d31HdIr+z1H7RC+gOvxz7b41GIQi/pBO0HcED9b4NK6J'\
    b'nVw3QPhn4\nRbC3aP8+Hp/3f6zAAT+igyJePYT+PlwbA9RnYhw/8JXmU0QF/T5GRO'\
    b'gx6E4gv9F49uEsCOG+IP3D\nyrK6tT6DFQ6wtP3YPkiT/q5q8kg/WRewT+ExOL91a'\
    b'5onBuZGP+/wAb+HCOC+wcpo5AO6Mb9G8+Fo\ngLRgv+86iuZjDWc/sxctZeF2IL9u'\
    b'bMOR2TVmPxP2fWGz5Fq/jdFclVTAQr+KNc9RDtxEv97MvJpJ\nFEg/8sz9hFDaQr+'\
    b'L2ROhdBKSP5Ee2MEha4o/1iZTAW7ygb+aHqfGwgVOv9SeoauvT3Q/hfinX5mR\nab'\
    b'/6413wMh11P6pNRx3Ir3o/r0E3sE/BsL8xHdVgBPtpv+HAj2usy1C/n9fbXozwgT'\
    b'8RynQnL6eG\nv1otLLdwG5S/oWvMASE9cT+S44K1xySdP4Rnk/1R95w/0GwJPrb94'\
    b'j+mes/hXfeRvzM8rM8+qLU/\noSaxtBVS07813NAoA+KjvzMELx8S2pY/Od8XWNK2'\
    b'rL+i5evTGTCUP4HQzdvaHMK/eDgd5ROsxb9I\nXvb+ETb6vg7wQVwhZs0+txr5jjc'\
    b'gHT8AS6fB1jj9vkTuzDqR2ek+icZy/pvp+z7LDe7MN4TmvtXC\njAlhHvG+jejYAj'\
    b'gh4b4IkBKgW3Ydv2jaaYuwITw/wg/UXSoPUL/edO23Co36Pjwkgomg5yK/3UWy\nW'\
    b'3LmMz/f+A9m3J4oP2DyHzMYnzy/PE/nx1+OQD+JXBjEo9Mav/eqlUVvXgc/gJ20S'\
    b'A7yDr9V6rbi\nRg4FvwaXa1tdo/e+npLQTxziGT/PctoIYukHP0QmBYeur9m+imGU'\
    b'c5NuAj8Op9H5RNvdvlwnZ6UK\nehM/BYoc7vu4DD/wtyu8oz4VP2LiJ0QAOwM/Hg1'\
    b'/5HKYAb+R7A0kZo4LP6xYVB3dRMO+Oe6pig1C\n+b7a6oVePY06P/mv1/zAlS+/4l'\
    b'5uTMk9EL9U2J1/5vQpP0n3qdj6sFA/fpQCq9U5RL/zzxHJ5vth\nv+/F+l854jk/E'\
    b'cMU29cUNT+6cnwxUWMyP6c4e5Vllhg/i0QK5NFGCb/kZBYYwBBtvyLXblpNt2i/\n'\
    b'8k9x3lgQaD/gR/i+GoU4P/RHaRd9kVW/TSsEQM8wTz9f9Ka1RUJlP3eZ0gmGZIC/'\
    b'3JqGloXtmD91\nvFpND+5Cv3GNZuRIKFY/J/NYTVCcY7/54fFB+pA1v4enG32hKlY'\
    b'/j72fSFR7Xb92YeAn0lKov899\nIqXydns/pXrP4V33kb9ZBAQ2zzSoP0DLwUI7dW'\
    b'6/dfyg0GzxkT+TAobWsyFvP/IMjoRk83g/VZPF\nKK1rmr/MrO7aBNpNv8kAtVt9f'\
    b'0i/bSzmuWggYj8zHrUZAuEBv+ceJV4o+u6+nYjCywNs7z4kO5K0\nLfHXvp3i9CbZ'\
    b'W+Q+58kwe/Ax2z4phazxlb0BP98QafVQENi+vKyXoJY03L7S7+8phvj0Pk1pPyWx'\
    b'\nDRw/sB9TNVCzIL/F5X/J/R2zPjEZe71nHBO/an/OkYIPAj/MpG1sAhK+Pix2Dcg'\
    b'Bfhq/GtMJbT5M\nIz8v+YXu2ETyPil3do+YkB8/Inhnf4fiIr8/JMlESzoGv+86sp'\
    b'uQd+u+mNJMAmrWDT/JWj1eMAvg\nPlnn0WgQh/y+yMLzef7mAj91XR76CkjOPp5Px'\
    b'G3ROMC+Ei5dH6IH8j7HnqsPgszsPh/UP6fOqNG+\nTWWxZlcVtb53vBlY5zv3vuEX'\
    b'5kDSVSS/iKNG3YXLCj+afDrtJkUvv/Vg04+Gufw+80uzVZeTEz8+\n7iqrN35CPyY'\
    b'OIKtGajq/HGYZAVeHKj/pUb5wxX46v04Y4fX5Njs/LBVbdEOeJ78TPzayTJYbPwa'\
    b'3\nQBMU9Sy/N/BPnf8EHT/KRlYpzLNrv94gwAfC1Ee/P+NwpJKpWj+ExpeDC841Pz'\
    b'LGSXwOBVG/iV8d\nSXrlRD9/FVh36gpsP0f70gNXMHu/iKhRTK7Gkz8IbigIJIszP'\
    b'8OOYznDckI/nTy/tMINWb/UIvTm\nB18Tv7RDvOXOWiY/rREzO50OZ7+T8h60P1l5'\
    b'P0oKbitMfrO/Mjyszz6otT9Dy8FCO3Vuv+5bkyue\nlLM/7zA3slOUt7+LC5EWT0J'\
    b'ivwpotN2fo3C/lJyI0gOsij8RGGQ+DhBOv5qLy3jYAog/W7JSJAyG\nl7/r02ibKi'\
    b'X6PnP4sx620Xc+StQfOWt/Ab+iW/IZI7TvPn3eCMxuaOK+obmFU7yb9L5TvDCHfQ'\
    b'rl\nPqV6KLK3FNs+0DI8b12r075gZwS6WV/2Prv/hlfrmSy/ouBKVWEgOT+5vpJuI'\
    b'vfgvv4GZp/OrRU/\nt2COWgJKHr9+bz7353wQv35fEPdZ2ys/8g4E2C0hNL8yrKoC'\
    b'kqPWPtLxDIKziBS/7hOE3RZ6KD8p\nGTdsFdkVPw/fl7QzcQM/3FfoUdFcD7/g6V5'\
    b'XN4z8vtLdfzDgz/w+fN2rcBj0CL8Dxd/xzErKvgnC\nxnardQu/lMgL+hn88L6EqS'\
    b'6CFwnbvpSTurw1PKa+yTqq383W5b6fw4pV3m52PrJy2aaaExA/8rDU\nXA32Dr8WC'\
    b'KD+6sggP2yZloM0zxG/Tb9OGC0j8b4KDKIZhqUhvw3gMqK0sSk/GC08Hdt7Fj8Rn'\
    b'cYj\nmIEkPxuv26Belym/kfdo1J8GFD8zEJH4wnMEv9tdizjytiM/IqQ8k5SeD7/H'\
    b'8tz/TUY7P41OXua3\nZzQ/2BuUJ3MFRT/FwCUSboj9vlL1NRYnJCw/aeVYA6CnHj9'\
    b'NZI4cxAtXv2u6DtW52BU/iPH12pO2\ncr/z3yZckAhMPwxsMmhVAEm/N7GEuPe1RD'\
    b'/2SrmrvzU0v0iuCcGH5yW/cOfC92cDTz+xeIMB0DuR\nv5N1G2Mlt7U/oiaxtBVS0'\
    b'790/KDQbPGRP+4wN7JTlLe/9qamYdhm1D+0PFR+7aZMP7bSbrz0ME+/\nFReYwQCl'\
    b'Ub8vVuwiNspTv466YiCCPoE/NNc87NvNi79OOkx2I5DVvrZABzlI+/Q+6fqlky7p'\
    b'5D5o\n3aBT4jfXvoTw7UO+LcM+VoujVKR/4T5gdpHC2EDXPtMpiJ4w3tm+kOQBEUw'\
    b'rw77EGxDNg6G4vttb\n3CH/he6+SJDj8TitFr+2W5yxjhjpPuSaZEaq6vM+6P2OWQ'\
    b'V2wz4XLBDsbKnqvqdvitDed+a+7x/I\nvjqP+T6lTs4WB3LrPscn2KQnVu2+Q84EG'\
    b'l3b5b7A8ArzuPvYPorl+/4Ndd++3VOpHakA6z67eTdS\nRNjJvgnJhPdSLcQ+ZLyf'\
    b'Pdxp4z4/sWe6ECCwvvM2omW3vK6+jcisshds4j4WbrcTIlEBP6gW87xg\nxQA/ndA'\
    b'N+d0J+74aeUS+ijsFv/C74tfwHyg/U8nNpxXEDb/FcwtBcwM6P6gndl8grSM/RyE'\
    b'aSXhT\nMr8G1mQ0Cf9Wv/tSYTLH7lo/+AIwyD5YU7/LAMNXu7Q9v12IoDB+E1S/eW'\
    b'2RGvwIUz9gUc+M4Hnv\nvvmqx/Wl8Bg/99WCkTWY+r5MAc3XRZZiP1QblChy8km/u'\
    b'TXGOUznMz88SjK2Uxsov8U/co9nzT0/\naUU7VXB4NL+s+rg2MpiVv9dEqTSXFoo/'\
    b'trsTE2T3b7/51FazTNlwvxgMvQOp5GW/F0ojvM4eUD9u\nSQIC4lloP+XytCtJ+0G'\
    b'/ZNplx/coQj/p8t0fwM/Pv5LzCVEhOsA/NNzQKAPio7+SAobWsyFvP4sL\nkRZPQm'\
    b'K/vTxUfu2mTD+MtBByRebQP9fAZS9y/cC/mwNiYn/dpD/6F8gd5vdiP+5wS/8IUF'\
    b'u/JOyw\n97WX+T7Pe2A81S8FP80lmp4YWOU+3ZKLhphEBz8Vt8tx3Ez/vsuHu5swg'\
    b'vC+gpdW2rf13T5aythg\nadMSv2D01BDo3uQ+6Goc5rWwCL/V5kepONG2vutpwi3r'\
    b'OT8/uzALUizrOb/eWg5GeoL/vkcanHqu\n6jC/ZoUSPiEVFT9go6jnXaDrPjMWaU3'\
    b'9UDu/MorLOeTBHD9hs0v/NVwPP/XwCZ+TCwE/eTl+TAHZ\nLj/KShDhsR4AP68grZ'\
    b'ePKvs+DmZb3Tav8j4hBhDd0A/2vhoNmeLInu++vslGjCAFDL9h3Ijt77HZ\nvmgZn'\
    b'6KNlQu/RJWrJVJ3A79NyHiQDQMMP+hbaGsvqvE+zAkwq8AXCr8LHOQHlbgxP0Exz'\
    b'DDj2TE/\nc2797jJgIr8PAJt+m0JFP/DWeRKUVke/j9lPvVZjJz85a8czNE8+P84K'\
    b'luJYpFu/0eNp9ACQKb9Y\n2CcGcklhv8egQpJkgC8/dsQWr6WdUj9YcwBx4wAzPyy'\
    b'SRJeQjUY/4kXvkk5DKr8mMVSzAfFhP2s3\ntcTv91a/VLHGN30tPj/yeS5RFGkSvy'\
    b'4AuGd+UiE/qF9jclaVNr+gf4ubWMORvwrOoMK+zIw/UKoQ\n9QV/er+nNG/zB2R8v'\
    b'zcF4nEykmG/22VfE5cvXT+9JfGlxTlZP9Bjs/+mbFm/UGE7AHRdTj+Z0/Jk\nyIvA'\
    b'P7SUtgTLV8C/MwQvHxLalj/wDI6EZPN4PwlotN2fo3C/t9JuvPQwT7/XwGUvcv3A'\
    b'vw5zlkYU\n4cA/yI4PEIl6mb9R7om28g2TP3BrjgyKuYC/odV5vzppej8EC/Kdmu7'\
    b'yvofFs/SFMhE/XWDOgGQe\nwb7b8XWsjigCv+rrrq4P9/Y+fAxGAwL14j7uGDs9WH'\
    b'cCP2m96Aio+r4+SAFGuhZM7z454FoimoYV\nP4HSaIOs3kO/l3+9EfwbJD/mg6Zsy'\
    b'AT4PkZLXSRsazM/RcWdwZxw/j6fOSlKhRg2v+3a9BLuQjY/\nLrTObQh2DT9jCpnC'\
    b'TlLrPsjrbsHymPQ+ftmwxkK2OL8RI5IO1xf5Pjgk+3wAh+W+yTiFKLAa+z5w\nxQd'\
    b'b4vbuvpwEp47QUfc+1/6Wzw7FDj/3EJm1p+qbPrFXgbjaQ/u+o1zKc2N88T4NoeW'\
    b'oiZ8Dv6fR\nuKhX3gC/1lf2caSSEj+NZKEcrAcgv/ZFowtKTje/9eRR1MxmLT/pu9'\
    b'SCFHRNvzgSoP0vxj8/1pDs\n2j1z7z4FGt855vctv3cT2/MkrSm/PzSe7AP6Rb/PJ'\
    b'JZvYX9fP6q7bc2IikQ/yt68RASVSr9jYF1H\nQa4mv0x1/+3/9j6/SLFzD+IXNT/u'\
    b'nRRKlPVTvx6sJYR8XUU/SW+oNjpZSb/jNhc30GkOvyslbFHg\nWSi/vppmfTAE6j5'\
    b'q4aupgbRyP50G0T6oAXK/+au96/yAJD/CUg7Q9hpkP1XZBTu301U/CYMyogeP\nSD'\
    b'9zVXLss64wv2XZgjZW5SK/ofZNF5eqRb81zR+QhXikv0qyGG1VUZc/ON8XWNK2rL'\
    b'9Vk8UorWua\nv5OciNIDrIo/FBeYwQClUb+cA2Jif92kP8mODxCJepm/rFytmO7Lq'\
    b'z+4ZiX0dIyTP/5HtfysBoG/\nBvKRorYlcj/sDCbGEWobPxHDHJADXcs+DKH6jR4/'\
    b'G7/xVEqOuf3/PsYxLtKULQS/HE72+gwxB79d\ngTtVZPgPv+I+EB44Ivs+EXeURYO'\
    b'C6D5M1HX09C/xPvIprclrbUO/+ESxkyUtUT+auYITss7YPh7H\nJXDVlDA/mcUUbU'\
    b'F3Mr8Xgvnv2kIpv/nqWnGrgEQ/GunK+AtoR78ak62AZg71Plej/GBf2y+/+7hK\nx'\
    b'yhEJT//HNeBLGQbP6E1RcX1tAk/Y4jB3RfAJL+ZzO+dKS8Bv2RQl0L76As/yRUTB'\
    b'LttEr9FOKFX\nCMfevnKJ43SIrgm/ycclaYVgB7+TE9dwlNjivjZJ6O7INuE+NQZs'\
    b'XRTz5L4aUO17hgcMP2t8gDI7\niRY/gydV63guEb8FQ/qIvaQBP9BqocjdgBK/RUQ'\
    b'3k2T4Bj8kxV5l4LUAv1wUjTleOw4/yfLCxw9B\nQr8AQMORtCQiv9InSeiHbTg/tf'\
    b'MJbtsQAr8rX+Y+GikLv2RVUuwAfx0/GwJEgPfB/b4Qhqx1hQ41\nvxEbvCt+L1A/v'\
    b'0vtp4fMSD9LrWuyDkUeP+NLW1DjdA4/w7qsfSmlDD/Z7WeVSvhiP5/v2oy/+JK/\n'\
    b'7x6602I4kr9ckjwrkF43P/Vcmubl42A/KHnGDV/TWz/RrqgGF5duv0YPWtXYYHm/'\
    b'O9RJM8YkN7/J\nXSoII3Cpvzdl0cjqwJo/o+Xr0xkwlD/WrO7aBNpNvyAYZD4OEE6'\
    b'/LlbsIjbKU7/3F8gd5vdiP1Hu\nibbyDZM/uWYl9HSMkz/J1F0NX3ipP0vT5vN4Wp'\
    b'e/fFW4+UHElb/oV+vov6ugvk4yYr3Eu+k+g/tb\nWKp6+L5gZCO7jtLKvoqfnx+hZ'\
    b'cU+9blHbucV4770FpwjIjLUPj2UIcJrErY+PENnRKsG7T4a3Ekg\nkXPEPiJzg4a0'\
    b'UzS/PXjq9CuxMD8EwUr+2vb5PrjrCEqnASQ/jk5ChMo/Cb9r1Pl7y/4Zvxtd8nLx'\
    b'\nbTA/9YojWv94Ir/LQW/HeuX2vq4hRsqciRy/IdFFkZyICL92K+gIS2YIP+hfC8F'\
    b'ik+U+l+NS0IVk\nCr9tqWAkVwrdvgeMXeXsJAA/19KrLGqD1L4leG0nYCvFvnN2n/'\
    b'JkHtC+8Tgy1dLLwz5KWRux/Fjy\nvnBVcuQQZ+6+aev03GWs8b6w0xKjbw8Bv3BsL'\
    b'CL2ZgY/RSvjJsSXBb/DrytR+UABPyUTsMT7cLe+\nps0zWeWA3j5FjQx4ifMSv18D'\
    b'7E5wCy8/FVvxVgDIPz9SZtwd92g3PzVL/otNuT+/Ej9244pMHb/8\nq6COrggFvwM'\
    b'x2odsqgY/HNy7sonyIL9dvDZlHeQYP/7KERU3EVw/pdlMYFmgVD/pmu0q+QUPPwZ'\
    b'b\nUA6jBxY/Vk3ps0LrJj+tQVu7Y71xP/LHPJEXwo6/kjji324wi7/J7wrq+hZHP7'\
    b'6ocsMqOWE/Cli+\nyYyFWT+J4nfrnqNvv8wEV4KV82a/etzEzheLIL87qdA12jqYP'\
    b'9wOtxvlJsu/gdDN29ocwr/JALVb\nfX9Iv5mLy3jYAog/jbpiIII+gT/PcEv/CFBb'\
    b'v29rjgyKuYC//ke1/KwGgb9N0+bzeFqXv6l2yrqg\ngMw/p7kEyIhnwz/Q3eim49b'\
    b'yvuUj8cFyOuW+0GCiCeTS9D41XPEDgFe4PvR2WqKfB80+GDCZ0WhF\ntD52dgDGhd'\
    b'noPgCJmlLSHrK+9OzjYz3J1L4v7f0spZn7vrhuNywtePA+SNkfj01iGr97N/U+GJ'\
    b'nb\nPhMWoY6ngug+abEKjej8Dj8uGgSe3WwJP9rjdr9g7fW+LIXM/6OeCz8hy1g/3'\
    b'vwIvzi4XgDNHAO/\nH2QI4YVc9b5Rct4bTSDiPmbIhIbZ19k+qXw4rsRIvT5altwM'\
    b'uIn3PtmEuZak9s4+b7Zuf8kBij5V\nRYevCPzQvvNq1hOA3/E+Af95o8B85j7yLyw'\
    b'sNpPUviGIx/FQS9c+XBX8SaqxvD4jXVILc/3Xvil5\nHYRFmvm+kKytGSwJ/D4qxi'\
    b'u+hXsFv23IjD7VRQU/4pQk5NEz8L68dnWXaHkEv0oZ4/rwfzC/sHLO\n9VaeI78fe'\
    b'sIQmCkmPyYjg4z7txg/yyEVFZ+GGj/1IsMCspkOPyVIhWGd8/++KkWSwvsbAT8xy'\
    b'UWB\no1Ajv+U5pRb9fCA/3mCpowMgQr+SW2s2iJDkvl/UcZ4WAQC/p08hzkRWDT98'\
    b'8w4igH4nP1UkTl6+\n0HM/4B8HA+IWYD+HseLl9PVAvxaG2WQ2xUQ/sugyZSJEHL9'\
    b'GZFd+WQZWP/VY5osGimA/lKB1nBnu\nXz8vukiHF32SP8ut7ROcPcK/eDgd5ROsxb'\
    b'9uLOa5aCBiP1uyUiQMhpe/Mdc87NvNi79y57D3tZf5\nPqPVeb86aXo/BfKRorYlc'\
    b'j97Vbj5QcSVv6i5BMiIZ8M/9ZVDWBhqxj9+7n7eqAPZPqDd0K0QS/K+\nyXYz/eQH'\
    b'yr4mT+ktvhiwPsj2BosY09s+RbJuLJF6uL7O2gL00iTgvjB9HXN8Wa++M6ulNF/j'\
    b'3D6f\nEqUoKdTDviO/kE6hUQM/R0X+TVpx1z6zSEJLK9jovlswURWcn9m+iRlAqFg'\
    b'/AD/xLzTQFdPpPurf\nV6Gk+QC/tfYj3VewBj/FbYB5mQPovj8i1B/Ug/8+nj7g9R'\
    b'AnEL8b9Td0bgr7viRlcrVY1eO+VSRt\nl/ap7T6WwBihpojmPlxZMnVIkuO+dF4Sm'\
    b'NyE5T6f+EM+QQTgPgSG66NRMe8+KTQwrRw7yr75orwo\n6wpIP1gt2LszTxQ/1MBs'\
    b'4QcyU7+id0+a+HOPv5Yb8Hz0nDY/7tDmhak0fz9eQpPA9G/Jv2lvp5i5\nL4s/MPk'\
    b'6Nmw4qj8Iw1FZPV6Xv4nVoUmVCoQ/uYMdvKvKZj9A5t9mz/RVPwfXxcao21w/sNo'\
    b'b05gX\nVL9wwdHYEVI5PzwR66SBdjm/WmyKKKrGZD8HFXFtbvYqvw0meQJ4FkO/Il'\
    b'Hvc3kgQT+OutW9fU0Y\nv5Y9CirUTRI/YHe+zm48/74CFTl5iio3PxoJXumCcjQ/0'\
    b'UXdLf6yMD+N1Immr2v6PqTYm2nQzCq/\nLQFe//eZDb/xxhJG5VT6vtlnfhbitiK/'\
    b'XMkhY1tlFL9Q6D1cSA79Pnwjtb/Dof6+P3j2/hE2+r6k\nHbUZAuEBvzLkaJsqJfo'\
    b'+sLtLdiOQ1b5SfmA81S8FPwYH8p2a7vK+aw8mxhFqGz9ESerov6ugvo7x\n6Kbj1v'\
    b'K+Vjh/3qgD2T4AaC1e6APgP2A/mqrywVq/sA3Hijs2pj+U13prXFK7vyGFZT0ptr'\
    b'c/MGXm\nrZWwsL/cyKtH1Za8v0SHK4smfL2/GBebUbEglr8suv/+tPK0vzkxriMJZ'\
    b'2u/rUGTbKjIjr++ZeTp\nh02UP9UIlEFYl3A/6X0Y3hFcRD8KOcSWvg+CP9don/7Y'\
    b'F1K/ggf85ibsaL8qCIdi4g15P70yL2nO\noDU/x9tvgYGcW7/7jWW8fkYjP8nygXE'\
    b'r3xK/syo6MbTRJb8eTsCd0qM/v9KrjpLZBTs/NHpY7XUI\nUb9+YoYDpKY0vw5zzR'\
    b'NoLTe/1xlHEWygVT/FosMH5ThLv3TSmNak8U6/xIdFDermYL9OpPHUpDiP\nv1GW3'\
    b'Yx333k/ne3ohbrkdj8ZcYFk2PJ7Pzcd2A+K8bO/9tRrih1RdT/ZjAbaewCgP+MG4'\
    b'8a3ZXS/\nMtFwL4eUgr/90jYl3ARbP3APaRd98XS/Rt4RRldkYT9WWfJ9fApmv9+R'\
    b'lTOhJWI/IdSAFEdtcz94\nht3PqTpbP7+9TRwIzhM/9LVP9YB0Mb8CPJC/t/Q0Pzh'\
    b'MynzJqDe/9exJvg9MND/DwtKFiGkTPwmt\nlI0QfDY/c3mVcHu4Q78XA+It3cUSP4'\
    b'LMShYjPCG/q6uNUDybID9DodcuC3LrviGTLlLHh+C+ACCq\n7Bj/HD+QA0HjvdDcP'\
    b'm3bklq/NQu/I4dBXCFmzT7gJSVeKPruvph2tB620Xc+EksHOUj79D5NDpqe\nGFjl'\
    b'PsLEs/SFMhE/MJAckANdyz4NM2K9xLvpPmkV8cFyOuW+CvDQrRBL8r71Ppqq8sFa'\
    b'v2zpDkEy\nW+Q/aqnyjqD4q78Ce6TCQN+2P014w3C+Jcq/Q89T2QfWuT+WeUa2ZgS'\
    b'9v4ZlZq+FI9C/woXq0XmL\nmr+v3MumopBavycDXSBJKbS/lulONHfnlb8ShRH2g2'\
    b'dnP5GMfGFvdpS/sRVY5fhdpb/iuxTyplBg\nP3CyGnPCXGk/7kjBZaieoD9n9sBBK'\
    b'ZBKv20rWnVF+W6/8tDwIBhwYj9ZOs8hDUUlPy3dM2lWyDQ/\nJrhBXvNmVz9NfVGD'\
    b'eeUov9Hk74CfguY+we4dsb04PT80lYrLTO4kPxHpnl6tXhu/SyyTLscqBr+P\nEFT'\
    b'Q6lxRv7MgiudzEl6/Q3ce2/3IZb9SpgSkosmRP9r+/G96AzS/SaYJIlhJUr9o7aK'\
    b'28jGkP/iP\nkkmxdHk/1Y+n0w9qtb8dapOuCm+Cv2tHuHnVyGO/7qAXD3EibT9i6q'\
    b'6Fd3dTv8kdPfixLW8/MY1l\ntnNqSz+eDF9Eu51vP9J+YiwY42g/Dnc2TSZhcD9mo'\
    b'7YwS5JJv2voC8zZcD+/6FNV9sBlVj9HHhua\nBWIjv5lOUwsdIzk/gHqyIwUtA7+C'\
    b'qqW5QSglP3RArR3+0DE/vLeymuqETL9uqfTHPI0aP0nCiLmf\nwCi/fzv9LM9yJz9'\
    b'v6c+yBdcbvzqcu3a9pNC+lwl2a5tkMz/Egnz4NRrpPsAZZsOUX+y+Ahn5jjcg\nHT'\
    b'8KiMLLA2zvPrXUHzlrfwG/twamky7p5D43kouGmEQHP/JkzoBkHsG+iqD6jR4/G7'\
    b'9u/VtYqnr4\nvkdeognk0vQ+i6Az/eQHyr6xDceKOzamP2ip8o6g+Ku/mGPPVAYD3'\
    b'j/11+AX13Gwv1Xs4yTjvLo/\n+2DGbFm4vb9S8LaWJYGWv1Fu6UWQIZ6/0Pjgj0n+'\
    b'qr/JRFkfOR1tvy+89ESUTZG/g+/Tb4g/xr9I\nl6CJlwQ8v1OsvLTm7ZK//3qFl7Q'\
    b'Ui7/JbyKiPNNgv4QIptSu13I/08axxmMHo78OjsEyb1IvP75R\nqWw2HGQ/DZZBaI'\
    b'cpbD8xXZTFI7YBvw0sPrd6UCS/LmJLM9nqSb9o01PwxF/kvmmdhUBuWSe/49g1\nZ'\
    b'Z+sIr/oCugMYyr/vqvdLx0fGyi/3UCwWsQQNr8D1uqmGSU4vwd2i13Apj6/CHdCB'\
    b'XBYNz+YBCZ8\n1ZhWP2wXHy5C4zQ/v56sYA96I7+iABIJnLGHv2YzD8cF/pc/3+Rc'\
    b'nzjukL9mRqHeqnV4v5/SA0QP\nwxW/YYUXcLdaXT/RvUdySdwePyJYJTd3alg/UsG'\
    b'drabEOb/pWQf4+XcKv57ShnnbLDO/P3clzIHf\nJz+gzXoo2VZAvwhznHNM7yC/1b'\
    b'i6FOQmND+/xy8H2YkXvzI6fe7mPxg/KNGkqI+KFr/jhslDFuwZ\nPwp++oQhvhm/r'\
    b'rdibnYAJz9rpy4JSL7qvviM5dYb4v0+7npwTeRt/L6Z4qikjwD1viShxAm1Iec+\n'\
    b'4kBq821CDL/AsHePWNPFPu83RMrCPfA+l1CnwdY4/b4yOpK0LfHXviNa8hkjtO8+'\
    b'rOGgU+I3175s\ntstx3Ez/vt7xdayOKAK/n1VKjrn9/z70WiO7jtLKvgJG8QOAV7g'\
    b'+BUfpLb4YsD6U13prXFK7vwN7\npMJA37Y/9dfgF9dxsL/2aVB7/Tm9P3vd60e2/r'\
    b'm/aTwuVj3HsD8B3nrMG0ODPwNfxtNJN4e/V97t\nOvfYhD/gyhtb3I5Bv48yuOxg3'\
    b'/a+jsvAbv0EWD9sujLgrjMLPyqGGReBP0e/IfwpCrk4Nz/1XjwC\nzQk/PyYJff/A'\
    b'oxy/2lMvRTV9TT98YXTIzDsmPxaYBV5cYQM/ySeU85CxJT8PanQJhz+6PmKNqEJY'\
    b'\nmPG+diIVvSFe/z6SyGWRdxfQPvuUwrhwdfA+cqTS8iNi477S+Ken7j7mPsApHmr'\
    b'k3QU/YwsIxcJA\nBj88t7S0emIiP+9yAspd9Do/A6aA0OQ7Kz+JNEfOhZpOv3eWcJ'\
    b'4XpSU/BMZEJoAcE7/Hc1GmrLs3\nv8hgiYJxJmO/pv0d5SwsMz+NOrGRdVRVP5r6R'\
    b'HSg/Uo/z2+Na7tCMT/8gEMJKdgov82pKiMOBxw/\nJNHwpD5QDT/DIj239S0KP3PF'\
    b'qhuvAN0+z4zGQmctNL+4VhnZOl8jP6Fvuar1CP6+/phgAQnA5r4b\n/yiv2Z/9vm3'\
    b'wfYS6Zfu+qZfuFeOC/z6jyf5M2SohvzOXUsw2FgU/T9vie+6pHr+YaMsv7JISPwO'\
    b'w\nocULZwu/0+7B0nzQ2T535gLGAbLhPqhpt3rIV/k+hvs30h32Cj9Oc2IB0/vAPt'\
    b'Qg+2YVzAM/PurM\nOpHZ6T6n4vQm2VvkPmneCMxuaOK+7uTtQ74twz7ehrubMILwv'\
    b'i3rrq4P9/Y+zDEu0pQtBL//pJ8f\noWXFPs14WqKfB80+ovMGixjT2z4hhWU9Kba3'\
    b'P014w3C+Jcq/VezjJOO8uj973etHtv65v98rT3aJ\nT8s/irHCAXpWvL8mpJkl5qG'\
    b'LP5LNWV4h6Y2/nC3UyHlFij+2r66EuRt2v+oRv4G26XE/N0DIQkra\ndL9wt7C9Q8'\
    b'EDP2UNxqzQrU8/NzRUn9IBQL+CqzcADgE/P2Z/RvkyIUQ/7WOOxCybUL/AmmAI1p'\
    b'4q\nP0LQy5nhWjY/GFkF9RyFOr/XNJ7/Fme0PseVeEZjBeA+1hDrzVlBH79e83WG+'\
    b'L73PijLBZ6zZME+\n3DATgximDD9Q4huhcJztPtClzXkC5Pg+ZvZ8Se75Gj88bnwG'\
    b'/KMNP5iMRQlvdRQ/xhnhk2o8Hj+w\nR/XaNpQxP5HgK5G8jjG/6NIw1nu6+r7jW4d'\
    b'ntap0P+2+u5Ivw4W/xvY9v0dpeD/p5Jr6PaBpP6z5\nZRBb6hW/3E2kThcJGz92Ag'\
    b'jEWfM+P8/Q4hfUwUK/WrZBcw/0ND+grVHUbYkuv8YO3rmhoAQ/HEX4\niW92Gr+vB'\
    b'rgW8hEsPxXeGrqSPwS/t7nEkm5pHr+/IMfGCt0HPxGXf6j5kRG//UijAQza6D7Rm'\
    b'Xe1\nZn3nvtbFTVI+cxA/6MwCI7z2ML+e4eqtddv5PoF+DunuIPK+zTFZEQwTDD8N'\
    b'LxaP8373vs4ZLg4d\nuAM/9+79QIOmHD9/Try7Obu6vn60ILyw0vW+j8Vy/pvp+z4'\
    b'MyjB78DHbPuS5hVO8m/S+WIqjVKR/\n4T6CmFbat/XdPjUMRgMC9eI+Ok72+gwxB7'\
    b'/VuUdu5xXjvsA4mdFoRbQ+GrRuLJF6uL4vZeatlbCw\nv0XPU9kH1rk/+WDGbFm4v'\
    b'b9pPC5WPcewP4qxwgF6Vry/V9DRViBYwD9lEq6s5E1wP7eL96zxEnK/\ndYQ/3EaZ'\
    b'ZD/5HLXrTL+PvzkYVhBfmZs/1Guf0sGrj79VsL4ahytgP/EMr/9kSUG/d0BDGgJD'\
    b'bD+R\n62dikc46vx6Hmb81K1y/0whCuxuIfr8/7aFMt8xGP/AzWnYRR1E/zKi+6H+'\
    b'KAT/a+PHyoUPmvqnK\nG2vm6xW/CPrC+FqdP79II1SYGpDMvtLufTkZk+w+x+NTXy'\
    b'uwNb948LWO5XEQv/yAVdMSWvS+rzjv\nrqAb976BNcQaDU0hvzThJzZkMya/IzoVW'\
    b'/RDM7+OZsWElgRhv4prJTVCDl4/a56MQHCkKr92cVoQ\nfPqTv6f+NHfG/py/YZ+u'\
    b'+3a8Wb++Y4HUWiJwP7ZnwJPa9Tm/TLpG6E4/LD8ObNferZ1Qvzf/Xfs9\nY1C/Cy+'\
    b'Hf1u3Rz9GnbR4zyBMv+CQ4MK2HVY/lHUn7ZCjTT+ZVG/WIcEvPyzNuBd1qzY/WIp'\
    b'9HCL8\nIb8g7vxVcrUGP99IG2oNaxa/hbryhkir8z6K1rN5i8Ivv/vXR1OAkzC/tr'\
    b'Ux6iO9KT/VKYamrXkS\nP2Ac/19WpxM/qVpla9PYDL+WrcoVN+ETP/Tggs2p8xI/Z'\
    b'xIGWepVEr9Kp52zbHDxvngncTtkcQM/\nqg7uzDeE5r7Khazxlb0BP368MId9CuU+'\
    b'4nyRwthA1z7HythgadMSv3UZOz1YdwI/A4E7VWT4D78Z\nIZwjIjLUPqFvAMaF2eg'\
    b'+EdQC9NIk4L7dyKtH1Za8v5d5RrZmBL2/UvC2liWBlr8C3nrMG0ODPyek\nmSXmoY'\
    b's/ZxKurORNcD/WceQjWh+/P+iX8BiL+L8/Thpex0gBkT8RSuFA6Fsxv0hFMOK38j'\
    b'c/msrF\nS4cGOD+UIUaLGQA0P9m5hhBVcUA/EequIOebUz8Ldh0nJckiP/vuzaj+R'\
    b'zk/1NVlS5GpND/Ep05+\nsnLPvopei4msYCI/CVUUOn1MLr+lbpfTpqfdPrxrUp89'\
    b'Pvm+9sfZ1Xey9b62imJmXrnoPkbLvhtp\nSuM+SodogBHTET9KufR6qp/sPs4H62w'\
    b'9xgQ/ePrN0nYoBT+gF/lr6zYfP4Dgqow6RhE/WXgctzdL\nFb+eLMm4/Aj7Pnbx51'\
    b'meKRI/i/yAE6p/PL/91eXzHzRiP7z0N7Bvc3U//J7peFvBYD+QJe25mDhJ\nP3y9q'\
    b'kPGxCc/fzL6W7cVDb8qTCmVG+Msv518DI8XkDQ/etc8NLoUNz9+f0lSUtlCPyVwY'\
    b'DIT+sK+\nzTpGPw//P79pwm5ddHYrv0jZLDfrpyQ/dbcIOPBWIb8XH58wT9a2Pi5B'\
    b'KAJExt8+nEAvX+FE+749\n8/nStSj5PklGHtf/tg6/6LdmxGT9ID8vgXQlak8Tv3K'\
    b'jMwr42Bg/IuF51WH97b4wT0mJleT4Ph+q\nxZTXtwq/NTni7WmFE79RZbW3HTLxvp'\
    b'XSl1XUfgK/lsWMCWEe8b46EGn1UBDYvtJ4KLK3FNs+DiyI\nnjDe2b5W99QQ6N7kP'\
    b'g+56Aio+r4+sz4QHjgi+z7lkyHCaxK2PjSGmlLSHrK+tH0dc3xZr75EhyuL\nJny9'\
    b'v4VlZq+FI9C/UW7pRZAhnr8CX8bTSTeHv5PNWV4h6Y2/t4v3rPEScr/ol/AYi/i/'\
    b'P+mRlyEP\nIdE/43aMmKsjoT9wnfnDQXJgvwenJ2LVaYK/3wsNroRUNj/5YFnet5U'\
    b'3v7lPf+UkNkw/nozIg85F\nZ7+2vxKXgAU7P6pn3VU00FA/5YwN2DovXz8xbMyFIV'\
    b'84PyYeXa8QJC8/g7en+PPePr/cwhQaOdP3\nvqRh0X5qYfE+CIRmBdmZ+z70lO5/L'\
    b'FUSv8R4u5bsA9m+eRcA0jhk4b7Fu6fJ/rMMPy1fbFTWWLy+\nTualuIvRKD9KFzjJ'\
    b'jvUQP+aWEAzTJiA/V84NtT1cPj9PBoZxBOtQP8UpW+tj8ka/lskwfSjkLz9+\np19'\
    b'tJvCAP/jqzxRTjYc/XuzIQCcrST9SkthirFpRv0YGY1sw3yQ/vcf3UpK9Kz8HPDh'\
    b'WvMPtviT8\nqVihhyc//OrS/N/YJr9JtB5jf+ZCP0rjbbXqsju/g2U/7xMeSr/NEB'\
    b'NTYTjwPq1r47hQUhm/blPk\nL1Lt+T5t/14CDr2MPmJnnmzIlgY/ca4Yi+6n475TI'\
    b'hziNX0Xv/xF5BnWtfS+i8qWSrKhCb8c8Gw6\noRoMP12KR/2xvOS+fA8iWf5b4D73'\
    b'7KAWtHz9PnDBoUZQ9vs+f0/64KCP/T6jKpfaVfPGPtCLj0sA\n5gE/COjYAjgh4b6'\
    b'frJegljTcvlYyPG9dq9O+6t4BEUwrw77xahzmtbAIv7YARroWTO8+wneURYOC\n6D'\
    b'6PQ2dEqwbtPuPt42M9ydS+HqulNF/j3D4YF5tRsSCWv8OF6tF5i5q/0vjgj0n+qr'\
    b'9W3u0699iE\nP5wt1Mh5RYo/c4Q/3EaZZD9NGl7HSAGRP+N2jJirI6E/8ajYzRirr'\
    b'D/ZcJpuJ6mSv085F/NfTZ+/\ntZiw1JGoST9xVxnRgXlkPx4UVJS6B2O/oEZ328Qo'\
    b'gL/daWpaTYxXP9ef7qTHA1w/0qlINF/6Ij8q\nDqD/2uoyP0YTgVm65xa/UZNdjTV'\
    b'B1D7s4o3ymTTtPu53VsILzQk/FGXxUWH5Gj9XrrfyAbsVv125\nxYUq+QQ/SVMbkj'\
    b'KfID9G9VwH9GsJP7IcqoLfnve+lVEk5l5aLD9w2FBBpR0lPzQf48q+vDU//pp0\nN'\
    b'ZoAQj8/HVykhPRZv4JfX2DlO04/ppaElgDlLr8AxEuSID9ov7btf24Gb3k/RxRfm'\
    b'4+roz/x1hEE\nGftKP6VrhsiazTY/POcn8WSYTL+F1crz0XBhv2AhqFCG8Ui/BL4K'\
    b'Q+1aUD9n5Qwkc7Ncv5yOi1K4\nMyk/D88X/P9dVr9K12r93LUvP5CrmrXG2RA/mVC'\
    b'asnemP7+qUc9ro0nMvpALZ1DDWhS/G2uY363E\n+T5Z+Lm92184v2C00KDFOUU/hx'\
    b'8F/CckTj94kZTUnmokP2T/OQd+GCO/rnyMujYULL9FgoGemwE4\nP9sMg8F5iUG/o'\
    b'14/az93N78Q9lt1tbv1vuQutD/zrQa/Q5ASoFt2Hb8m8e8phvj0PqtkBLpZX/Y+\n'\
    b'IroPzYOhuL53+EepONG2vm3gWiKahhU/IdN19PQv8T5+40kgkXPEPvrr/Sylmfu+'\
    b'exWlKCnUw74s\nuv/+tPK0v5zcy6aikFq/yERZHzkdbb8Iyxtb3I5Bv6+vroS5G3a'\
    b'/9hy160y/j7+bSOFA6Fsxv3Wd\n+cNBcmC/13Cabiepkr+Gk5uGSkfNP9eFvEOfA5'\
    b'W/V5CtfgCdhz/5IYSTbIqzv7UdkSEiRYk/A0t0\nmHVUcD9y0jVR/p2wv3bf1PtiP'\
    b'nk/iVOEueTegr9jv+Ddm4M7P9CcDpRMiGg/d9SC9AqIeL8UAC3Z\nSqjcvlXiyuv1'\
    b'xSW/NOVma/WuK7+KAuRbqLpMP5JG1LGJsEG/4Q/UC6DlWz8va4ygpkxNPxTtNf9f'\
    b'\nBlk/msOhjs6dXb/5Uuxg9kIqPyTK1UppXzk/6Z6tI0vMVj/LyxMdEWlkv7FTR7s'\
    b'kr1O/xHzfwwnd\nPz+ue6EyhrFzP3uQE+XKJVC/MMtEx8C8eb8VYhuWlN1Tv1O53t'\
    b'Gm2uc+lqXve8FjYD9bhjc1VIti\nP7ioFsedhUc/Ivo7kdspQr8IV2yBdlRgv8xja'\
    b'e2A3Vm/hCE//kCoQj8WhGazgCQqv4tnAea3Q0S/\n1Pz3TguCUT8E54oEXtcOP3zj'\
    b'3OR76RC/DfDej0QN+D6MtRr/D+1kP0wUCkfK80Y/wNaL9bN6Yb/a\nDsUugaVRv31'\
    b'uyPTLsgE/7XHS4mLlQD8rfJmBTcpUv+H5YhFjLTG/o0zvss7NST8Aiy9/od8ZP02'\
    b'r\nMqH0cjW/0dppi7AhPD+4aD8lsQ0cP2L/hlfrmSy/EmDcIf+F7r4tasIt6zk/P5'\
    b'nSaIOs3kO/CCqt\nyWttQ79kc4OGtFM0v45xNywtePA+C76QTqFRAz8aMa4jCWdrv'\
    b'ywDXSBJKbS/L7z0RJRNkb+kLrjs\nYN/2vv4Rv4G26XE/OhhWEF+Zmz8WRjDit/I3'\
    b'PwGnJ2LVaYK/UDkX819Nn7/XhbxDnwOVvzEG0Jew\nQ+0/TKcfRtItwz/dgvn0kPG'\
    b'FP9oLH7AZR+S/e25o93cgy7/N3G9+CE2HPwK549mQ2cm/ol22nD+P\ntD9R2nirp0'\
    b'FSP62JzRZTuII/MkflHKeFnT/bdBrdxwg4v2FmOGhGYFe/qHxAGWg1e7960yZByd'\
    b'hl\nPxHsjypr+io/2jSn7xFugr+90UY7FmBwv4/PADXUmUI//uti1Zyogb9Im0y3o'\
    b'G4bP925+LU2BUI/\nMO9qRNR+Kj9zvradPIRlvzei4A+zt0a/IhEaShN+Tz9D6vpT'\
    b'u3GeP7jwOCSuGH6/TljOhnVpkb+h\nii7MOeh2vyYP9IvDnGc/kb619OmNQb8Sudp'\
    b'bZP5MP8L3n31aUWO/lnTIweuvVD8G5bMjM2lnv+jo\njwPk2lc/UuujQ74JC78fuq'\
    b'x5AONcv7GEx2GriWo/YExZfy6vZ79vds5xzLj+Pp2Y/BTXIve+nbSW\nadLAM7/9s'\
    b'H9WFIdLvzTYuLoAq1e/ARcPGl8Ufj8HWFCAOsMvv8kuLuNnxkg/jfu80XS6V7/ZX'\
    b'YCb\ngCtRP/U3QgaNnk2/r12DMe1aZ7/cGwXpafYjv5htapaX4iw/5A/UXSoPUL+o'\
    b'H1M1ULMgv5fgSlVh\nIDk/co/j8TitFr/DMAtSLOs5v5t/vRH8GyQ//kSxkyUtUT9'\
    b'AeOr0K7EwP7HZH49NYhq/0Ez+TVpx\n1z6oQZNsqMiOv5jpTjR355W/gu/Tb4g/xr'\
    b'+Ry8Bu/QRYPzRAyEJK2nS/1Wuf0sGrj7+AysVLhwY4\nPw8MDa6EVDY/w5iw1JGoS'\
    b'T9WkK1+AJ2HP0ynH0bSLcM/DW5pLnRg5D/sYGG+BWVWP/YUl13aR8y/\nAwn5lhB1'\
    b'yb8VIVy1TkSIvyA0id4Fwa8/vq7RhzZMx7943XbY+vtzv0yrWr/lxaI/KRsHuyKN'\
    b'qr8J\nURtkDglLvzq/E84f02m/6UlWQag7c7++3E5xPOY3v1N3SH7lpmQ/y8rmw5L'\
    b'yYD+bX/EOmRBJP5jF\nTyql2WM/93AL6szoYz/Hje5vqm0av8SeicIE8w2/Angttp'\
    b'euM78ABsN1bcw7P5+BL6Swzio/YDH7\n5KDsFj9MYFG9fXM7v9CYEyK4rFi/L2ewv'\
    b'A7/bL913X9RDQ8yP1U17aiNwCy/jR6ZqD8qHD/08MuN\nDXA8P7TJeHkGbZe+ehPi'\
    b'GlYOKr+FuWAvn+hCP0gzfOhQ95i+tR1VYKAn8z46NLKZ6ucLvwDgRdvL\ntiA/GV2'\
    b'X8/DbDz9xUB3YqGQIP+/3Q/GW68e+M91o8OpuA7/wdkLDUcsDP9mGqd2D4UK/nYH'\
    b'o7VZ2\nML80MPP7t64LvxPrdHK6qSE/ZweTEeDzDD9xnBW4l5UVv4YjXmaqrDk/ZL'\
    b'WvDfmTIj8ibMZWq7Dt\nvt7kIIz4uQU/sF3ttwqN+j5u+X/J/R2zPq6+km4i9+C+j'\
    b'1KcsY4Y6T6BXA5GeoL/vgZ+pmzIBPg+\nCtuCE7LO2D5qukr+2vb5PhnK9T4Ymds+'\
    b'PFZCSyvY6L63ZeTph02UPw2FEfaDZ2c/xJagiZcEPL+4\nuTLgrjMLP8W3sL1DwQM'\
    b'/WLC+GocrYD9EIkaLGQA0P8RhWd63lTe/a1cZ0YF5ZD/3IYSTbIqzv9qC\n+fSQ8Y'\
    b'U/92BhvgVlVj/WamdwtfOfP9aCOyeDz5G/QXOAS+7kY7+ZvWT4JGudP0PhlnwhxX'\
    b'c/zBTd\nVdAAUb/ChWvbOiV4v/5LwSD+qy0/F99FF2A9Xj9xJ3wlvF9HPwR7AsYNh'\
    b'yY/GNEQIC/lJT8PxPgX\nz+PVvhU4umJm4C2/UFZP4egsEb8BrgX9sdA1P80isHZM'\
    b'fTm/Xcgk6+wGEr9B8Gkrp9cYv1lE0178\npDG/x6FU888uSb9uXfx4I3JKP+M+KSB'\
    b'8j0I/NI4u/u/vIj/thxn8qeR3v38qybOubDG/DzzScDNm\nPj/y9VBHF01cP/ePEJ'\
    b'XONzW/4GDheusrUL9yHrxHsBdQvxVwuJm0diE/HIvMAK7nEz9Onw6AvINd\nP9MAn'\
    b'RxcejY/k69m14ZyLr/1p3WlhEI2Pz7nZSRg5xG/3YYn17U6JL9P4OAgKLQYv4o0Y'\
    b'RcZ6hI/\nPBycX0XxBD+lmykQzQpXv75RxdwCHzC/00gcUR6hKz/ptdBcYXZHPxL+'\
    b'ULO5Dya/lqOegMxQHL8m\nigeK009CP2rSnWRfrDM/sGf5fEN6Db/YZwrWaC36vmH'\
    b'aSXc+yjA/6yOCiaDnIr+SGXu9ZxwTv0YH\nZp/OrRU//ZpkRqrq8z7PGZx6ruowvz'\
    b'9LXSRsazM/K8clcNWUMD/h6whKpwEkPzcQoY6ngug+QzlR\nFZyf2b7iCJRBWJdwP'\
    b'3+MfGFvdpS/Uay8tObtkr9ShhkXgT9Hv+MNxqzQrU8/6wyv/2RJQb8BuoYQ\nVXFA'\
    b'P2JPf+UkNkw/FRRUlLoHY7+1HZEhIkWJP9oLH7AZR+S/+BSXXdpHzL/Tgjsng8+R'\
    b'vwpJtAwr\nMOg/TmeMFZ+ZzT+o2euzHmx1P8EmLmOiwbq/nu+mwtbTlz8372IC7c9'\
    b'hv5Y03lllelq/Ro0DRI1r\nk78aiun+UYkeP1VhHnWlSkw/4EOw+ZJ8Uj8Fxkf7lT'\
    b'BLv1BZAow7ZSA/r12lNi6Kbj8gQPNnHx1Z\nP5NmDQmC1jg/YyfdxkSKaT9aA+4N2'\
    b'LQXv7+P7Ee+ijO/4kJGd/X9EL+vvQDVFbgqPwActXJBBhY/\nxzfO8OlULj9NZGam'\
    b'l9N7v9VFWDL8s0G/s0yJGsQ88r6wXNSAnDJjP6DGC7h7fUS/SC0ONCg7J79G\nrix'\
    b'59XFHv2SVu0Yzz1A/nDH1a0DuNb+LyDYBFLlUP737nwa7skC/ZpUJz4qKHj/38D1'\
    b'j5nNKPxcI\ntGYU/lS/UR7U2/vCUT/y29cH6GQXv7kDRKxWAce+rBwLEbFOIj+35X'\
    b'N/HwAev4s4N5kT90E/pjiP\ndl5BZr8jGs8EaMkwP0NWnyBLwjO/mPxHdC12QT8q7'\
    b'U1sGm8mv4P4aYZsvzY/jlt+zaLjUD/F7++M\n4Hn5Ppqo2/kdNNE+lkWyW3LmMz/h'\
    b'fs6Rgg8CP55fjloCSh6/E6+OWQV2wz6MhRI+IRUVP53EncGc\ncP4+qcUUbUF3Mr+'\
    b'ET0KEyj8JvyuzCo3o/A4//RtAqFg/AD/jfhjeEVxEP7cVWOX4XaW/BXuFl7QU\ni7'\
    b'/S+ykKuTg3P5YzVJ/SAUC/Z0BDGgJDbD856q4g55tTP8CMyIPORWe/oEZ328QogL'\
    b'8BS3SYdVRw\nP3luaPd3IMu/Agn5lhB1yb9Cc4BL7uRjv09njBWfmc0/uEKARYW3x'\
    b'z8jI53aIoxfvzUpwHYzo5g/\nmDuGSbWymz/dHZkgEC9WP82q77M5N2A/cq6ngYK1'\
    b'Kj8y4nPg/+X9PiNOwrHDD2E/44+UQLKUYT9R\ngpzqojdav/auZ+qUpAm/shyN3RY'\
    b'bUj/ga0Lei0RXP9nF/m2FzEi/RNputUG6TD8zba2dp64Bv4hG\nE8BcYRm/OLrQTF'\
    b'6HAj8kSl51I446P5JqABUj6yi/f8NQ7x7D5r48Df76Ma5UvwNTP9OFT1A/0raL\nk'\
    b'vKOXL/K+8D8U3pJv2tuHh+dTiY/uAKN76AxYD+94UYUuthRPwYFpNihFTA/dG5pQ'\
    b'INRQb9VrTfJ\n2vAgv4Pr0On/LUO/BvyUqohqLz88sI4r9JHzvjttieA4CwK/+ZuP'\
    b'ZGnJNj/BtVMIicLsPjr9unzZ\nxu2+dq5lVkX31z4D4jy4lZVFP0DMP31/uUe/mjV'\
    b'+Go82Ub9ykWfTuyAyv5fqrneqmi4/r5Dktj6Q\nLj/a9sTiYOw9vx158DkEizw/r+'\
    b'Qi5kabPT80YHXPZjUdP9nJF4OvFu8+a/oPZtyeKD9fJ25sAhK+\nPiVvPvfnfBC/U'\
    b'iMQ7Gyp6r7yvqjnXaDrPng4KUqFGDa/vYL579pCKb9o1Pl7y/4Zv/kZBJ7dbAk/\n'\
    b'xR800BXT6T4MOcSWvg+CP867FPKmUGA/428iojzTYL8RXjwCzQk/P5irNwAOAT8/'\
    b'E+tnYpHOOr9h\neR0nJckiPx/AEpeABTs/8mlqWk2MVz9w0jVR/p2wv87cb34ITYc'\
    b'/GSFctU5EiL+evWT4JGudP6jZ\n67MebHU/RSOd2iKMX79eLVyE1Ey1PzFTXBkdlY'\
    b'm/+Mp3TXd6pD9wYa6Z+W+xvxKP1/15NXK/5yPq\n5Dzjir+2KUZlXZghv1y/XbffQ'\
    b'EC/ASWqGhpkab8HWBuYEXGEP1W7Z7gLikM/HOXu74lAmT/nzhUU\nmw5lP33H/nlo'\
    b'f2a/HoOjPHgxob+cbvIWCX3xPhxZVk2HJhs/Ida0OLLGNr9QOA98QmfSPrli3uUr'\
    b'\nJ0I/1RJW/BzFET939RW+/X1tP89WJ5YhcVi/PLkhT1iFQD9Qnt5Q6tM0PxvaB2p'\
    b'vSEQ/gJyJN/8F\nWL+ePikV2LlMvwat5MEkXlC/fbk0Z1TPRz/7uf/mnO0wP02iAi'\
    b'wOQUw/OIpDYMStRL8RCNS+p18v\nvxp0NJSTjVI/iIhuZVwnUr8ajVb2u4PxvlpsX'\
    b'ew3Mvo+d1Tt0FRqF7/21Ir5H4xav7eUND+pm1C/\no09UUgNJZD+olZppb4RAPw9K'\
    b'5BGndDI/6yTOd5EUQb9ALaY+UzpMP+vi/iewsBY/Lpa/1pMRTb/S\nOrkYURohv4u'\
    b'd+70ZgzA/nfIfMxifPL+Cdg3IAX4av3RfEPdZ2ys/MW6K0N535r4eFmlN/VA7v7/'\
    b'a\n9BLuQjY/COtacauARD8kXfJy8W0wP1Tmdr9g7fW+Nd1XoaT5AL+KaJ/+2BdSv0'\
    b'CyGnPCXGk/ggim\n1K7Xcj/rCH3/wKMcv5t/RvkyIUQ/E4eZvzUrXL/p7s2o/kc5P'\
    b'75n3VU00FA/2J/upMcDXD+A39T7\nYj55PwG549mQ2cm/IjSJ3gXBrz9F4ZZ8IcV3'\
    b'P8cmLmOiwbq/MynAdjOjmD8nU1wZHZWJv0qmpAk4\n2tg/zBYQ37OGsr/FxLRfidw'\
    b'/v5KC5qmqPLi/YDhvWWlvkb+9iXFsFnNwPyHYi1NggoU/deHnrpCB\noT9YTAYh3B'\
    b'F4v0eoLIjGvFy/8zl50hYSkr9XzsSUCBhVP/05EMF+QFq/gLuLsVMXkb94Ft08sd'\
    b'Ii\nv3m5rvwpuTm/7SzcpYxB+b4ICh4Tyf9OPxDTDTpKjj2/SDVsRXFzGb/drlY0/'\
    b'x98v0A/b1cu4mA/\nEhiUtpyrFb9ng5IuOrNlP6fg6wePAUu/rLGNMIVPSr8HuGn1'\
    b'LMxRvxhYySvDHWA/Gbl96JRNTL87\ngxStfTZWP+3v70LvxD+/lvd/GmwwQD+9/z2'\
    b'BnjlSPwoe3NRg6GK/0/HJGhUUYD/8tWosgIkRv0vs\ntTuYgvu+7/mo9FitKD9egV'\
    b'wuNxA9P9oxsoYOjlE/3zunwTbcc7+lZdDii/4nP+gKXQGWo0C/KW2a\nXvOhTz8Mk'\
    b'FbGSspDv3XiLk9Nvz0/Sqfle1ZEYD+NvXbXmbX9PhIuhH7Up/C+vU/nx1+OQD870'\
    b'wlt\nPkwjPw4PBNgtITS/lyLIvjqP+T7Gics55MEcPx+zzm0Idg0/xejK+AtoR7/K'\
    b'iiNa/3giv715zP+j\nngs/V/Ij3VewBj9LB/zmJuxov+tIwWWonqA/1caxxmMHo7/'\
    b'nUy9FNX1NP/NjjsQsm1C/1ghCuxuI\nfr9l1WVLkak0Pz2NDdg6L18/JalINF/6Ij'\
    b'+MU4S55N6Cv6Rdtpw/j7Q/vq7RhzZMx78aFd1V0ABR\nv5XvpsLW05c/ljuGSbWym'\
    b'z/5yndNd3qkP8oWEN+zhrK/lhkv8hLo3T8QRiWmfCuPv1y9/NiW86a/\nQm0tzdf/'\
    b'xr82huAHdStsv4TxKbpHV3q/ftd+Ag8goL+l083igb5Fv2xyXiDOfX2/Y5ePblGj'\
    b'n79a\nbgcGQkp7vxmHJvlSanq/ms5tWqwZmb/1ESnaDbfcPt5i+EHkxJ4+YeN+l9B'\
    b'//b5d2LpEp2n7vtRN\nF82/VxK/vs5gUdWT4z7Cr89eBXI5P5+wMzrpWTY/0HgJFz'\
    b'KOXL+qUy/fhHA/v7+cyp3m1hC/k+iH\nNL6eET9tRyi1o3IyP1xjPOQ4hRy/W2klv'\
    b'cM3Az8UX5xv69YEP+5ochSFNPa+IQ1cibBlHD8xxF69\nS7IlvzlM00DZYxE/y1Co'\
    b'/irTML+pd5lukUm6PowAKinlewE/xiqTknnPc76XaNwTdWG3vrNsJSuC\n0Uo/rpc'\
    b'z+tP8Tj851ZoVIwUTPw3sQ7EPLy6/RYYBLTWaKb8M8HNS2xgxP5YBOjlKTES/CIz'\
    b'EFuWR\nQb/tPWfNDELIPhGNTOWxGiC/c1sYxKPTGr+B9oXu2ETyPgSeqgKSo9Y+Ll'\
    b'3OFgdy6z7NsUv/NVwP\nPzkNmcJOUus+u5CtgGYO9T4RPm/HeuX2vuzEWD/e/Ai/q'\
    b'WSAeZkD6L4lCIdi4g15PyX2wEEpkEq/\n143BMm9SLz9xYHTIzDsmP/SaYAjWnio/'\
    b'Hu2hTLfMRj85OE9+snLPvqdszIUhXzg/Dw6g/9rqMj8j\nwODdm4M7Py7aeKunQVI'\
    b'/dt122Pr7c7+zhWvbOiV4vzPvYgLtz2G/yB2ZIBAvVj9wYa6Z+W+xv5XE\ntF+J3D'\
    b'+/DkYlpnwrj78Yhyw1mkrjPxhBTVqA83C/ZE7fmd89lL9AS+gfDR2sv4EVqJ14VX'\
    b'4/A7wo\n0qOCXr853SI+IoXMv54YIK8Mqr0/E8zMHPQ4sD8PidDQZofQv4Or6fv5C'\
    b'r6/Wfn7fy8Wl78wiiLx\nGrL6PrGYz4+U1Pk++YPcJF63Az8FGDcBVCEYv4CXAWpp'\
    b'VR2/R449rM4EDD8wmExGcdA1PzxFyTee\nqxy/hRSq75XxDT853qxsrU0cPzcvZKM'\
    b'QbRi/qvFrIl5pFL93Uf3P1Twqv63M7hkTPyM/l/gUi38+\nED8Eyg7ZTfcRv7ZL+U'\
    b'g8uPk+MHqQamNW875pXycLIPYCvy1fxTnjFSa/TG4yZ7JlAD+4XM01mOu/\nvmNru'\
    b'5BStOW+0HVyIr/f9D7PUkl10YM9P21NpepGYk4/agNx1/quID+Ha+dWkDkmv6unJ'\
    b'Kho0yi/\n4GH02cwm6r7UuVuFJvgQv8BS+/JvdkS/Z9mxKIktLb+WdD7B57gAvx4z'\
    b'QSTNlxC/+aiVRW9eBz/W\ndnaPmJAfPyLyDIKziBS/Ix7YpCdW7b458gmfkwsBP2v'\
    b'wbsHymPQ+8qL8YF/bL7/JIUbKnIkcv/XC\nXgDNHAO/ECvUH9SD/z5NMy9pzqA1P2'\
    b'srWnVF+W6/wVGpbDYcZD/QlwVeXGEDPyfQy5nhWjY/9zNa\ndhFHUT+9XouJrGAiP'\
    b'+EfXa8QJC8/rhOBWbrnFr/GnA6UTIhoP7yJzRZTuII/S6tav+XFoj/VScEg\n/qst'\
    b'P9I03lllelq/36rvszk3YD8Hj9f9eTVyv5GC5qmqPLi/Xb382Jbzpr8aQU1agPNw'\
    b'v9FmJ9UA\nJ+Q/ZsPNjpMgpb8jWFrHg714Pw6X2S9BddK/CpgQ6NIwuT9kxYgEWn6'\
    b'9P2mJFe/+yMC/LKzveEub\npb9JjNQeWd+9v7+w0fgsmb6/aE2HSKGuiL+KXDZQlM'\
    b'LRPjsCAlCSahE/f+YGJLlkA7+P1rM+y84D\nvwfn1+FOviU/jt2FjKsdDT+YFitJt'\
    b'8lcPxKKnbTsx+A+TQCBEIJdR79KYgrnmpNZvyk5t0V4ZTI/\n7uKV9cj0UD9JaA1I'\
    b'kXRcP4+zFcY5Mk6/Qss7K3pb+j61wOLQJcpHv1IL6IjdbSW/waLtKe9rET9j\nbbs'\
    b'WEdJHv+gnqGmWgVI/gv7J9LG6T7/j6PslBBMEPzB7cbrRPvA+Vjg/4TKHGb9vQfY'\
    b'ogRVAPzjK\nZ7JY/Di/pkIz1v1UYz+AeX+DQQQ6v9E667zttTI/ZFwclvuKPb+gRb'\
    b'eKW2kFv00y2Pq04DO/0H0I\nEhonU796lseBWYwgP8EPHd9R3SK/s6W0SA7yDr8ue'\
    b'Gd/h+Iiv/QThN0Weig/h8kEGl3b5b5JOX5M\nAdkuP3XZsMZCtji/kblKxyhEJT/Q'\
    b'zkWRnIgIvwNeCOGFXPW+qDng9RAnEL/H22+BgZxbv/HQ8CAY\ncGI/Q5ZBaIcpbD+'\
    b'wJ5TzkLElP2JZBfUchTq/FKm+6H+KAT/MVBQ6fUwuv5C3p/jz3j6/m3VdjTVB\n1D'\
    b'5v1IL0Coh4vzdH5RynhZ0/LRsHuyKNqr8h30UXYD1eP0mNA0SNa5O/IbCngYK1Kj'\
    b'/sI+rkPOOK\nv184b1lpb5G/QG0tzdf/xr9oTt+Z3z2Uv2TDzY6TIKW/pO41wMl33'\
    b'D/w/juj8pZav82ufsH+Ybw/\nENCcB3Ixt7+lOvG3leGxP27wGhTlIai/4aJSLQob'\
    b's79eftQT7t2ev1dyM37r6JC/J8HqSC+Nqr9a\nGRFHbEbNvoXr6j4w36E+guO46oX'\
    b'I0r4eDZllCnHsvlSRvR1oUvs+6VcpuVp0wD5oGV2BpMH+PhYj\nfFFNpfe++i+pwK'\
    b'NYE7+AGccH+Rj8vowmVtCQT/8+fCfrpz6eB78Hud7pgK8IP3r7N+9anBy/IJkc\nt'\
    b'QrgBj/JhGUNx7D3Pi91RwQTHf0+gQxFmjHN8r7MP0BCRyDyvrEx5nLgcRQ/Qgcjn'\
    b'K0JJL9+CXOz\n6HHpvmwFrbKvIOo+F6uvm3Tnsr5RastVBmQtv4+II8J4GxE/Nmqy'\
    b'rN05Oz+LgkIEWGIXP5BhSZ5L\nXgG/lOkEIZXBFL/E7E6dwrQlP/8meA413/W+pXt'\
    b'y9ZCyLr9tlFOUaU7BvtrVftEL6A6//eu24kYO\nBb+JJMlESzoGv1IZN2wV2RU/3u'\
    b'4K87j72D7gShDhsR4APwUjkg7XF/k+6xzXgSxkGz9IK+gIS2YI\nP1px3htNIOI+7'\
    b'/Y3dG4K+76rjWW8fkYjP2A6zyENRSU/Wl2UxSO2Ab8mc3QJhz+6PgMynv8WZ7Q+\n'\
    b'H/nx8qFD5r7jbZfTpqfdPhvDFBo50/e+NeKN8pk07T4v/SzZSqjcvul0Gt3HCDi/'\
    b'AlEbZA4JS7+R\nJ3wlvF9HP9WK6f5RiR4/h+Fz4P/l/T7zKUZlXZghv8GJcWwWc3A'\
    b'/NYbgB3UrbL8/S+gfDR2svyRY\nWseDvXg/6/47o/KWWr8PZjXYXomoP1UlemMuNX'\
    b'+/ZGOOm0BLdT8GLAUwfllfP7UHM0/XtZY//USy\nYpU9iL/FJFaUTW5xP4/AR7wmN'\
    b'Zm/K1N3d1jdiT9YFaG7UCbXvj/QiwPHh90+wRzl1AK88L7XIbRD\ntYXxPhqlMrvW'\
    b'EQc/xgYNr4Ljtr5639ShxNEdP67Ir4aF/bw+0geAi3Kyxr6QHvXBimwiv1urud4K'\
    b'\niQg/AGhQJ/bDBz+PMNoOASUjP5wNniJ5eyC/pYTf2jIa+T5ur01LD0brvpEP3dT'\
    b'FG9G+bcJqx5e8\n/r55pgM5uR4Tv4GhrOOUFCI/nG1t3RIsIL8gEUZ9/+rLPqSmdz'\
    b'buprk+0yBcf4Ll7b7/DbY5w4EC\nv52HP9EeXhC/XfEEPK6aND/Qe39YHzABv23Ps'\
    b'0xePBI/rhitGjz1Bb/zAy+Yd+gRP8JZlK95wua+\ntJTIo2NmJ7+bn9Eesh7PPoT7'\
    b'b41GIQi/lZxrW12j9760OrKbkHfrvobel7QzcQM/YuT7/g11374U\nIa2Xjyr7PhI'\
    b'k+3wAh+W+GjZFxfW0CT8ZYQvBYpPlPqvBhIbZ19k+amZytVjV475p8oFxK98Sv+H'\
    b'c\nM2lWyDQ/8Cs+t3pQJL+2jahCWJjxvkOWeEZjBeA+lcoba+brFb8pa1KfPT75vq'\
    b'Ji0X5qYfE+63dW\nwgvNCT9w4srr9cUlv1JmOGhGYFe/Ob8Tzh/Tab/regLGDYcmP'\
    b'1dhHnWlSkw/IU7CscMPYT9mv123\n30BAvx/Yi1NggoU/ifEpukdXer+AFaideFV+'\
    b'Pw2X2S9BddK/za5+wf5hvD9VJXpjLjV/v1n5uH3u\nM9M/NVP75RdJvb9rA5z8/vx'\
    b'iv4QNGEA+c4e/2BMvu4FbeT/snZz8NRBpPxstKxjZoIW/woEoGoKk\ndD/jijAkSv'\
    b'65vkJBW4wKpv2+7ZXx8Pomyb5IGU9+2qf9Pln46Qn49A+/JIVpyr0n0r4ecu59f9'\
    b'45\nv4cmgtEDaCM/k6W8JH8wCz9pObC03JwgP8nZMTGODxa/MQHJP7sE9T7F3HO+i'\
    b'tAQv/eRMubrnic/\ni8/PyJ31F79caK22RRgWP9VAZ8QAdQy/xZNOJuZ8Cz8NIeF8'\
    b'7lYQP0qfoTqkgiu/g6oVFeyWKD9n\ntIaXnkPPPtwQ4rYRouS+adyIRJtE7z7h2V6'\
    b'pbcgwP1tVALp8cis/oJKqIGlOPL+zP2MKQsEYv7x2\nyoK2Pgy/ZgkR6bfmFT8qzO'\
    b'9c3EYhvxtMMp5RVQG/gOCeHP9rJT+GgFBXvOfmvhIUtB3BA/W+I5LQ\nTxziGT+T0'\
    b'kwCatYNP3BY6FHRXA+/IVapHakA6z49ZlvdNq/yPqc4hSiwGvs+VYjB3RfAJL+Z4'\
    b'1LQ\nhWQKv/5tOK7ESL0+KiJtl/ap7T64KjoxtNElvzO4QV7zZlc/NGJLM9nqSb8l'\
    b'IhW9IV7/PjAR681Z\nQR+/CPrC+FqdP78iyNnVd7L1vn2EZgXZmfs+DWXxUWH5Gj9'\
    b'i5WZr9a4rv6p8QBloNXu/5klWQag7\nc7/U0BAgL+UlP91DsPmSfFI/4o+UQLKUYT'\
    b'8EJaoaGmRpv3Xh566QgaE/e9d+Ag8goL/4uyjSo4Je\nvweYEOjSMLk/ENCcB3Ixt'\
    b'79hY46bQEt1PzVT++UXSb2/3ICrfxt3vT/1W7CR0bdVv5BFaU8D8IO/\nqjjpY2ir'\
    b'fD/BvlabF/dWP2HtqbTzVnG/PUgK1NpXbj+JZ5pyvACsPpTcktFJsm+++sQQOSOX'\
    b'1D4Q\n12t5xB3sPsxuToAQ+f2+eLbXg5VY8D67lt53cIQZv2UDVWegr/Y+SAfrzZM'\
    b'6Bj/Y3X8/qF4BP5FB\nAWKa3nk+hWwI6/7sDT/Lgs5cTPXqPquR/aCalBc/Q5rZU0'\
    b'0pEb8wVXbEjGPhPsgwrMTqcAa/ef/g\nXKpyzL6l/yVXLqzwPhQBqVesov2+l5/gL'\
    b'Cd/ID+O1Lh8FnPtPuL1iBBA5Ou+ZelI/Xh3uL6N26Hx\nNwEiPxcfNXb7dzK/NCNX'\
    b'YkO6Or9mTC761foXv22sQIxnkB0/D342xqfbFz8p8Qv3pqojv9clxtdP\nvCg/C8u'\
    b'bqPfcKz+bHjxctbKHvoIqomdXDdA+7nXaCGLpBz/4Wj1eMAvgPrnpXlc3jPy+wZE'\
    b'3UkTY\nyb6vBRDd0A/2vvfFB1vi9u6+mczvnSkvAb/Pq2AkVwrdvsWZ3Ay4ifc+lr'\
    b'8YoaaI5j5ATsCd0qM/\nvxh9UYN55Si/BdRT8MRf5L6Ey2WRdxfQPq/zdYb4vvc+E'\
    b'SxUmBqQzL4NimJmXrnoPvuU7n8sVRK/\nbK638gG7Fb+OAuRbqLpMP3nTJkHJ2GU/'\
    b'wdxOcTzmN7/PxvgXz+PVvhrGR/uVMEu/UoKc6qI3Wr8H\nWBuYEXGEP1hMBiHcEXi'\
    b'/pdPN4oG+Rb873SI+IoXMv2PFiARafr0/pjrxt5XhsT8ULAUwfllfP2UD\nnPz+/G'\
    b'K//FuwkdG3Vb+4HW7KRUbNPzo+zerj97+/fw2p94X4sr89i2DW24CSv54SyN8opp'\
    b'A/gUlg\nt8PSgD9g57gu7JrJvoq8lzeRBsS+MtLwUxcE0T5QwO/BttGCPsX9qF9h3'\
    b'ew+KSe54Z688L5gFbsi\nfeciv+LOLIlwrRW/ePnA3tAx8j5885kjzMwgP0mhNh6F'\
    b'Zfu+8qHxvtIkE7/Nzc+6jH0cv2vw8l2P\nyPs+gigrZ2G01T6Wi6VvJgYIP0Lt5A5'\
    b'VbgA/NAHR1bfbwT4WFmp0cOwQPxl6SsNULwm/51/FhhtT\n+z5t4L24AMbhvlpWS2'\
    b'nkJuM+4bQkWv640T7x6+/dMs4kvwC/qDiwTRm/7tKSr/wWFr9U0Vp+IacU\nP8+pQ'\
    b'I024ei+zDAdMVqF5D7yxlRbkLT/PqcvRQCsjBs/sscakB/dED/X7C8ytGKovur4R'\
    b'bC3aP8+\n8TIFh66v2b4z59FoEIf8vujdfzDgz/w+KMSE91ItxD4pD5niyJ7vvksE'\
    b'p47QUfc+gVCXQvvoCz9w\njF3l7CQAP9qEuZak9s4+q1oydUiS4768q46S2QU7P53'\
    b's74CfguY+Z52FQG5ZJ7/WlMK4cHXwPgfQ\nBZ6zZME+Xe99ORmT7D6fy74baUrjPg'\
    b'l6u5bsA9m+N7nFhSr5BD+dRtSxibBBv0rsjypr+io/U3dI\nfuWmZD9+OLpiZuAtv'\
    b'w1ZAow7ZSA/Qq5n6pSkCb9au2e4C4pDP1aoLIjGvFy/bnJeIM59fb+cGCCv\nDKq9'\
    b'P2mJFe/+yMC/bfAaFOUhqL+1BzNP17WWP3kNGEA+c4e/lkVpTwPwg788Ps3q4/e/'\
    b'vxKXe6rX\nJsE/BZwI0oRCrD8RHxkbUFqLv4xNNrVle4Q/Qo/lKTK7ej965vzu0Ne'\
    b'zvhkWoRoKSNm+NPTx82Q0\nw75OebHcKhz6PhgVv/Pr4uu+Y1t8ngzv6L7P5Y4ykg'\
    b'4xv+aFqixBlxS/R3ILQKAaOT/FZCCclksz\nPyEeV1KMIgm/B3ADt6/8JL8iG42Vf'\
    b'Qc0v3UoPEr0XiU/VETz+ttW5L7A3zAi2XAdP6WCmTkUQwo/\nL97suq6F/76CgsT6'\
    b'eQIeP1m3ncvZFyi/BZpeQNNWJz9csK1AFd7GPgfiYfp9XN6+DVKueoF67T6W\nB05'\
    b'QAmj2Pu95IM2vSdm+slWUWfrGP79SDgYe/ZLyPi1/raEK5/i+3nwpJGWiFj+/qVW'\
    b'S+EcUv9TS\nhQIuSBc/uPb0XfcuMD/ptlo+nRP3vi6f93+swAE/j2GUc5NuAj/Hwv'\
    b'N5/uYCP3zdq3AY9Ai/ALyf\nPdxp4z6yyUaMIAUMv9H+ls8OxQ4/yhUTBLttEr8d0'\
    b'6ssaoPUvhnYbn/JAYo+Hl4SmNyE5T5Feljt\ndQhRv9TuHbG9OD0/htg1ZZ+sIr8l'\
    b'pdLyI2LjvhIxE4MYpgw/1eNTXyuwNb9Sh2iAEdMRPwoYANI4\nZOG+RFMbkjKfID/'\
    b'kD9QLoOVbP9s0p+8RboK/0srmw5LyYD9fVk/h6CwRv7BdpTYuim4/rRyN3RYb\nUj'\
    b'8e5e7viUCZP/Q5edIWEpK/ZJePblGjn78UzMwc9DiwPyus73hLm6W/4aJSLQobs7'\
    b'/5RLJilT2I\nv9ETL7uBW3k/qzjpY2irfD9/Dan3hfiyvwScCNKEQqw/q8EJgquIt'\
    b'z/jUyvgt6hlv3rpxpDptGo/\nuNg+wa1PbT8PD1IzLVyuvruGLIsZP9m+vKoCiznM'\
    b'4b6qlyMm+Rj6vkvoch9DO+C+cOcbOiMN3r7o\nFvnDKhzAvkU9DNMc0uQ+3k6/C9t'\
    b'CFD9YMjpBw2sQP2NOQsAdVe8+1g3inCyZD78nYDJgW5ISv8GE\n5lG9W/e+PplGUP'\
    b'GH+j5tX/haq7LUPogpX7GBGgM/bi66OvGpmr6qF9H47wD0Pn+SfWLzgOS+WByk\nI'\
    b'12g9L7BrYJH3BzpvvXR7okc8cM+6wzH9NZUwz5a3CYN6Crmvm/edyrgjAs/HtNb2'\
    b'ARXBz9gDa3d\nDcwBP0cHchf4iA2/WlyFnQqF/L4pZ7vDKKjuPgsusi3JfQq/iwbV'\
    b'cCjEAb//f11JEI+nvhWFIl49\nhP4+PpnR+UTb3b60Wh76CkjOPjnA3/HMSsq+0KF'\
    b'nuhAgsL7A24jt77HZvr4NmbWn6ps+BzuhVwjH\n3r6veW0nYCvFvohFh68I/NC+Nv'\
    b'RDPkEE4D6fYoYDpKY0v1uVistM7iQ/hQvoDGMq/77s96en7j7m\nPnvjG6FwnO0+c'\
    b'PC1juVxEL8zt/R6qp/sPiW7p8n+sww/SvVcB/RrCT8ma4ygpkxNP73RRjsWYHC/\n'\
    b'ol/xDpkQST//rQX9sdA1PxtA82cfHVk/5mtC3otEVz/kzhUUmw5lP1POxJQIGFU/'\
    b'XW4HBkJKe78Q\nidDQZofQv0qM1B5Z372/X37UE+7dnr/RJFaUTW5xP9+dnPw1EGk'\
    b'/zb5Wmxf3Vj85i2DW24CSvxwf\nGRtQWou/6lMr4LeoZb+1NpdyJDDRP27/rS6gYM'\
    b'A/RnDpnHRkoj879ImEYLDbPgvjdCS0MeO+8kP1\nfOnM8D6dTnD8zI/yPrF+iC/az'\
    b'Aq/2Dv7M5LM477RujVPZc0vv02AoclILRA/4QYaTIL9Hz/p+5rN\njL4lPxU4AD+U'\
    b'3QO/XxnHc7mtDz91ZJkEdU8dv2bsyyfTGik/oO1IR+ONGb92wLc0A8DYvvbqpgjd'\
    b'\nXwK/Ooj6yMlt/D5yxAe0A5caP6IuNOxjDCC/TTijUuQEMj9z0Toy9WLkPoYeD0A'\
    b'lE+a+eWONw0A8\n2D7z4zrqYl8dP5NNPdfIY0G/zECsRyacS7+hev4A7ssPv9WSr6'\
    b'P6XR8/aim4mKgtJD83FOeaMywz\nv5gdTqt0mjk/XzyI8utMQT9p5M4PARnIPgYaA'\
    b'9RnYhw/eSdnpQp6Ez88T8Rt0TjAvsPBxnardQu/\nfBOiZbe8rr7nGZ+ijZULv/dW'\
    b'gbjaQ/u+qYnjdIiuCb8Gdp/yZB7Qvphs1hOA3/E+T47ro1Ex7z4N\nc80TaC03v97'\
    b'onl6tXhu/sN0vHR8bKL9sKR5q5N0FPyOmzXkC5Pg+EYBV0xJa9L6qB+tsPcYEP+V'\
    b'Y\nbFTWWLy+8Byqgt+e974Y7TX/XwZZP5/PADXUmUI/mMVPKqXZYz/HIrB2TH05v6'\
    b'dmDQmC1jg/4MX+\nbYXMSL97x/55aH9mvwQ6EMF+QFq/GIcm+VJqer+Cq+n7+Qq+v'\
    b'7+w0fgsmb6/V3IzfuvokL+UwEe8\nJjWZvxYtKxjZoIW/Z+2ptPNWcb+gEsjfKKaQ'\
    b'P39NNrVle4Q/cenGkOm0aj9v/60uoGDAP0IYsN9T\n7b4/IiIrFyTAlj/xLxE+qvX'\
    b'jvrkk7kSA3/S+4jpmofRW777URDPXW/0VP77GjqCy1gC/q3c0jrTx\n1j7Pzz/olJ'\
    b'09v8CbpcVa3Ss/kxkM5AsFMr8lx3RgPLIlP1bVtiwvMxC/T1bcUikmFz/rgUnfCP'\
    b'L3\nvpt+rUrgUSA/ERJOKtLbGb/sztsUEYAkP3P1Z7NMMPi+rq1AHK8xxL6ILkWBQ'\
    b'aoQP1HQjyOSjRu/\nNRMyYV99Iz9cVoE8pQvnPqlIz19qo82+ljZfQFLPzj7acLYS'\
    b'Nv8NP915zUFiny6/Q4mw+is0Pr8v\nesAJUisBv+pL9tPvBQ0/Vxayh/ZWGD99gHw'\
    b'lNSUlv5B+dV24JjE/p6VQ3/nSKz/hk02wPHnhvtWV\n5lNEBf0+X4wc7vu4DD/6LV'\
    b'0fogfyPs/IC/oZ/PC+wMmsshds4j5ZlaslUncDv6NcynNjfPE+xscl\naYVgB7/IO'\
    b'TLV0svDPpT8eaPAfOY+nTgwrRw7yr7gGUcRbKBVP9wsky7HKga/5kCwWsQQNr+DC'\
    b'wjF\nwkAGP8f2fEnu+Ro/EjnvrqAb974h+83SdigFP03mpbiL0Sg/llEk5l5aLD+a'\
    b'w6GOzp1dv/3rYtWc\nqIG/+XAL6szoYz+lyCTr7AYSv2An3cZEimk/TtputUG6TD8'\
    b'fg6M8eDGhv4C7i7FTF5G/ms5tWqwZ\nmb9W+ft/LxaXv2tNh0ihroi/J8HqSC+Nqr'\
    b'8pU3d3WN2JP8KBKBqCpHQ/M0gK1NpXbj9+SWC3w9KA\nP0SP5Skyu3o/pdg+wa1Pb'\
    b'T9EcOmcdGSiPyMiKxckwJY/H/unTgnzsD8=\n'

if __name__ == "__main__":
    main()
