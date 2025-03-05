"""
besmarts.mechanics.force_pairwise
"""

import math
import subprocess
import tempfile
import os
import shutil

from besmarts.core import topology
from besmarts.core import assignments
# from besmarts.core import trees
from besmarts.core import codecs
# from besmarts.core import hierarchies
from besmarts.core import primitives
# from besmarts.core import graphs
# from besmarts.core import perception

from besmarts.mechanics import molecular_models as mm

# import math
INF = math.inf

# electrostatics
def energy_function_coulomb_mix(*, eps, c, s, qq, x) -> float:
    return [[
        INF if xi == 0.0 else
        s[0]*eps[0]*q/xi if xi < c[0] else
        0.0 for q in qq
    ] for x0 in x for xi in x0]

def force_function_coulomb_mix(*, eps, c, s, qq, x) -> float:
    return [[
        INF if xi == 0.0 else
        s[0]*eps[0]*q/(xi*xi) if xi < c[0] else
        0.0 for q in qq
    ] for x0 in x for xi in x0]

def force_gradient_function_coulomb_mix(*, eps, c, s, qq, x) -> float:
    return [[
        INF if xi == 0.0 else
        (q*-s[0]*eps[0]*2/(xi*xi*xi)) if xi < c[0] else
        0.0 for q in qq
    ] for x0 in x for xi in x0]

def force_system_coulomb_mix(*, eps, c, s, qq, x) -> float:
    return [[
        (q, INF) if xi == 0.0 else
        (q, s[0]*eps[0]/(xi*xi)) if xi < c[0] else
        0.0 for q in qq
    ] for x0 in x for xi in x0]

def force_gradient_system_coulomb_mix(*, eps, c, s, qq, x) -> float:
    return [
        [
            (q, INF) if xi == 0.0 else
            (q, -s[0]*eps[0]*2/(xi*xi*xi)) if xi < c[0] else
            0.0 for q in qq
        ]
        for xi in x[0]
    ]

# vdW
def energy_function_lennard_jones_combined(*, s, c, ee, rr, x) -> float:
    xx = [
        [
            INF if xi == 0.0 else
            math.pow(ri/xi, 6) if xi < c[0] else
            0.0 for ri in rr
        ] for x0 in x for xi in x0
    ]
    return [
        [4.0*s[0]*ei*(xi*xi - xi) for ei, xi in zip(ee, xxi)]
        for xxi in xx
    ]

def force_function_lennard_jones_combined(*, s, c, ee, rr, x) -> float:
    xx = [
        [
            INF if (ri > 0.0 and xi == 0.0) else
            math.pow(ri/xi, 6) if xi < c[0]
            else 0.0 for ri in rr
        ] for x0 in x for xi in x0
    ]
    return [
        [
            0.0 if ri == 0.0 else
            24.0*s[0]*ei*(2.0*xi*xi - xi)/x0
            for ei, xi, ri, x0 in zip(ee, xxi, rr, x[0])
        ]
        for xxi in xx
    ]

def force_gradient_function_lennard_jones_combined(*, s, c, ee, rr, x) -> float:
    xx = [
        [
            INF if (ri > 0.0 and xi == 0.0) else
            math.pow(ri/xi, 6) if xi < c[0] else
            0.0 for ri in rr
        ] for x0 in x for xi in x0
    ]
    return [
        [
            0.0 if ri == 0.0 else
            (-24.0*ei*s[0]*(26.0*xi*xi - 7.0*xi)/(x0*x0))
            for ei, xi, ri, x0 in zip(ee, xxi, rr, x[0])
        ]
        for xxi in xx
    ]

def force_system_lennard_jones_combined(*, s, c, ee, rr, x) -> float:
    xx = [
        [
            INF if (ri > 0.0 and xi == 0.0) else
            math.pow(ri/xi, 6) if xi < c[0] else
            0.0 for ri in rr
        ] for x0 in x for xi in x0
    ]
    return [
        [
            (ei, 0.0) if ri == 0.0 else
            (ei, 24.0*s[0]*(2.0*xi*xi - xi)/x0)
            for ei, xi, ri, x0 in zip(ee, xxi, rr, x[0])
        ]
        for xxi in xx
    ]

def force_gradient_system_lennard_jones_combined(*, s, c, ee, rr, x) -> float:
    xx = [
        [
            INF if (ri > 0.0 and xi == 0.0) else
            math.pow(ri/xi, 6) if xi < c[0] else
            0.0 for ri in rr
        ] for x0 in x for xi in x0
    ]
    return [
        [
            (ei, 0.0) if ri == 0.0 else
            (ei, -24.0*s[0]*(26.0*xi*xi - 7.0*xi)/(x0*x0))
            for ei, xi, ri, x0 in zip(ee, xxi, rr, x[0])
        ]
        for xxi in xx
    ]


class chemical_model_procedure_antechamber(mm.chemical_model_procedure):
    """
    """

    def __init__(self, topology_terms):
        self.name = ""
        self.topology_terms = topology_terms
        self.procedure_parameters = {}
        self.minimize = False

    def assign(self, cm: mm.chemical_model, pm: mm.physical_model, overrides=None) -> mm.physical_model:
        """
        """
        symbol = "qq"
        pm.values = []
        cdc = codecs.primitive_codec_formal_charge()
        tmpfolder = tempfile.mkdtemp()
        for pi, pos in enumerate(pm.positions):
            
            charges = {}
            labels = {}

            q = int(cdc.count_charge_smiles(pos.graph.nodes))
            nconfs = min((len(x) for x in pos.selections.values()))
            # average over the confs
            # if nconfs > 0:
                # print(f"NConfs for charging: {nconfs}. Using first in gas phase")
            for ci in range(1, nconfs+1):
                
                conf_charges = {}
                if ci > 1:
                    break
                if len(pos.graph.nodes) == 1:
                    conf_charges = {next(iter(pos.graph.nodes)): q}
                else:
                    with open(os.path.join(tmpfolder, "mdin"), "w") as f:
                        steps = "10000" if self.minimize else "0"
                        f.write(f"\n&qmmm\n")
                        f.write(f"qm_theory='AM1', maxcyc={steps}, grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700, qmcharge={q:d},\n")

                        if "sqm" in self.procedure_parameters:
                            sqm_opts = self.procedure_parameters['sqm']
                            f.write(f"{sqm_opts}\n")

                        f.write(" /\n")

                        for j, (n, xyz) in enumerate(pos.selections.items(), 1): 
                            x,y,z = xyz[ci-1]
                            elem = pos.graph.nodes[n[0]].primitives["element"].on()[0]
                            name = primitives.element_tr[str(elem)] + str(j)
                            # print(f"{elem} {name} {x:12.9f} {y:12.9f} {z:12.9f}")
                            f.write(f"{elem} {name} {x:12.9f} {y:12.9f} {z:12.9f}\n")
                        f.write("\n")

                    # print("RUNNING ANTECHAMBER FOR CHARGES")
                    subprocess.run([
                        "sqm", "-O",
                        "-i", "mdin",
                        "-o", "mdout"
                        ], 
                        cwd=tmpfolder,
                    )

                    subprocess.run([
                        "antechamber",
                        "-c", "bcc", 
                        "-nc", f"{q}",
                        "-pf", "y",
                        "-dr", "n",
                        "-fi", "sqmout",
                        "-i", "mdout",
                        "-fo", "mol2",
                        "-o", "out.mol2"
                        ],
                        cwd=tmpfolder,
                        capture_output=True
                    )
                    subprocess.run([
                        "antechamber",
                        "-c", "wc", 
                        "-cf", "q.dat", 
                        "-nc", f"{q}",
                        "-pf", "y",
                        "-dr", "n",
                        "-fi", "mol2",
                        "-i", "out.mol2",
                        "-fo", "mol2",
                        "-o", "out.mol2"
                        ],
                        cwd=tmpfolder,
                        capture_output=True
                    )
                    try:
                        with open(os.path.join(tmpfolder, "q.dat")) as f:
                            qdat = f.read().split()
                            qm_q =  [*map(float, qdat)]
                            # call this am1 error
                            print(pos.smiles)
                            print("Charges:", pos.smiles, qm_q)
                            delta = abs(sum(qm_q) - q)
                            if delta < .5 and delta > 0.0:
                                if sum(qm_q) > 0.0:
                                    p = [x for x in qm_q if x > 0]
                                else:
                                    p = [x for x in qm_q if x < 0]
                                s = (sum(p) - sum(qm_q))/sum(p)
                                qm_q = [s*x if x > 0 else x for x in qm_q]
                            else:
                                assert q == sum(qm_q)
                            conf_charges = {
                                i: x for i, x in zip(pos.graph.nodes, qm_q)
                            }

                    except Exception as e:
                        print("Error:", e)
                        print("Could not generate charges, input of error is below:")
                        mdin = os.path.join(tmpfolder, "mdin")
                        if os.path.exists(mdin):
                            with open(mdin) as f:
                                for line in f:
                                    print(line, end='')
                        print("Could not generate charges, output of error is below:")
                        mdout = os.path.join(tmpfolder, "mdout")
                        if os.path.exists(mdout):
                            with open(mdout) as f:
                                for line in f:
                                    print(line, end='')
                        raise e
                for i, qi in conf_charges.items():
                    key = (pi, i)
                    if i not in charges:
                        charges[key] = {"q": [0.0]}
                        # don't add labels because we want to hide these from
                        # functions that flatten the parameter list for fitting.
                        # Clearly we cannot optimize a black-box method such as
                        # this one.
                        # If we do want this, make sure to add the values to
                        # cm.topology_terms['q'].values
                        # labels[i,] = {"q": "am1bcc"}
                    charges[key]["q"][0] = ((ci-1)*charges[key]["q"][0] + qi)/ci
            pm.values.append(charges)
            pm.labels.append(labels)
        shutil.rmtree(tmpfolder)
        return pm

class chemical_model_procedure_combine_coulomb(mm.chemical_model_procedure):
    """
    """

    def __init__(self, top_parm):
        self.name = ""
        assert "qq" in top_parm

    def assign(self, cm, pm, overrides=None):
        params = pm.values

        pos = pm.positions

        indices = assignments.graph_assignment_matrix_pair_indices(pos)
        indices.extend(assignments.graph_assignment_matrix_bond_indices(pos))

        for ki, kj in indices:
        # for pi, posi in enumerate(pm.positions):

            # mixed = {}
            # pairs = assignments.smiles_assignment_geometry_distances(
            #     posi,
            #     graphs.graph_pairs(posi.graph)
            # ).selections

            # bonds = assignments.smiles_assignment_geometry_distances(
            #     posi,
            #     graphs.graph_bonds(posi.graph)
            # ).selections

            # pairs.update(bonds)

            if len(ki) == 2:
                pi, i = ki
            elif len(ki) == 3:
                pi, i, ci = ki
            else:
                assert False
            # i = i,
            parami = params[pi]

            if len(kj) == 2:
                pj, j = kj
            elif len(kj) == 3:
                pj, j, cj = kj
            else:
                assert False
            # j = j,
            paramj = params[pj]

            if (pi, i) not in parami or (pj, j) not in paramj:
                continue

            key = ki, kj
            if key not in parami:
                parami[key] = {}

            ei = parami[(pi, i)].get("q")
            ej = paramj[(pj, j)].get("q")

            parami[key]["qq"] = [(e1*e2) for e1,e2 in zip(ei, ej)]

            # parami.update(mixed)

            # for pj, posj in enumerate(pm.positions):
            #     if pj < pi:
            #         continue
            #     paramj = params[pj]
            #     pairs = [((pi, i), (pj, j)) for i in posi.graph.nodes for j in posj.graph.nodes]
            #     for (pi, i),(pj, j) in pairs:
            #         if (pi, i) not in parami or (pj, j) not in paramj:
            #             continue
            #         key = (pi, i),(pj, j)
            #         if key not in mixed:
            #             mixed[key] = {}
            #         ei = parami[(pi, i)].get("q")
            #         ej = paramj[(pj, j)].get("q")
            #         mixed[key]["qq"] = [(e1*e2) for e1,e2 in zip(ei, ej)]
            #     parami.update(mixed)

        return pm

class chemical_model_procedure_combine_lj_lorentz_berthelot(mm.chemical_model_procedure):
    """
    """

    def __init__(self, top_parm):
        assert "ee" in top_parm
        assert "rr" in top_parm

    def assign(self, cm, pm, overrides=None):
        params = pm.values
        pos = pm.positions

        indices = assignments.graph_assignment_matrix_pair_indices(pos)
        indices.extend(assignments.graph_assignment_matrix_bond_indices(pos))

        for ki, kj in indices:

            # pi, i = ki[:2]
            # # i = i,
            # parami = params[pi]

            # pj, j = kj[:2]
            # # j = j,
            # paramj = params[pj]

            # if (pi, i) not in parami or (pj, j) not in paramj:
            #     continue

            # # key = ki, kj
            # key = (pi, i), (pj, j)
            if len(ki) == 2:
                pi, i = ki
            elif len(ki) == 3:
                pi, i, ci = ki
            else:
                assert False
            # i = i,
            parami = params[pi]

            if len(kj) == 2:
                pj, j = kj
            elif len(kj) == 3:
                pj, j, cj = kj
            else:
                assert False
            # j = j,
            paramj = params[pj]

            if (pi, i) not in parami or (pj, j) not in paramj:
                continue

            key = ki, kj
            if key not in parami:
                parami[key] = {}

            ei = parami[(pi, i)].get("e")
            ej = paramj[(pj, j)].get("e")

            parami[key]["ee"] = [(e1*e2)**.5 for e1,e2 in zip(ei, ej)]

            ri = parami[(pi, i)].get("r")
            rj = paramj[(pj, j)].get("r")

            parami[key]["rr"] = [(r1+r2)/2.0 for r1,r2 in zip(ri, rj)]

            # parami = params[pi]
            # for i,j in pairs:
            #     if (pi, i) not in parami or (pi, j) not in parami:
            #         continue
            #     key = (pi, i),(pi, j)
            #     if key not in mixed:
            #         mixed[key] = {}
            #     ei = parami[i,].get("e")
            #     ej = parami[j,].get("e")
            #     if ei is None or ej is None:
            #         continue
            #     # mixed[key]["ee"] = [(e1*e2) for e1,e2 in zip(ei, ej)]
            #     mixed[key]["ee"] = [(e1*e2)**.5 for e1,e2 in zip(ei, ej)]
            #     ri = parami[i,].get("r")
            #     rj = parami[j,].get("r")
            #     if ri is None or rj is None:
            #         continue
            #     mixed[key]["rr"] = [(r1+r2)/2.0 for r1,r2 in zip(ri, rj)]

            # parami.update(mixed)

            # for pj, posj in enumerate(pm.positions):
            #     if pj < pi:
            #         continue
            #     paramj = params[pj]
            #     pairs = [(i, j) for i in posi.graph.nodes for j in posj.graph.nodes]
            #     for i,j in pairs:
            #         if (pi, i) not in parami or (pj, j) not in paramj:
            #             continue

            #         key = (pi, i),(pi, j)
            #         if key not in mixed:
            #             mixed[key] = {}
            #         ei = parami[i,].get("e")
            #         ej = parami[j,].get("e")
            #         if ei is None or ej is None:
            #             continue
            #         # mixed[key]["ee"] = [(e1*e2) for e1,e2 in zip(ei, ej)]
            #         mixed[key]["ee"] = [(e1*e2)**.5 for e1,e2 in zip(ei, ej)]
            #         ri = parami[i,].get("r")
            #         rj = parami[j,].get("r")
            #         if ri is None or rj is None:
            #             continue
            #         mixed[key]["rr"] = [(r1+r2)/2.0 for r1,r2 in zip(ri, rj)]

            #     parami.update(mixed)

        return pm

def chemical_model_coulomb(perception):

    cm = mm.chemical_model("Q", "Electrostatics", topology.pair)

    cm.energy_function = energy_function_coulomb_mix
    cm.force_function = force_function_coulomb_mix
    cm.force_gradient_function = force_gradient_function_coulomb_mix

    cm.force_system = force_system_coulomb_mix
    cm.force_gradient_system = force_gradient_system_coulomb_mix

    cm.internal_function = assignments.graph_assignment_geometry_pair_matrix
    cm.derivative_function = assignments.graph_assignment_jacobian_pair_matrix

    cm.system_terms = {
        "c": mm.system_term("cutoff", "c", "float", "A", [10.0], ""),
        "eps": mm.system_term("dielectric", "eps", "float", "kcal/mol*A/e/e", [332.06371329919216], ""),
        # originally using 332.063711 
        #OpenMM ONE_4PI_EPS0 gives me 138.93545764438198 (in kJ and nm presumably)
        # ONE_4PI_EPS0 is in openmm/platforms/reference/include/SimTKOpenMMRealType.h
        #converting to kcal and ang gives me 332.06371329919216 and this makes me match
        #every digit down to 1e-16
    }

    cm.topology_terms = {
        "qq": mm.topology_term("qq", "partial_charge", "float", "e*e", {}, "", {}),
        "q": mm.topology_term("q", "partial_charge", "float", "e", {}, "", {}),
    }
    cm.topology_terms["s"] = mm.topology_term("s", "scale", "float", "", {}, "", {})

    return cm

def chemical_model_lennard_jones(perception) -> mm.chemical_model:

    cm = mm.chemical_model("N", "vdW", topology.pair)

    cm.energy_function = energy_function_lennard_jones_combined
    cm.force_function = force_function_lennard_jones_combined
    cm.force_gradient_function = force_gradient_function_lennard_jones_combined

    cm.force_system = force_function_lennard_jones_combined
    cm.force_gradient_system = force_gradient_system_lennard_jones_combined

    cm.internal_function = assignments.graph_assignment_geometry_pair_matrix
    cm.derivative_function = assignments.graph_assignment_jacobian_pair_matrix

    cm.topology_terms = {
        "e": mm.topology_term("e", "depth", "float", "kcal/mol", {}, "", {}),
        "r": mm.topology_term("r", "sigma", "float", "A", {}, "", {}),
        "ee": mm.topology_term("ee", "depth_combine", "float", "kcal/mol", {}, "", {}),
        "rr": mm.topology_term("rr", "sigma_combine", "float", "A", {}, "", {}),
        "s": mm.topology_term("s", "scale", "float", "", {}, "", {})
    }

    cm.system_terms = {
        "c": mm.system_term("cutoff", "c", "float", "A", [9.0], ""),
    }

    return cm
