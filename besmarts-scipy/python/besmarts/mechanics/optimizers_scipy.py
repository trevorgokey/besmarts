
"""
besmarts.mechanics.optimizers_scipy
"""

import copy
import scipy.optimize
import datetime

from besmarts.core import geometry
from besmarts.core import assignments
from besmarts.core import compute
from besmarts.core import arrays
from besmarts.core import logs
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import objectives
from besmarts.mechanics import fits


def optimize_positions_scipy(csys, psys: mm.physical_system):

    pos = copy.deepcopy(psys.models[0].positions[0])
    
    args, keys = objectives.array_flatten_assignment(pos.selections)

    # jac = objectives.array_geom_gradient
    # jac = None
    jac = True

    result = scipy.optimize.minimize(
        objectives.array_geom_energy_gradient,
        args,
        jac=jac,
        args=(keys, csys, psys),
        options={'disp': False, 'gtol': 1e-3, 'ftol': 1e-4, 'maxiter': 1000},
        method='L-BFGS-B',
    )

    n_confs = len(list(pos.selections.values())[0])
    for (c, n, i), v in zip(keys, result.x):
        pos.selections[n][c][i] = v

    return pos

def ff_obj(x, keys, csys, refpsys, ref, refcsys):
    # print("input array:")
    # print(x)
    pos = copy.deepcopy(ref)
    recalc = set([x[0] for x in keys])

    for k, v in zip(keys, x):
        mm.chemical_system_set_value(csys, k, v)

    reuse = set(range(len(csys.models))).intersection(recalc)
    psys = mm.chemical_system_to_physical_system(csys, pos, ref=refpsys, reuse=reuse)
    optpos = optimize_positions_scipy(csys, psys)

    obj = 0
    for ic, rxyz in ref[0].selections.items():
        oxyz = optpos.selections[ic]
        for a, b in zip(rxyz, oxyz):
            dx = geometry.array_distance(a, b)
            dx2 = dx*dx
            obj += dx2
            # print(obj, dx2)
    # for k in keys:
    #     newv = mm.chemical_system_get_value(csys, k)
    #     refv = mm.chemical_system_get_value(refcsys, k)
    #     print(k, f"| New: {newv} Ref {refv} Diff {newv-refv}")
    # print("Total obj:", obj)
    return obj

def objective_gdb(csys, gdb, obj, psysref=None, reuse=None, ws=None):

    tasks = {}
    for i, x in obj.items():
        psys = None
        if psysref:
            psys = psysref[x.addr.eid[0]]
        tasks[i] = x.get_task(gdb, csys, psys=psys, reuse=reuse)

    # print("Starting tasks")
    if ws:
        # wq = compute.workqueue_local('127.0.0.1', 55555)
        # ws = compute.workqueue_new_workspace(wq, shm={})

        results = compute.workspace_submit_and_flush(
            ws,
            fits.objective_run_distributed,
            {i: ([x], {}) for i, x in tasks.items()}
        )
        # print("############################################")
        
        # results = {i:v for i, v in enumerate(
        #     compute.workspace_submit_and_flush(
        #         ws,
        #         fits.objective_run,
        #         tasks.values()
        #     )
        # )}
        # ws.close()

    else:
        results = {}
        for i, task in tasks.items():
            # print(f"RUNNING {i}")
            results[i] = task.run()
        # results = {i: task.run() for i, task in tasks.items()}
    # print("Tasks done. Calculating objective")
    X = 0
    for i, x in obj.items():
        # print(f"Objective {i}", x)
        xi = x.compute_total(assignments.graph_db_get_entries(gdb, x.addr.eid), results[i])
        X += xi*x.scale
    return X

def objective_gradient_gdb(args, keys, csys, gdb, obj, history=None, psysref=None, reuse=None, ws=None, verbose=False):

    if history is None:
        history = []
    n = len(history)
    # csys = copy.deepcopy(csys)
    tasks = {}
    h = []
    # print("REUSE IS", reuse)
    logs.dprint(f"{logs.timestamp()} Generating {len(obj)} objectives", on=verbose)
    for i, x in obj.items():
        psys = None
        if psysref:
            psys = psysref[x.addr.eid[0]]
            # logs.dprint(f"{logs.timestamp()} Reparameterizing with reuse={reuse}", on=verbose)
            
            # get a copy
            # psys = mm.chemical_system_to_physical_system(
            #     csys,
            #     psys.models[0].positions,
            #     ref=psys,
            #     reuse=list(range(len(psys.models)))
            # )
            # now change values
            for k, v in zip(keys, args):
                mm.physical_system_set_value(psys, k, v)
        else:
            print("WARNING: No parameterized system given. This will recharge the molecules and is likely not intended.")
            csys = copy.deepcopy(csys)
            for k, v in zip(keys, args):
                mm.chemical_system_set_value(csys, k, v)
            psys = mm.chemical_system_to_physical_system(
                csys,
                psys.models[0].positions,
                ref=None,
                # reuse=[0,2,3,4,5]
            )
            reuse=list(range(len(psys.models)))

            #     logs.dprint(f"PSYS: setting {k}={v}", on=verbose)
        # dcsys = copy.deepcopy(csys)
        dcsys = csys
        dreuse=list(range(len(psys.models)))
        # dreuse=[2,3,4,5]
        tasks[(n, (i, 0))] = x.get_task(gdb, dcsys, keys={}, psys=psys, reuse=dreuse)
        for j, (k, v) in enumerate(zip(keys, args), 1):

            # hi = v*1e-4
            hi = 1e-6
            h.append(hi)
            # dreuse = [i for i in range(len(csys.models)) if i != k[0]]
            # dreuse = [2,3,4,5]
            # dreuse = reuse
            # dcsys = copy.deepcopy(csys)
            # dcsys = csys
            # mm.chemical_system_set_value(dcsys, k, v+hi)
            tasks[(n, (i,j))] = x.get_task(gdb, dcsys, keys={k: v+hi}, psys=psys, reuse=dreuse)
            # dcsys = copy.deepcopy(csys)
            # dcsys = csys
            # mm.chemical_system_set_value(dcsys, k, v-hi)
            tasks[(n, (i,-j))] = x.get_task(gdb, dcsys, keys={k: v-hi}, psys=psys, reuse=dreuse)
            # mm.chemical_system_set_value(csys, k, v)

    logs.dprint(f"{logs.timestamp()} Starting {len(tasks)} tasks", on=verbose)
    if ws:
        # wq = compute.workqueue_local('127.0.0.1', 55555)
        # ws = compute.workqueue_new_workspace(wq, shm={})

        # print(f"Starting {len(tasks)} tasks")
        results = compute.workspace_submit_and_flush(
            ws,
            fits.objective_run_distributed,
            {i: ([x], {}) for i, x in tasks.items()},
            chunksize=50,
            verbose=verbose
        )

    else:
        results = {}
        for i, task in tasks.items():
            print(f"RUNNING TASK {i}/{len(tasks)}")
            results[i] = task.run()
        # results = {i: task.run() for i, task in tasks.items()}
    # print("Tasks done. Calculating objective")
    logs.dprint(f"{logs.timestamp()} Calculating {len(tasks)} objectives", on=verbose)
    X = 0
    grad = list([0.0] * len(keys))
    for i, x in obj.items():
        # print(f"Objective {i}", x)
        gradx = list([0.0] * len(keys))
        ref = assignments.graph_db_get_entries(gdb, x.addr.eid)
        dx = x.compute_diff(ref, results[(n, (i,0))])
        # if verbose:
        #     print("Grad Differences:")
        #     print(dx)
        x2 = x.scale*arrays.array_inner_product(dx,dx)
        X += x2
        # if verbose:
        #     print(f"MM Grad for obj {i}:")
        # print("dx", dx)
        for j in range(1, len(keys)+1):
            dxb = results[(n, (i,j))]
            dxa = results[(n, (i,-j))]
            # print("DXA", dxa[0][0].graphs[0].rows[0].columns[0].selections)
            # print("DXB", dxb[0][0].graphs[0].rows[0].columns[0].selections)
            dxdp = x.scale*x.compute_gradient_2pt(dx, dxa, dxb, h[j-1])
            # print("dxdp", dxdp)
            grad[j-1] += dxdp
            gradx[j-1] += dxdp
            # if verbose:
            #     print(keys[j-1], gradx[j-1], h[j-1])
        
        gnormx = arrays.array_inner_product(gradx, gradx)**.5
        if verbose:
            print(f"  {i:04d} | X2= {x2:10.5g} |g|= {gnormx:10.5g}")

        
    gnorm = arrays.array_inner_product(grad, grad)**.5
    if verbose:
        print(f">>> X2= {X:10.5g} |g|={gnorm:10.5g}")
    logs.dprint(f"{logs.timestamp()} Done. {len(tasks)} tasks complete", on=verbose)
    history.append(X)
    return X, grad


def fit_gdb(args, keys, csys, gdb, obj, psysref=None, reuse=None, ws=None):

    for k, v in zip(keys, args):
        mm.chemical_system_set_value(csys, k, v)

    X = objective_gdb(csys, gdb, obj, psysref=psysref, reuse=reuse, ws=ws)

    return X

def fit_grad_gdb(args, keys, csys, gdb, obj, history=None, psysref=None, reuse=None, ws=None, verbose=False):

    if history is None:
        history = []

    # csys = copy.deepcopy(csys)
    # if reuse is None:
    #     reuse = set(range(len(csys.models)))
    for k, v in zip(keys, args):
        v0 = mm.chemical_system_get_value(csys, k)
        # mm.chemical_system_set_value(csys, k, v)
        # if k[0] in reuse:
        #     reuse.remove(k[0])
        if verbose:
            dv = v-v0
            # if abs(dv) < 1e-6:
            #     dv = 0.0
            # mm.chemical_system_set_value(csys, k, v0)
            print(f"Setting {k} from {v0:.10g} to {v:.10g} d={v-v0:.10g}")

    # reuse = list(reuse)
    X, g = objective_gradient_gdb(args, keys, csys, gdb, obj, history=history, psysref=psysref, reuse=reuse, ws=ws, verbose=verbose)
    # print(f"RETURN IS {X}")

    return X, g

def singlepoint_forcefield_gdb_scipy(x0, args, bounds=None):
    keys, csys, gdb, obj, psysref, reuse, ws = args
    y0, g0 = objective_gradient_gdb(args, keys, csys, gdb, obj, psysref=psysref, reuse=reuse, ws=ws)
    return y0, g0

def optimize_forcefield_gdb_scipy(x0, args, bounds=None, step_limit=None):

    keys, csys, gdb, obj, history, psysref, reuse, ws, verbose = args
    if verbose:
        print(datetime.datetime.now(), "Calculating initial obj")
        # for eid, psys in psysref.items():
        #     print(f"Positions EID {eid}")
        #     for ic, sel in psys.models[0].positions[0].selections.items():
        #         print(ic, sel)

        #     print(f"Parameterization EID {eid}")
        #     for m, pm in enumerate(psys.models[:2]):
        #         print(csys.models[m].name)
        #         for lbls, vals in zip(pm.labels, pm.values):
        #             for ic in lbls:
        #                 v = vals.get(ic)
        #                 print(" ", ic, lbls[ic], ":", v) 
    # y0 = None
    # y0, g0 = objective_gradient_gdb(x0, keys, csys, gdb, obj, history=history, psysref=psysref, reuse=reuse, ws=ws, verbose=verbose)

    opts = {'disp': verbose, 'maxls': 100, 'iprint': 1000 if verbose else 0, 'ftol': 1e-4, 'gtol': 1e-3}
    if step_limit:
        opts['maxiter'] = int(step_limit)
    if verbose:
        print(datetime.datetime.now(), "Calculating fit")
    result = scipy.optimize.minimize(
        fit_grad_gdb,
        x0,
        args=args,
        bounds=bounds,
        jac=True,
        options=opts,
        method='L-BFGS-B',
    )
    y0 = history[0]
    return result.x, y0, result.fun, result.jac


def optimize_forcefield_scipy(csys, pos):
    """
    Take a csys and a psys, then fit the params according to the flattened list
    The flattened list is the psys values.
    I need to append all psys values into the same list.
    Then, I need to take all modified values and reassign to psys
    """
    kv = mm.chemical_system_iter_keys(csys)
    kv = {k:v for k, v in kv.items() if k[0] in [0,1,2] and k[1] in "kl"}
    keys = list(kv.keys())
    x0 = list(kv.values())
    psys = mm.chemical_system_to_physical_system(csys, pos)
    # x0 = [x*.5 for x in x0]

    result = scipy.optimize.minimize(
        ff_obj,
        x0,
        args=(keys, csys, psys, pos, copy.deepcopy(csys)),
        options={'disp': True, 'gtol': .001, 'ftol': .001},
        method='L-BFGS-B'
    )
    return {k:v for k,v in zip(keys, result.x)}
