
"""
besmarts.mechanics.optimizers_scipy
"""

import copy
import scipy.optimize
import numpy as np
import datetime

from besmarts.core import geometry
from besmarts.core import assignments
from besmarts.core import compute
from besmarts.core import arrays
from besmarts.core import logs
from besmarts.core import returns
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import objectives
from besmarts.mechanics import fits


def optimize_positions_scipy(csys, psys: mm.physical_system, step_limit=1000, tol=1e-10):

    pos = copy.deepcopy(psys.models[0].positions[0])
    
    args, keys = objectives.array_flatten_assignment(pos.selections)

    # jac = objectives.array_geom_gradient
    # jac = None
    jac = True
    hess = objectives.array_geom_hessian

    opts = {'disp': False, 'maxls': 10, 'gtol': tol, 'ftol': 1e-9, 'maxiter': step_limit}
    # opts = {'iprint': 11000}
    # opts = {'disp': True}
    method = 'L-BFGS-B'
    opts = {'disp': False, 'ftol': tol, 'gtol': tol, 'maxiter': step_limit, 'maxls': 1000, 'maxcor': len(args)**2}
    hess = None
    # method = 'Newton-CG'
    # opts = {'disp': False, 'maxiter': step_limit}
    # method = 'Nelder-Mead'
    # print("SciPy Minimizer Starting")
    # method = 'CG'
    # method = 'BFGS'
    # opts = {'disp': False, 'gtol': tol, 'maxiter': step_limit}
    # method = 'trust-exact'
    # opts = {'disp': False, 'gtol': tol, 'maxiter': step_limit, 'max_trust_radius': 1e-1, 'initial_trust_radius': 1e-3 }
    result = scipy.optimize.minimize(
        objectives.array_geom_energy_gradient,
        args,
        jac=jac,
        hess=hess,
        args=(keys, csys, psys),
        options=opts,
        tol=tol,
        method=method,
    )

    n_confs = len(list(pos.selections.values())[0])
    for (c, n, i), v in zip(keys, arrays.array_round(result.x, 12)):
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

def objective_hessp(x, p, *args):
    history = args[6]
    # return arrays.array_inner_product(list(history[-1][3]), list(p))
    return np.array(arrays.array_multiply(list(history[-1][3]), list(p)))

def objective_hess(x, *args):

    history = args[6]

    diag = history[-1][3]

    hess = []
    for i in range(len(diag)):
        hess.append([0 for j in range(len(diag))])
    for i in range(len(diag)):
        hess[i][i] = diag[i]

    return hess

def objective_gradient_gdb(args, keys, csys, gdb, obj, priors, penalties=None, history=None, psysref=None, reuse=None, ws=None, verbose=False):

    if history is None:
        history = []
    if penalties is None:
        penalties = []
    n = len(history)

    out = []
    args = arrays.array_round(args, 12)
    # csys = copy.deepcopy(csys)
    tasks = {}
    h = []
    # print("REUSE IS", reuse)
    line = f"{logs.timestamp()} Generating {len(obj)} objectives"
    out.append(line)
    if verbose:
        print(line)
    dcsys = csys
    if ws:
        dcsys = None
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
            reapply = set() 
            kv = {}
            for k, v, p in zip(keys, args, priors):
                v = p[1]*v + p[0]
                kv[k] = v
                # print(f"Reassigning {k}")
                # print(f"Setting pval to {k}={v}")
                mm.physical_system_set_value(psys, k, v)
                reapply.add(k[0])
            for m in reapply:
                procs = csys.models[m].procedures
                if len(procs) > 1:
                    # print(f"Reapplying {m}:")
                    # print(psys.models[m].labels)
                    # print(psys.models[m].values)
                    for _ in range(1, len(psys.models[m].values)):
                        psys.models[m].values.pop()
                        psys.models[m].labels.pop()
                    for proc in procs[1:]:
                        psys.models[m] = proc.assign(csys.models[m], psys.models[m], overrides={k[1:]: v for k, v in kv.items() if k[0] == m})
                    # print(f"New values {m}:")
                    # print(psys.models[m].labels)
                    # print(psys.models[m].values)
        else:
            line = "WARNING: No parameterized system given. This will recharge the molecules and is likely not intended."
            out.append(line)
            if verbose:
                print(line)
            csys = copy.deepcopy(csys)
            for k, v in zip(keys, args):
                v = p[1]*v + p[0]
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
        dreuse=list(range(len(psys.models)))
        # dreuse=[2,3,4,5]
        # tasks[(n, (i, 0))] = x.get_task(gdb, dcsys, keys=[{}], psys=psys, reuse=dreuse)
        # print(psys.models[2].labels)
        # print(f"Opt torsions:")
        # print(psys.models[2].labels)
        # print(f"Opt torsion values:")
        # print(psys.models[2].values)

        # next part assumes the f_0 is first in the list
        # so just force it here for simplicity
        grad_keys = [{}]
        tasks[(n, (i,(0,)))] = x.get_task(gdb, dcsys, keys=grad_keys, psys=psys, reuse=dreuse)
        
        grad_keys = []
        if x.batch_size is None:
            batch_size = len(keys)
        else:
            batch_size = x.batch_size

        for batch in arrays.batched(enumerate(zip(keys, args, priors), 1), batch_size):

            kbatch = tuple([b[0] for b in batch])
            for j, (k, v, p) in batch:
                v = round(v*p[1] + p[0], 12)
                # hi = v*1e-4
                hi = 1e-6
                h.append(hi)
                grad_keys.extend([{k: v-hi}, {k:v+hi}])
                # dreuse = [i for i in range(len(csys.models)) if i != k[0]]
                # dreuse = [2,3,4,5]
                # dreuse = reuse
                # dcsys = copy.deepcopy(csys)
                # dcsys = csys
                # mm.chemical_system_set_value(dcsys, k, v+hi)
                # tasks[(n, (i,j))] = x.get_task(gdb, dcsys, keys=[{k: v-hi}, {k:v+hi}] , psys=psys, reuse=dreuse)
                # dcsys = copy.deepcopy(csys)
                # dcsys = csys
                # mm.chemical_system_set_value(dcsys, k, v-hi)
                # tasks[(n, (i,-j))] = x.get_task(gdb, dcsys, keys={k: v-hi}, psys=psys, reuse=dreuse)
                # mm.chemical_system_set_value(csys, k, v)
            tasks[(n, (i,kbatch))] = x.get_task(gdb, dcsys, keys=grad_keys, psys=psys, reuse=dreuse)
            grad_keys = []


    line = f"{logs.timestamp()} Running {len(tasks)} tasks"
    out.append(line)
    if verbose:
        print(line)

    if ws:
        # wq = compute.workqueue_local('127.0.0.1', 55555)
        # ws = compute.workqueue_new_workspace(wq, shm={})

        # print(f"Starting {len(tasks)} tasks")
        chunksize = 100
        if len(tasks)//40 < ws.nproc:
            chunksize = len(tasks)
        results = compute.workspace_submit_and_flush(
            ws,
            fits.objective_run_distributed,
            {i: ([x], {}) for i, x in tasks.items()},
            chunksize=chunksize,
            verbose=verbose,
            batchsize=10000,
            timeout=0.00
        )

    else:
        results = {}
        for i, task in tasks.items():
            # print(f"RUNNING TASK {i}/{len(tasks)}")
            results[i] = task.run()
        # results = {i: task.run() for i, task in tasks.items()}
    # print("Tasks done. Calculating objective")
    X = 0
    Y = 0
    grad = list([0.0] * len(keys))
    hess = list([0.0] * len(keys))
    grady = list([0.0] * len(keys))
    obj_tasks = {}

    h = tuple(1e-6/j[1] for j in priors)

    i = -1

    for i, x in obj.items():
        # print(f"Objective {i}", x)
        if not x.include and not verbose:
            continue
        if verbose:
            gradx = list([0.0] * len(keys))

        ref = assignments.graph_db_get(gdb, x.addr)
        my_results = {}
        for (n, (ii, jlist)), work in results.items():

            if ii != i:
                continue

            if (n, (i, 0)) not in my_results:
                # just pull the first directly
                my_results[(n, (i, 0))] = results[(n, (i,(0,)))]

            # then skip if that is the key currently have
            if jlist == (0,):
                continue

            my_results[(n, (i, 0))].extend(work)
        obj_tasks[("obj", n, i)] = [(x, ref, my_results[(n, (i, 0))], h), {"verbose": verbose}]


    lines = [
        f"{logs.timestamp()} Calculating {len(obj_tasks)} objectives",
        f"=== Objectives ==="
    ]
    if verbose:
        print(lines[0])
        print(lines[1])
    out.extend(lines)

    if ws:
        obj_results = compute.workspace_submit_and_flush(
            ws,
            run_objective,
            obj_tasks,
            chunksize=None,
            verbose=verbose,
            batchsize=5000,
            timeout=.00
        )
    else:
        obj_results = {i : run_objective(*args, **kwds)
            for i, (args, kwds) in obj_tasks.items()
        }
    for idx in sorted(obj_results, key=lambda x: x[2]):
        if idx[0] != "obj":
            continue
        jid, n, i = idx
        ret = obj_results[idx]
        (x2, gradxi, hessxi, y2, gradyi) = ret.value
        out.extend(ret.out)
        if verbose:
            for line in ret.out:
                print(line)
        X += x2
        Y += y2
        grad = arrays.array_add(grad, gradxi)
        hess = arrays.array_add(hess, hessxi)
        grady = arrays.array_add(grady, gradyi)

        # if verbose:
        #     print(f"  {i:04d} | "+ output)

    kv = {k: v*p[1] + p[0] for k,v,p in zip(keys,args,priors)}
    # kv = dict(zip(keys,args))

    for i, pen in enumerate(penalties, i+1):
        pen: fits.objective_config_penalty


        # generate the reference target values and then run to get deltas
        dx = pen.get_task().run(kv)
        ret = pen.compute_diff(dx, verbose=verbose)
        out.extend(ret.out)

        if verbose:
            for line in ret.out:
                print(line)
        x2 = pen.scale*ret.value
        X += x2
        # if verbose:
        #     print(f"MM Grad for obj {i}:")
        # print("dx", dx)
        gradx = list([0.0] * len(keys))
        for j, k in enumerate(keys, 0):
            if k not in dx:
                continue
            kvi = {k: dx[k]}
            dx2dp = pen.compute_gradient(kvi)
            dx2dp *= pen.scale
            dx2dp = round(dx2dp, 12)
            # print("dx2dp", dx2dp)
            grad[j] += dx2dp
            gradx[j] += dx2dp

            d2x2dp2 = pen.compute_hessian(kvi)
            d2x2dp2 *= pen.scale
            d2x2dp2 = round(d2x2dp2, 12)
            # print("dxdp", dxdp)
            hess[j] += d2x2dp2

            gnormx = arrays.array_inner_product(gradx, gradx)**.5
            line = f">> S= {pen.scale: 14.6f} R2= {x2: 14.6f} |g|= {gnormx: 14.6f}"
            if verbose:
                print(line)
            out.append(line)

    for k, g in zip(keys, grad):
        line = f"{k} {g}"
        out.append(line)
        if verbose:
            print(line)

    gnorm = arrays.array_inner_product(grad, grad)**.5
    hnorm = arrays.array_inner_product(hess, hess)**.5
    if verbose:
        gnormy = arrays.array_inner_product(grady, grady)**.5
        dX = 0.0
        if len(history):
            dX = X - min((x[0] for x in history))
        out.append(f">>> {datetime.datetime.now()} Step= {len(history)+1:3d} X2= {X: 14.6e} |g|= {gnorm: 14.6e} |h|= {hnorm: 14.6e} DX2= {dX: 14.6e} Y2= {Y: 14.6e} |gy|={gnormy: 14.6e}")
        if verbose:
            print(out[-1])
    out.append("Total Parameter Grad:")
    if verbose:
        print(out[-1])
    out.append(f"{logs.timestamp()} Done. {len(tasks)} tasks complete")
    if verbose:
        print(out[-1])
    history.append((X, grad, args, hess, out))
    return X, grad

def run_objective(x, ref, results, h, verbose=False, shm=None):
    X = 0
    Y = 0
    N = len(h)
    gradx = list([0.0] * N)
    hessx = list([0.0] * N)
    grady = list([0.0] * N)

    ret = x.compute_diff(ref, results[0], verbose=verbose)
    dx = ret.value
    x2 = x.scale*arrays.array_inner_product(dx,dx)
    
    # if ret.out and verbose:
    #     logs.dprint(ret.out, on=verbose)
    # if ret.err:
    #     logs.dprint(ret.err, on=True)

    # dx = arrays.array_sum(dx)

    for j in range(N):
        # dxb = results[(n, (i,j))]
        dxa = results[2*j+1]
        dxb = results[2*j+2]
        # print("DXA", dxa[0][0].graphs[0].rows[0].columns[0].selections)
        # print("DXB", dxb[0][0].graphs[0].rows[0].columns[0].selections)

        # this is the objective grad wrt parameter j; nobj numbers

        dxdp, d2xdp2 = x.compute_gradient_2pt(ref, dx, dxa, dxb, h[j])

        # this is the objective2 grad wrt parameter j (summed over objectives)
        dx2dp = 2 * arrays.array_inner_product(dx, dxdp)

        # d2xdp2 = x.compute_gradient_2pt(ref, None, dxa, dxb, h*h)
        d2x2dp2 = 2 * (arrays.array_inner_product(dx, d2xdp2) + sum(d2xdp2))

        dx2dp *= x.scale
        d2x2dp2 *= x.scale

        # print("dxdp", dxdp)
        if x.include:
            gradx[j] += dx2dp
            hessx[j] += d2x2dp2
        else:
            grady[j] += dx2dp



    sym = ""
    if x.include:
        sym = "X2"
        X += x2
        gnorm = arrays.array_inner_product(gradx, gradx)**.5
        hnorm = arrays.array_inner_product(hessx, hessx)**.5
    else:
        sym = "Y2"
        Y += x2
        gnorm = arrays.array_inner_product(grady, grady)**.5
        hnorm = 0

    output = f">> S= {x.scale: 14.6f} {sym}= {x2: 14.6e} |g|= {gnorm: 14.6e} |h|= {hnorm:14.6e}"
    ret.out.append(output)
    return returns.success((X, gradx, hessx, Y, grady), out=ret.out)


def fit_gdb(args, keys, csys, gdb, obj, psysref=None, reuse=None, ws=None):

    for k, v in zip(keys, args):
        mm.chemical_system_set_value(csys, k, v)

    X = objective_gdb(csys, gdb, obj, psysref=psysref, reuse=reuse, ws=ws)

    return X

def fit_grad_gdb(args, keys, csys, gdb, obj, priors, penalties=None, history=None, psysref=None, reuse=None, ws=None, verbose=False):

    if history is None:
        history = []

    if penalties is None:
        penalties = []

    args = arrays.array_round(args, 12)

    out = []
    # csys = copy.deepcopy(csys)
    # if reuse is None:
    #     reuse = set(range(len(csys.models)))
    for k, v, p in zip(keys, args, priors):
        v0 = mm.chemical_system_get_value(csys, k)
        # mm.chemical_system_set_value(csys, k, v)
        # if k[0] in reuse:
        #     reuse.remove(k[0])
        v1 = p[0] + v*p[1]

        out.append(f"Setting {str(k):20s} from {v0:15.7g} to {v1:15.7g} d={v1-v0:15.7g}")
        if verbose:
            dv = v-v0
            # if abs(dv) < 1e-6:
            #     dv = 0.0
            # mm.chemical_system_set_value(csys, k, v0)
            print(out[-1])

    # reuse = list(reuse)
    X, g = objective_gradient_gdb(args, keys, csys, gdb, obj, priors, penalties=penalties, history=history, psysref=psysref, reuse=reuse, ws=ws, verbose=verbose)
    # print(f"RETURN IS {X}")
    return X, g

def singlepoint_forcefield_gdb_scipy(x0, args, bounds=None):
    keys, csys, gdb, obj, psysref, reuse, ws = args
    y0, g0 = objective_gradient_gdb(args, keys, csys, gdb, obj, psysref=psysref, reuse=reuse, ws=ws)
    return y0, g0

def print_parameterization(csys, psysref):
    for eid, psys in psysref.items():
        print(f"Positions EID {eid}")
        for ic, sel in psys.models[0].positions[0].selections.items():
            print(ic, sel)

        print(f"Parameterization EID {eid}")
        for m, pm in enumerate(psys.models[:2]):
            print(csys.models[m].name)
            for lbls, vals in zip(pm.labels, pm.values):
                for ic in lbls:
                    v = vals.get(ic)
                    print(" ", ic, lbls[ic], ":", v) 

def optimize_forcefield_gdb_scipy(x0, args, bounds=None, step_limit=None, maxls=20, anneal=False):

    keys, csys, gdb, obj, priors, penalties, history, psysref, reuse, ws, verbose = args

    # if verbose:
    #     print_parameterization(csys, psysref)

    out = []

    hessp = objective_hessp
    hess = None
    if 1:
        method = 'L-BFGS-B'
        opts = {
            'maxls': maxls,
            'maxcor': len(args)**2,
            'iprint': 1000 if verbose else 0,
            'ftol': 1e-9,
            'gtol': 1e-9
        }
        hessp = None

    # slow
    elif 0:
        method = 'TNC'
        opts = {
            'xtol': 1e-6,
            'disp': False,
            'maxCGit': maxls,
            # 'iprint': 1000 if verbose else 0,
            'xtol': 1e-3,
            'ftol': 1e-4,
            'gtol': 1e-3,

        }
    elif 0:
        method = 'Newton-CG'
        opts = {
            'xtol': 1e-6,
            'disp': verbose,
            # 'iprint': 1000 if verbose else 0,
            # 'tol': 1e-4,
            # 'gtol': 1e-3,

        }
        # bounds=None

    elif 0:
        # decent first step but seems to hit local minima
        method = 'TNC'
        opts = {
            'xtol': 1e-6,
            'disp': True,
            'maxCGit': maxls,
            # 'iprint': 1000 if verbose else 0,
            'xtol': 1e-4,
            'ftol': 1e-4,
            'gtol': 1e-4,
            'scale': np.full(len(x0), 1.0),
            'offset': np.full(len(x0), 0.0),
            'maxfun': step_limit,
            'eta': .50
        }
    elif 0:
        # slow
        method = 'trust-krylov'
        opts = {
            'disp': verbose,
            'gtol': 1e-6,
            'inexact': True,

        }
    elif 0:
        # slow
        method = 'trust-exact'
        opts = {
            'disp': verbose,
            'gtol': 1e-3,

        }
        hess=objective_hess
        hessp=None
    elif 0:
        # slow
        method = 'trust-ncg'
        opts = {
            'disp': verbose,
            'gtol': 1e-4,
        }


    if step_limit:
        opts['maxiter'] = int(step_limit)


    if verbose:
        print(datetime.datetime.now(), "Starting physical parameter optimization")

    if anneal is True:
        anneal = 100
    else:
        anneal = int(anneal)
    if anneal:
        # all x must be bound, if unbounded set to += 50 priors
        print(datetime.datetime.now(), f"Annealing enabled n={anneal}")
        B = scipy.optimize.Bounds()
        B.lb = []
        B.ub = []
        for b in bounds:

            if b[0] is None:
                B.lb.append(-50)
            else:
                B.lb.append(b[0])

            if b[1] is None:
                B.ub.append(50)
            else:
                B.ub.append(b[1])




        kwds = {
            "jac": True,
            "hessp": hessp,
            "hess": hess,
            "tol": 1e-4,
            "method": method,
            "options": opts,
            "bounds": bounds,
            "args": args,
        }
        # args = tuple((*args, False))
        # result = scipy.optimize.dual_annealing(
        #     fit_grad_gdb,
        #     B,
        #     args=args,
        #     minimizer_kwargs=kwds,
        #     x0=x0,
        # )
        result = scipy.optimize.basinhopping(
            fit_grad_gdb,
            x0,
            disp=True,
            niter=anneal,
            minimizer_kwargs=kwds,
        )

    else:
        result = scipy.optimize.minimize(
            fit_grad_gdb,
            x0,
            args=args,
            bounds=bounds,
            jac=True,
            hessp=hessp,
            hess=hess,
            tol=1e-4,
            options=opts,
            method=method,
        )
    y0 = history[0][0]

    min_step = history[0]
    for step in history:
        if step[0] < min_step[0]:
            min_step = step

    y1, g1, x1, h1, _ = min_step

    out = [line for step in history for line in step[4]]
    
    return returns.success((x1, y0, y1, g1), out=out, err=[])


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
