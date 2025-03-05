
"""
besmarts.mechanics.optimizers_scipy
"""

import copy
import scipy.optimize
import numpy as np
import datetime
import math
import functools

from besmarts.core import geometry
from besmarts.core import assignments
from besmarts.core import compute
from besmarts.core import arrays
from besmarts.core import logs
from besmarts.core import configs
from besmarts.core import returns
from besmarts.mechanics import molecular_models as mm
from besmarts.mechanics import objectives
from besmarts.mechanics import fits

PRECISION = configs.precision

def optimize_positions_scipy(
    csys,
    psys: mm.physical_system,
    step_limit=1000,
    tol=1e-10
):

    pos = copy.deepcopy(psys.models[0].positions)
    args, keys = objectives.array_flatten_matrix_assignment(pos)

    # jac = objectives.array_geom_gradient
    # jac = None
    jac = True

    hess = objectives.array_geom_hessian

    method = 'L-BFGS-B'
    opts = {
        'disp': False,
        'ftol': tol,
        'gtol': tol,
        'maxiter': step_limit,
        'maxls': 1000,
        'maxcor': len(keys)**2
    }
    hess = None

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

    for (pi, c, n, i), v in zip(keys, arrays.array_round(result.x, PRECISION)):
        pos[pi].selections[n][c][i] = v
    
    # print("Final pos", pos.selections)
    return pos


def ff_obj(x, keys, csys, refpsys, ref, refcsys):

    pos = copy.deepcopy(ref)
    recalc = set([x[0] for x in keys])

    for k, v in zip(keys, x):
        mm.chemical_system_set_value(csys, k, v)

    reuse = set(range(len(csys.models))).intersection(recalc)
    psys = mm.chemical_system_to_physical_system(
        csys,
        pos,
        ref=refpsys,
        reuse=reuse
    )
    optpos = optimize_positions_scipy(csys, psys)

    obj = 0
    for ic, rxyz in ref[0].selections.items():
        oxyz = optpos.selections[ic]
        for a, b in zip(rxyz, oxyz):
            dx = geometry.array_distance(a, b)
            dx2 = dx*dx
            obj += dx2
    return obj


def objective_hessp(x, p, *args):
    history = args[6]
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


def collect(gdb, obj, full_results, new_results, batch_map, hp, verbose=False):
    ready = {}
    for (n, (ii, jlist)), work in new_results.items():

        if (n, (ii, 0)) not in full_results:
            full_results[(n, (ii, 0))] = [*([None]*(len(batch_map[ii])))]

        if jlist == (0,):
            full_results[(n, (ii, 0))][0] = work
        else:
            result_idx = batch_map[ii].index(tuple(jlist))
            full_results[(n, (ii, 0))][result_idx] = work

        if all([_ is not None for _ in full_results[(n, (ii, 0))]]):
            full_result = [
                x for y in full_results.pop((n, (ii, 0))) for x in y
            ]
            x = obj[ii]
            ref = assignments.graph_db_get(gdb, x.addr)
            ready[("obj", n, ii)] = [
                (ii, x, ref, full_result, hp), {"verbose": verbose}
            ]
    return ready, full_results


def process(finished, X, Y, grad, hess, grady, verbose=False):
    out = []
    for (_, _, i), ret in finished.items():
        (x2, gradxi, hessxi, y2, gradyi) = ret.value
        if verbose:
            for line in ret.out:
                print(line)
        X[i] = x2
        Y[i] = y2
        grad[i] = gradxi
        hess[i] = hessxi
        grady[i] = gradyi
        out.extend(ret.out)

    return X, Y, grad, hess, grady, out


def objective_gradient_gdb(
    args,
    keys,
    csys,
    gdb,
    objbatches,
    priors,
    penalties=None,
    history=None,
    psysref=None,
    reuse=None,
    ws=None,
    verbose=False,
    minstep=10**(-PRECISION),
    return_gradient=True
):

    # if ws:
    #     ws.close()
    # ws = None

    if history is None:
        history = []
    if penalties is None:
        penalties = []
    n = len(history)

    X = 0
    Y = 0
    Xi = {}
    Yi = {}
    gradi = {}
    hessi = {}
    gradyi = {}

    X_t = {}
    Y_t = {}
    grad_t = {}
    grady_t = {}
    hessi_t = {}

    grad = list([0.0] * len(keys))
    hess = list([0.0] * len(keys))
    grady = list([0.0] * len(keys))

    out = []
    args = arrays.array_round(args, PRECISION)

    h = []

    dX = 0.0
    beststep = None
    t = arrays.array_magnitude(args)
    dt0 = 0
    dt = 0
    dargs = arrays.array_scale(args, 0)
    dargs0 = arrays.array_scale(args, 0)

    if len(history):

        beststep = arrays.argmin([x[0] for x in history])
        dargs0 = arrays.array_difference(args, history[beststep][2])
        dt0 = arrays.array_magnitude(dargs0)

        dargs = arrays.array_difference(args, history[-1][2])
        dt = arrays.array_magnitude(dargs)

        if dt < minstep:
            if dt > 0:
                out.append(f"Step size too small: {dt:.15e} ({minstep:.15e}). Skipping evaluation")
                if verbose:
                    print(out[-1])
            X, grad, args, hess, out, X_t = history[-1]
            history.append(history[-1])

            if return_gradient:
                return X, grad
            else:
                return X

    # big job, try to start computing while tasks are being generated
    # also, reap objective as it comes due to memory consumption
    # currently doesn't work (fixed?)
    # async_compute = True and (ws is not None) and len(args)*len(objlst) > 1000
    # async_compute = True and (ws is not None)
    async_compute = False

    verbose = verbose and not async_compute

    dcsys = csys
    if ws:
        dcsys = None
        ws.reset()

    z = [v*p[1] + p[0] for v, p in zip(args, priors)]
    z = arrays.array_round(z, PRECISION)
    hi = 1e-2
    h = tuple(([hi] * len(z)))

    hp = tuple()
    if return_gradient:
        hp = tuple(hi/j[1] for j in priors)

    n_finished = 0

    chunksize = configs.compute_runtime.get("task_chunksize", 100)

    # targetbatch = 1000*len(keys)
    # if not return_gradient:
    #     targetbatch *= 3
    # objbatches = []
    # cur_batch_score = 0
    # cur_batch = []

    # total_cost = 0
    # for i, obj in objlst.items():
    #     bsz = obj.batch_size
    #     if bsz is None:
    #         total_cost += 1
    #     else:
    #         total_cost += len(keys) // bsz + bool(len(keys) % bsz)

    # total_cost = 0

    # for i, obj in objlst.items():
    #     bsz = obj.batch_size
    #     if bsz is None:
    #         cost = len(keys)
    #     else:
    #         cost = len(keys) // bsz + bool(len(keys) % bsz)
    #     total_cost += cost

    # targetbatch = max(1, total_cost//configs.compute_runtime.get("task_batches", 1))

    # for i, obj in objlst.items():
    #     bsz = obj.batch_size
    #     if bsz is None:
    #         cost = len(keys)
    #     else:
    #         cost = len(keys) // bsz + bool(len(keys) % bsz)

    #     if (cur_batch_score + cost <= targetbatch or len(cur_batch) == 0):
    #         cur_batch.append((i, obj))
    #         cur_batch_score += cost
    #     else:
    #         objbatches.append(cur_batch)
    #         cur_batch = []
    #         cur_batch_score = 0

    # if cur_batch:
    #     objbatches.append(cur_batch)
    #     # cur_batch = []
    #     # cur_batch_score = 0

    # objbatchsize = 500 if len(objlst) > 500 else len(objlst)

    n_objbatches = len(objbatches)

    # n_objbatches = len(objlst)//objbatchsize + bool(len(objlst)%objbatchsize)
    # for obsz, obj in enumerate(arrays.batched(objlst.items(), objbatchsize), 1):
    batches_cum_sum = 0
    for obsz, obj in enumerate(objbatches, 1):
        tasks = {}
        full_results = {}
        batch_map = {}
        chunk = {}
        batches_cum_sum += len(obj)
        # for idx, (i, x) in enumerate(obj, 1 + objbatchsize*(obsz-1)):
        line = f"{logs.timestamp()} Objective batch {obsz}/{n_objbatches}"
        logs.append(line, out, verbose)

        objlst = dict(obj)

        for idx, (i, x) in enumerate(obj, batches_cum_sum - len(obj) + 1):

            # print(idx, i, x)
            # line = f"{logs.timestamp()} PSystem {idx}/{len(obj)}"
            # logs.append(line, out, verbose)

            eids = x.addr.eid
            psys = {}
            reapply = set()
            if psysref:
                for eid in eids:
                    psys[eid] = psysref[eid]

                    # kv = {}
                    # for k, v, p in zip(keys, args, priors):
                    #     v = p[1]*v + p[0]
                    #     kv[k] = v
                    #     # print(f"Reassigning {k}")
                    #     # print(f"Setting pval to {k}={v}")
                    #     mm.physical_system_set_value(psys[eid], k, v)
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
                    #                     k[1:]: v
                    #                     for k, v in kv.items() if k[0] == m
                    #                 }
                    #             )
            # else:
                # line = (
                    # "WARNING: No parameterized system given. This will "
                    # "recharge the molecules and is likely not intended."
                # )
                # logs.append(line, out, verbose)
                # csys = copy.deepcopy(csys)
                # for k, v in zip(keys, args):
                    # v = p[1]*v + p[0]
                    # mm.chemical_system_set_value(csys, k, v)
                # psys = fits.gdb_to_physical_systems(gdb, csys)
                # reuse = list(range(len(csys.models)))

            dreuse = [x for x in list(range(len(csys.models))) if x not in reapply]
            grad_keys = [{}]
            task = x.get_task(
                gdb,
                dcsys,
                keys=grad_keys,
                psys=psys,
                reuse=reuse
            )
            tasks[(n, (i, (0,)))] = task

            if async_compute:
                chunk[(n, (i, (0,)))] = task

            batch_map[i] = [(0,)]

            if return_gradient:

                # line = f"{logs.timestamp()} PSystem grads"
                # logs.append(line, out, verbose)

                grad_keys = []
                if x.batch_size is None:
                    batch_size = len(keys)
                else:
                    batch_size = x.batch_size

                batches = [
                    tuple(x)
                    for x in arrays.batched(range(1, 1+len(keys)), batch_size)
                ]
                batch_map[i].extend(batches)

                for batch in arrays.batched(enumerate(zip(keys, z, h), 1), batch_size):

                    kbatch = tuple([b[0] for b in batch])
                    for j, (k, v, hi) in batch:
                        if x.grad_mode == "c2":
                            grad_keys.extend([{k: v-hi}, {k:v+hi}])
                        elif x.grad_mode == "f1":
                            grad_keys.append({k:v+hi})

                    task = x.get_task(
                        gdb,
                        dcsys,
                        keys=grad_keys,
                        psys=psys,
                        reuse=dreuse
                    )
                    tasks[(n, (i, kbatch))] = task
                    if async_compute:
                        chunk[(n, (i, kbatch))] = task
                        if len(chunk) >= chunksize:
                            compute.workspace_local_submit(
                                ws,
                                {key: (fits.objective_run_distributed, [t], {})
                                 for key, t in chunk.items()}
                            )
                            chunk = {}

                        task_keys = list(tasks)

                        new_results = compute.workspace_flush(
                            ws,
                            set(task_keys),
                            0.0,
                            maxwait=0.0,
                            verbose=False
                        )
                        for key in new_results:
                            if key in tasks:
                                tasks.pop(key)
                                n_finished += 1

                        ready, full_results = collect(
                            gdb,
                            objlst,
                            full_results,
                            new_results,
                            batch_map,
                            hp,
                            verbose=verbose
                        )
                        new_results.clear()

                        obj_results = compute.workspace_submit_and_flush(
                            ws,
                            run_objective,
                            ready,
                            chunksize=chunksize,
                            verbose=False,
                            batchsize=chunksize*1000,
                            timeout=0.0,
                            clear=False
                        )
                        ready.clear()

                        Xi, Yi, gradi, hessi, gradyi, retouti = process(
                            obj_results,
                            Xi,
                            Yi,
                            gradi,
                            hessi,
                            gradyi,
                            verbose=verbose
                        )
                        obj_results.clear()
                        out.extend(retouti)

                    grad_keys = []
            if verbose and async_compute:
                print(
                    f"\r{logs.timestamp()}",
                    f"{obsz:10d}/{n_objbatches:<10d}",
                    f"{idx:10d}/{len(objlst)}",
                    f"Async: {async_compute}",
                    f"Complete: {n_finished:10d}/{len(tasks):<10d}",
                    end=''
                )

        if async_compute and len(chunk):
            compute.workspace_local_submit(
                ws,
                {key: (fits.objective_run_distributed, [t], {})
                 for key, t in chunk.items()}
            )
            chunk.clear()

        if verbose:
            print()

        line = f"{logs.timestamp()} Running {len(tasks)} tasks"
        logs.append(line, out, verbose)

        if ws:
            if async_compute:
                new_results = compute.workspace_flush(
                    ws,
                    set(tasks),
                    0.0,
                    verbose=True
                )
                for key in new_results:
                    if key in tasks:
                        tasks.pop(key)

                ready, full_results = collect(
                    gdb,
                    objlst,
                    full_results,
                    new_results,
                    batch_map,
                    hp,
                    verbose=verbose
                )
                new_results.clear()
                obj_results = compute.workspace_submit_and_flush(
                    ws,
                    run_objective,
                    ready,
                    chunksize=chunksize,
                    verbose=False,
                    batchsize=chunksize*1000,
                    timeout=0.0,
                    clear=False
                )
                ready.clear()

                # X, Y, grad, hess, grady, retout = process(
                #     obj_results,
                #     Xi,
                #     Yi,
                #     grad,
                #     hess,
                #     grady,
                #     verbose=verbose
                # )
                Xi, Yi, gradi, hessi, gradyi, retout = process(
                    obj_results,
                    Xi,
                    Yi,
                    gradi,
                    hessi,
                    gradyi,
                    verbose=verbose
                )
                for (i, x) in obj:
                    if i in Xi:
                        t = type(x)
                        v = Xi[i]
                        X_t[t] = X_t.get(t, 0) + v
                        grad_t[t] = arrays.array_add(grad_t.get(t, [0.0]*len(keys)), gradi[i])
                        hessi_t[t] = arrays.array_add(hessi_t.get(t, [0.0]*len(keys)), hessi[i])
                        Xi.pop(i)
                    elif i in Yi:
                        t = type(x)
                        v = Yi[i]
                        Y_t[t] = Y_t.get(t, 0) + v
                        grady_t[t] = arrays.array_add(grady_t.get(t, [0.0]*len(keys)), gradyi[i])
                        Yi.pop(i)
                # for xi, v in Xi.items():
                #     t = type(obj[xi][1])
                #     X_t[t] = X_t.get(t, 0) + v
                #     grad_t[t] = arrays.array_add(grad_t.get(t, [0.0]*len(keys)), gradi[xi])
                #     hessi_t[t] = arrays.array_add(hessi_t.get(t, [0.0]*len(keys)), hessi[xi])
                # for xi, v in Yi.items():
                #     t = type(obj[xi][1])
                #     grady_t[t] = arrays.array_add(grady_t.get(t, [0.0]*len(keys)), gradyi[xi])
                #     Y_t[t] = Y_t.get(t, 0) + v
                obj_results.clear()
                out.extend(retout)

            if tasks:

                line = f"{logs.timestamp()} Physical prop compute"
                logs.append(line, out, verbose)

                new_results = compute.workspace_submit_and_flush(
                    ws,
                    fits.objective_run_distributed,
                    {i: ([x], {}) for i, x in tasks.items()},
                    chunksize=chunksize,
                    verbose=verbose,
                    batchsize=chunksize*1000,
                    timeout=0.1, #max(0, math.log(len(tasks), 1000)-1),
                    clear=False
                )
                tasks.clear()

                line = f"{logs.timestamp()} Collecting"
                logs.append(line, out, verbose)

                ready, full_results = collect(
                    gdb,
                    objlst,
                    full_results,
                    new_results,
                    batch_map,
                    hp,
                    verbose=verbose
                )
                new_results.clear()

                line = f"{logs.timestamp()} Objective compute"
                logs.append(line, out, verbose)

                obj_results = compute.workspace_submit_and_flush(
                    ws,
                    run_objective,
                    ready,
                    chunksize=chunksize,
                    verbose=verbose,
                    batchsize=chunksize*1000,
                    timeout=0.0
                )
                ready.clear()

                line = f"{logs.timestamp()} Objective collect"
                logs.append(line, out, verbose)

                Xi, Yi, gradi, hessi, gradyi, retout = process(
                    obj_results,
                    Xi,
                    Yi,
                    gradi,
                    hessi,
                    gradyi,
                    verbose=verbose
                )
                for (i, x) in obj:
                    if i in Xi:
                        t = type(x)
                        v = Xi[i]
                        X_t[t] = X_t.get(t, 0) + v
                        grad_t[t] = arrays.array_add(grad_t.get(t, [0.0]*len(keys)), gradi[i])
                        hessi_t[t] = arrays.array_add(hessi_t.get(t, [0.0]*len(keys)), hessi[i])
                        Xi.pop(i)
                    elif i in Yi:
                        t = type(x)
                        v = Yi[i]
                        Y_t[t] = Y_t.get(t, 0) + v
                        grady_t[t] = arrays.array_add(grady_t.get(t, [0.0]*len(keys)), gradyi[i])
                        Yi.pop(i)
                obj_results.clear()
                out.extend(retout)

                assert not full_results

        else:
            for i, task in tasks.items():

                new_results= {i: task.run()}

                ready, full_results = collect(
                    gdb,
                    objlst,
                    full_results,
                    new_results,
                    batch_map,
                    hp,
                    verbose=verbose
                )
                new_results.clear()
                for oit, obj_task in ready.items():
                    ret = run_objective(*obj_task[0], **obj_task[1])
                    Xi, Yi, gradi, hessi, gradyi, retouti = process(
                        {oit: ret},
                        Xi,
                        Yi,
                        gradi,
                        hessi,
                        gradyi,
                        verbose=verbose
                    )
                    out.extend(retouti)
                    for (i, x) in obj:
                        if i in Xi:
                            t = type(x)
                            v = Xi[i]
                            X_t[t] = X_t.get(t, 0) + v
                            grad_t[t] = arrays.array_add(grad_t.get(t, [0.0]*len(keys)), gradi[i])
                            hessi_t[t] = arrays.array_add(hessi_t.get(t, [0.0]*len(keys)), hessi[i])
                            Xi.pop(i)
                        elif i in Yi:
                            t = type(x)
                            v = Yi[i]
                            Y_t[t] = Y_t.get(t, 0) + v
                            grady_t[t] = arrays.array_add(grady_t.get(t, [0.0]*len(keys)), gradyi[i])
                            Yi.pop(i)
                    ret = None
                ready.clear()

            full_results.clear()
            tasks.clear()

    line = f"{logs.timestamp()} Restraints"
    logs.append(line, out, verbose)

    kv = {k: v * p[1] + p[0] for k, v, p in zip(keys,args,priors)}

    for i, pen in enumerate(penalties, len(objlst)):
        pen: fits.objective_config_penalty

        # generate the reference target values and then run to get deltas
        dx = pen.get_task().run(kv)
        ret = pen.compute_diff(dx, verbose=verbose)
        out.extend(ret.out)

        if verbose:
            for line in ret.out:
                print(line)
        x2 = pen.scale*ret.value
        t = type(pen)
        X_t[t] = X_t.get(t, 0) + x2

        gradx = list([0.0] * len(keys))
        hessx = list([0.0] * len(keys))
        if return_gradient:
            for j, k in enumerate(keys, 0):
                if k not in dx:
                    continue
                kvi = {k: dx[k]}
                dx2dp = pen.compute_gradient(kvi)
                dx2dp *= pen.scale
                dx2dp = round(dx2dp, 12)
                # print("dx2dp", dx2dp)
                #grad[j] += dx2dp
                gradx[j] += dx2dp

                d2x2dp2 = pen.compute_hessian(kvi)
                d2x2dp2 *= pen.scale**2
                d2x2dp2 = round(d2x2dp2, 12)
                # print("dxdp", dxdp)
                hessx[j] += d2x2dp2
        grad_t[t] = arrays.array_add(grad_t.get(t, list([0.0] * len(keys))), gradx)
        hessi_t[t] = arrays.array_add(hessi_t.get(t, list([0.0] * len(keys))), hessx)

        gnormx = arrays.array_inner_product(gradx, gradx)**.5
        line = (
            f">> RID={i:02d} S= {pen.scale: 14.6f} "
            f"R2= {x2: 14.6f} "
            f"|g|= {gnormx: 14.6f}"
        )
        if verbose:
            print(line)
        out.append(line)

    line = f"{logs.timestamp()} Finalizing"
    logs.append(line, out, verbose)

    X = sum(X_t.values())
    if hessi_t:
        hess = functools.reduce(arrays.array_add, hessi_t.values())
        hnorm = arrays.array_inner_product(hess, hess)**.5
    else:
        hnorm = 0.0

    if grad_t:
        grad = functools.reduce(arrays.array_add, grad_t.values())
        gnorm = arrays.array_inner_product(grad, grad)**.5
    else:
        gnorm = 0.0

    Y = sum(Y_t.values())
    if grady_t:
        grady = functools.reduce(arrays.array_add, grady_t.values())
        gnormy = arrays.array_inner_product(grady, grady)**.5
    else:
        gnormy = 0.0

    t = arrays.array_magnitude(args)
    if len(history):
        beststep = arrays.argmin([x[0] for x in history])
        dX = X - history[beststep][0]
    out.append(
        f">>> {datetime.datetime.now()} Totals "
        f"Step= {len(history)+1:3d} "
        f"|t| = {t:10.4e} "
        f"|dt0| = {dt0:10.4e} "
        f"|dt| = {dt:10.4e} "
        f"X2= {X: 12.5e} "
        f"|g|= {gnorm: 12.5e} "
        f"|h|= {hnorm: 12.5e} "
        f"DX2= {dX: 12.5e} "
        f"Y2= {Y: 12.5e} "
        f"|gy|={gnormy: 12.5e}"
    )
    if verbose or async_compute:
        print(out[-1])
    for t in set(list(X_t) + list(Y_t)):
        gradi = grad_t.get(t, [0])
        gnormi = arrays.array_inner_product(gradi, gradi)**.5

        hessi = hessi_t.get(t, [0])
        hnormi = arrays.array_inner_product(hessi, hessi)**.5

        gradyi = grady_t.get(t, [0])
        gnormyi = arrays.array_inner_product(gradyi, gradyi)**.5
        
        dX = 0.0
        if beststep is not None:
            dX = X_t.get(t, 0) - history[beststep][5].get(t, 0)
        out.append(
            f"==> {str(t)} "
            f"Step= {len(history)+1:3d} "
            f"X2= {X_t.get(t, 0): 12.5e} "
            f"|g|= {gnormi: 12.5e} "
            f"|h|= {hnormi: 12.5e} "
            f"DX2= {dX: 12.5e} "
            f"Y2= {Y_t.get(t, 0): 12.5e} "
            f"|gy|={gnormyi: 12.5e}"
        )
        if verbose or async_compute:
            print(out[-1])
    out.append(f"{logs.timestamp()} Done.")
    if verbose:
        print(out[-1])
    out.append("Total Parameter Step and Grad:")
    if verbose:
        print(out[-1])
    for ii, (k, g, t, dti0, dti) in enumerate(zip(keys, grad, args, dargs0, dargs)):
        line = f"{ii:4d} {str(k):20s} t= {t:11.4e} dt0= {dti0:11.4e} dt= {dti:11.4e} g= {g:11.4e}"
        out.append(line)
        if verbose:
            print(line)

    if history and dt < 1e-6:
        X, grad, args, hess, out, X_t = history[-1]
        history.append(history[-1])
    else:
        history.append((X, grad, args, hess, out, X_t))
    if return_gradient:
        return X, grad
    else:
        return X


def run_objective(oid, x, ref, results, h, verbose=False, shm=None):
    X = 0
    Y = 0


    N = len(h)
    gradx = list([0.0] * N)
    hessx = list([0.0] * N)
    grady = list([0.0] * N)

    # if it is 1 then we are not doing gradients
    ret = x.compute_diff(ref, results[0], verbose=verbose)
    dx = ret.value
    if x.weights is not None:
        dx = arrays.array_multiply(dx, x.weights)
    x2 = x.scale*arrays.array_inner_product(dx, dx)
    if len(results) > 1:
        for j in range(N):
            if x.grad_mode == "c2":
                dxa = results[2*j+1]
                dxb = results[2*j+2]
                # print("DXA", dxa[0][0].graphs[0].rows[0].columns[0].selections)
                # print("DXB", dxb[0][0].graphs[0].rows[0].columns[0].selections)

                dxdp, d2xdp2 = x.compute_gradient(ref, [dxa, results[0], dxb], h[j])
            elif x.grad_mode == "f1":
                # dxa = results[0]
                dxb = results[j+1]
                # print("DXA", dxa[0][0].graphs[0].rows[0].columns[0].selections)
                # print("DXB", dxb[0][0].graphs[0].rows[0].columns[0].selections)

                dxdp, d2xdp2 = x.compute_gradient(ref, [results[0], dxb], h[j])

            if x.weights is not None:
                dxdp = arrays.array_multiply(dxdp, x.weights)
                d2xdp2 = arrays.array_multiply(d2xdp2, x.weights)
            dx2dp = 2 * arrays.array_inner_product(dx, dxdp)
            d2x2dp2 = 2 * (
                arrays.array_inner_product(dx, d2xdp2) + arrays.array_inner_product(dxdp, dxdp)
            )
            dx2dp *= x.scale
            d2x2dp2 *= x.scale**2
            # ret.out.append("dx2dp: " + str(dx2dp))

            if x.include:
                gradx[j] += dx2dp
                hessx[j] += d2x2dp2
            else:
                grady[j] += dx2dp
        if x.verbose > 2:
            ret.out.append(
                f"OID {oid:05d} Gradient Decomposition scale: {x.scale} (sorted ascending magnitude):"
            )
            for j, gx in enumerate(gradx):
                ret.out.append(
                    f"GD Param: {j:4d} "
                    f"d2dpx: {gx:14.6e} "
                    f"d2dpx2: {hessx[j]:14.6e}"
                )

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

    # addr = x.addr.eid[0]
    output = (
        f"\n>> OID={oid:05d} S= {x.scale: 14.6f} " +
        f"{sym}= {x2: 14.6e} " +
        f"|g|= {gnorm: 14.6e} " +
        f"|h|= {hnorm:14.6e}"
    )
    ret.out.append(output)
    return returns.success((X, gradx, hessx, Y, grady), out=ret.out)


def singlepoint_forcefield_gdb_scipy(
    args,
    keys,
    csys,
    gdb,
    obj,
    priors,
    penalties=None,
    history=None,
    psysref=None,
    reuse=None,
    ws=None,
    verbose=False,
    minstep=10**(-PRECISION),
):

    if history is None:
        history = []

    if penalties is None:
        penalties = []

    args = arrays.array_round(args, PRECISION)

    out = []
    if len(history):

        beststep = arrays.argmin([x[0] for x in history])

        dargs = arrays.array_difference(args, history[-1][2])
        dt = arrays.array_magnitude(dargs)

        if dt < minstep:
            X, grad, args, hess, out, X_t = history[-1]
            line = f"Step size too small: {dt:.15e} ({minstep:.15e}). Skipping evaluation"
            if verbose and line not in out:
                print(line)
            # history.append((X, grad, args, hess, out, X_t))
            return X
    # csys = copy.deepcopy(csys)
    # if reuse is None:
    #     reuse = set(range(len(csys.models)))
    for k, v, p in zip(keys, args, priors):
        v0 = mm.chemical_system_get_value(csys, k)
        # if k[0] in reuse:
        #     reuse.remove(k[0])
        v1 = p[0] + v*p[1]
        mm.chemical_system_set_value(csys, k, v)

        out.append(f"Setting {str(k):20s} from {v0:15.7g} to {v1:15.7g} d={v1-v0:15.7g}")
        if verbose:
            # dv = v-v0
            # if abs(dv) < 1e-6:
            #     dv = 0.0
            # mm.chemical_system_set_value(csys, k, v0)
            print(out[-1])

    # reuse = list(reuse)
    X = objective_gradient_gdb(args, keys, csys, gdb, obj, priors, penalties=penalties, history=history, psysref=psysref, reuse=reuse, ws=ws, verbose=verbose, return_gradient=False)
    # print(f"RETURN IS {X}")
    return X


def fit_grad_gdb(
    args,
    keys,
    csys,
    gdb,
    obj,
    priors,
    penalties=None,
    history=None,
    psysref=None,
    reuse=None,
    ws=None,
    verbose=False,
    minstep=10**(-PRECISION),
):

    if history is None:
        history = []

    if penalties is None:
        penalties = []

    args = arrays.array_round(args, PRECISION)
    out = []

    if len(history):

        dargs = arrays.array_difference(args, history[-1][2])
        dt = arrays.array_magnitude(dargs)

        if dt < minstep:
            X, grad, args, hess, out, X_t = history[-1]
            line = f"Step size too small: {dt:.15e} ({minstep:.15e}). Skipping evaluation"
            if verbose and line not in out:
                print(line)
            history.append((X, grad, args, hess, out, X_t))
            return X, grad

    for i, (k, v, p) in enumerate(zip(keys, args, priors)):
        v0 = mm.chemical_system_get_value(csys, k)
        v1 = p[0] + v*p[1]

        # if abs(v1-v0) < 1e-6:
        #     args[i] = (v0 - p[0])/p[1]
        #     v1 = v0

        out.append(
            f"Setting {str(k):20s} "
            f"from {v0:15.7g} to {v1:15.7g} d={v1-v0:15.7g}"
        )
        if verbose:
            print(out[-1])

        mm.chemical_system_set_value(csys, k, v1)

    X, g = objective_gradient_gdb(
        args,
        keys,
        csys,
        gdb,
        obj,
        priors,
        penalties=penalties,
        history=history,
        psysref=psysref,
        reuse=reuse,
        ws=ws,
        verbose=verbose,
        minstep=minstep
    )

    # return X
    return X, g


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


def callback(xk):
    d = arrays.array_magnitude(xk)

    print(f"CALLBACK: distance is {d}")
    if d > 25:
        raise StopIteration

def optimize_forcefield_gdb_scipy(
    x0,
    args,
    bounds=None,
    step_limit=None,
    maxls=20,
    ftol=1e-5,
    gtol=1e-5,
    anneal=False
):

    (
        keys,
        csys,
        gdb,
        obj,
        priors,
        penalties,
        history,
        psysref,
        reuse,
        ws,
        verbose,
        minstep
    ) = args

    out = []

    method = obj.method

    hessp = None
    hess = None

    configtab = {
        "L-BFGS-B": {
            "options": {
                'maxls': maxls,
                #'maxcor': len(x0)**2,
                'iprint': 0 if verbose else 0,
                'ftol': ftol,
                'gtol': gtol
            },
            "hessp": None,
            "hess": None,
        },
        "TNC": {
            "options": {
                'xtol': 1e-7,
                'disp': False,
                'maxCGit': maxls,
                'ftol': ftol,
                'gtol': gtol,
            },
            "hessp": None,
            "hess": None,
        },
        "CG": {
            "options": {
                'gtol': 1e-5,
                'disp': verbose,
            },
            "hessp": objective_hessp,
            "hess": None,
        },
        "trust-krylov": {
            "options": {
                'gtol': 1e-5,
                'disp': verbose,
                'inexact': True,
            },
            "hessp": objective_hessp,
            "hess": None,
        },
        "trust-ncg": {
            "options": {
                'gtol': 1e-5,
                'disp': verbose,
                'inexact': True,
            },
            "hessp": objective_hessp,
            "hess": None,
        },
        "trust-exact": {
            "options": {
                'gtol': 1e-5,
                'disp': verbose,
                'inexact': True,
            },
            "hessp": None,
            "hess": objective_hess,
        }
    }
    config = configtab[method]
    config.update(obj.method_config)
    opts = config['options']
    hess = config['hess']
    hessp = config['hessp']

    for xi, x in obj.objectives.items():
        if method == "L-BFGS-B":
            x.grad_mode = "f1"
        else:
            x.grad_mode = "c2"

    objlst = obj.objectives
    chunksize = configs.compute_runtime.get("task_chunksize", 100)

    objbatches = []
    cur_batch_score = 0
    cur_batch = []
    total_cost = 0
    line = (
        f"{logs.timestamp()} Generating {len(objlst)} objectives. "
    )
    out.append(line)
    if verbose:
        print(line)

    for i, obj in objlst.items():
        bsz = obj.batch_size
        if bsz is None:
            cost = len(keys)
        else:
            cost = len(keys) // bsz + bool(len(keys) % bsz)
        total_cost += cost

    targetbatch = max(1, total_cost//configs.compute_runtime.get("task_batches", 1))

    for i, obj in objlst.items():
        bsz = obj.batch_size
        if bsz is None:
            cost = len(keys)
        else:
            cost = len(keys) // bsz + bool(len(keys) % bsz)

        if (cur_batch_score + cost <= targetbatch or len(cur_batch) == 0):
            cur_batch.append((i, obj))
            cur_batch_score += cost
        else:
            objbatches.append(cur_batch)
            cur_batch = []
            cur_batch_score = 0

    if cur_batch:
        objbatches.append(cur_batch)

    if verbose:

        line = (
            f"{logs.timestamp()} Processing objectives in "
            f"{len(objbatches)} batches (target cost {targetbatch} total {total_cost})"
        )
        print(line)

    args = (
        keys,
        csys,
        gdb,
        objbatches,
        priors,
        penalties,
        history,
        psysref,
        reuse,
        ws,
        verbose,
        minstep
    )

    # if 0:
    #     method = 'L-BFGS-B'
    #     opts = {
    #         'maxls': maxls,
    #         'maxcor': len(args)**2,
    #         'iprint': 1000 if verbose else 0,
    #         'ftol': 1e-7,
    #         'gtol': 1e-7
    #     }
    #     # hessp = None

    # # slow
    # elif 0:
    #     method = 'TNC'
    #     opts = {
    #         'xtol': 1e-6,
    #         'disp': False,
    #         'maxCGit': maxls,
    #         # 'iprint': 1000 if verbose else 0,
    #         'xtol': 1e-3,
    #         'ftol': 1e-4,
    #         'gtol': 1e-3,

    #     }
    # elif 0:
    #     method = 'CG'
    #     opts = {
    #         'gtol': 1e-7,
    #         'disp': verbose,
    #         # 'iprint': 1000 if verbose else 0,
    #         # 'tol': 1e-4,
    #         # 'gtol': 1e-3,

    #     }
    #     # bounds=None
    #     hessp = objective_hessp
    #     hess = None

    # elif 0:
    #     # decent first step but seems to hit local minima
    #     method = 'TNC'
    #     opts = {
    #         'xtol': 1e-7,
    #         'disp': True,
    #         'maxCGit': maxls,
    #         # 'iprint': 1000 if verbose else 0,
    #         'xtol': 1e-7,
    #         'ftol': 1e-7,
    #         'gtol': 1e-7,
    #         'scale': np.full(len(x0), 1.0),
    #         'offset': np.full(len(x0), 0.0),
    #         'maxfun': step_limit,
    #         'eta': .50
    #     }
    # elif 0:
    #     # slow
    #     method = 'trust-krylov'
    #     opts = {
    #         'disp': verbose,
    #         'gtol': 1e-7,
    #         'inexact': True,
    #     }
    #     hessp = objective_hessp
    #     hess = None
    # elif 0:
    #     # slow
    #     method = 'trust-exact'
    #     opts = {
    #         'disp': verbose,
    #         'gtol': 1e-7,

    #     }
    #     hess = objective_hess
    #     hessp = None
    # elif 1:
    #     # more hessp but less iterations than newton-cg
    #     method = 'trust-ncg'
    #     opts = {
    #         'disp': verbose,
    #         'gtol': 1e-7,
    #     }
    #     hessp = objective_hessp
    #     hess = None
    # elif 0:
    #     method = 'Newton-CG'
    #     opts = {
    #         'xtol': 1e-7,
    #         'disp': verbose,
    #         # 'iprint': 1000 if verbose else 0,
    #         # 'tol': 1e-4,
    #         # 'gtol': 1e-3,

    #     }
    #     # bounds=None
    #     hessp = objective_hessp
    #     hess = None

    if step_limit:
        opts['maxiter'] = int(step_limit)

    if verbose:
        print(datetime.datetime.now(), "Starting physical parameter optimization")

    if anneal is True:
        anneal = 100
    else:
        anneal = int(anneal)
    if anneal and step_limit != 0:
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
            "tol": 1e-5,
            "method": method,
            "options": opts,
            "bounds": bounds,
            "args": args,
        }
        scipy.optimize.basinhopping(
            fit_grad_gdb,
            x0,
            disp=True,
            niter=anneal,
            stepsize=.01,
            T=len(args[3].objectives) * 10,
            minimizer_kwargs=kwds,
        )

    elif step_limit == 0:
        singlepoint_forcefield_gdb_scipy(x0, *args)
    else:
        scipy.optimize.minimize(
            fit_grad_gdb,
            x0,
            args=args,
            bounds=bounds,
            jac=True,
            hessp=hessp,
            hess=hess,
            tol=1e-5,
            options=opts,
            method=method,
        )
    y0 = history[0][0]

    min_step = history[0]
    for step in history:
        if step[0] < min_step[0]:
            min_step = step

    y1, g1, x1, h1, _, _ = min_step

    out = [line for step in history for line in step[4]]

    return returns.success((x1, y0, y1, g1), out=out, err=[])

