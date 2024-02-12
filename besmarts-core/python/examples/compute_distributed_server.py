

# from multiprocessing import Process
from threading import Thread
import time

from besmarts.core import compute
from besmarts.core import codecs
from besmarts.core import graphs
from besmarts.core import mapper
from besmarts.codecs import codec_rdkit
from besmarts.core.arrays import batched

def start_workqueue():

    # starts the server to host remote connections
    wq = compute.workqueue_local('', 55555)
    # t = Process(target=m.mgr.get_server().serve_forever)
    # ctx = multiprocessing.get_context('fork')
    
    # t = Thread(target=hello_printer, args=(wq,))
    # t.start()
    # m.threads['server'] = t
    # t.start()
    # time.sleep(2)
    # t.join()
    return wq

def construct(smi, shm=None):
    gcd = codec_rdkit.graph_codec_rdkit()
    icd = codecs.intvec_codec(gcd.primitive_codecs, gcd.atom_primitives, gcd.bond_primitives)
    g = gcd.smiles_decode(smi)
    ic = graphs.graph_to_structure_angles(g)
    if ic:
        u = mapper.union_list(ic)
        # i = icd.structure_encode(u)
        return u

def construct_idx(idx, shm=None):
    gcd = codec_rdkit.graph_codec_rdkit()
    icd = codecs.intvec_codec(gcd.primitive_codecs, gcd.atom_primitives, gcd.bond_primitives)

    smi = shm.smiles[idx]

    g = gcd.smiles_decode(smi)
    ic = graphs.graph_to_structure_angles(g)
    if ic:
        u = mapper.union_list(ic)
        # i = icd.structure_encode(u)
        return u

def hello_submit(ws, units, chunksize=10):
    n = len(units)
    i = 0
    for chunk in batched(units.items(), chunksize):
        # for i in range(20):
            # compute.workspace_submit(ws, print, (f"{i:05d} hello!",), {})
        
        i += len(chunk)
        # chunk = {idx: (construct, (smi,), {}) for idx, smi in chunk}
        chunk = {idx: (construct_idx, (smi,), {}) for idx, smi in chunk}
        # construct(smi)
        # compute.workspace_local_submit(ws, i, print, (f"{i:05d} hello!",), {})
        # 
        compute.workspace_local_submit(ws, chunk)
        # compute.workspace_local_submit(ws, idx, construct, (smi,), {})
        # print(f"Submitting: {i/n*100:5.2f}%", end='\r')
            # compute.workspace_submit(ws, time.sleep, (5,), {})
    # print()
    # ws.submitted.set()


def hello_printer(wq):
    # m = compute.manager('127.0.0.1', 55555)
    # gcd = graph_codec_rdkit()
    smiles = open('./smiles.dat').read().split()

    print("hello_printer start")
    if not smiles[-1]:
        smiles.pop()

    shm = compute.shm_construct(smiles)

    ws = compute.workqueue_new_workspace(wq, shm)
    # ws.start()
    print("Workspace started")
    # smiles = smiles[:100]
    # smiles = list(smiles * 10)
    # compute.workspace_submit(ws, "END", (), {})
    # smiles = smiles[:100]
    smiles_list = list(sorted(smiles, key=lambda x: len(x), reverse=False))
    smiles = dict(zip(range(1,1+len(smiles_list)), range(len(smiles_list))))
    # smiles = dict(zip(range(1,1+len(smiles_list)), smiles_list))
    for i, smi_i in smiles.items():
        smi = smiles_list[smi_i]
        # smi = smi_i
        print(f"{i:04d} {smi}")
    # smiles = dict(zip(range(1,1+len(smiles_list)), smiles_list))

    results = {}
    smiles_to_submit = smiles
    # t = None
    t = compute.workspace_local_run(ws)
    ws.reset()
    while len(results) < len(smiles):
        # print("Resetting workspace")
        # ws.submitted.clear()
        # print(f"\n** SUBMITTING {len(smiles_to_submit)} ({len(results)}/{len(smiles)} done) **")
        # if len(smiles_to_submit) < 100:
        #     print(list(smiles_to_submit))
        # ti = Thread(target=hello_submit, args=(ws, smiles_to_submit))
        # ti.start()
        # # time.sleep(1)
        # print("Joining submitter")
        # ti.join()
        hello_submit(ws, smiles_to_submit)


        these_results = compute.workspace_flush(ws, set(smiles_to_submit), timeout=10.)
        results.update(these_results)
        # print(f"\nReceived {len(these_results)} new results. Total: {len(results)/len(smiles)*100:5.2f}%")
        # print(f"{list(these_results)}")
        # else:
        #     time.sleep(10.0)
        # t = Thread(target=compute.workspace_join, args=(ws, len(smiles)))


        # print("Setting to done")
        # ws.done.set()
        # print("Joining runner")
        # t.join()
        smiles_to_submit = {i: j for i, j in smiles.items() if i not in results}
        # time.sleep(10)
    for idx, result in sorted(results.items(), key=lambda x: x[0]):
        print(idx, result)

    ws.done.set()
    # print("Joining runner")
    # t.join()

    # q = ws.get_oqueue()
    # while not q.empty():
    #     distfun = q.get()
    #     print(f"Got result {distfun}")
    #     q = ws.get_oqueue()
    print("Done. Closing")
    ws.close()
    

if __name__ == "__main__":
    wq = start_workqueue()
    hello_printer(wq)
    
# m.threads['server'].join()
# submit_work(m)

# while m.threads:
#     for addr, ws in list(m.threads.items()):
#         if ws.mgr.state().get('status', 0) == 0:
#             print("DONE!")
#             m.threads.pop(addr)
#     time.sleep(1)
# t.join()
