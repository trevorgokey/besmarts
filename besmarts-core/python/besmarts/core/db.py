"""
besmarts.core.db
"""
import dbm
import glob
import os

from typing import Sequence, Dict, List

from besmarts.core import arrays
from besmarts.core import codecs
from besmarts.core import compute

class db_dict:
    def __init__(self, icd, name=""):
        self.icd = icd
        self.name = name
        if name:
            self.kv = None
            assert db_intvec_create(name)
        else:
            self.kv = {}

    def write_intvec(self, kv, prefix=""):
        if self.name:
            return db_intvec_write(self.name, kv, prefix=prefix)
        else:
            self.kv.update(kv)
            return len(kv)

    def read_intvec(self, keys, prefix=""):
        if self.name:
            return db_intvec_read(self.name, kv, prefix=prefix)
        else:
            return self.kv[k]

    def read_intvec_list(self, keys, prefix=""):
        if self.name:
            return db_intvec_read_list(self.name, kv, prefix=prefix)
        else:
            return [self.kv[k] for k in keys]

    def delete_intvec(self, keys, prefix=""):
        if self.name:
            db_intvec_delete(self.name, keys, prefix=prefix)
        else:
            for k in keys:
                del self.kv[k]
        
    def write_subgraph(self, kv, prefix=""):
        if self.name:
            return db_intvec_write(
                self.name,
                {k: self.icd.subgraph_encode(v) for k,v in kv.items()},
                prefix=prefix
            )
        else:
            self.kv.update({k: self.icd.subgraph_encode(v) for k,v in kv.items()})
            return len(kv)

    def write_structure(self, kv, prefix=""):
        if self.name:
            return db_structure_write(
                self.icd,
                self.name,
                kv,
                prefix=prefix
            )
        else:
            self.kv.update({k: self.icd.structure_encode(v) for k,v in kv.items()})
            return len(kv)

    def write_graph(self, kv, prefix=""):
        if self.name:
            return db_intvec_write(
                self.name,
                {k: self.icd.graph_encode(v) for k,v in kv.items()},
                prefix=prefix
            )
        else:
            self.kv.update({k: self.icd.graph_encode(v) for k,v in kv.items()})
            return len(kv)

    def read_graph(self, key, prefix=""):
        if self.name:
            return self.icd.graph_decode(db_intvec_read(
                self.name,
                key,
                prefix=prefix
            ))
        else:
            return self.icd.graph_decode(self.kv[k])

    def read_graph_list(self, keys, prefix=""):
        if self.name:
            return [
                self.icd.graph_decode(x) for x in db_intvec_read_list(
                    self.name,
                    keys,
                    prefix=prefix
                )
            ]
        else:
            return [self.icd.graph_decode(self.kv[x]) for x in keys]

    def read_structure_list(self, keys, prefix=""):
        if self.name:
            return [
                self.icd.structure_decode(x) for x in db_intvec_read_list(
                    self.name,
                    keys,
                    prefix=prefix
                )
            ]
        else:
            return [self.icd.structure_decode(self.kv[x]) for x in keys]

    def read_structure(self, key, prefix=""):
        if self.name:
            return self.icd.structure_decode(db_intvec_read(self.name, key, prefix=prefix))
        else:
            return self.icd.structure_decode(self.kv[key])

    def remove(self):
        for fn in glob.glob(self.name+"*"):
            os.remove(fn)

def db_intvec_create(db_name) -> bool:
    try:
        with dbm.open(db_name, 'c') as db:
            return True
    except Exception:
        return False

def db_graph_write(icd: codecs.intvec_codec, db_name, pairs, prefix=""):

    if prefix:
        prefix = prefix + ":"

    with open(dbm.open(db_name), 'w') as db:
        for k,v in pairs.items():
            db[prefix+str(k)] = icd.graph_encode(v).tobytes()

def db_structure_write_distributed(pairs, shm=None):
    prefix = shm.prefix
    if prefix:
        prefix = prefix + ":"
    for k,v in pairs:
        shm.db[prefix+str(k)] = shm.icd.graph_encode(v).tobytes()


def db_structure_write(icd: codecs.intvec_codec, db_name, pairs, prefix=""):

    if prefix:
        prefix = prefix + ":"

    with open(dbm.open(db_name), 'wf') as db:
        wq = compute.workqueue_local('', 0)
        ws = compute.workqueue_new_workspace(
            wq,
            address=('127.0.0.1', 0),
            shm={"icd": icd, "db": db, "prefix": prefix}
        )
        
        compute.workspace_submit_and_flush(db_structure_write, arrays.batched(pairs, 10000), chunksize=10)
        for k,v in pairs.items():
            db[prefix+str(k)] = icd.graph_encode(v).tobytes()

        ws.close()
        wq.close()
        db.sync()

    return len(pairs)

def db_intvec_write(db_name, pairs, prefix=""):

    if prefix:
        prefix = prefix + ":"

    with dbm.open(db_name, 'wf') as db:
        for k,v in pairs.items():
            db[prefix+str(k)] = v.v.tobytes()
        db.sync()
    return len(pairs)

def db_intvec_delete(db_name, keys, prefix=""):

    if prefix:
        prefix = prefix + ":"

    with dbm.open(db_name, 'wf') as db:
        for k,v in pairs.items():
            del db[prefix+str(k)]
        db.sync()

## readers

def db_intvec_read(db_name, keys, prefix=""):

    if prefix:
        prefix = prefix + ":"

    with dbm.open(db_name, 'rfu') as db:
        v = arrays.intvec()
        v.v.frombytes(db[prefix+str(k)])
        return v

def db_intvec_read_list(db_name, keys, prefix=""):

    if prefix:
        prefix = prefix + ":"

    vals = []
    with dbm.open(db_name, 'rfu') as db:
        for k in keys:
            v = arrays.intvec()
            v.v.frombytes(db[prefix+str(k)])
            vals.append(v)

    return vals

def db_graph_read(icd: codecs.intvec_codec, db_name, keys, prefix=""):

    if prefix:
        prefix = prefix + ":"

    keyvals = {}
    with dbm.open(db_name, 'rfu') as db:
        for k in keys:
            v = arrays.intvec()
            v.v.frombytes(db[prefix+str(k)])
            keyvals[k] = icd.graph_decode(v)

    return keyvals


