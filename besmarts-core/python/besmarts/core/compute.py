"""
besmarts.core.compute

Responsible for setting up and distributing large compute jobs

Architecture-wise, it interfaces the multiprocessing module

I need something that has a remote queue that pushes work to outside, but also
a more local version that allows faster local versions, but hopefully these can
be the same interface

The main process shown span some threads to handle the server thread and then
the main process will throw work at it. I should have specialized functions for
certain tasks, which can be handled here

I think the way to go is to have each function call a manager instance which registers its
required methods. We have a connectionmanager which will launch new managers,
one for each function, and so the main process will control the conman
"""


import traceback
from typing import Tuple, Callable, Sequence, Mapping, Dict
import itertools
import io
import ctypes
import struct
import random
import time
import os
import sys
import pprint
import array
import functools
import logging
from datetime import datetime

# networking
import queue
import socket

# process management
import multiprocessing
from multiprocessing import (
    util,
    managers,
    sharedctypes,
    queues,
    context,
    connection,
    process,
)
import multiprocessing.pool

import threading
import pickle

from besmarts.core import configs
from besmarts.core import graphs
from besmarts.core import arrays

distributed_function = Tuple[Callable, Sequence, Mapping]

# search/replace this if using a single file
# remote_compute_enable = configs.remote_compute_enable

TIMEOUT = 20

SHM_GLOBAL = None

## from https://stackoverflow.com/questions/34361035/python-thread-name-doesnt-show-up-on-ps-or-htop

LIB = "libcap.so.2"
try:
    libcap = ctypes.CDLL(LIB)
except OSError:
    print("Library {} not found. Unable to set thread name.".format(LIB))
# else:
#     def _name_hack(self):
#         # PR_SET_NAME = 15
#         # try:
#         #     libcap.prctl(15, "thread".encode())
#         # except Exception as e:
#         #     print(f"Warning, exception {e}")
#         threading.Thread._bootstrap_original(self)

#     threading.Thread._bootstrap_original = threading.Thread._bootstrap
#     threading.Thread._bootstrap = _name_hack

###


def thread_name_set(name):
    libcap.prctl(15, name.encode())


class workspace_status:
    INVALID = -1
    EMPTY = 0
    INACTIVE = 1
    SUBMITTING = 2
    WAITING = 3
    RUNNING = 4
    DONE = 5


def dprint(*args, **kwds):
    if False:
        print(*args, **kwds)


def dispatch(c, id, methodname, args=(), kwds={}):
    """
    Send a message to manager using connection `c` and return response
    """
    if methodname in ["incref", "decref"]:
        return
    t0 = time.perf_counter_ns()
    c.send((id, methodname, args, kwds))
    dprint(f"SENDING time: {(time.perf_counter_ns() - t0)*1e-9}")
    dprint(f"RECEIVING...")
    kind, result = c.recv()
    dprint(f"TOTAL time: {(time.perf_counter_ns() - t0)*1e-9}")
    if kind == "#RETURN":
        dprint(f"RETURNING...")
        return result
    raise managers.convert_to_error(kind, result)


managers.dispatch = dispatch


class Connection(connection.Connection):
    """
    Connection class based on an arbitrary file descriptor (Unix only), or
    a socket handle (Windows).
    """

    _write = connection.Connection._write
    _read = connection.Connection._read

    def _send(self, buf, write=_write):
        remaining = len(buf)
        while True:
            n = write(self._handle, buf)
            remaining -= n
            if remaining == 0:
                break
            buf = buf[n:]

    def _recv(self, size, read=_read):
        buf = io.BytesIO()
        handle = self._handle
        remaining = size
        while remaining > 0:
            chunk = read(handle, remaining)
            n = len(chunk)
            if n == 0:
                if remaining == size:
                    raise EOFError
                else:
                    raise OSError("got end of file during message")
            buf.write(chunk)
            remaining -= n
        return buf

    def _send_bytes(self, buf):
        n = len(buf)
        # For wire compatibility with 3.2 and lower
        # header = struct.pack("!i", n.to_bytes)

        # The main reason we subclass connections: sometimes shm is large
        # and so we send an insanely large buf size. Default is a 4 byte signed
        # int which is about 2GB. This should remove the problem
        header = n.to_bytes(8, byteorder="big", signed=False)
        # print("SHIP IT:", n,  header)
        if n > 16384:
            # The payload is large so Nagle's algorithm won't be triggered
            # and we'd better avoid the cost of concatenation.
            self._send(header)
            self._send(buf)
        else:
            # Issue #20540: concatenate before sending, to avoid delays due
            # to Nagle's algorithm on a TCP socket.
            # Also note we want to avoid sending a 0-length buffer separately,
            # to avoid "broken pipe" errors if the other end closed the pipe.
            self._send(header + buf)

    def _recv_bytes(self, maxsize=None):
        # buf = bytes(self._recv(8).getbuffer())
        # buf = [self._recv(4).getvalue() for x in range(128//4)]
        buf = self._recv(8).getvalue()
        size = int.from_bytes(buf, byteorder="big", signed=False)
        # size, = struct.unpack("!i", buf.getvalue())

        if maxsize is not None and size > maxsize:
            return None
        return self._recv(size)


connection.Connection = Connection


def SocketClient(address, timeout=TIMEOUT):
    """
    Return a connection object connected to the socket given by `address`
    """
    family = connection.address_type(address)
    with socket.socket(getattr(socket, family)) as s:
        dprint(f"Connecting socket to {address}")
        try:
            # s.settimeout(timeout)
            s.setblocking(True)
            t0 = time.perf_counter()
            s.connect(address)
            dprint(
                f"Socket connect to {address} time: {time.perf_counter() - t0}"
            )
        except Exception as e:
            # print(f"Failed to connect: {e}")
            raise e
        # print(f"Connected")
        sd = s.detach()
        # print(f"Detached")
        # traceback.print_stack()
        return Connection(sd)


# disable passwords
# connection.answer_challenge = lambda x, y: True
# connection.deliver_challenge = lambda x, y: True


def Client(address, family=None, authkey=None, timeout=TIMEOUT):
    """
    Returns a connection to the address of a `Listener`
    """
    t0 = time.perf_counter()
    family = family or connection.address_type(address)
    connection._validate_family(family)
    # print("Client TB:")
    # traceback.print_stack()
    if family == "AF_PIPE":
        c = connection.PipeClient(address)
    else:
        c = SocketClient(address, timeout=timeout)
    if authkey is not None and not isinstance(authkey, bytes):
        raise TypeError("authkey should be a byte string")

    if authkey is not None:
        connection.answer_challenge(c, authkey)
        connection.deliver_challenge(c, authkey)
    dprint(f"STARTING Client init time: {time.perf_counter() - t0}")

    return c


class Listener(connection.Listener):
    def __init__(self, *args, **kwargs):
        kwargs["backlog"] = 64000
        super().__init__(*args, **kwargs)
        dprint("INIT LISTENER WITH LARGE BACKLOG")

    def accept(self):
        """
        Accept a connection on the bound socket or named pipe of `self`.

        Returns a `Connection` object.
        """
        # t0 = time.perf_counter()
        if self._listener is None:
            raise OSError("listener is closed")
        c = self._listener.accept()
        if self._authkey:
            deliver_challenge(c, self._authkey)
            answer_challenge(c, self._authkey)
        # t = time.perf_counter()
        # dprint(f"SERVER ACCEPT TIME: {t - t0:.6f}")
        return c


class Server(managers.Server):
    def __init__(self, *args, **kwargs):
        dprint("HELLO FROM SERVER")
        super().__init__(*args, **kwargs)
        # self.shm_pool = multiprocessing.Pool(32, server_shm_init, )
        self.shm_pool = None
        self.shm_msg_cache = None
        self.shm_msg_cache_lock = threading.Lock()

    def serve_client(self, conn):
        """
        Handle requests from the proxies in a particular process/thread
        """
        util.debug(
            "starting server thread to service %r",
            threading.current_thread().name,
        )

        recv = conn.recv
        send = conn.send
        id_to_obj = self.id_to_obj

        while not self.stop_event.is_set():
            name = "none"
            methodname = "none"
            try:
                methodname = obj = None
                thread_name_set("recv")
                request = recv()
                ident, methodname, args, kwds = request

                try:
                    obj, exposed, gettypeid = id_to_obj[ident]
                except KeyError as ke:
                    try:
                        obj, exposed, gettypeid = self.id_to_local_proxy_obj[
                            ident
                        ]
                    except KeyError:
                        raise ke

                if methodname not in exposed:
                    raise AttributeError(
                        "method %r of %r object is not in exposed=%r"
                        % (methodname, type(obj), exposed)
                    )

                function = getattr(obj, methodname)
                name = type(obj).__name__ + "_"

                try:
                    dprint(f"T{name} RUNNING ON OBJ {function} {args} {kwds}")
                    thread_name_set(name + methodname)
                    t0 = time.perf_counter()
                    # these are always readonly so we can fork and skip GIL
                    if False and name.startswith("shm"):
                        res = self.shm_pool.apply(function, args, kwds)
                    else:
                        res = function(*args, **kwds)
                    t = time.perf_counter()
                    dprint(
                        f"T{name} RUNNING ON OBJ DONE t={t-t0:.6f} {function} {args} {kwds}"
                    )
                except Exception as e:
                    msg = ("#ERROR", e)
                else:
                    typeid = gettypeid and gettypeid.get(methodname, None)
                    if typeid:
                        rident, rexposed = self.create(conn, typeid, res)
                        token = managers.Token(typeid, self.address, rident)
                        msg = ("#PROXY", (rexposed, token))
                    else:
                        msg = ("#RETURN", res)

            except AttributeError:
                if methodname is None:
                    msg = ("#TRACEBACK", format_exc())
                else:
                    try:
                        fallback_func = self.fallback_mapping[methodname]
                        result = fallback_func(
                            self, conn, ident, obj, *args, **kwds
                        )
                        msg = ("#RETURN", result)
                    except Exception:
                        msg = ("#TRACEBACK", traceback.format_exc())

            except EOFError:
                util.debug(
                    "got EOF -- exiting thread serving %r",
                    threading.current_thread().name,
                )
                sys.exit(0)

            except Exception:
                msg = ("#TRACEBACK", traceback.format_exc())

            try:
                try:
                    if (
                        name.startswith("shm")
                        and methodname == "get"
                        and msg[0] == "#RETURN"
                    ):
                        with self.shm_msg_cache_lock:
                            if self.shm_msg_cache is None:
                                self.shm_msg_cache = (
                                    connection._ForkingPickler.dumps(msg)
                                )
                        conn.send_bytes(self.shm_msg_cache)
                    else:
                        send(msg)
                except Exception:
                    send(("#UNSERIALIZABLE", traceback.format_exc()))
            except Exception as e:
                util.info(
                    "exception in thread serving %r",
                    threading.current_thread().name,
                )
                util.info(" ... message was %r", msg)
                util.info(" ... exception was %r", e)
                conn.close()
                sys.exit(1)
            break

    def accept_connection(self, c, name):
        """
        Spawn a new thread to serve this connection
        """
        threading.current_thread().name = name
        c.send(("#RETURN", None))
        self.serve_client(c)

    def serve_forever(self):
        """
        Run the server forever
        """
        self.stop_event = threading.Event()
        process.current_process()._manager_server = self
        try:
            # self.accepter()
            accepter = threading.Thread(target=self.accepter, name="listener")
            accepter.daemon = True
            accepter.start()
            try:
                while not self.stop_event.is_set():
                    self.stop_event.wait(1)
            except (KeyboardInterrupt, SystemExit):
                pass
        finally:
            if sys.stdout != sys.__stdout__:  # what about stderr?
                util.debug("resetting stdout, stderr")
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
            sys.exit(0)

    def accepter(self):
        dprint("ACCEPTING CONNECTIONS")
        # pool = ThreadPool(32)
        # with ThreadPool(32) as pool:
        thread_name_set("listener")
        recvs = {}
        sends = {}
        procs = {}
        idx = 0
        while True:
            # threads = [x for x in threads if x.is_alive()]
            if True:
                try:
                    # print("LISTENING")
                    # t0 = time.perf_counter()
                    c = self.listener.accept()
                    # t = time.perf_counter()
                    # print(f"CONNECTION RECV time {t-t0:.6f}")
                except OSError as e:
                    # except Exception as e:
                    # print(f"acceptor exception {type(e)} {e}")
                    # print(
                    #     f"Warning, there are {len(threads)} active connections, "
                    #     f"waiting 10 seconds before accepting new connections."
                    # )
                    # time.sleep(10.0)
                    continue
                # if len(threads) > 600:
                #     print(f"there are {len(threads)} active connections, waiting 2")
                #     time.sleep(2)
                # if len(threads) > 700:
                #     print(f"there are {len(threads)} active connections, waiting 3")
                #     time.sleep(3)
                # if len(threads) > 800:
                #     print(f"there are {len(threads)} active connections, waiting 4")
                #     time.sleep(4)
                # if len(threads) > 900:
                #     print(f"there are {len(threads)} active connections, waiting 5")
                #     time.sleep(5)
                # self.handle_request(c)
                t = threading.Thread(target=self.handle_request, args=(c,))
                # pool.apply_async(self.handle_request, (c,))
                t.daemon = True
                t.start()

            else:
                if len(recvs) < 32:
                    new_out = []
                    new_err = []
                    new_conn = []
                    new_t = threading.Thread(
                        target=self.recv_request,
                        args=(new_conn, new_out, new_err),
                        name=f"recv_{idx}",
                    )
                    new_t.daemon = True
                    new_t.start()
                    recvs[idx] = (new_conn, new_t, new_out, new_err)
                    idx += 1

                for i, (c, tr, out, err) in list(recvs.items()):
                    if tr.is_alive():
                        continue

                    if not c:
                        tr.join()
                        del recvs[i]
                        continue

                    if err and not out:
                        t = threading.Thread(
                            target=self.send_request,
                            args=(c[0], err[0]),
                            name=f"senderr_{i}",
                        )
                        t.daemon = True
                        t.start()
                        sends[i] = (c, t, None, None)
                    elif out:
                        msg = []
                        fn = out[0][0], out[0][1]

                        new_pt = threading.Thread(
                            target=self.run_request,
                            args=(c[0], out[0], msg),
                            name=f"run_{i}",
                        )
                        new_pt.daemon = True

                        if fn not in procs:
                            procs[fn] = []

                        if len(procs[fn]) == 0:
                            # print(f"STARTED RUN THREAD {i}")
                            new_pt.start()
                        else:
                            new_pt.start()

                        # print(f"APPEND RUN THREAD {i}")
                        procs[fn].append((i, c, new_pt, msg))

                        # t = threading.Thread(target=self.send_request, args=(c[0], out, []))
                        # t.daemon = True
                        # t.start()
                    tr.join()
                    del recvs[i]

                for fn in procs:
                    if len(procs[fn]) > 0:
                        j, pc, pt, pmsg = procs[fn][0]
                        # print(f"CHECKING RUN THREAD {j}, this fn has len={len(procs[fn])}")
                        if not pt.is_alive():
                            pt.join()
                            procs[fn].pop(0)
                            if pmsg:
                                st = threading.Thread(
                                    target=self.send_request,
                                    args=(pc[0], pmsg[0]),
                                    name=f"send_{j}",
                                )
                                st.daemon = True
                                st.start()
                                sends[j] = (pc, st, None, None)

                            # if len(procs[fn]) > 0:
                            #     procs[fn][0][2].start()

                for i, (c, t, _, _) in list(sends.items()):
                    if not t.is_alive():
                        t.join()
                        del sends[i]

            # t = threading.Thread(target=self.handle_request, args=(c, out, err))
            # pool.apply_async(self.handle_request, (c,))
            # t.daemon = True
            # print("TRYING TO JOIN....")
            # t.start()
            # threads.append(t)
            # t.join(timeout=None)
            # print("JOINED")

    def recv_request(self, conn, out, err):
        try:
            c = self.listener.accept()
        except OSError:
            return
        conn.append(c)
        request = None
        # print("HANDLING REQUEST")
        t0 = time.perf_counter()
        try:
            # print("HANDLING CHALLENGE")
            connection.deliver_challenge(c, self.authkey)
            connection.answer_challenge(c, self.authkey)
            # print("HANDLING RECV")
            request = c.recv()
            ignore, funcname, args, kwds = request
            assert funcname in self.public, "%r unrecognized" % funcname
            func = getattr(self, funcname)
            out.append((func, args, kwds))
            # print("RECV DONE")
        except Exception:
            msg = ("#TRACEBACK", traceback.format_exc())
            # tf = time.perf_counter()
            err.append(msg)

    def run_request(self, c, request, msg):
        func, args, kwds = request
        try:
            name = threading.current_thread().name
            # print(f"T{name} RUNNING {func} {args} {kwds}")
            t0 = time.perf_counter()
            result = func(c, *args, **kwds)
            t = time.perf_counter()
            # print(
            #     f"T{name} RUNNING DONE t={t-t0:.6f} {func} {args} {kwds}"
            # )
        except Exception:
            msg.append(("#TRACEBACK", traceback.format_exc()))
        else:
            msg.append(("#RETURN", result))

    def send_request(self, c, msg):
        # msg = msg[0]

        try:
            c.send(msg)
            c.close()
        except Exception as e:
            try:
                c.send(("#TRACEBACK", traceback.format_exc()))
                c.close()
            except Exception:
                pass
            util.info("Failure to send message: %r", msg)
            # util.info(" ... request was %r", request)
            util.info(" ... exception was %r", e)

    def handle_request(self, conn):
        """
        Handle a new connection
        """
        try:
            self._handle_request(conn)
        except SystemExit:
            # Server.serve_client() calls sys.exit(0) on EOF
            pass
        finally:
            conn.close()

    def _handle_request(self, c):
        thread_name_set("request")
        request = None
        dprint("HANDLING REQUEST")
        t0 = time.perf_counter()
        try:
            dprint("HANDLING CHALLENGE")
            connection.deliver_challenge(c, self.authkey)
            connection.answer_challenge(c, self.authkey)
            dprint("HANDLING RECV")
            request = c.recv()
            ignore, funcname, args, kwds = request
            assert funcname in self.public, "%r unrecognized" % funcname
            func = getattr(self, funcname)
        except Exception:
            msg = ("#TRACEBACK", traceback.format_exc())
            tf = time.perf_counter()
        else:
            tf = time.perf_counter()
            try:
                name = threading.current_thread().name
                dprint(f"T{name} RUNNING {func} {args} {kwds}")
                t0 = time.perf_counter()
                result = func(c, *args, **kwds)
                t = time.perf_counter()
                dprint(
                    f"T{name} RUNNING DONE t={t-t0:.6f} {func} {args} {kwds}"
                )
            except Exception:
                msg = ("#TRACEBACK", traceback.format_exc())
            else:
                msg = ("#RETURN", result)

        t1 = time.perf_counter()

        try:
            c.send(msg)
        except Exception as e:
            try:
                c.send(("#TRACEBACK", traceback.format_exc()))
            except Exception:
                pass
            util.info("Failure to send message: %r", msg)
            util.info(" ... request was %r", request)
            util.info(" ... exception was %r", e)
        t2 = time.perf_counter()
        dprint(
            f"REQUEST TOTAL {t2-t0:.6f} RECV {tf-t0:.6f} PROC {t1 - tf:.6f} SEND {t2-t1:.6f}"
        )

    def shutdown(self, c):
        if self.shm_pool:
            self.shm_pool.terminate()
        super().shutdown(c)

    def incref(self, c, ident):
        with self.mutex:
            try:
                self.id_to_refcount[ident] += 1
            except KeyError as ke:
                # If no external references exist but an internal (to the
                # manager) still does and a new external reference is created
                # from it, restore the manager's tracking of it from the
                # previously stashed internal ref.
                if ident in self.id_to_local_proxy_obj:
                    self.id_to_refcount[ident] = 1
                    self.id_to_obj[ident] = self.id_to_local_proxy_obj[ident]
                    obj, exposed, gettypeid = self.id_to_obj[ident]
                    util.debug("Server re-enabled tracking & INCREF %r", ident)
                else:
                    raise ke

    def decref(self, c, ident):
        if (
            ident not in self.id_to_refcount
            and ident in self.id_to_local_proxy_obj
        ):
            util.debug("Server DECREF skipping %r", ident)
            return

        with self.mutex:
            if self.id_to_refcount[ident] <= 0:
                raise AssertionError(
                    "Id {0!s} ({1!r}) has refcount {2:n}, not 1+".format(
                        ident,
                        self.id_to_obj[ident],
                        self.id_to_refcount[ident],
                    )
                )
            self.id_to_refcount[ident] -= 1
            if self.id_to_refcount[ident] == 0:
                del self.id_to_refcount[ident]

        if ident not in self.id_to_refcount:
            # Two-step process in case the object turns out to contain other
            # proxy objects (e.g. a managed list of managed lists).
            # Otherwise, deleting self.id_to_obj[ident] would trigger the
            # deleting of the stored value (another managed object) which would
            # in turn attempt to acquire the mutex that is already held here.
            self.id_to_obj[ident] = (None, (), None)  # thread-safe
            util.debug("disposing of obj with id %r", ident)
            with self.mutex:
                del self.id_to_obj[ident]


class BaseProxy(managers.BaseProxy):
    def __init__(self, *args, **kwds):
        kwds["incref"] = False
        super().__init__(*args, **kwds)

    def _after_fork(self):
        return


def MakeProxyType(name, exposed, _cache={}):
    """
    Return a proxy type whose methods are given by `exposed`
    """
    exposed = tuple(exposed)
    try:
        return _cache[(name, exposed)]
    except KeyError:
        pass

    dic = {}

    for meth in exposed:
        exec(
            """def %s(self, *args, **kwds):
        return self._callmethod(%r, args, kwds)"""
            % (meth, meth),
            dic,
        )

    ProxyType = type(name, (BaseProxy,), dic)
    ProxyType._exposed_ = exposed
    _cache[(name, exposed)] = ProxyType
    return ProxyType


managers.MakeProxyType = MakeProxyType
managers.listener_client["pickle"] = (Listener, Client)


def AutoProxy(
    token,
    serializer,
    manager=None,
    authkey=None,
    exposed=None,
    incref=True,
    manager_owned=False,
):
    """
    Return an auto-proxy for `token`
    """
    _Client = Client

    if exposed is None:
        conn = _Client(token.address, authkey=authkey)
        try:
            exposed = dispatch(conn, None, "get_methods", (token,))
        finally:
            conn.close()

    if authkey is None and manager is not None:
        authkey = manager._authkey
    if authkey is None:
        authkey = process.current_process().authkey

    ProxyType = managers.MakeProxyType("AutoProxy[%s]" % token.typeid, exposed)
    proxy = ProxyType(
        token,
        serializer,
        manager=manager,
        authkey=authkey,
        incref=incref,
        manager_owned=manager_owned,
    )
    proxy._isauto = True
    return proxy


managers.AutoProxy = AutoProxy

# util.get_logger()
# util.log_to_stderr(level=logging.DEBUG)

# managers.BaseManager._Server = Server


def _repopulate_pool(self):
    """Bring the number of pool processes up to the specified number,
    for use after reaping workers which have exited.
    """
    for i in range(self._processes - len(self._pool)):
        w = self.Process(
            target=multiprocessing.pool.worker,
            args=(
                self._inqueue,
                self._outqueue,
                self._initializer,
                self._initargs,
                self._maxtasksperchild,
                self._wrap_exception,
            ),
        )
        self._pool.append(w)
        w.name = w.name.replace("Process", "PoolWorker")
        w.daemon = False
        w.start()
        # util.debug('added worker')


multiprocessing.pool.Pool._repopulate_pool = _repopulate_pool


def Process(obj, *args, **kwds):
    if type(obj) is multiprocessing.pool.Pool:
        return obj._ctx.Process(*args, **kwds)
    else:
        return obj.Process(*args, **kwds)


multiprocessing.pool.Pool.Process = Process


class workspace_manager(managers.SyncManager):
    _Server = Server

    def get_workspaces(self) -> Dict:
        pass

    def get_state(self) -> Dict:
        pass

    def get_status(self) -> workspace_status:
        pass

    def get_table(self) -> Dict:
        pass

    def get_result(self) -> Dict:
        pass

    def create(self, *args, **kwargs):
        breakpoint()
        return super().create(*args, **kwargs)


class workqueue_manager(managers.SyncManager):
    def get_iqueue(self) -> queue.Queue:
        pass

    def get_rqueue(self) -> queue.Queue:
        pass

    def get_oqueue(self) -> queue.Queue:
        pass

    def get_state(self) -> Dict:
        pass


class workqueue:
    def __init__(self, addr, port):
        self.mgr = workqueue_manager(address=(addr, port), authkey=b"0")
        self.mgr._Client = Client

    def get_workspaces(self):
        pass

    def get_state(self):
        pass

    def get_status(self):
        return self.get_state().get("status", workspace_status.INVALID)


class workqueue_local(workqueue):
    def __init__(self, addr, port):
        super().__init__(addr, port)
        self.threads = {}
        self.workspaces = {}
        # print(f"Workspace local addr is {id(self.workspaces)}")
        self.state = dict({"status": workspace_status.INACTIVE})
        # self.status = multiprocessing.Value('d', -1)

        self.mgr._Client = functools.partial(Client, timeout=TIMEOUT)

        self.mgr.register("get_workspaces", self.get_workspaces)
        self.mgr.register("get_state", self.get_state)
        self.mgr.register("get_status", self.get_status)

        self.mgr.start()
        self.remote_workspaces = self.mgr.get_workspaces()
        # self.mgr.get_server().serve_forever()
        # self.mgr.workspaces = self.mgr.dict()
        # print(f"Workspace remote addr is {id(self.workspaces)}")

    def get_threads(self):
        return self.threads

    def get_workspaces(self):
        # print(f"GET LOCAL Workspace addr is {id(self.workspaces)} values {self.workspaces}")
        # return dict(self.mgr.get_workspaces().items())
        return self.workspaces

    def put_workspaces(self, wss):
        # workspaces = self.mgr.get_workspaces()
        self.remote_workspaces.update(wss)
        self.workspaces.update(wss)
        # print(f"PUT LOCAL Workspace addr is {id(wss)} values {dict(wss.items())}")
        return True

    def remove_workspace(self, ws):
        pass

    def close(self):
        # self.pool.close()
        # self.pool.terminate()
        self.mgr.shutdown()

    # def get_state(self):
    #     return self.state
    def get_state(self):
        return self.state


def manager_connect(mgr: managers.BaseManager, success: threading.Event):
    success.clear()
    try:
        print("Manager connect thread hello...")
        mgr.connect()
        # print("Connected...")
        success.set()
        # print("Manager connect setting success...")
    # except ConnectionError:
    #     print("manager_connect: ConnectionError")
    except EOFError:
        print("manager_connect: EOFError")
    except TimeoutError:
        print("manager_connect: TimeoutError")
    except AssertionError as e:
        print(f"manager_connect: AssertionError {e}")
    except Exception as e:
        print(f"manager_connect: {type(e)} {e}")


def remote_connect(mgr, timeout=TIMEOUT):
    success = threading.Event()

    print("Attempting to connect to manager...")
    t = threading.Thread(
        target=manager_connect,
        args=(
            mgr,
            success,
        ),
    )
    t.start()

    connected = True

    t.join(timeout=timeout)
    print("Waited on thread")
    if t.is_alive() or not success.is_set():
        connected = False

    print(f"Connection success is {connected}")
    return connected


def manager_remote_get_status_thread(ws, out):
    thread_name_set("get_status")
    try:
        # print("manager_remote_get_state_thread: calling get_state")
        # state = mgr.get_state()
        status = ws.remote_state.get("status", workspace_status.INVALID)
        # print("manager_remote_get_state_thread: calling get_state items")
        # state = dict(state.items())
        out.append(status)
        print(f"manager_remote_get_status_thread: success {out}")
        # success.set()
    except ConnectionError:
        print("manager_remote_get_status_thread: ConnectionError")
    except EOFError:
        print("manager_remote_get_status_thread: EOFError")
    except TimeoutError:
        print("manager_remote_get_status_thread: TimeoutError")
    except AssertionError as e:
        print(f"manager_remote_get_status_thread: AssertionError {e}")
    except queue.Empty:
        # print("queue.Empty")
        pass


def manager_remote_get_state_thread(ws, out):
    thread_name_set("get_state")
    try:
        # print("manager_remote_get_state_thread: calling get_state")
        # state = mgr.get_state()
        state = ws.mgr.get_state()
        print("manager_remote_get_state_thread: calling get_state items")
        # state = dict(state.items())
        out.append(state)
        print(f"manager_remote_get_state_thread: success {out}")
        # success.set()
    except ConnectionError:
        print("manager_remote_get_state_thread: ConnectionError")
    except EOFError:
        print("manager_remote_get_state_thread: EOFError")
    except TimeoutError:
        print("manager_remote_get_state_thread: TimeoutError")
    except AssertionError as e:
        print(f"manager_remote_get_state_thread: AssertionError {e}")
    except queue.Empty:
        # print("queue.Empty")
        pass


def manager_remote_get_status(wq, timeout=5.0):
    out = []
    state = {}

    t = threading.Thread(
        target=manager_remote_get_status_thread,
        args=(wq, out),
    )
    print("manager_remote_get_status: starting thread")
    t.start()

    t.join(timeout=timeout)
    if len(out):
        return out[0]
    else:
        return workspace_status.INVALID


def manager_remote_get_state(wq, timeout=TIMEOUT):
    out = []
    state = {}

    t = threading.Thread(
        target=manager_remote_get_state_thread,
        args=(wq, out),
    )
    print("manager_remote_get_state: starting thread")
    t.start()

    t.join(timeout=timeout)
    if len(out):
        return out[0]
    else:
        return None


def workqueue_remote_get_workspaces_thread(wss, out, success):
    thread_name_set("get_ws")
    success.clear()
    try:
        ret = dict(wss.items())
        out.update(ret.items())
        success.set()
    except ConnectionError:
        print("workqueue_remote_get_workspaces_thread: ConnectionError")
    except EOFError:
        print("workqueue_remote_get_workspaces_thread: EOFError")
    except TimeoutError:
        print("workqueue_remote_get_workspaces_thread: TimeoutError")
    except AssertionError as e:
        print(f"workqueue_remote_get_workspaces_thread: AssertionError {e}")
    except queue.Empty:
        # print("queue.Empty")
        pass


def workqueue_remote_put_workspaces_thread(mgr, inp, success):
    thread_name_set("put_ws")
    try:
        mgr.get_workspaces().update(inp)
        # wss.update(inp)
        success.set()
    except ConnectionError:
        print("workqueue_remote_put_workspaces_thread: ConnectionError")
    except EOFError:
        print("workqueue_remote_put_workspaces_thread: EOFError")
    except TimeoutError:
        print("workqueue_remote_put_workspaces_thread: TimeoutError")
    except AssertionError as e:
        print(f"workqueue_remote_put_workspaces_thread: AssertionError {e}")
    except queue.Empty:
        # print("queue.Empty")
        pass


def workqueue_remote_get_workspaces(wq, timeout=TIMEOUT):
    out = {}
    wss = {}
    success = threading.Event()

    t = threading.Thread(
        target=workqueue_remote_get_workspaces_thread,
        args=(wq.remote_workspaces, out, success),
    )
    t.start()

    t.join(timeout=timeout)
    if (not t.is_alive()) or success.is_set():
        wss.update(out)

    print(f"returning workspaces {wss}")
    return wss


def workqueue_remote_put_workspaces(wq, wss, timeout=TIMEOUT):
    success = threading.Event()

    t = threading.Thread(
        target=workqueue_remote_put_workspaces_thread,
        args=(wq.mgr, wss, success),
    )
    t.start()
    t.join(timeout=timeout)

    return success.is_set()


def manager_remote_get_iqueue_thread(mgr, out):
    thread_name_set("get_iq")
    try:
        iq = mgr.get_iqueue()
        out.append(iq)
    except ConnectionError:
        print("manager_remote_get_iqueue_thread: ConnectionError")
    except EOFError:
        print("manager_remote_get_iqueue_thread: EOFError")
    except TimeoutError:
        print("manager_remote_get_iqueue_thread: TimeoutError")
    except AssertionError as e:
        print(f"manager_remote_get_iqueue_thread: AssertionError {e}")
    except queue.Empty:
        print("manager_remote_get_oqueue_thread: queue.Empty")
        pass
    except Exception as e:
        print("manager_remote_get_oqueue_thread: Exception {e}")


def manager_remote_get_oqueue_thread(mgr, out):
    thread_name_set("get_oq")
    try:
        # print("manager_remote_get_oqueue_thread: STARTING mgr.get_oqueue")
        oq = mgr.get_oqueue()
        out.append(oq)
    except ConnectionError:
        print("manager_remote_get_oqueue_thread: ConnectionError")
    except EOFError:
        print("manager_remote_get_oqueue_thread: EOFError")
    except TimeoutError as e:
        print(f"manager_remote_get_oqueue_thread: TimeoutError {e}")
    except AssertionError as e:
        print(f"manager_remote_get_oqueue_thread: AssertionError {e}")
    except queue.Empty:
        print("manager_remote_get_oqueue_thread: queue.Empty")
        pass
    except Exception as e:
        print("manager_remote_get_oqueue_thread: Exception {e}")


def manager_remote_queue_put_thread(oq, obj, n, success):
    thread_name_set("q_put")
    success.clear()
    try:
        oq.put(obj, block=False, n=n)
        success.set()
    except ConnectionError:
        print("manager_remote_oqueue_put_thread: ConnectionError")
    except EOFError:
        print("manager_remote_oqueue_put_thread: EOFError")
    except TimeoutError:
        print("manager_remote_oqueue_put_thread: TimeoutError")
    except AssertionError as e:
        print(f"manager_remote_oqueue_put_thread: AssertionError {e}")
    except queue.Empty:
        # print("queue.Empty")
        pass
    except struct.error as e:
        # some madness if we are trying to ship something too large.
        # break the list into halves, and if it is a single return try to
        # send each task individually
        if n > 1:
            print(
                f"\nWarning, tried to send too much data. Breaking into pieces and retrying (current n={n})."
            )

            first = obj[: len(obj) // 2]
            if len(first) == 1:
                manager_remote_oqueue_put_thread(oq, first[0], 1, success)
            else:
                manager_remote_oqueue_put_thread(
                    oq, first[0], len(first), success
                )

            second = obj[len(obj) // 2 :]
            if len(second) == 1:
                manager_remote_oqueue_put_thread(oq, second[0], 1, success)
            else:
                manager_remote_oqueue_put_thread(
                    oq, second, len(second), success
                )
        elif n == 1 and len(obj) > 1:
            for k, v in obj.items():
                print(
                    f"\nWarning, tried to send too much data. Breaking into pieces and retrying (current tasks={len(obj)} key={k})."
                )
                manager_remote_oqueue_put_thread(oq, {k: v}, 1, success)
        else:
            success.clear()
            raise e


def queue_get_nowait_thread(q, out, n, block, timeout):
    thread_name_set("q_get")
    try:
        result = q.get(block=block, timeout=timeout, n=n)
        # print(f"I GOT RESULT from q {q} : {result}")
        out.append(result)
    except ConnectionError:
        # print("ConnectionError")
        pass
    except EOFError:
        pass
        # print("queue_get_nowait_thread: EOFError")
    except TimeoutError:
        pass
        # print("queue_get_nowait_thread: TimeoutError")
    except AssertionError as e:
        pass
        # print(f"queue_get_nowait_thread: AssertionError {e}")
    except AttributeError as e:
        pass
        # print(f"queue_get_nowait_thread: AttributeError {e}")
    except queue.Empty:
        if n == 1:
            out.append({})
        else:
            out.append([])
        # print("queue_get_nowait_thread: queue.Empty")
    except Exception as e:
        pass
        # print(f"queue_get_nowait_thread: Exception {e}")


def manager_remote_queue_qsize_thread(q, out):
    thread_name_set("q_qsize")
    try:
        # result = q.get(block=True, timeout=TIMEOUT)
        result = q.qsize()
        # print(f"I GOT RESULT from q {q} : {result}")
        out.append(result)
        # success.set()
    except ConnectionError:
        # print("ConnectionError")
        pass
    except EOFError:
        pass
        # print("queue_get_nowait_thread: EOFError")
    except TimeoutError:
        pass
        # print("queue_get_nowait_thread: TimeoutError")
    except AssertionError as e:
        pass
        # print(f"queue_get_nowait_thread: AssertionError {e}")
    except AttributeError as e:
        pass
        # print(f"queue_get_nowait_thread: AttributeError {e}")
    except Exception as e:
        pass
        # print(f"queue_get_nowait_thread: Exception {e}")


def queue_get_nowait(q, block=False, timeout=TIMEOUT, n=1):
    # print(f"queue_get_nowait")
    # success = threading.Event()
    out = []
    thread_timeout = timeout
    if timeout is not None and timeout > 1.0:
        # give some space since the timeout is really for the thread timeout
        timeout -= 0.5
    t = threading.Thread(
        target=queue_get_nowait_thread, args=(q, out, n, block, timeout)
    )
    t.start()
    # print(f"queue_get_nowait starting thread")
    t.join(timeout=thread_timeout)
    # print(f"queue_get_nowait joined")

    # alive = t.is_alive()
    # print(f"queue_get_nowait alive {alive}")
    # is_success = success.is_set()
    # print(f"queue_get_nowait success {is_success}")
    result = None
    # print(f"out is {out}")
    if out:
        result = out[0]

    # print(f"queue_get_nowait returning")
    return result


def manager_remote_queue_qsize(q, timeout=TIMEOUT):
    out = []
    sz = None

    t = threading.Thread(
        target=manager_remote_queue_qsize_thread, args=(q, out)
    )
    # print(f"manager_remote_queue_qsize starting thread")
    t.start()

    t.join(timeout=timeout)
    # print(f"manager_remote_get_iqueue joined")
    if out:
        sz = out[0]

    # print(f"Get iqueue status: alive: {t.is_alive()} success {success.is_set()}")

    return sz


def manager_remote_get_iqueue(mgr, timeout=TIMEOUT):
    out = []
    iq = None

    t = threading.Thread(
        target=manager_remote_get_iqueue_thread, args=(mgr, out)
    )
    print(f"manager_remote_get_iqueue starting thread")
    t.start()

    t.join(timeout=timeout)
    # print(f"manager_remote_get_iqueue joined")
    if out:
        iq = out[0]

    # print(f"Get iqueue status: alive: {t.is_alive()} success {success.is_set()}")

    return iq


def manager_remote_queue_put(oq, obj, timeout=TIMEOUT, n=1):
    success = threading.Event()

    t = threading.Thread(
        target=manager_remote_queue_put_thread, args=(oq, obj, n, success)
    )
    t.start()
    t.join(timeout=timeout)

    return (not t.is_alive()) or success.is_set()


def manager_remote_get_oqueue(mgr, timeout=TIMEOUT):
    out = []
    oq = None

    print(f"manager_remote_get_oqueue starting thread")
    t = threading.Thread(
        target=manager_remote_get_oqueue_thread, args=(mgr, out)
    )
    t.start()

    t.join(timeout=timeout)
    if out:
        oq = out[0]

    return oq


class workqueue_remote(workqueue):
    def __init__(self, addr, port):
        super().__init__(addr, port)
        self.mgr.register("get_workspaces")
        self.mgr.register("get_state")
        self.is_connected = self.connect()
        self.remote_state = None
        if self.is_connected:
            self.remote_state = self.mgr.get_state()
            self.remote_workspaces = self.mgr.get_workspaces()

    def connect(self):
        return remote_connect(self.mgr)

    def get_workspaces(self):
        print(f"GET REMOTE Workspace")
        return workqueue_remote_get_workspaces(self)

    def put_workspaces(self, wss):
        print(f"PUT REMOTE Workspace")
        return workqueue_remote_put_workspaces(self, wss)

    def get_state(self):
        if self.remote_state is None:
            self.remote_state = manager_remote_get_state(self)
        return state

    def get_status(self):
        return manager_remote_get_status(self)


class workspace:
    def __init__(self, addr, port):
        self.mgr = workspace_manager(address=(addr, port), authkey=b"0")
        # self.state = None
        self.holding = set()

    def get_iqueue(self):
        pass

    def get_rqueue(self):
        pass

    def get_oqueue(self):
        pass

    def get_state(self):
        pass

    def get_status(self):
        pass


_ForkingPickler = context.reduction.ForkingPickler


class myqueue(queues.Queue):
    def __init__(self, maxsize=0):
        super().__init__(maxsize=maxsize, ctx=context._default_context)

    def get(self, block=True, timeout=None, n=1):
        ret = []
        if self._closed:
            raise ValueError(f"Queue {self!r} is closed")
        if block and timeout is None:
            with self._rlock:
                for i in range(n):
                    res = self._recv_bytes()
                    ret.append(res)
                    self._sem.release()
        else:
            if block:
                deadline = time.monotonic() + timeout
            if not self._rlock.acquire(block, timeout):
                raise queue.Empty
            try:
                for i in range(n):
                    if block:
                        timeout = deadline - time.monotonic()
                        if not self._poll(timeout):
                            if ret:
                                break
                            else:
                                raise queue.Empty
                    elif not self._poll():
                        if ret:
                            break
                        else:
                            raise queue.Empty
                    res = self._recv_bytes()
                    ret.append(res)
                    self._sem.release()
            finally:
                self._rlock.release()
        # unserialize the data after having released the lock
        ret = [_ForkingPickler.loads(res) for res in ret]

        if n == 1:
            return ret[0]
        else:
            return ret

    def put(self, obj, block=True, timeout=None, n=1):
        if self._closed:
            raise ValueError(f"Queue {self!r} is closed")
        if not self._sem.acquire(block, timeout):
            raise Full

        with self._notempty:
            if self._thread is None:
                self._start_thread()
            if n == 1:
                obj = [obj]
            self._buffer.extend(obj)
            self._notempty.notify()

    # def put(self, obj):

    #     # print("HELO", self)
    #     # print("WLOCK PUT", self._wlock)
    #     # serialize the data before acquiring the lock
    #     # print(f"PUTTING {obj}")
    #     obj = _ForkingPickler.dumps(obj)
    #     if self._wlock is None:
    #         # writes to a message oriented win32 pipe are atomic
    #         self._writer.send_bytes(obj)
    #     else:
    #         with self._wlock:
    #             # print("LOCKED PUT", self._wlock)
    #             self._writer.send_bytes(obj)
    #             # print("LOCKED PUT SENT", self._wlock)
    #     # print("WUNLOCK PUT", self._wlock)


class myqueue(queue.Queue):
    def put(self, item, block=True, timeout=None, n=1):
        """Put an item into the queue.

        If optional args 'block' is true and 'timeout' is None (the default),
        block if necessary until a free slot is available. If 'timeout' is
        a non-negative number, it blocks at most 'timeout' seconds and raises
        the Full exception if no free slot was available within that time.
        Otherwise ('block' is false), put an item on the queue if a free slot
        is immediately available, else raise the Full exception ('timeout'
        is ignored in that case).
        """

        remain = []
        with self.not_full:
            if self.maxsize > 0:
                if not block:
                    if self._qsize() >= self.maxsize - (n - 1):
                        raise queue.Full
                elif timeout is None:
                    while self._qsize() >= self.maxsize - (n - 1):
                        self.not_full.wait()
                elif timeout < 0:
                    raise ValueError("'timeout' must be a non-negative number")
                else:
                    endtime = time.time() + timeout
                    while self._qsize() >= self.maxsize - (n - 1):
                        remaining = endtime - time.time()
                        if remaining <= 0.0:
                            raise queue.Full
                        self.not_full.wait(remaining)
            if n == 1:
                item = [item]
            for o in item:
                self._put(o)

            self.not_empty.notify_all()
            # with self.not_empty:

    def get(self, block=True, timeout=None, n=1):
        """Remove and return an item from the queue.

        If optional args 'block' is true and 'timeout' is None (the default),
        block if necessary until an item is available. If 'timeout' is
        a non-negative number, it blocks at most 'timeout' seconds and raises
        the Empty exception if no item was available within that time.
        Otherwise ('block' is false), return an item if one is immediately
        available, else raise the Empty exception ('timeout' is ignored
        in that case).
        """
        ret = []
        with self.not_empty:
            for i in range(n):
                if not block:
                    if not self._qsize():
                        if ret:
                            break
                        else:
                            raise queue.Empty
                elif timeout is None:
                    while not self._qsize():
                        self.not_empty.wait()
                elif timeout < 0:
                    raise ValueError("'timeout' must be a non-negative number")
                else:
                    endtime = time.time() + timeout
                    while not self._qsize():
                        remaining = endtime - time.time()
                        if remaining <= 0.0:
                            if ret:
                                break
                            else:
                                raise queue.Empty
                        self.not_empty.wait(remaining)
                item = self._get()
                ret.append(item)

            self.not_full.notify_all()

        if n == 1:
            ret = ret[0]
        return ret


# myqueue = multiprocessing.Queue
# myqueue = queue.Queue


class myiqueue(queue.Queue):
    def get(self, block=True, timeout=None):
        raise Exception()
        print(f"GETTING IQ from {id(self)}")
        return super().get(block=block, timeout=timeout)


def shm_init(proxy):
    """ """
    print(f"{datetime.now()} shm_init: building shm")
    shm = shm_local()
    data = dict(proxy.get())
    shm.__dict__.update(data)
    print(f"{datetime.now()} shm_init: shm has members {list(data.keys())}")
    return shm


class shm_local:
    def __init__(self, procs_per_task=1, data=None):
        self.procs_per_task = procs_per_task

        if data is not None:
            self.__dict__.update(data)

    def get(self):
        return self.__dict__

    def remote_init(self):
        return shm_init


class workspace_local(workspace):
    """
    Assumes that we are process-local to all needed resources and do not need
    to use the proxy interface i.e. no connection needed
    """

    def __init__(self, addr, port, shm: shm_local = None, nproc=-1):
        super().__init__(addr, port)

        self.mgr._Client = functools.partial(Client, timeout=TIMEOUT)
        self.state = dict({"status": workspace_status.EMPTY})

        self.iqueue = myqueue()
        self.oqueue = myqueue()
        self.holding_remote = {}
        self.holding_remote_lock = threading.Lock()
        self.finished = 0
        self.finished_remote = 0
        self.mgr.register("get_iqueue", lambda: self.iqueue)
        self.mgr.register("get_oqueue", lambda: self.oqueue)

        self.mgr.register("get_state", lambda: self.state)

        # this will load whatever interface that shm has
        self.mgr.register("get_shm", lambda: self.shm)
        # self.mgr.register('init_shm', functools.partial(shm_init, self.shm))

        if shm is None:
            self.shm: shm_local = shm_local()
        elif type(shm) is dict:
            self.shm: shm_local = shm_local()
            self.shm.__dict__.update(shm)
        else:
            self.shm = shm

        if nproc == -1:
            self.nproc: int = max(
                1,
                configs.processors
                if configs.processors
                else os.cpu_count() - 1,
            )
        elif nproc is None:
            self.nproc = os.cpu_count() - 1
        else:
            self.nproc = nproc
        # self.shared = shared_mem()
        # self.mgr.register('get_shared', lambda: self.shared)

        # self.nproc = min(os.cpu_count(), 2)

        ntasks = 1
        if self.shm.procs_per_task > 0:
            ntasks = max(1, self.nproc // self.shm.procs_per_task)
        else:
            self.shm.procs_per_task = self.nproc

        global SHM_GLOBAL
        SHM_GLOBAL = self.shm
        # workspace_run_init(self.shm)
        self.pool = multiprocessing.pool.Pool(
            ntasks,
            workspace_run_init,
            (self.shm.procs_per_task,),
            context=multiprocessing.get_context("fork"),
        )

        self.done = threading.Event()
        self.done.clear()

        self.loadbalance_stop = threading.Event()
        self.loadbalance_stop.clear()

        self.gather_stop = threading.Event()
        self.gather_stop.clear()

        self.run_stop = threading.Event()
        self.run_stop.clear()

        self.gather_thread = None  #
        self.loadbalance_thread = None
        self.run_thread = None

        self.mgr.start()
        self.remote_iqueue = None
        self.remote_oqueue = None
        self.remote_oqueue_size = 0
        self.remote_iqueue_size = 0
        self.remote_state = None
        self.run_thread = None

        # self.pool = None
        self.reset()
        print(f"Started local workspace on {self.mgr.address}")

    def start(self):
        self.mgr.start()
        self.reset()

    def reset(self):
        self.done.clear()
        self.gather_stop.clear()
        self.loadbalance_stop.clear()
        self.run_stop.clear()

        self.iqueue = myqueue()
        self.oqueue = myqueue()

        # print("Starting iqueue reference")
        self.remote_iqueue = self.mgr.get_iqueue()
        # self.remote_iqueue._Client = Client
        # print("Starting oqueue reference")
        self.remote_oqueue = self.mgr.get_oqueue()
        # self.remote_oqueue._Client = Client

        self.remote_oqueue_size = 0
        self.remote_oqueue_size_lock = threading.Lock()
        self.remote_iqueue_size = 0
        self.remote_iqueue_size_lock = threading.Lock()

        self.remote_state = self.mgr.get_state()

        if self.pool:
            self.pool.close()
            self.pool.terminate()
            ntasks = 1
            if self.shm.procs_per_task > 0:
                ntasks = max(1, self.nproc // self.shm.procs_per_task)
            global shm_global
            shm_global = self.shm
            self.pool = multiprocessing.pool.Pool(
                ntasks,
                workspace_run_init,
                (self.shm.procs_per_task,),
                context=multiprocessing.get_context("fork"),
            )
        # self.close()
        # self.pool = None
        # if self.nproc > 0:
        #     self.pool = multiprocessing.Pool(self.nproc)
        self.holding.clear()
        self.holding_remote.clear()
        # self.set_status(workspace_status.EMPTY)

        # starts all threads
        self.run_thread = workspace_local_run(self)

    def close(self):
        try:
            self.set_status(workspace_status.DONE)
            self.gather_stop.set()
            self.loadbalance_stop.set()
            self.run_stop.set()
            self.done.set()
            if self.run_thread:
                self.run_thread.join()
            if self.pool is not None:
                self.pool.close()
                self.pool.terminate()
                self.pool = None
            self.mgr.shutdown()
            self.remote_iqueue = None
            self.remote_oqueue = None
        except BrokenPipeError:
            pass
        # self.pool.close()
        # print("Setting to done")
        # self.done.set()
        # print("Joining gather thread")
        # if self.gather_thread and self.gather_thread.is_alive():
        #     self.gather_thread.join(timeout=1.0)

    def set_status(self, status):
        self.remote_state.update({"status": status})
        # self.mgr.get_state()
        self.state["status"] = status

    def get_iqueue(self):
        return self.iqueue

    def get_oqueue(self):
        return self.oqueue

    def get_state(self):
        return self.state

    def get_status(self):
        return self.state.get("status", workspace_status.INVALID)

    def get_table(self):
        print("LOCAL GET TABLE")
        return self.mgr.get_table()

    def get_result(self):
        print("LOCAL GET RESULT")
        return self.mgr.get_result()


class workspace_remote(workspace):
    """
    Assumes that we need to go through the proxy interface to access resources
    that are not local i.e. need to go through a connection
    """

    def __init__(self, addr, port, nproc=1):
        super().__init__(addr, port)

        self.error_count: int = 0
        self.error_limit: int = 1

        self.mgr._Client = functools.partial(Client, timeout=TIMEOUT)

        print("Registering remote functions")
        self.mgr.register("get_iqueue")
        self.mgr.register("get_rqueue")
        self.mgr.register("get_oqueue")
        self.mgr.register("get_status")
        self.mgr.register("get_state")
        self.mgr.register("get_workers")
        self.mgr.register("get_table")
        self.mgr.register("get_result")
        self.mgr.register("get_shm")

        print(f"workspace_remote_init: Connecting to {addr}:{port}")
        self.is_connected = self.connect()

        if nproc is None:
            self.nproc: int = max(1, os.cpu_count() - 1)
        else:
            self.nproc = int(nproc)

        self.pool = None

        self.iqueue = None
        self.oqueue = None

        # stop everything
        self.done = threading.Event()

        # stop receiving data, but allow flushing remaining
        self.stop = threading.Event()

        self.gather_thread = None
        self.pusher_thread = None

        self.remote_state = None

    def start(self):
        if self.is_connected:
            print("Connecting queues...")
            self.remote_iqueue = self.get_iqueue()
            self.remote_oqueue = self.get_oqueue()
            print("Connecting state...")
            self.remote_state = self.get_state()
            self.iqueue = myqueue(maxsize=1)
            self.oqueue = myqueue()

            if not (
                self.remote_iqueue and self.remote_oqueue and self.remote_state
            ):
                self.is_connected = False
                print("Connecting queues or state failed.")
                return False
            else:
                print(f"Starting pool with {self.nproc} processes")
                ntasks = 1
                if self.shm.procs_per_task > 0:
                    ntasks = max(1, self.nproc // self.shm.procs_per_task)
                global shm_global
                shm_global = self.shm
                self.pool = multiprocessing.Pool(
                    ntasks, workspace_run_init, (self.shm.procs_per_task,)
                )
                print("Launching threads")
                self.gather_thread = threading.Thread(
                    target=workspace_remote_local_gather_thread, args=(self,)
                )  #
                self.gather_thread.start()
                self.pusher_thread = threading.Thread(
                    target=workspace_remote_local_pusher_thread, args=(self,)
                )  #
                self.pusher_thread.start()
                return True
        else:
            return False

    def close(self):
        if self.pool is not None:
            # print("Closing pool")
            self.pool.terminate()
            self.pool.close()
        self.stop.set()
        self.done.set()

        self.remote_iqueue = None
        self.remote_oqueue = None
        self.remote_state = None

    def get_iqueue(self):
        q = manager_remote_get_iqueue(self.mgr)
        if q is None:
            print("get_iqueue: error")
            self.error_count += 1
        return q

    def get_oqueue(self):
        q = manager_remote_get_oqueue(self.mgr)
        if q is None:
            print("get_oqueue: error")
            self.error_count += 1
        return q

    def oqueue_put(self, obj, n=1):
        success = manager_remote_queue_put(self.remote_oqueue, obj, n=n)
        if not success:
            print("oqueue_put: error")
            self.error_count += 1
        return success

    def iqueue_get(self, n=1):
        item = queue_get_nowait(self.remote_iqueue, n=n)
        if item is None:
            print("iqueue_get: error")
            self.error_count += 1

        return item

    def iqueue_put(self, obj, n=1):
        success = manager_remote_queue_put(self.remote_iqueue, obj, n=n)
        if not success:
            print("iqueue_put: error")
            self.error_count += 1
        return success

    def connect(self):
        success = remote_connect(self.mgr)
        if not success:
            self.error_limit += 1
        return success

    def get_state(self):
        if self.remote_state is None:
            self.remote_state = manager_remote_get_state(self)
        if self.remote_state is None:
            print("get_state: error")
            self.error_count += 1
        return self.remote_state

    def get_status(self):
        status = manager_remote_get_status(self)
        if status == workspace_status.INVALID:
            print("get_status: error")
            self.error_count += 1
        return status


def workspace_remote_local_pusher_thread(ws: workspace_remote):
    """
    pull results from the remote workers and put them in the local queue
    """
    thread_name_set("pusher")
    # print("GATHER THREAD ACTIVE")
    # remote_q = ws.get_oqueue()
    # remote_q = None
    # while remote_q is None or ws.done.is_set():
    #     # there are two timeouts; this one is for the thread
    #     # there is another timeout for the socket so we loop here...
    #     # we need to block on a local socket, but not block on remote..
    #     remote_q = manager_remote_get_oqueue(ws.mgr, timeout=None)
    # q = ws.get_oqueue()
    while ws.oqueue is None:
        time.sleep(0.01)
    good = 0
    bad = 0

    packets = []
    sleepiness = 0.0
    while not ws.done.is_set():
        # print(f"REMOTE GET FUNCTIONS from q {iq}" )

        # print("workspace_remote_get_functions: Getting item from iq...")

        # time.sleep(5.0)
        t0 = time.perf_counter()
        n = max(2, 100 - len(packets))
        if n <= 100:
            try:
                items = ws.oqueue.get(block=False, n=n)
                packets.extend(items)
                # for item in items:
                #     packet.update(item)
            except queue.Empty:
                sleepiness += 2.0
                time.sleep(sleepiness)
        lt = time.perf_counter() - t0
        if packets:
            # need to make sure packet is not too large
            # try to guess 10MB chunks
            # sections = len(pickle.dumps(packet)) // int(10e6) + 1
            # unsent = {}
            # failed = True
            # t0 = time.perf_counter()
            tosend = packets
            if len(packets) > 1:
                n = len(packets)
            else:
                n = 1
                tosend = packets[0]

            # print(f"{datetime.now()} STARTING PUSH {len(packets)} \nVALUES ARE\n{pprint.pformat(packets)}")
            failed = not ws.oqueue_put(tosend, n=n)
            rt = time.perf_counter() - t0
            if failed:
                bad += 1

            else:
                good += n
                packets.clear()
                sleepiness = 0.0
            # print(
            #     f"{datetime.now()} REMOTE PUSH success {not failed} good={good} bad={bad} this={len(packets)} rtime: {rt:6.3f} ltime: {lt:6.3f}"
            # )
            # if ws.done.is_set():
            #     break
            if failed:
                sleepiness += 2.0
                time.sleep(sleepiness)
            # packets.clear()
            # packet.update(unsent)
            # time.sleep(1.0)
        # print(f"REMOTE GATHER THREAD PUSHED {item[0]}")


def workspace_remote_local_gather_thread(ws: workspace_remote):
    """
    pull results from the remote workers and put them in the local queue
    """
    thread_name_set("gather")
    # print("GATHER THREAD ACTIVE")
    # remote_q = ws.get_oqueue()
    # remote_q = None
    # while remote_q is None or ws.done.is_set():
    #     # there are two timeouts; this one is for the thread
    #     # there is another timeout for the socket so we loop here...
    #     # we need to block on a local socket, but not block on remote..
    #     remote_q = manager_remote_get_oqueue(ws.mgr, timeout=None)

    # might be None be meh
    iq = ws.remote_iqueue

    sleepiness = 0.0
    while not (ws.done.is_set() or ws.stop.is_set()):
        # print(f"REMOTE GET FUNCTIONS from q {iq}" )

        # print("workspace_remote_get_functions: Getting item from iq...")

        t0 = time.perf_counter()
        item = ws.iqueue_get(n=1)
        # print(f"workspace_remote_get_functions: received {result} {result is not None}")
        rt = time.perf_counter() - t0
        if item is not None and len(item):
            if type(item) is not dict:
                print(
                    f"Warning, gather thread received malformed taskset:\n{item}"
                )
            t0 = time.perf_counter()
            # this has a maxsize and will block
            # if we timeout
            try:
                ws.iqueue.put(item, n=1, timeout=60)
                sleepiness = 0.0
            except queue.Full:
                # this is an attempt to not hold on to jobs and idle other
                # workers if there are few jobs
                iqs = manager_remote_queue_qsize(iq)
                if iqs is not None and iqs == 0:
                    ws.iqueue_put(item, n=1)

            lt = time.perf_counter() - t0
            # print(
            #     f"{datetime.now()} REMOTE GATHER THREAD RECEIVED {len(item)} rtime: {rt:6.3f} ltime: {lt:6.3f}"
            # )
        else:
            sleepiness += 2.0
            time.sleep(sleepiness)
        # else:
        #     time.sleep(.5)


def workspace_local_remote_loadbalance_thread(ws: workspace_local):
    """
    pull results from the remote workers and put them in the local queue
    """
    thread_name_set("loadbal")
    # print(f"LOADBALANCE THREAD ACTIVE on {id(ws.iqueue)}")
    # remote_q = ws.get_oqueue()
    # remote_q = None
    # while remote_q is None or ws.done.is_set():
    #     # there are two timeouts; this one is for the thread
    #     # there is another timeout for the socket so we loop here...
    #     # we need to block on a local socket, but not block on remote..
    #     remote_q = manager_remote_get_oqueue(ws.mgr, timeout=None)
    # remote_q = ws.mgr.get_oqueue()
    remote_q = ws.remote_iqueue
    # print("GATHER THREAD HAVE Q")
    packet = []
    i = 0
    sleepiness = 0.0
    while (
        not ws.done.is_set()
        and not ws.loadbalance_stop.is_set()
        and ws.remote_iqueue is not None
    ):
        t0 = time.perf_counter()
        try:
            # rt = time.perf_counter() - t0
            # t0 = time.perf_counter()
            rqs = None
            while rqs is None:
                rqs = manager_remote_queue_qsize(remote_q)
                if ws.done.is_set() or ws.loadbalance_stop.is_set():
                    break
                time.sleep(min(10.0, sleepiness))

            if rqs is None:
                continue

            with ws.remote_iqueue_size_lock:
                ws.remote_iqueue_size = rqs
            iqs = ws.iqueue.qsize()

            if rqs > 10000 and iqs > 10000:
                time.sleep(5 * 60.0)
                continue

            elif rqs > 1000 and iqs > 1000:
                time.sleep(1 * 60.0)
                continue

            elif rqs > 500 and iqs > 500:
                time.sleep(1 * 10.0)
                continue

            if rqs < iqs:
                n = max(1, (iqs - rqs) // 2)
                n = min(100, n)
                obj = ws.iqueue.get(block=False, n=n)
                if obj is not None:
                    if n == 1 or (n > 1 and len(obj) == 1):
                        if n > 1 and len(obj) == 1:
                            n = 1
                            obj = obj[0]
                        # print(f"\nPULLED {len(obj)} from LOCAL rqs {rqs} iqs {iqs}")
                        if not manager_remote_queue_put(remote_q, obj, n=n):
                            ws.iqueue.put(obj, n=n)
                        else:
                            if n == 1:
                                with ws.holding_remote_lock:
                                    ws.holding_remote.update(obj)
                            else:
                                with ws.holding_remote_lock:
                                    for o in obj:
                                        ws.holding_remote.update(o)
                            with ws.remote_iqueue_size_lock:
                                ws.remote_iqueue_size = (
                                    manager_remote_queue_qsize(remote_q)
                                )
                    elif n > 1:
                        n = len(obj)
                        # print(f"\nPULLED single from LOCAL rqs {rqs} iqs {iqs}")
                        # remote_q.put(obj)
                        # ws.iqueue.put(obj)
                        if not manager_remote_queue_put(remote_q, obj, n=n):
                            ws.iqueue.put(obj, n=n)
                        else:
                            with ws.holding_remote_lock:
                                for o in obj:
                                    ws.holding_remote.update(o)
                            with ws.remote_iqueue_size_lock:
                                ws.remote_iqueue_size = (
                                    manager_remote_queue_qsize(remote_q)
                                )
            elif iqs < rqs:
                n = max(1, (rqs - iqs) // 2)
                n = min(100, n)
                obj = queue_get_nowait(remote_q, n=n)
                if obj is not None and len(obj):
                    if n > 1:
                        n = len(obj)
                        # with ws.holding_remote_lock:
                        #     for o in obj:
                        #         for k in o:
                        #             if k in ws.holding_remote:
                        #                 ws.holding_remote.pop(k)
                        if n == 1:
                            obj = obj[0]
                    # else:
                    #     with ws.holding_remote_lock:
                    #         for k in obj:
                    #             if k in ws.holding_remote:
                    #                 ws.holding_remote.pop(k)
                    ws.iqueue.put(obj, n=n)
                    with ws.remote_iqueue_size_lock:
                        ws.remote_iqueue_size = manager_remote_queue_qsize(
                            remote_q
                        )

            else:
                sleepiness += 5.0
                time.sleep(min(60.0, sleepiness))

        except queue.Empty:
            pass
        t1 = time.perf_counter()
        dt = t1 - t0
        if dt < 0.0:
            time.sleep(0.0 - dt)
        # time.sleep(1.0)


def workspace_local_remote_gather_thread(ws: workspace_local):
    """
    pull results from the remote workers and put them in the local queue
    """
    thread_name_set("gather")
    # print("GATHER THREAD ACTIVE")
    # remote_q = ws.get_oqueue()
    # remote_q = None
    # while remote_q is None or ws.done.is_set():
    #     # there are two timeouts; this one is for the thread
    #     # there is another timeout for the socket so we loop here...
    #     # we need to block on a local socket, but not block on remote..
    #     remote_q = manager_remote_get_oqueue(ws.mgr, timeout=None)
    # remote_q = ws.mgr.get_oqueue()
    remote_q = ws.remote_oqueue
    # print("GATHER THREAD HAVE Q")

    i = 0
    sleepiness = 0.0
    while (
        not ws.done.is_set()
        and not ws.gather_stop.is_set()
        and ws.remote_oqueue is not None
    ):
        t0 = time.perf_counter()
        try:
            # t0 = time.perf_counter()
            # rt = time.perf_counter() - t0
            # t0 = time.perf_counter()
            # we want this to block
            obj = queue_get_nowait(
                remote_q, block=False, timeout=None, n=10000
            )
            if obj is not None and len(obj):
                with ws.remote_oqueue_size_lock:
                    ws.remote_oqueue_size = manager_remote_queue_qsize(
                        remote_q
                    )
                ws.finished_remote += sum((len(x) for x in obj))
                n = len(obj)
                if n == 1:
                    obj = obj[0]
                ws.oqueue.put(obj, block=False, n=n)
                sleepiness = 0.0
            # lt = time.perf_counter() - t0
            # print(
            #     f"\n{id(ws.iqueue)} LOCAL GATHER THREAD RECEIVED {list(obj.keys())} rtime: {rt:6.3f} ltime: {lt:6.3f}"
            # )

            i += 1
            # time.sleep(1.0)
        except queue.Empty:
            sleepiness += 1.0
            time.sleep(min(sleepiness, 60.0))
            # pass
        t1 = time.perf_counter()
        dt = t1 - t0
        if dt < 0.0:
            time.sleep(0.0 - dt)
        # else:
        #     time.sleep(1.0)
    # print(f"GATHER THREAD DONE; transferred {i} items")


def workqueue_push_workspace(wq: workqueue_local, ws: workspace):
    addr = ws.mgr.address
    # print("Getting workspaces")
    wss = wq.get_workspaces()
    v = wss.get(addr, 1)
    # v = wq.mgr.get_workspaces().get(addr, 0)
    v -= 1
    # print("Putting workspaces")
    if not wq.put_workspaces({addr: v}):  # .update({addr: v})
        print(f"Could not push workspace.")


def workspace_local_run_thread(ws: workspace_local):
    thread_name_set("run")
    # print(f"PROCESSING QUEUE {ws}")
    # print(f"iq {ws.iqueue}")
    # print(f"oq {ws.oqueue}")
    started = False
    work = []
    # print("local_run_thread: getting mgr.iqueue")
    iq = ws.remote_iqueue
    # remote_workers = ws.mgr.get_workers()
    # remote_workers = 1
    # print(dir(remote_workers))
    # print("Starting local pool")
    # self.pool = None
    # if self.nproc > 0:
    #     self.pool = multiprocessing.Pool(self.nproc)

    ntasks = 1
    if ws.shm.procs_per_task:
        ntasks = max(1, ws.nproc // ws.shm.procs_per_task)
    local_n = 0
    pool = ws.pool

    # if ws.nproc and (ws.pool is None):
    #     print("workspace_local_run_thread: Starting pool")
    #     pool = multiprocessing.pool.Pool(ws.nproc, workspace_run_init, (ws.shm,), context=multiprocessing.get_context("fork"))
    #     # ws.pool = pool

    # elif ws.pool is not None:
    #     print("workspace_local_run_thread: Using pool")
    #     pool = ws.pool
    # else:
    if pool is None:
        print("WARNING POOL IS NONE")
        return

    if ws.mgr.address[0] != "127.0.0.1":
        ws.gather_stop.clear()
        ws.gather_thread = threading.Thread(
            target=workspace_local_remote_gather_thread, args=(ws,)
        )
        # print("Starting gather thread...")
        ws.gather_thread.start()

        ws.loadbalance_stop.clear()
        ws.loadbalance_thread = threading.Thread(
            target=workspace_local_remote_loadbalance_thread, args=(ws,)
        )
        # print("Starting loadbalance thread...")
        ws.loadbalance_thread.start()

    # print("workspace_local_run_thread: Entering run loop")
    try:
        while (
            not ws.done.is_set() and not ws.run_stop.is_set()
        ) or ws.iqueue.qsize():
            t0 = time.perf_counter()
            # time.sleep(.1)
            # remote_p = int(remote_workers)
            idx = None
            local_n = len(work)
            # workspace_local_remote_loadbalance_thread(ws)
            # print("workspace_local_run_thread: In loop")
            functions = {}
            put_back = {}
            stole = False
            try:
                if local_n < ntasks:
                    # print(f"\nworkspace_local_run_thread: Getting local functions {id(ws.iqueue)}")
                    functions = ws.iqueue.get(block=False)
                # if ws.iqueue.qsize():
                # else:
                #     time.sleep(1.0)
                # riqsize = iq.qsize()
            except queue.Empty:
                # print(f"\nworkspace_local_run_thread: LOCAL EMPTY")
                pass
            except Exception:
                # qsize can fail
                print(f"\nworkspace_local_run_thread: Exception {type(e)} {e}")
                continue
            # print(f"\nworkspace_local_run_thread: received {len(functions)} tasks")
            try:
                # this locks with the loadbalancer, so avoid for now
                if not functions:
                    if (
                        False
                        and ntasks > 0
                        and (len(work) < ntasks or len(ws.holding) < ntasks)
                    ):
                        # steal from the remote queue
                        time.sleep(1.0)
                        riqsize = (
                            ws.remote_iqueue_size
                        )  # manager_remote_queue_qsize(iq)
                        while riqsize and not functions:
                            # print("workspace_local_run_thread: Trying to steal remote functions")

                            put_back.clear()
                            functions = queue_get_nowait(iq)
                            # functions = iq.get(block=False)
                            if functions is None or len:
                                continue

                            for idx in list(functions):
                                if idx in ws.holding:
                                    put_back[idx] = functions.pop(idx)

                            if put_back:
                                # print("putting back existing functions")
                                manager_remote_queue_put(
                                    iq, put_back, block=False
                                )
                                # iq.put(put_back, block=False)
                        # while idx in ws.holding:
                        #     put_back.append((idx, distfun))
                        #     idx, distfun = iq.get(block=False)
                        # for oldidx, olddistfun in put_back:
                        # print(f"Put back {oldidx} to IQ")
                        # iq.put((oldidx, olddistfun))
                        # print(f"Stole {idx} from IQ")

                    # replicate remote running jobs
                    elif len(ws.holding_remote) < ntasks * 2 and (
                        len(work) < ntasks or len(ws.holding) < ntasks
                    ):
                        # throttle calls to qsize
                        riqsize = ws.remote_iqueue_size
                        if riqsize == 0:
                            with ws.holding_remote_lock:
                                functions = dict(
                                    list(ws.holding_remote.items())[
                                        -ntasks * 2 :
                                    ]
                                )
                            for idx in list(functions):
                                if idx in ws.holding:
                                    functions.pop(idx)
                    # if functions:
                    #     print("Taking an assumed dead function")
                    # print(f"Stole {idx} from RIP")
                    # ws.holding_remote.pop(idx)

            except queue.Empty:
                pass
            # check here since the above q ops are slow
            if ws.done.is_set():
                break

            if functions:
                # print(f"There are {len(functions)} to run")

                local_submit = True

                # if ws.nproc > 0:
                #     if local_n >= ws.nproc:
                #         # case 1: we have local and remote workers
                #         local_submit = False
                # else:
                #     local_submit = False

                if local_submit:
                    for idx, distfun in functions.items():
                        # print("local_n = ", local_n, "remote_n:", remote_n)
                        # print("LOCAL SUBMIT?", local_submit)
                        if local_submit and idx not in ws.holding:
                            # print(f"Pushing job {idx} to local q {iq}")
                            work.append(
                                (
                                    idx,
                                    pool.apply_async(
                                        workspace_run, (distfun,), {}
                                    ),
                                )
                            )
                            ws.holding.add(idx)
                else:
                    remote_put = {
                        idx: unit
                        for idx, unit in functions.items()
                        if idx not in ws.holding_remote
                    }
                    remote_puts = []
                    if remote_put:
                        remote_puts.append(remote_put)
                    # put 50% of the work on the queue
                    # if iq.qsize() == 0:
                    #     for i in range(int(.5*ws.iqueue.qsize())):
                    #         try:
                    #             functions = ws.iqueue.get(block=False)
                    #             remote_puts.append(functions)
                    #         except queue.Empty:
                    #             break

                    # print(f"\nCollected {len(remote_puts)}/{int(0.5 * ws.iqueue.qsize())} items for remote push")
                    for remote_put in remote_puts:
                        if remote_put:
                            # print(f"Pushing job {remote_put.keys()} to remote q {iq}")
                            manager_remote_queue_put(
                                iq, remote_put, block=False
                            )
                            # iq.put(remote_put, block=False)
                            # print(f"Pushed")
                            with ws.holding_remote_lock:
                                ws.holding_remote.update(remote_put)
                            if ws.nproc:
                                for idx, distfun in remote_put.items():
                                    if (
                                        local_n < ntasks
                                        and idx not in ws.holding
                                    ):
                                        work.append(
                                            (
                                                idx,
                                                pool.apply_async(
                                                    workspace_run,
                                                    (distfun,),
                                                    {},
                                                ),
                                            )
                                        )
                                        ws.holding.add(idx)

            if work:
                drop = set()
                # print(f"Scanning {len(work)} work units")
                for i in range(len(work)):
                    idx, unit = work[i]
                    if unit.ready():
                        drop.add(i)
                        result = unit.get()
                        # print(f"Unit {idx} is ready")
                        # put_functions(ws.oqueue, [(idx, result)])
                        # if idx == 32:
                        #     breakpoint()
                        if idx in ws.holding:
                            ws.holding.remove(idx)
                        with ws.holding_remote_lock:
                            if idx in ws.holding_remote:
                                ws.holding_remote.pop(idx)
                        ws.oqueue.put({idx: result}, block=False)
                        # print(f"Unit {idx} is done")

                working = [x for i, x in enumerate(work) if i not in drop]
                # print(f"Just finished:\n{drop}")
                # print(f"Still have:\n{working}")
                ws.finished += len(drop)

                work.clear()
                work.extend(working)

                working.clear()
            else:
                # if we have no work then we are not holding anything
                ws.holding.clear()
            t1 = time.perf_counter()
            dt = t1 - t0
            if dt < 0.0:
                time.sleep(0.0 - dt)

    except BrokenPipeError as e:
        print(f"Warning, BrokenPipeError {e}")

    # if ws.pool is None:
    #     pool.terminate()
    #     pool.close()
    ws.loadbalance_stop.set()
    if ws.loadbalance_thread:
        ws.loadbalance_thread.join(timeout=None)

    ws.gather_stop.set()
    if ws.gather_thread:
        ws.gather_thread.join(timeout=None)
    # print("PROCESSING THREAD DONE")
    return


def workspace_local_run(ws: workspace_local):
    """
    Take jobs from the input queue and distribute to the processing queues,
    which can be a (low latency) local pool or a (high latency) manager
    """
    ws.set_status(workspace_status.RUNNING)
    # global shm
    # shm = ws.shm
    # print(f"Run thread pool is {ws.pool}")
    t = threading.Thread(target=workspace_local_run_thread, args=(ws,))
    t.start()
    return t


def workspace_submit_and_flush(
    ws, fn, iterable, chunksize=1, timeout=1.0, batchsize=0
):
    results = {}
    j = len(results)

    n = len(iterable)

    if n == 0:
        return results

    if batchsize == 0:
        batchsize = n

    while len(results) < n:
        todo = {
            idx: unit for idx, unit in iterable.items() if idx not in results
        }
        for batch in arrays.batched(todo.items(), batchsize):
            for chunk in arrays.batched(batch, chunksize):
                tasks = {}
                for idx, (args, kwds) in chunk:
                    tasks[idx] = (fn, args, kwds)

                workspace_local_submit(ws, tasks)
            results.update(
                workspace_flush(
                    ws, set((x[0] for x in batch)), timeout=timeout
                )
            )
        j = len(results)
    print(f"Batch: {j/n*100:5.2f}%  {j:8d}/{n}")
    return results


def workspace_flush(ws: workspace_local, indices, timeout: float = TIMEOUT):
    # print(f"JOINING 1000 functions")
    if len(indices) == 0:
        return {}
    results = {}
    ws.finished = 0
    ws.finished_remote = 0
    oq = ws.oqueue
    iq = ws.iqueue

    # riq = ws.mgr.get_iqueue()
    # riq = manager_remote_get_iqueue(ws.mgr)
    riq = ws.remote_iqueue
    if not riq:
        return {}
    # roq = manager_remote_get_oqueue(ws.mgr)
    roq = ws.remote_oqueue
    if not roq:
        return {}
    n = len(indices)
    at_least_one = False
    waited = 0.0
    waittime = 5.0

    totalwait = None
    if timeout is not None:
        waittime = min(timeout, waittime)
        totalwait = timeout
    # ws.holding_remote.clear()
    for idx in list(indices):
        if idx in ws.holding_remote:
            ws.holding_remote.pop(idx)

    ttp = 10  # times to print; 100 is every 1%
    update = set()
    # force_update = set()
    force_update = True
    remote_finished = 0
    iqsize0 = 0
    oqsize0 = 0
    i = 0
    first = True
    ti = time.monotonic()
    t0 = time.monotonic()
    sleepiness = 0.0

    while (totalwait is not None and waited < totalwait) or (
        ws.holding or ws.iqueue.qsize() or ws.oqueue.qsize()
    ):
        t00 = time.perf_counter()
        # oq = ws.mgr.get_oqueue()
        # print("JOINING GET", oq)
        # try:
        # if len(ws.holding_remote):
        #     time.sleep(.5)
        count = len(indices.intersection(results))
        dt = time.monotonic() - ti
        if i != count:
            force_update = True
        else:
            force_update = False
        # i = len(results)
        i = count
        # print("\nGetting sizes")
        iqsize = ws.iqueue.qsize()
        oqsize = ws.oqueue.qsize()

        if dt < 10.0:
            force_update = False
        else:
            force_update = True
        # force_update = True

        progress = int(i / n * ttp) % ttp
        if i == n:
            progress = ttp
        if force_update or progress not in update:
            update.add(progress)
            ti = time.monotonic()
            riqsize = ws.remote_iqueue_size
            if riqsize is None:
                riqsize = -1
            roqsize = ws.remote_oqueue_size
            if roqsize is None:
                roqsize = -1
            erc = 0
            if ws.finished > 0:
                erc = ws.finished_remote * ws.nproc / ws.finished
            print(
                f"\r{datetime.now()} P: {i/n*100:6.2f}%  {n-i:4d}/{n} "
                f"IQ: {iqsize:4d} OQ: {oqsize:4d} "
                f"IP: {len(ws.holding):4d} "
                f"LF: {ws.finished:4d} "
                f"RF: {ws.finished_remote:4d} "
                f"RIQ: {riqsize:4d} ROQ: {roqsize:4d} "
                f"RIP: {len(ws.holding_remote):4d} ",
                f"ERC: {erc:6.1f} ",
                end="",
            )
            # show the first entry for timing comparsions
            if first:
                print()
                first = False
            # if waited:
            #     print(" (waited) ", end="")
            force_update = False
        if len(indices.difference(results)) == 0:
            break

        # idx, result = oq.get(block=True, timeout=timeout)
        # packets = oq.get(block=False, n=1000)
        # print("\nTrying  to get...")
        packets = queue_get_nowait(oq, timeout=5.0, n=1000)
        if (packets is None or len(packets) == 0) and not (
            ws.holding or iqsize or oqsize
        ):
            if waited >= totalwait:
                print(f"Done waiting")
                break
            sleepiness = 0.0
            time.sleep(waittime)
            waited += waittime
            # print(f"Waited {int(waited)}/{totalwait}")
        elif packets:
            # print("\nReceived packet")
            sleepiness = 0.0
            force_update = True
            # print(f"RECEIVED type {type(packets)} values are \n{pprint.pformat(packets)}")
            for packet in packets:
                # print(f"    packet is type {type(packet)}")
                results.update(packet)
                for idx in packet:
                    at_least_one = True
                    waited = False
                    # print(f"\nUnit {idx} is ready: {result}")

                    # if idx in ws.holding:
                    #     ws.holding.remove(idx)

                    with ws.holding_remote_lock:
                        if idx in ws.holding_remote:
                            ws.holding_remote.pop(idx)
                # ws.holding.remove(idx)
                # except queue.Empty:
                #     pass

        t01 = time.perf_counter()
        dt = t01 - t00
        if dt < 0.0:
            time.sleep(0.0 - dt)
        # print("JOINING GOT", i)

    # print("\nDone flushing.")
    print()
    return results
    print()
    print("DONE: PRINTING")
    for ret in sorted(results, key=lambda x: x[0]):
        print(ret)


def workqueue_new_workspace(
    wq: workqueue_local, address=None, shm=None, nproc=-1
):
    if address is None:
        if configs.remote_compute_enable:
            ip = ""
        else:
            ip = "127.0.0.1"
        port = 0  # wq.mgr.address[1] + 1
    else:
        ip, port = address
        if not configs.remote_compute_enable:
            if ip == "":
                ip = "127.0.0.1"
            else:
                assert ip in ["127.0.0.1", "localhost", "::1"]
    # port = wq.mgr.address[1] + 1
    ws = workspace_local(ip, port, shm=shm, nproc=nproc)

    address = ws.mgr.address
    if address[0] == "0.0.0.0":
        address = ("127.0.0.1", address[1])

    wq.threads[address] = ws

    if configs.remote_compute_enable and ip != "127.0.0.1":
        workqueue_push_workspace(wq, ws)
        # print(f"pushed workspace at {address} values {wq.workspaces}")
    else:
        print(
            f"workspace listening on local host. Remote connections prohibited."
        )

    return ws


def workqueue_remove_workspace(wq: workqueue_local, ws: workspace_local):
    address = ws.mgr.address

    print(f"Removing workspace {address}")
    for addr, t in list(wq.threads.items()):
        if len(addr) == 1:
            continue
        if address[1] == addr[1] and addr in wq.threads:
            wq.threads.pop(addr)
    for addr, status in list(wq.workspaces.items()):
        if len(addr) == 1:
            continue
        if address[1] == addr[1] and addr in wq.workspaces:
            wq.workspaces.pop(addr)
    for addr, status in list(wq.mgr.get_workspaces().items()):
        # print(f"checking {addr}")
        if len(addr) == 1:
            continue
        if address[1] == addr[1]:
            wq.mgr.get_workspaces().pop(addr)
            # print(f"Found workspace. Removed {address}")


def workqueue_get_workspace(wq: workqueue_remote, addr, wq_port) -> workspace:
    print("Getting workspace...")
    ws = None
    wss = wq.get_workspaces()
    wss = list(
        [x for x in wss.items() if type(x[0]) is tuple and len(x[0]) == 2]
    )
    if wss:
        for (address, port), v in sorted(wss, key=lambda x: x[1]):
            if v != 0:
                continue
            ws = workspace_remote(addr, port)
            if ws.is_connected:
                wq.put_workspaces({addr: v + 1})
                print(f"Received workspace {addr}:{port}")
                # now I will have a fully copied shm in the remote ws
                success = workspace_remote_shm_init(ws, timeout=None)
                print(f"Initializing shared memory success {success}")
                if not success:
                    ws = None
                else:
                    print(f"Shared memory is {ws.shm}")
                    print(
                        f"Shared memory has members {list(ws.shm.__dict__.keys())}"
                    )
                    print(
                        f"Requested processors per task is {ws.shm.procs_per_task}"
                    )
            else:
                print(f"Unable to connect workspace")
                ws = None
    else:
        print("No workspaces")

    return ws


def workspace_local_submit(ws, work):
    # iq = ws.iqueue
    ws.iqueue.put(work, block=True)


def workspace_is_active(ws: workspace_remote):
    # print("workspace_is_active: Getting status...")
    s = ws.get_status()
    is_active = s not in [workspace_status.INVALID, workspace_status.DONE]
    # print(f"workspace status: {is_active}")
    return is_active


def workqueue_is_active(wq):
    print("workqueue_is_active: Checking status...")
    s = wq.get_status()
    is_active = s not in [workspace_status.INVALID]
    print(f"Status: {is_active}")
    return is_active


def workqueue_remote_is_active(wq):
    print("workspace_remote_is_active: Checking status")
    s = wq.mgr.get_status()
    return s not in [workspace_status.INVALID]


def workspace_remote_compute(wq: workqueue_remote, ws: workspace_remote):
    print("workspace_remote_compute: starting")
    work = []
    success = True
    # launch a thread that constantly pulls data from iqueue (which is remote)
    # and then we just pull from iqueue
    waits = 0
    waits_max = 6 * 5
    timeout = 10.0
    completed = 0
    processes = 1
    if ws.shm.procs_per_task > 0:
        processes = max(1, ws.nproc // ws.shm.procs_per_task)
    sleepiness = 0.0

    pool = ws.pool
    try:
        print(f"workspace_remote_compute: Starting compute")
        while success:
            t0 = time.perf_counter()
            iqsize = ws.iqueue.qsize()
            oqsize = ws.oqueue.qsize()
            print(
                f"{datetime.now()} Finished: {completed:4d} IQ: {iqsize:4d} OQ: {oqsize:4d} IP: {len(ws.holding):4d} E: {ws.error_count}/{ws.error_limit}",
                end="\n",
            )

            if ws.error_count >= ws.error_limit:
                print("Too many errors, exiting.")
                break

            force_update = False
            functions = {}
            if len(ws.holding) <= processes:
                try:
                    for f in ws.iqueue.get(block=False, n=2):
                        force_update = True
                        if type(f) is dict:
                            functions.update(f)
                        else:
                            print(
                                f"Warning, received a malformed taskset:\n{f}"
                            )
                            assert type(f) is list
                            assert type(f[0]) is dict
                            functions.update(f[0])
                    # functions[idx] = distfun
                except queue.Empty:
                    pass
                if functions is None:
                    functions = {}
                # functions = workspace_remote_get_functions(ws)
                # if functions:
                #     print(f"workspace_remote_compute: {len(functions)} tasks")

            if not (ws.holding or work or iqsize or oqsize):
                if waits < waits_max:
                    waits += 1
                    # print(
                    #     f"Waiting {waits} {len(ws.holding)} {len(work)} {iqsize} {oqsize}"
                    # )
                    # if waits == waits_max - 1:
                    #     ws.stop.set()
                    sleepiness += 2.0
                    time.sleep(sleepiness)
                    if waits % 6 == 0:
                        if not workspace_is_active(ws):
                            break
                else:
                    print("Waited but no work. Bailing.")
                    break
            else:
                waits = 0
                sleepiness = 0.0
                # sleepiness += 2.0
                # time.sleep(sleepiness)
            # if (not success and not work):
            #     break
            success = True
            try:
                # if functions:
                #     print(f"workspace_remote_compute: starting {len(functions)} tasks")
                for idx, distfun in list(functions.items()):
                    # print(f"workspace_remote_compute: starting task {idx}")

                    if idx not in ws.holding:
                        work.append(
                            (
                                idx,
                                pool.apply_async(
                                    workspace_run, (distfun,), {}
                                ),
                            )
                        )
                        ws.holding.add(idx)
                    # print("compute_remote: putting task result to oqueue")
                if work:
                    drop = set()
                    # print(f"Scanning {len(work)} work units")
                    while not drop or len(work) >= processes:
                        drop = set()
                        for i in range(len(work)):
                            idx, unit = work[i]
                            if unit.ready():
                                drop.add(i)
                                result = unit.get()
                                force_update = True
                                # print(f"\nUnit {idx} is ready: {result}")
                                # print(f"\nUnit {idx} is ready")
                                # put_functions(ws.oqueue, [(idx, result)])
                                # if idx == 32:
                                #     breakpoint()
                                if idx in ws.holding:
                                    ws.holding.remove(idx)
                                ws.oqueue.put({idx: result}, block=False)
                                success = True
                                completed += 1
                                # success |= ws.oqueue_put(
                                #     (idx, result)
                                # )  # ws.mgr.get_oqueue()
                                # print(f"Pushing result {idx} success {success}")
                                # ws.oqueue.put((idx, result))
                        working = [
                            x for i, x in enumerate(work) if i not in drop
                        ]

                        work.clear()
                        work.extend(working)

                        working.clear()

                t1 = time.perf_counter()
                if t1 - t0 < 1.0 and not force_update:
                    time.sleep(1 - (t1 - t0))
            except Exception as e:
                print(f"workspace_remote_compute exception: {e}")
                break

    except BrokenPipeError:
        success = False
    print("Leaving workspace...")

    return success


def workspace_compute_local(wq, ws):
    print("workspace_compute_local: starting")
    while workspace_is_active(ws):
        # a function that works inside the workspace, likely pulling shared mem
        # from the workspace
        distfun: distributed_function = workspace_get_function(ws)

        # run the function
        if type(distfun[0]) is not str:
            work = ws.pool.apply_async(*distfun)
            # workspace_run(distfun)
        workspace_push_function(ws, distfun, work)

    workqueue_push_workspace(wq, ws)


def workspace_remote_shm_init_thread(mgr, out):
    try:
        shm_proxy = mgr.get_shm()
        remote_init = shm_proxy.remote_init()
        out.extend([shm_proxy, remote_init])
    except ConnectionError:
        print("workspace_remote_shm_init_thread: ConnectionError")
    except EOFError:
        print("workspace_remote_shm_init_thread: EOFError")
    except TimeoutError:
        print("workspace_remote_shm_init_thread: TimeoutError")
    except AssertionError as e:
        print(f"workspace_remote_shm_init_thread: AssertionError {e}")


def remote_init_thread(remote_init, shm_proxy, out):
    try:
        out.append(remote_init(shm_proxy))
    except ConnectionError:
        print("remote_init_thread: ConnectionError")
    except EOFError:
        print("remote_init_thread: EOFError")
    except TimeoutError:
        print("remote_init_thread: TimeoutError")
    except AssertionError as e:
        print(f"remote_init_thread: AssertionError {e}")


def workspace_remote_shm_init(ws, timeout=TIMEOUT):
    out = []

    t = threading.Thread(
        target=workspace_remote_shm_init_thread, args=(ws.mgr, out)
    )
    t.start()
    t.join(timeout=timeout)

    if out:
        shm_proxy = out[0]
        remote_init = out[1]
        out.clear()
        t = threading.Thread(
            target=remote_init_thread, args=(remote_init, shm_proxy, out)
        )
        t.start()
        t.join(timeout=timeout)
        if out:
            ws.shm = out[0]
            return True

    ws.error_count += 1
    return False


def compute_remote(addr, port, processes=1):
    retry = 0
    # connect to the work workqueue, which will serve workspaces
    retry_n = 240
    wq = None
    success = False
    configs.processes = min(processes, os.cpu_count())
    process.current_process().name = "WS_" + str(int(random.random() * 1e6))
    while retry < retry_n:
        print(f"Connecting to workqueue {addr}:{port}")
        wq = workqueue_remote(addr, port)
        success = wq.is_connected
        if success:
            print(f"Connected to workqueue {addr}:{port}")
        else:
            print(
                f"compute_remote: Failed to connect to workqueue {addr}:{port}. Waiting"
            )
            wq = None

        try:
            if success:
                success = False
                # workspaces hold global/shared space for the functions to operate
                # In this case, I will need to setup the shared memory here first, but
                # perhaps likely do this in the main process
                ws: workspace_remote = workqueue_get_workspace(wq, addr, port)
                # now this is supposed to copy the proxy objects appropriately

                if ws is not None:
                    ws.nproc = processes
                    success = ws.start()

                    # this needs to copy everything
                    if success:
                        print(
                            f"Starting remote compute from {addr} {port} with {processes} processes"
                        )
                        success = workspace_remote_compute(wq, ws)
                        print(f"compute_remote: success is {success}")
                    # if not success:
                    #     retry += 1
                    # else:
                    #     success = workqueue_is_active(wq)
                    ws.close()
        except EOFError:
            print("Disconnected from workspace.")
        except Exception as e:
            print(f"compute_remote: Exception. {e}")

        if success:
            retry = 0
        else:
            retry += 1
            # wq = None
            # ws = None
            print(f"Failed. Retries {retry}/{retry_n}")
            time.sleep(5)

    return retry < retry_n


def workspace_run(distfun: distributed_function):
    global SHM_GLOBAL
    fn = distfun[0]
    args = distfun[1]
    kwargs = distfun[2]
    return fn(*args, **kwargs, shm=SHM_GLOBAL)


def workspace_run_init(procs_per_task, t0=None):
    if procs_per_task > 0:
        configs.processors = min(configs.processors, procs_per_task)
