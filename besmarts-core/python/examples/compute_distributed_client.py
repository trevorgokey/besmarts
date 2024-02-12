from besmarts.core import compute

import sys
import os

from compute_distributed_server import construct
from compute_distributed_server import construct_idx

ip = '127.0.0.1'
pt = 55555
processes = 14

if len(sys.argv) > 1:
    ip = sys.argv[1]
if len(sys.argv) > 2:
    pt = int(sys.argv[2])
if len(sys.argv) > 3:
    processes=int(sys.argv[3])

if __name__ == "__main__":
    try:
        compute.compute_remote(ip, pt, processes=processes)
        print("Exiting normally.")
        sys.exit(0)
    except EOFError:
        sys.exit("Disconnected")
    except ConnectionRefusedError:
        sys.exit("Connection refused.")
    finally:
        os._exit(0)

