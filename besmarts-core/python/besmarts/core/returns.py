"""
besmarts.core.returns

should be using this with compute!
"""

from enum import Enum
from typing import Generic, TypeVar
import time

T = TypeVar("T")


class return_status(Enum):
    NONE = None
    SUCCESS = 0
    FAIL = 1


class return_value(Generic[T]):
    def __init__(self, value: T, out: str, err: str, status: return_status):
        self.value: T = value
        self.out: str = out
        self.err: str = err
        self.time_ns: int = 0
        self.status: return_status = status


def success(value: T = None, out="", err="", t0=None):
    """
    Indicate a successful function call

    Parameters:
    value : T
        The return value of the underlying function call
    out : str
        The output of the function
    err : str
        The errors of the function

    Returns
    -------
    return_value
        The return object with a status set to success
    """

    ret = return_value(value, out, err, return_status.SUCCESS)

    if t0 is not None:
        ret.time_ns = time.perf_counter_ns() - t0

    return ret


def fail(value: T = None, out="", err="", t0=None):
    """
    Indicate a failed function call

    Parameters:
    value : T
        The return value of the underlying function call
    out : str
        The output of the function
    err : str
        The errors of the function

    Returns
    -------
    return_value
        The return object with a status set to fail
    """
    ret = return_value(value, out, err, return_status.FAIL)

    if t0 is not None:
        ret.time_ns = time.perf_counter_ns() - t0

    return ret
