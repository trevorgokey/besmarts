"""
besmarts.core.returns

should be using this with compute!
"""

from enum import Enum
from typing import Generic, TypeVar

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
        self.status: return_status = status


def success(value: T = None, out="", err=""):
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

    return return_value(value, out, err, return_status.SUCCESS)


def fail(value: T = None, out="", err=""):
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
    return return_value(value, out, err, return_status.FAIL)
