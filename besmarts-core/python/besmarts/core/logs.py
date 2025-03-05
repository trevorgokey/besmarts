"""
besmarts.core.logs

Configure logging and console output
"""

import logging
import datetime

def cout_logger(name):
    return logging.getLogger(name)

def file_logger(name):
    return logging.getLogger(name)

def timestamp():
    return datetime.datetime.now()

def append(line, out, verbose=False):
    out.append(line)
    if verbose:
        print(line)

def dprint(*args, **kwargs):
    on = kwargs.get("on", False)
    if "on" in kwargs:
        kwargs.pop("on")
    if on:
        print(timestamp(), *args, **kwargs)
