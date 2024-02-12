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

def dprint(*args, **kwargs):
    on = kwargs.get("on", False)
    if "on" in kwargs:
        kwargs.pop("on")
    on = False
    if on:
        print(timestamp(), *args, **kwargs)
