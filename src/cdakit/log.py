import logging
import logging.handlers
from functools import wraps
from multiprocessing import Queue
from typing import Optional


def listener_configurer():
    root = logging.getLogger()
    h = logging.StreamHandler()
    f = logging.Formatter('|%(asctime)s|%(process)d|%(module)s|%(levelname)-8s| %(message)s', '%y/%m/%d %H:%M:%S')
    h.setFormatter(f)
    root.addHandler(h)


def loglistener(queue: Queue):
    listener_configurer()
    while True:
        try:
            record = queue.get()
            if record is None:  # We send this as a sentinel to tell the listener to quit.
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)  # No level or filter logic applied - just do it!
        except Exception:
            import sys, traceback
            print('Whoops! Problem:', file=sys.stderr)
            traceback.print_exc(file=sys.stderr)


def worker_configurer(queue: Optional[Queue] = None, level=20):
    if queue is None:
        h = logging.StreamHandler()
    else:
        h = logging.handlers.QueueHandler(queue)  # Just the one handler needed
    root = logging.getLogger()
    root.addHandler(h)
    root.setLevel(level)


class logit:
    def __init__(self, level=None):
        self.level = level

    def __call__(self, func):
        @wraps(func)
        def wrapper(logqueue: Optional[Queue] = None, logconf: dict = {}, *args, **kwargs):
            if self.level is not None:
                logconf["level"] = self.level
            worker_configurer(logqueue, **logconf)
            return func(*args, **kwargs)
        return wrapper

