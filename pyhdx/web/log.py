import datetime
import logging
import sys


def get_default_handler(stream=None):
    sh = logging.StreamHandler(stream)
    formatter = logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s")
    sh.setFormatter(formatter)

    return sh


# https://stackoverflow.com/questions/7621897/python-logging-module-globally
def setup_custom_logger(name):
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - %(module)s - %(message)s"
    )

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger


def setup_md_log(name, md_window, log_level=logging.DEBUG):
    logger = logging.getLogger(name)
    sh = logging.StreamHandler(md_window)

    formatter = logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s")
    sh.setFormatter(formatter)
    sh.setLevel(log_level)
    logger.addHandler(sh)


class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    https://stackoverflow.com/questions/19425736/how-to-redirect-stdout-and-stderr-to-logger-in-python
    """

    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self.linebuf = ""

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass


def logger(root_name):
    def decorator(function):
        def wrapper(*args, **kwargs):
            wrapper.calls += 1
            dt = datetime.datetime.now().strftime("%Y%m%d")
            logger = logging.getLogger(
                f"{root_name}.{function.__name__}.{dt}_{wrapper.calls}"
            )
            logger.setLevel(logging.DEBUG)
            sys.stderr = StreamToLogger(logger, logging.DEBUG)
            wrapper.logger = logger
            return function(*args, **kwargs)

        wrapper.calls = 0
        return wrapper

    return decorator
