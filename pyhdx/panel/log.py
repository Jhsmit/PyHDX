import logging

#https://stackoverflow.com/questions/7621897/python-logging-module-globally
def setup_custom_logger(name):
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger


def setup_md_log(name, md_window, log_level=logging.DEBUG):
    logger = logging.getLogger(name)
    sh = logging.StreamHandler(md_window)

    formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s')
    sh.setFormatter(formatter)
    sh.setLevel(log_level)
    logger.addHandler(sh)