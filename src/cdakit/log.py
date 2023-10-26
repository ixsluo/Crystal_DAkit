import logging


def init(level):
    logging.basicConfig(
        format='|%(asctime)s|%(module)s|%(levelname)s| %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=level,
    )

