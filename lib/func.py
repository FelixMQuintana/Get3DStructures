import json
import logging
import os
from pathlib import Path
from typing import Dict


def change_directory(directory: Path):
    """

    :param directory:
    :param skip:
    :return:
    """
    if not directory.exists():
        os.mkdir(directory)
    logging.debug("Changing directory to path %s" % directory)
    os.chdir(directory)
    return True


def load_json(filename: Path) -> Dict:
    """

    Parameters
    ----------
    filename

    Returns
    -------

    """
    try:
        with open(filename, "r") as f:
            data = json.load(f)
    except FileNotFoundError as ex:
        raise ex
    return data