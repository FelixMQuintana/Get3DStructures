import json
import logging
import os
from pathlib import Path
from typing import Dict, List


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


def get_files(file_path: Path, regular_expression: str) -> List[Path]:
    return list(file_path.rglob(regular_expression))
 #   file_path.glob('**/*')
 #   files: List[Path] = []
 #   for file in os.listdir(file_path):

#    return [Path(str(file)) for file in os.listdir(file_path)]
