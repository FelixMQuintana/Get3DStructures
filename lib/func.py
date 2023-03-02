import os
from pathlib import Path


def change_directory(directory: Path, skip=True):
    """

    :param directory:
    :param skip:
    :return:
    """
    if skip:
        return False
    if not directory.exists():
        os.mkdir(directory)
    os.chdir(directory)
    return True

def count_structures(directory: Path) -> tuple[int, int]:
    os.chdir(directory)
    crystal_structures = 0
    homology_modelling = 0
    for file in os.listdir(directory):
        if str(file).startswith("AF"):
            