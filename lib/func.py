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
