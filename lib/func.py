import os
from pathlib import Path


def change_directory(directory: Path):
    """

    :param directory:
    :param skip:
    :return:
    """
    if not directory.exists():
        os.mkdir(directory)
    os.chdir(directory)
    return True
