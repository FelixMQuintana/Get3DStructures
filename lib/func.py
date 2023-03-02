import os
from pathlib import Path
from typing import List, Optional

from Commands.Structure import CrystalStructure, HomologyStructure
from lib.const import ALLOWED_EXT


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


def get_structure_files(directory: Path) -> tuple[List[CrystalStructure], List[HomologyStructure]]:
    os.chdir(directory)
    crystal_structures: List[CrystalStructure] = []
    homology_modelling: List[HomologyStructure] = []
    for file in os.listdir(directory):
        if str(file).startswith("AF"):
            homology_modelling.append(HomologyStructure(Path(str(file))))
        elif str(file).startswith(ALLOWED_EXT.PDB.value):
            crystal_structures.append(CrystalStructure(Path(str(file))))
    return crystal_structures, homology_modelling
