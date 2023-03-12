import os
from abc import ABC
from pathlib import Path
from typing import List, Optional

from Commands.Structure import CrystalStructure, HomologyStructure, StructureFile
from Commands.command import Command, UniProtID
from lib.const import ALLOWED_EXT


def read_af_piddt(alpha_fold_structure_path: HomologyStructure) -> List:
    read_alpha_fold_file_obj = open(alpha_fold_structure_path.path, "r")
    current_res = 0
    plddt = []
    for line in read_alpha_fold_file_obj:
        if line.startswith("ATOM") and line.split()[2] == "C":
            if int(line.split()[8]) > current_res:
                plddt.append(float(line.split()[14]))
                current_res = int(line.split()[8])
    return plddt


def get_structure_files(directory: Path) -> tuple[List[CrystalStructure], List[HomologyStructure]]:
    os.chdir(directory)
    crystal_structures: List[CrystalStructure] = []
    homology_modelling: List[HomologyStructure] = []
    for file in os.listdir(directory):
        if str(file).startswith("AF-" + directory.name) and str(file).endswith(ALLOWED_EXT.CIF.value) or \
                str(file).startswith("sp") and str(file).endswith(ALLOWED_EXT.CIF.value):
            homology_structure = HomologyStructure(Path(str(file)))
            homology_structure.piddt = read_af_piddt(homology_structure)
            homology_modelling.append(homology_structure)
        elif str(file).endswith(ALLOWED_EXT.CIF.value):
            crystal_structures.append(CrystalStructure(Path(str(file))))
    return crystal_structures, homology_modelling


class StructureResults:

    def __init__(self, uniprot: UniProtID, structures: tuple[List[CrystalStructure], List[HomologyStructure]]) -> None:
        self._accession: UniProtID = uniprot
        self._crystal_structures: List[CrystalStructure] = structures[0]
        self._homology_structures: List[HomologyStructure] = structures[1]

    @property
    def id(self) -> str:
        return self._accession.id

    @property
    def crystal_structures(self) -> List[CrystalStructure]:
        return self._crystal_structures

    @property
    def homology_structures(self) -> List[HomologyStructure]:
        return self._homology_structures

    @property
    def all_structures(self) -> List[StructureFile]:
        return [*self._crystal_structures, *self._homology_structures]

    @property
    def crystal_structure_count(self) -> int:
        return len(self._crystal_structures)

    @property
    def homology_structure_count(self) -> int:
        return len(self._homology_structures)


class PostProcessing(Command, ABC):

    def __init__(self, working_directory: str, all_files: bool, specific_file: str) -> None:
        """

        Args:
            working_directory:
            all_files:
            specific_file:
        """
        super().__init__(working_directory)
        my_tuple: Optional[tuple[List[CrystalStructure], List[HomologyStructure]]] = None
        if str(specific_file).startswith("AF") and str(specific_file).endswith(ALLOWED_EXT.CIF.value) or \
                str(specific_file).startswith("sp") and str(specific_file).endswith(ALLOWED_EXT.CIF.value):
            structure_file: HomologyStructure = HomologyStructure(Path(specific_file))
            my_tuple = ([], [structure_file])
        elif str(specific_file).endswith(ALLOWED_EXT.CIF.value):
            structure_file: CrystalStructure = CrystalStructure(Path(specific_file))
            my_tuple = ([structure_file], [])
        if all_files:
            self._structure_results: List[StructureResults] = \
                [StructureResults(UniProtID(str(directories)),
                                  get_structure_files(self.working_directory.joinpath(str(directories))))
                 for directories in os.listdir(self.working_directory)]
        else:
            self._structure_results: List[StructureResults] = [
                StructureResults(UniProtID(str(Path(specific_file).parent)), my_tuple)]
