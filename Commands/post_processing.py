import os
from abc import ABC
from pathlib import Path
from typing import List, Optional

from Commands.Structure import CrystalStructure, HomologyStructure, StructureFile
from Commands.command import Command, UniProtID


def get_structure_files(directory: Path, structure_type: str) -> tuple[List[CrystalStructure], List[HomologyStructure]]:
    os.chdir(directory)
    crystal_structures: List[CrystalStructure] = []
    homology_modelling: List[HomologyStructure] = []
    for file in os.listdir(directory):
        if str(file).startswith("AF-" + directory.name) and str(file).endswith("." + structure_type) or \
                str(file).startswith("sp") and str(file).endswith("." + structure_type):
            homology_structure = HomologyStructure(directory.joinpath(Path(str(file))))
            homology_modelling.append(homology_structure)
        elif str(file).endswith("." + structure_type):
            crystal_structures.append(CrystalStructure(directory.joinpath(Path(str(file)))))
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

    def __init__(self, specific_file: Path) -> None:
        """

        Args:
            working_directory:
            all_files:
            specific_file:
        """
        super().__init__()
        my_tuple: Optional[tuple[List[CrystalStructure], List[HomologyStructure]]] = None
        if str(specific_file).startswith("AF") and str(specific_file).endswith("." + self.structure_type) or \
                str(specific_file).startswith("sp") and str(specific_file).endswith("." + self.structure_type):
            structure_file: HomologyStructure = HomologyStructure(Path(specific_file))
            my_tuple = ([], [structure_file])
            self._structure_results: List[StructureResults] = [
                StructureResults(UniProtID(str(Path(specific_file).parent), self.working_directory), my_tuple)]
        elif str(specific_file).endswith("." + self.structure_type):
            structure_file: CrystalStructure = CrystalStructure(Path(specific_file))
            my_tuple = ([structure_file], [])
            self._structure_results: List[StructureResults] = [
                StructureResults(UniProtID(str(Path(specific_file).parent), self.working_directory), my_tuple)]
        else:
            self._structure_results: List[StructureResults] = \
                [StructureResults(UniProtID(str(directories), self.working_directory),
                                  get_structure_files(self.working_directory.joinpath(str(directories)),
                                                      self.structure_type))
                 for directories in os.listdir(self.working_directory)]
