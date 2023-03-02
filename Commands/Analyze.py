import os
from builtins import function
from pathlib import Path
from typing import List

import matplotlib as matplotlib

from Commands.Structure import CrystalStructure, HomologyStructure
from Commands.command import UniProtID, Command
from lib.func import get_structure_files


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
    def crystal_structure_count(self) -> int:
        return len(self._crystal_structures)

    @property
    def homology_structure_count(self) -> int:
        return len(self._homology_structures)


class Analyze(Command):

    def __init__(self, working_directory: str, all: bool) -> None:
        super().__init__(working_directory)
        if all:
            self._mode: function = self.__check_structures

    def run(self) -> None:
        self._mode()

    def __check_structures(self) -> None:
        structure_results: List[StructureResults] = [
            StructureResults(UniProtID(str(directories)), get_structure_files(Path(str(directories))))
            for directories in os.listdir(self.working_directory)]
        number_of_uniprot_ids = len(structure_results)
        number_of_crystal = sum([structure.crystal_structures for structure in structure_results])
        number_of_homology = sum([structure.homology_structures for structure in structure_results])
        print(f"Number of uniprotIDs {number_of_uniprot_ids}")
        print(f"Number of crystals {number_of_crystal}")
        print(f"Number of homology {number_of_homology}")
        #explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

        #fig, ax = plt.subplots()
        #ax.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
            #   shadow=True, startangle=90)
        #plt.show()
