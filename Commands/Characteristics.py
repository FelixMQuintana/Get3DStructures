import os

from Commands.Structure import StructureFile
from Commands.post_processing import PostProcessing
from pathlib import Path

from lib.const import StructureCharacteristicsMode
from lib.func import change_directory


class Characteristics(PostProcessing):

    def __init__(self, working_directory: Path, specific_file: Path, mode: StructureCharacteristicsMode):
        super().__init__(working_directory, specific_file)
        if mode == StructureCharacteristicsMode.AUTOSITE.value:
            self.mode = self.finding_ligand_binding_pockets
    def run(self) -> None:
        [self.thread_pool.map(self.mode, structures.all_structures) #self.mode(structure) for structure in structures.all_structures]
         for structures in self._structure_results]

    def finding_ligand_binding_pockets(self, structure_file: StructureFile):
        os.system("prepare_ligand -l %s -o %s" % (structure_file.path, structure_file.path.with_suffix("pdbqt")) )