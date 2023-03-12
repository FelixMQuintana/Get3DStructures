import os
from pathlib import Path
from typing import Optional

from Commands.Structure import StructureFile
from Commands.command import Command
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from Commands.post_processing import PostProcessing
from lib.const import ALLOWED_EXT
from lib.func import change_directory


class RepairPDB(PostProcessing):

    def __init__(self, working_directory: str, all_files: bool, specific_pdb_file: Optional[str], dataset_directory: str):
        super().__init__(working_directory, all_files, specific_pdb_file)
        self.dataset_directory: Path = Path(dataset_directory)
        if not self.dataset_directory.exists():
            os.mkdir(self.dataset_directory)

    def run(self) -> None:
        [[self.repair_pdb(structure, uniprot_id) for structure in uniprot_id.all_structures] for uniprot_id in
         self._structure_results if change_directory(self.working_directory.joinpath(uniprot_id.id), skip=False)]

    def repair_pdb(self, pdb_structure: StructureFile, uniprot_id) -> None:
        working_dir: Path = self.dataset_directory.joinpath(uniprot_id.id)
        if not working_dir.exists():
            os.mkdir(working_dir)
        print(f"Repairing {pdb_structure.path} as {working_dir.joinpath(pdb_structure.path)}")
        try:
            fixer = PDBFixer(filename=str(pdb_structure.path))
        except IndexError as IE:
            print(f"Index error {IE}, skipping. Copying over original file")
            os.system(f"cp {pdb_structure.path} {working_dir.joinpath(pdb_structure.path)}")
            return None
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        #        fixer.removeHeterogens(True)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        PDBFile.writeFile(fixer.topology, fixer.positions, open(working_dir.joinpath(pdb_structure.path), "w"))
