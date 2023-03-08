from pathlib import Path
from typing import Optional

from Commands.Structure import StructureFile
from Commands.command import Command
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from Commands.post_processing import PostProcessing


class RepairPDB(PostProcessing):

    def __init__(self, working_directory: str, all_files: bool, specific_pdb_file: Optional[str]):
        super().__init__(working_directory, all_files, specific_pdb_file)

    def run(self) -> None:
        [[self.repair_pdb(structure) for structure in uniprot_id.all_structures] for uniprot_id in
         self._structure_results]

    def repair_pdb(self, pdb_structure: StructureFile) -> None:
        fixer = PDBFixer(filename=pdb_structure.path.as_uri())
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        #        fixer.removeHeterogens(True)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_structure.path.as_uri(), "w"))
