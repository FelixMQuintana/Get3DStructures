import logging
import os
import threading
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

    def __init__(self, specific_pdb_file: Optional[Path], dataset_directory: Path):
        super().__init__(specific_pdb_file)
        self.dataset_directory: Path = Path(dataset_directory)
        if not self.dataset_directory.exists():
            logging.info("Database doesn't appear to exist. Building it now!")
            logging.info("Building directory: %s" % self.dataset_directory)
            os.mkdir(self.dataset_directory)
        [(logging.info("Building database for %s." % self.dataset_directory.joinpath(structure_result.id)),
          os.mkdir(self.dataset_directory.joinpath(structure_result.id)))
         for structure_result in self._structure_results
         if not self.dataset_directory.joinpath(structure_result.id).exists()]
        self.active_threads =[]
    def run(self) -> None:
     #   threads = [[self.thread_pool.starmap(self.repair_pdb, [structure, uniprot_id]) for structure in
      #  threads = [[self.thread_pool.apply_async(self.repair_pdb, [structure, uniprot_id]) for structure in
        threads = [[self.repair_pdb(structure, uniprot_id) for structure in
      #  threads = [[threading.Thread(target= self.repair_pdb, args=[structure, uniprot_id]) for structure in
                    uniprot_id.all_structures] for uniprot_id in
                   self._structure_results]

        [[thread.wait() for thread in list_of_threads] for list_of_threads in self.active_threads]
    #    [[thread.wait() for thread in list_of_threads] for list_of_threads in threads]
    #    [[thread.start() for thread in list_of_threads] for list_of_threads in threads]
    #    [[thread.join() for thread in list_of_threads] for list_of_threads in threads]
    def repair_pdb(self, pdb_structure: StructureFile, uniprot_id) :
        working_dir: Path = self.dataset_directory.joinpath(uniprot_id.id)
        logging.info(f"Repairing {pdb_structure.path} as {working_dir.joinpath(pdb_structure.path.name)}")
        try:
            fixer = PDBFixer(filename=str(pdb_structure.path))
        except IndexError as IE:
            print(f"Index error {IE}, skipping. Copying over original file")
            os.system(f"cp {pdb_structure.path} {working_dir.joinpath(pdb_structure.path)}")
            return None
        logging.debug("Finding Missing Residues %s" % pdb_structure.path.name)
        fixer.findMissingResidues()
        logging.debug("Finding nonstandard residues %s" % pdb_structure.path.name)
        fixer.findNonstandardResidues()
        logging.debug("Fixing Nonstandard Residues %s" % pdb_structure.path.name)
        fixer.replaceNonstandardResidues()
        #        fixer.removeHeterogens(True)
        logging.debug("Finding Missing Atoms %s" % pdb_structure.path.name)
        fixer.findMissingAtoms()
        self.active_threads.append(self.thread_pool.apply_async(self.fuckme, [fixer, pdb_structure, working_dir]))
    def fuckme(self, fixer, pdb_structure, working_dir):
        logging.debug("Adding Missing Atoms %s" % pdb_structure.path.name)
        fixer.addMissingAtoms()
        logging.info("Writing file %s" % working_dir.joinpath(pdb_structure.path.name))
        PDBFile.writeFile(fixer.topology, fixer.positions, open(working_dir.joinpath(pdb_structure.path.name), "w"))

