import logging
import os
import threading
from pathlib import Path
from typing import Optional, List

from Commands.Structure import StructureFile
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from Commands.post_processing import PostProcessing

from lib.const import ALLOWED_EXT


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
        self.active_threads = []

    def exception_handler(self, args):
        print(f'Thread failed: {args.exc_value}')

    def run(self) -> None:
        threading.excepthook = self.exception_handler

        threads = [[self.thread_pool.apply_async(self.repair_pdb, [structure, uniprot_id]) for structure in
                    uniprot_id.all_structures] for uniprot_id in
                   self._structure_results]
      #  [self.repair_pdb(structure, uniprot_id)
      #                                    for uniprot_id in self._structure_results
      #                                    for structure in uniprot_id.all_structures]
      #  threads: List[threading.Thread]= [threading.Thread (target=self.repair_pdb, args=[structure, uniprot_id],)
           #                               for uniprot_id in self._structure_results
          #                                for structure in uniprot_id.all_structures]
     #   [thread.start() for thread in threads]
     #   [thread.join() for thread in threads]
        [[thread.wait() for thread in list_of_threads] for list_of_threads in threads]

    def repair_pdb(self, pdb_structure: StructureFile, uniprot_id):
        working_dir: Path = self.dataset_directory.joinpath(uniprot_id.id)
        if working_dir.joinpath(pdb_structure.path.name).exists():
            logging.info("File exists %s skipping" % working_dir.joinpath(pdb_structure.path.name))
            return None
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
        if fixer.nonstandardResidues != {}:
            logging.debug("Fixing Nonstandard Residues %s" % pdb_structure.path.name)
            fixer.replaceNonstandardResidues()
        logging.debug("Finding Missing Atoms %s" % pdb_structure.path.name)
        fixer.findMissingAtoms()
        if fixer.missingAtoms != {}:
            chains = list(fixer.topology.chains())
            keys = list(fixer.missingResidues.keys())
            for key in keys:
                chain = chains[key[0]]
                if key[1] == 0 or key[1] == len(list(chain.residues())):
                    fixer.missingResidues[key] = []
                else:
                    if len(fixer.missingResidues[key]) > 10:
                        fixer.missingResidues[key] = []
            logging.debug("Adding Missing Atoms %s" % pdb_structure.path.name)
            print(fixer.missingAtoms)
            fixer.addMissingAtoms()
        fixer.addMissingHydrogens()
        fixer.removeHeterogens(keepWater=False)
        logging.info("Writing file %s" % working_dir.joinpath(pdb_structure.path.name))
        pdb_file = open(working_dir.joinpath(pdb_structure.path.name), "w")
        PDBFile.writeFile(fixer.topology, fixer.positions, pdb_file)
        os.system(f"cp {pdb_structure.path.parent.joinpath(uniprot_id.id).with_suffix(ALLOWED_EXT.FASTA.value)} "
                  f"{working_dir}")
        pdb_file.close()
        del fixer
