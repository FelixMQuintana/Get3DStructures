import logging
import os
from pathlib import Path

from pdbfixer import PDBFixer
from openmm.app import PDBFile

from Commands.command import Command
from dataStructure.collections import Collection, HomologyStructureFetcher, ExperimentalStructureFetcher, \
    UniProtAcessionFetcher
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure import StructureFile
from lib.const import AllowedExt


class RepairStructures(Command):

    def __init__(self, repaired_dataset: Path):
        super().__init__()
        self.repaired_dataset: Path = repaired_dataset
        self.collection = Collection(self.working_directory, HomologyStructureFetcher(75),
                                     UniProtAcessionFetcher())
        if not self.repaired_dataset.exists():
            logging.info("Database doesn't appear to exist. Building it now!")
            logging.info("Building directory: %s" % self.repaired_dataset)
            os.mkdir(self.repaired_dataset)
        [(logging.info("Building database for %s." % self.repaired_dataset.joinpath(protein_structures.id)),
          os.mkdir(self.repaired_dataset.joinpath(protein_structures.id)))
         for protein_structures in self.collection.protein_structure_results.values()
         if not self.repaired_dataset.joinpath(protein_structures.id).exists()]

    def run(self) -> None:
        threads = [self.thread_pool.apply_async(self.command, [structure_file, protein_structures]) for
                   protein_structures in
                   self.collection.protein_structure_results.values() for structure_file in
                   protein_structures.all_structures]
        [thread.wait() for thread in threads]

    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures):
        working_dir: Path = self.repaired_dataset.joinpath(structure_file.id)
        if working_dir.joinpath(structure_file.path.name).exists():
            logging.info("File exists %s skipping" % working_dir.joinpath(structure_file.path.name))
            return None
        logging.info(f"Repairing {structure_file.path} as {working_dir.joinpath(structure_file .path.name)}")
        try:
            fixer = PDBFixer(filename=str(structure_file.path))
        except IndexError as IE:
            print(f"Index error {IE}, skipping. Copying over original file")
            os.system(f"cp {structure_file.path} {working_dir.joinpath(structure_file.path)}")
            return None
        logging.debug("Finding Missing Residues %s" % structure_file.path.name)
        fixer.findMissingResidues()
        logging.debug("Finding nonstandard residues %s" % structure_file.path.name)
        fixer.findNonstandardResidues()
        if fixer.nonstandardResidues != {}:
            logging.debug("Fixing Nonstandard Residues %s" % structure_file.path.name)
            fixer.replaceNonstandardResidues()
        logging.debug("Finding Missing Atoms %s" % structure_file.path.name)
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
            logging.debug("Adding Missing Atoms %s" % structure_file.path.name)
            fixer.addMissingAtoms()
        fixer.addMissingHydrogens()
        fixer.removeHeterogens(keepWater=False)
        logging.info("Writing file %s" % working_dir.joinpath(structure_file.path.name))
        pdb_file = open(working_dir.joinpath(structure_file.path.name), "w")
        PDBFile.writeFile(fixer.topology, fixer.positions, pdb_file)
        os.system(f"cp {structure_file.path.parent.joinpath(structure_file.id).with_suffix(AllowedExt.FASTA.value)} "
                  f"{working_dir}")
        os.system(f"cp {structure_file.path.parent.joinpath(structure_file.id).with_suffix(AllowedExt.JSON.value)} "
                  f"{working_dir}")
