import logging
import os
import subprocess
import threading
from typing import List

from Commands.Structure import StructureFile
from Commands.post_processing import PostProcessing
from pathlib import Path

from lib.const import StructureCharacteristicsMode


class Characteristics(PostProcessing):

    def __init__(self, specific_file: Path, mode: StructureCharacteristicsMode, binding_site_database: Path):
        super().__init__(specific_file)
        if mode == StructureCharacteristicsMode.AUTOSITE:
            self.mode = self.finding_ligand_binding_pockets
        self.autosite_location = self.args[StructureCharacteristicsMode.AUTOSITE.value]
        self.binding_site_database: Path = binding_site_database
        if not self.binding_site_database.exists():
            logging.info("Database doesn't appear to exist. Building it now!")
            logging.info("Building directory: %s" % self.binding_site_database)
            os.mkdir(self.binding_site_database)
        [(logging.info("Building database for %s." % self.binding_site_database.joinpath(structure_result.id)),
          os.mkdir(self.binding_site_database.joinpath(structure_result.id)))
         for structure_result in self._structure_results
         if not self.binding_site_database.joinpath(structure_result.id).exists()]

    def run(self) -> None:
        threads: List = [
            [self.thread_pool.apply_async(self.mode, [structures, structure_result.id])
             for structures in structure_result.all_structures]
            for structure_result in self._structure_results]
        [[thread.wait() for thread in thread_list] for thread_list in threads]

    def finding_ligand_binding_pockets(self, structure_file: StructureFile, id: str):

        with threading.Lock():
            logging.info("Running the command:  %s/prepare_receptor -r %s -o %s"
                         % (self.autosite_location, structure_file.path, structure_file.path.name.split(".")[0] + ".pdbqt"))
            p1 = subprocess.Popen([self.autosite_location + "prepare_receptor -r" + structure_file.path.as_posix() + " -o " + structure_file.path.name.split(".")[0] + ".pdbqt" ],
                                  shell=True)
            p1.wait()
            logging.info("Running the command: %sautosite -r %s -o %s" % (self.autosite_location,
                                                                          structure_file.path.name.split(".")[
                                                                              0] + ".pdbqt",
                                                                          self.binding_site_database.joinpath(
                                                                              id).joinpath(
                                                                              structure_file.path.name.split(".")[0])))
            if not self.binding_site_database.joinpath(id).joinpath(structure_file.path.name.split(".")[0]).exists():
                os.mkdir(self.binding_site_database.joinpath(id).joinpath(structure_file.path.name.split(".")[0]))
            os.system("%s/autosite -r %s -o %s" % (self.autosite_location,
                                                   structure_file.path.name.split(".")[0] + ".pdbqt",
                                                   self.binding_site_database.joinpath(id).joinpath(
                                                       structure_file.path.name.split(".")[0])))
