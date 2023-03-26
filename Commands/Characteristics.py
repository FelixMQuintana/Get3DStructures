import logging
import os
import subprocess
import threading
from typing import List

from Commands.Structure import StructureFile
from Commands.post_processing import PostProcessing, StructureResults
from pathlib import Path

from lib.const import StructureCharacteristicsMode
from lib.func import change_directory


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
            [self.mode(structures, structure_result.id)
             # [self.thread_pool.apply_async(self.mode, [structures, structure_result.id])
             for structures in structure_result.all_structures]
            # self.thread_pool.apply_async(self.mode, structures.all_structures)
            # self.mode(structure) for structure in structures.all_structures]
            for structure_result in self._structure_results]
        [[thread.wait() for thread in thread_list] for thread_list in threads]
        # [result.wait() for result in results]

    def finding_ligand_binding_pockets(self, structure_file: StructureFile, id: str):



        return self.thread_pool.apply_async(self.get_binding_site, [structure_file, id])

        # p1 = subprocess.Popen("%s/prepare_ligand -l %s" % (self.autosite_location, structure_file.path.name))

    def get_binding_site(self, structure_file: StructureFile, id: str ):
       # p1 = self.thread_pool.apply_async(lambda autosite_location, structure_file_name:
       #                              os.system("%sprepare_ligand -l %s" % (autosite_location, structure_file_name)),
       #                              [self.autosite_location, structure_file.path.name])
        with threading.Lock():
          #  logging.info(
          #      "Changing working directory %s to structure_file's to run autosite program." % structure_file.path.parent)
          #  os.chdir(structure_file.path.parent)

          #  if not structure_file.path.parent.joinpath(structure_file.path.name.split(".")[0]).exists():
          #      logging.info(
          #          "Making temporary sub-directory for prepare ligand script called %s" % structure_file.path.parent.joinpath(
          #              structure_file.path.name.split(".")[0]))
          #      os.mkdir(structure_file.path.parent.joinpath(structure_file.path.name.split(".")[0]))
          #  logging.info("Changing working directory to %s" % structure_file.path.parent.joinpath(
          #      structure_file.path.name.split(".")[0]))
          #  os.chdir(structure_file.path.parent.joinpath(structure_file.path.name.split(".")[0]))
            #logging.info("Copying %s to working directory." % structure_file.path)
            #os.system("cp %s ./" % structure_file.path)
            logging.info("Running the command:  %s/prepare_ligand -l %s"
                    % (self.autosite_location, structure_file.path))
           # os.system("%sprepare_ligand -l %s" % (self.autosite_location, structure_file.path.name))
            p1 = subprocess.Popen([self.autosite_location+"prepare_ligand -l" + structure_file.path.as_posix()], shell=True)
            p1.wait()
            logging.info("Running the command: %sautosite -r %s -o %s" % (self.autosite_location,
                                                                      structure_file.path.name.split(".")[0] + ".pdbqt",
                                                                      self.binding_site_database.joinpath(id)))
        # p2 = subprocess.Popen("%s/autosite -r %s -o %s" % (self.autosite_location,
        #                                           structure_file.path.name.split(".")[0] + ".pdbqt",
        #                                           self.binding_site_database.joinpath(id)))
        #p2= self.thread_pool.apply_async(lambda autosite_location, structure_file_name, database: os.system(
        #    "%s/autosite -r %s -o %s" % (autosite_location,
        #                                 structure_file_name,
        #                                 database)), [self.autosite_location,
        #                                              structure_file.path.name.split(".")[0] + ".pdbqt",
        #                                              self.binding_site_database.joinpath(id)])
            os.system("%s/autosite -r %s -o %s" % (self.autosite_location,
                                               structure_file.path.name.split(".")[0] + ".pdbqt",
                                               self.binding_site_database.joinpath(id)))