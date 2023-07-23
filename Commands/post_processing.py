import abc
import logging
import os
import re
import statistics
from abc import ABC
from pathlib import Path
from typing import List, Optional, Type

import pandas

from lib.const import SupportedFileTypeRegex
from lib.func import get_files

log = logging.getLogger()


def get_structure_files(directory: Path, structure_type: str, id) -> tuple[
    List[CrystalStructure], List[HomologyStructure]]:
    os.chdir(directory)
    crystal_structures: List[CrystalStructure] = []
    homology_modelling: List[HomologyStructure] = []
    for file in os.listdir(directory):
        if str(file).startswith("AF-" + directory.name) and str(file).endswith("." + structure_type) or \
                str(file).startswith("sp") and str(file).endswith("." + structure_type):
            homology_structure = HomologyStructure(directory.joinpath(Path(str(file))), id)
            #        print(homology_structure.path)
            #        if statistics.mean(homology_structure.piddt) ==0:
            #            homology_structure = copy_plddt(homology_structure)
            #  if len(homology_structure.fasta) < 50:
            #      log.info(f"Skipping this structure {homology_structure.path} because length is below threshold of "
            #               f"50 with length {len(homology_structure.fasta)}")
            #      continue
            #        if statistics.mean(homology_structure.piddt) < 85:
            #            log.info(f"Skipping this structure {homology_structure.path} because length is below threshold of "
            #                     f"85 PLDDT with score {statistics.mean(homology_structure.piddt)}")
            #            continue
            homology_modelling.append(homology_structure)
    #   elif str(file).endswith("." + structure_type):
    #       crystal_structure = CrystalStructure(directory.joinpath(Path(str(file))),id)
    #    if len(crystal_structure.fasta) < 50:
    #        log.info(f"Skipping this structure {crystal_structure.path} because length is below threshold of "
    #                 f"50 with length {len(crystal_structure.fasta)}")
    #        continue
    #      crystal_structures.append(crystal_structure)
    return crystal_structures, homology_modelling


def copy_plddt(struct: HomologyStructure):
    log.info("Copying PLDDT")
    with open(str(struct.path).replace("done", "raw_structs"), "r") as ref:
        plddt = []
        for line in ref:
            if line.startswith("ATOM"):  # and line.split()[2] == "C":
                #         if int(line.split()[5]) > current_res or int(line.split()[5])==0:
                plddt.append(float(line.split()[-2]))
    struct._piddt = plddt
    return struct




#    def __init__(self, specific_file: Path) -> None:
#        """#

#        Args:
#            working_directory:
#            all_files:
#            specific_file:
#        """
#        super().__init__()

#        my_tuple: Optional[tuple[List[CrystalStructure], List[HomologyStructure]]] = None
#        if str(specific_file).startswith("AF") and str(specific_file).endswith("." + self.structure_type) or \
#                str(specific_file).startswith("sp") and str(specific_file).endswith("." + self.structure_type) or \
#                str(specific_file).startswith("tr") and str(specific_file).endswith("." + self.structure_type):
#            structure_file: HomologyStructure = HomologyStructure(Path(specific_file),
#                                                                  UniProtID(str(Path(specific_file).parent)))
#            my_tuple = ([], [structure_file])
#            self._structure_results: List[StructureResults] = [
#                StructureResults(UniProtID(str(Path(specific_file).parent), self.working_directory), my_tuple)]
#        elif str(specific_file).endswith("." + self.structure_type):
#            structure_file: CrystalStructure = CrystalStructure(Path(specific_file),
#                                                                UniProtID(str(Path(specific_file).parent)))
#            my_tuple = ([structure_file], [])
#            self._structure_results: List[StructureResults] = [
#                StructureResults(UniProtID(str(Path(specific_file).parent), self.working_directory), my_tuple)]
#        else:
#            self._structure_results: List[StructureResults] = \
#                [StructureResults(UniProtID(str(directories), self.working_directory),
#                                  get_structure_files(self.working_directory.joinpath(str(directories)),
#                                                      self.structure_type, str(directories)))
#                 for directories in os.listdir(self.working_directory)
#                 if self.working_directory.joinpath(str(directories)).is_dir()]
#        self.go_terms = pandas.read_csv("/home/felix/Data/GOTerms.tsv", delimiter="\t").columns.to_numpy()

        #  print(len(uniprot_id.all_structures))
# def remove_low_plddt_regions
