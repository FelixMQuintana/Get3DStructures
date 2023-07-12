import logging
import os
import statistics
from abc import ABC
from pathlib import Path
from typing import List, Optional

from Commands.Structure import CrystalStructure, HomologyStructure, StructureFile
from Commands.command import Command, UniProtID
import pandas

log = logging.getLogger()

def get_structure_files(directory: Path, structure_type: str, id) -> tuple[List[CrystalStructure], List[HomologyStructure]]:
    os.chdir(directory)
    crystal_structures: List[CrystalStructure] = []
    homology_modelling: List[HomologyStructure] = []
    for file in os.listdir(directory):
        if str(file).startswith("AF-" + directory.name) and str(file).endswith("." + structure_type) or \
                str(file).startswith("sp") and str(file).endswith("." + structure_type):
            homology_structure = HomologyStructure(directory.joinpath(Path(str(file))),id)
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
    with open(str( struct.path).replace("done","raw_structs"),"r") as ref:
        plddt = []
        for line in ref:
            if line.startswith("ATOM"):  # and line.split()[2] == "C":
            #         if int(line.split()[5]) > current_res or int(line.split()[5])==0:
                plddt.append(float(line.split()[-2]))
    struct._piddt = plddt
    return struct

class StructureResults:

    def __init__(self, uniprot: UniProtID, structures: tuple[List[CrystalStructure], List[HomologyStructure]]) -> None:
        self._accession: UniProtID = uniprot
        self._crystal_structures: List[CrystalStructure] = structures[0]
        self._homology_structures: List[HomologyStructure] = structures[1]

    @property
    def id(self) -> str:
        return self._accession.id

    @property
    def fasta(self) -> str:
        return self._accession.fasta

    @property
    def crystal_structures(self) -> List[CrystalStructure]:
        return self._crystal_structures

    @property
    def homology_structures(self) -> List[HomologyStructure]:
        return self._homology_structures

    @property
    def all_structures(self) -> List[StructureFile]:
        return [*self._crystal_structures, *self._homology_structures]

    @property
    def crystal_structure_count(self) -> int:
        return len(self._crystal_structures)

    @property
    def homology_structure_count(self) -> int:
        return len(self._homology_structures)

    @property
    def homology_fastas(self) -> List[str]:
        return [struct.fasta for struct in self.homology_structures]

    @property
    def crystal_fastas(self) -> List[str]:
        return [struct.fasta for struct in self.crystal_structures]

    @property
    def uniprotID(self, ) ->UniProtID:
        """

        Returns
        -------

        """
        return self._accession

class PostProcessing(Command, ABC):

    def __init__(self, specific_file: Path) -> None:
        """

        Args:
            working_directory:
            all_files:
            specific_file:
        """
        super().__init__()
        my_tuple: Optional[tuple[List[CrystalStructure], List[HomologyStructure]]] = None
        if str(specific_file).startswith("AF") and str(specific_file).endswith("." + self.structure_type) or \
                str(specific_file).startswith("sp") and str(specific_file).endswith("." + self.structure_type) or \
                str(specific_file).startswith("tr") and str(specific_file).endswith("." + self.structure_type) :
            structure_file: HomologyStructure = HomologyStructure(Path(specific_file), UniProtID(str(Path(specific_file).parent)))
            my_tuple = ([], [structure_file])
            self._structure_results: List[StructureResults] = [
                StructureResults(UniProtID(str(Path(specific_file).parent), self.working_directory), my_tuple)]
        elif str(specific_file).endswith("." + self.structure_type):
            structure_file: CrystalStructure = CrystalStructure(Path(specific_file), UniProtID(str(Path(specific_file).parent)))
            my_tuple = ([structure_file], [])
            self._structure_results: List[StructureResults] = [
                StructureResults(UniProtID(str(Path(specific_file).parent), self.working_directory), my_tuple)]
        else:
            self._structure_results: List[StructureResults] = \
                [StructureResults(UniProtID(str(directories), self.working_directory),
                                  get_structure_files(self.working_directory.joinpath(str(directories)),
                                                      self.structure_type, str(directories)))
                 for directories in os.listdir(self.working_directory)
                 if self.working_directory.joinpath(str(directories)).is_dir()]
        self.go_terms = pandas.read_csv("/home/felix/Data/GOTerms.tsv", delimiter="\t").columns.to_numpy()


                  #  print(len(uniprot_id.all_structures))
   # def remove_low_plddt_regions
