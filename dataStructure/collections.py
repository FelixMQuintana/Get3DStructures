"""

"""
import abc
import logging
import os
import re
from abc import ABC, abstractmethod
from pathlib import Path
from statistics import mean
from typing import List, Type, Dict

from tqdm import tqdm

from dataStructure.protein.accession import UniProtIDFastaFile, UniProtAcessionFile
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure.structure import HomologyStructure, GOLabelsFile, CrystalStructure, StructureFile, \
    SupportedFileType, MeshFile
from lib.const import SupportedFileTypeRegex, AllowedExt, METRIC

log = logging.getLogger()


class DataFetcher(ABC):

    def __init__(self, optional_custom_regex: str = None):
        self.optional_custom_regex = optional_custom_regex

    @property
    @abstractmethod
    def file_type_regex(self) -> SupportedFileTypeRegex:
        """

        Returns A Supported Structure Type to specify how to filter structures
        -------

        """
        raise NotImplementedError

    @property
    @abstractmethod
    def file_type(self) -> Type[SupportedFileType]:
        """

        Returns
        -------

        """
        raise NotImplementedError

    def get_data(self, file_path: Path) -> List[SupportedFileType]:
        """

        Parameters
        ----------
        file_path: Directory to search for StructureFiles.

        Returns A list of all StructureFiles found in given directory
        -------

        """
        if self.optional_custom_regex is not None:
            return [self.file_type(file) for file in file_path.rglob(self.optional_custom_regex)]
        else:
            return [self.file_type(file) for file in file_path.rglob(self.file_type_regex.value)]


class ExperimentalStructureFetcher(DataFetcher):

    @property
    def file_type(self) -> Type[CrystalStructure]:
        return CrystalStructure

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.EXPERIMENTAL_STRUCTURE


class MeshFetcher(DataFetcher):

    @property
    def file_type(self) -> Type[SupportedFileType]:
        return MeshFile

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.PLY_FILE


class AtomFetcher(DataFetcher):

    @property
    def file_type(self):
        return StructureFile

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.ATOM_STRUCTURE


class HomologyStructureFetcher(DataFetcher):

    def __init__(self, optional_custom_regex: str = None, quality_cut_off: int = 0, ):
        """

        Parameters
        ----------
        quality_cut_off
        """
        super().__init__(optional_custom_regex)
        self._quality_cut_off = quality_cut_off

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.HOMOLOGY_STRUCTURE

    @property
    def file_type(self) -> Type[HomologyStructure]:
        return HomologyStructure

    def get_data(self, file_path: Path) -> List[SupportedFileType]:
        if self.optional_custom_regex is not None:
            homology_files = [self.file_type(file) for file in file_path.rglob(self.optional_custom_regex)]
        else:
            homology_files = [self.file_type(file) for file in file_path.rglob(self.file_type_regex.value)]
        return_list = []
        for file in tqdm(homology_files,desc="Parsing through files"):
           # print(file.path)
           # if mean(file.piddt) == 0:
           #     file = self.copy_plddt(file)
           # if mean(file.piddt) < self._quality_cut_off:
           #     log.info(f"Skipping {file.path} since this file's mean plddt score {mean(file.piddt)} is less than 90.")
           #     continue
            return_list.append(file)
        return return_list

    @staticmethod
    def copy_plddt(struct: HomologyStructure):
        log.info("Copying PLDDT")
        with open(str(struct.path).replace("Repaired_structs", "raw_structs"), "r") as ref:
            plddt = []
            current_res_num = 0
            for line in ref:
                if line.startswith("ATOM"):
                    if current_res_num < 1000 and line.split()[5].isnumeric():
                        if int(line.split()[5]) > current_res_num:  # and line.split()[2] == "C":
                            plddt.append(float(line.split()[-2]))
                            current_res_num = int(line.split()[5])
                    else:
                        if int(line.split()[4][1:]) > current_res_num:
                            plddt.append(float(line.split()[-2]))
                            current_res_num = int(line.split()[4][1:])
            struct._piddt = plddt
        return struct


class GOTermFetcher(DataFetcher):

    @property
    def file_type(self) -> Type[SupportedFileType]:
        return GOLabelsFile

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.CSV_FILE


class UniProtFastaFetcher(DataFetcher):

    @property
    def file_type(self) -> Type[SupportedFileType]:
        return UniProtIDFastaFile

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.FASTA_FILE


class UniProtAcessionFetcher(DataFetcher):

    @property
    def file_type(self) -> Type[SupportedFileType]:
        return UniProtAcessionFile

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.JSON_FILE


class Collection:

    def __init__(self, working_directory: Path, *args: DataFetcher):
        self.working_directory: Path = working_directory
        self._results = {}
        self._protein_structures: Dict[ProteinStructures] = {}
        self.add_fetchers(*args)

    def add_fetchers(self, *args: DataFetcher):

        collection = []
        for data_fetcher in args:
            collection.extend(data_fetcher.get_data(self.working_directory))
        for file in collection:
            if isinstance(file, StructureFile):
                if isinstance(file, UniProtAcessionFile) and file.path.name != (file.id + AllowedExt.JSON.value):
                    continue
                if self._results.get(file.id) is None:
                    self._results[file.id] = [file]
                else:
                    self._results[file.id].append(file)
            else:
                self._results[file.id] = [file]
        for key, value in self._results.items():
            self._protein_structures[key] = ProteinStructures(value)

    @property
    def protein_structure_results(self):  # -> Dict[ProteinStructures]:
        return self._protein_structures

    def remove_redundant_proteins(self):
        unique_sequences = set()
        ids_to_remove = []
        for proteins in self.protein_structure_results.values():
           # proteins_with_non_empty_annotation = [protein for protein in proteins.all_structures if protein.annotation]
            # Add non-empty annotation proteins to the set of unique sequences
          #  for protein in proteins_with_non_empty_annotation:
            if proteins.annotations is not None:
                if proteins.fasta_file.fasta not in unique_sequences:
                    unique_sequences.add(proteins.fasta_file.fasta)
                else:
                    ids_to_remove.append(proteins.id)
        # Remove redundant proteins, prioritizing those with empty annotations
        for id in ids_to_remove:
            self._protein_structures.pop(id)
        ids_to_remove = []
        # Additionally, remove any proteins that were not previously removed but have a redundant sequence
        for proteins in self.protein_structure_results.values():
            if proteins.fasta_file.fasta not in unique_sequences:
                unique_sequences.add(proteins.fasta_file.fasta)
            else:
                ids_to_remove.append(proteins.id)
        # Remove redundant proteins
        for id in ids_to_remove:
            self._protein_structures.pop(id)

    #TODO this is not the way to do this properly but until I completely refactor the code, this will do
    def write_collection_to_db(self, db_path:Path) -> None:
        for proteins in self.protein_structure_results.values():
            os.mkdir(db_path.joinpath(proteins.id))
            for file in proteins.all_files:
                os.system(f"cp {file.path} {db_path.joinpath(proteins.id)}")


class AnalyzeCollection:

    @staticmethod
    def get_counts_of_annotation(collection, annotation) -> int:
        count = 0
        for proteins in collection.protein_structure_results.values():
            if proteins.annotations is not None and annotation in proteins.annotations:
                count += 1
        return count

#  def remove_redundant_proteins(self):
  #      unique_sequences = []
  #      ids_to_remove = []
  #      for proteins in self.protein_structure_results.values():
  #          for protein in proteins.all_structures:
  #              if protein.fasta not in unique_sequences:
  #                  unique_sequences.append(protein.fasta)
  #              elif protein.annotation is not None:
  #                  ids_to_remove.append(protein.id)
  #      for id in ids_to_remove:
  #          self._protein_structures.pop(id)