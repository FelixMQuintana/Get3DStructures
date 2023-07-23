"""

"""
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Type, Dict

from dataStructure.protein.accession import UniProtIDFastaFile, UniProtAcessionFile
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure import HomologyStructure, GOLabelsFile, CrystalStructure, StructureFile, \
    SupportedFileType
from lib.const import SupportedFileTypeRegex, ALLOWED_EXT


class DataFetcher(ABC):

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
        return [self.file_type(file) for file in file_path.rglob(self.file_type_regex.value)]


class ExperimentalStructureFetcher(DataFetcher):

    @property
    def file_type(self) -> Type[CrystalStructure]:
        return CrystalStructure

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.EXPERIMENTAL_STRUCTURE


class HomologyStructureFetcher(DataFetcher):

    @property
    def file_type_regex(self) -> SupportedFileTypeRegex:
        return SupportedFileTypeRegex.HOMOLOGY_STRUCTURE

    @property
    def file_type(self) -> Type[HomologyStructure]:
        return HomologyStructure


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
       #elf._collection = []
        self._results = {}
        self._protein_structures: Dict[ProteinStructures] = {}
        self.add_fetchers(*args)

    def add_fetchers(self, *args: DataFetcher):
        collection = []
        for data_fetcher in args:
            collection.extend(data_fetcher.get_data(self.working_directory))
        for file in collection:
            if isinstance(file, StructureFile):
                if isinstance(file, UniProtAcessionFile) and file.path.name != (file.id+ALLOWED_EXT.JSON.value):
                    continue
                if self._results.get(file.id) is None:
                    self._results[file.id] = [file]
                else:
                    self._results[file.id].append(file)
        for key, value in self._results.items():
            self._protein_structures[key] = ProteinStructures(value)

    @property
    def protein_structure_results(self):#-> Dict[ProteinStructures]:
        return self._protein_structures
