import abc
import logging
from abc import ABC
from pathlib import Path
from typing import Optional

import numpy
import numpy as np
import pandas
from tqdm import tqdm

from lib.const import Metrics, CoordinateType
from Commands.command import Command, FactoryBuilder
from dataStructure.protein.structure import StructureFile
from dataStructure.collections import ExperimentalStructureFetcher, HomologyStructureFetcher, Collection


class Calculate(Command, ABC):
    """

    """


fetcher_type = {
    Metrics.HOMOLOGY_STRUCTURES: HomologyStructureFetcher,
    Metrics.CRYSTAL_STRUCTURES: ExperimentalStructureFetcher

}


class SimilarityScore(Calculate, ABC):

    def __init__(self, structure_type: Metrics.HOMOLOGY_STRUCTURES, coord_type: CoordinateType, out_file_name: str):
        super().__init__()
        self.collection = Collection(self.working_directory, fetcher_type[structure_type]())
        self._distances: Optional[np.ndarray] = None
        self.coord_type = coord_type
        self._out_file_name: str = out_file_name

    @staticmethod
    @abc.abstractmethod
    def metric(structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType,
               *args) -> float:
        """

        Returns
        -------

        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def distances_matrix(self) -> numpy.ndarray:
        raise NotImplementedError

    def run(self) -> None:
        ordering = open(self.working_directory.joinpath("ordering_of_distances.txt"))

        for index, structures in tqdm(enumerate(self.collection.protein_structure_results.values())):
            for structure in structures.all_structures:
                for inner_index, structures_inner in enumerate(self.collection.protein_structure_results.values()):
                    for structure_inner in structures_inner.all_structures:
                        self.distances_matrix[index, inner_index] = self.metric(structure, structure_inner,
                                                                                self.coord_type)
                        ordering.write(structure.id + "," + structure_inner.id + "\n")
        np.save(self._out_file_name)
        ordering.close()


class CalculateClusters(Calculate, ABC):

    def __init__(self, data_file: Path, ordering_file: Optional[Path] = None):
        super().__init__()
        self.data = np.load(data_file.as_uri())
        if ordering_file is None:
            logging.info("Ordering file was not provided, trying to check working dir")
            if self.working_directory.joinpath("ordering_of_distances.txt").is_file():
                self.ordering_file = self.working_directory.joinpath("ordering_of_distances.txt")
        else:
            self.ordering_file = ordering_file

    @staticmethod
    @abc.abstractmethod
    def clustering_algorithm(data: np.ndarray, ordering_file):
        raise NotImplementedError

    def read_ordering_file(self) -> np.ndarray:
        return pandas.read_csv(self.ordering_file, header=None)

    def run(self) -> None:
        self.clustering_algorithm(self.data,self.ordering_file)

calculate_commands = {

}

class CalculateClassFactory(FactoryBuilder):
    @staticmethod
    def build(*args, **kwargs) -> Command:

