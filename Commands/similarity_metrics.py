import abc
from abc import ABC

import numpy as np
from Bio import PDB
from tqdm import tqdm
from sklearn.metrics import rand_score
from Commands.command import FactoryBuilder, Command
from dataStructure.collections import Collection, HomologyStructureFetcher
from lib.const import AminoAcids, SimilarityDistance, ClusterAnalysisType
from lib.func import granthem_distance, find_rigid_alignment, calculate_rmsd


class SimilarityMethod(Command, ABC):

    def __init__(self, distances_file_name: str):
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())
        self.distances_file_path = self.working_directory.joinpath(distances_file_name)

    def run(self) -> None:
        distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                    len(self.collection.protein_structure_results.values())))
        distances = self.distance_calculation(distances)
        np.save(file=str(self.distances_file_path), arr=distances)

    @abc.abstractmethod
    def distance_calculation(self, distances):
        raise NotImplementedError


class GranthamDistance(SimilarityMethod):
    """

    """

    def distance_calculation(self, distances):
        for index, structures in tqdm(enumerate(self.collection.protein_structure_results.values())):
            for structure in structures.all_structures:
                for index2, structures_inner in enumerate(self.collection.protein_structure_results.values()):
                    for structure_inner in structures_inner.all_structures:
                        distances[index][index2] = granthem_distance(sequence1=AminoAcids.get_rep(structure.fasta),
                                                                     sequence2=AminoAcids.get_rep(
                                                                         structure_inner.fasta))
        return distances


def get_coords(structure):
    pdb_parser = PDB.PDBParser()
    protein_b = pdb_parser.get_structure(structure.id, structure.path)
    coords_protein_b = []
    for x in protein_b.get_atoms():
        if x.name in ["N", "CA", "C", "O"]:
            coords_protein_b.append(list(x.coord))
    return np.array(coords_protein_b)


class BackboneGeometry(SimilarityMethod):

    def distance_calculation(self, distances):
        for index, structures in tqdm(enumerate(self.collection.protein_structure_results.values())):
            for structure in structures.all_structures:
                for index2, structures_inner in enumerate(self.collection.protein_structure_results.values()):
                    for structure_inner in structures_inner.all_structures:
                        structure_coords = get_coords(structure)
                        structure_inner_coords = get_coords(structure_inner)
                        R, t = find_rigid_alignment(structure_coords, structure_inner_coords)
                        aligned_structure_coords = (R.dot(structure_coords.T)).T + t
                        distances[index][index2] = calculate_rmsd(aligned_structure_coords, structure_inner_coords)
        return distances


class ClusterAnalysis:

    def __init__(self, ground_truth, labels, metric: ClusterAnalysisType):
        self.ground_truth = ground_truth
        self.labels = labels
        self.metric = metric

    def run(self) -> None:
        if self.metric == ClusterAnalysisType.RAND_INDEX:
            print(rand_score(self.ground_truth, self.labels))


supported_commands = {
    SimilarityDistance.GRANTHAM: GranthamDistance,
    SimilarityDistance.BACKBONE: BackboneGeometry
}


class SimilarityMetricBuilder(FactoryBuilder):

    @staticmethod
    def build(mode: SimilarityDistance, *args, **kwargs):
        return supported_commands[mode](*args, **kwargs)
