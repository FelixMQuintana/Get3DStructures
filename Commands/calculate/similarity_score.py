import numpy
import numpy as np

from Commands.calculate.calculate import SimilarityScore
from dataStructure.protein.structure import StructureFile
from lib.const import AminoAcids, grantham_distance_matrix_row_dict, grantham_distance_matrix, CoordinateType
from lib.func import find_rigid_alignment, get_coords


class CalculateSubstructureDistances(SimilarityScore):

    @property
    def distances_matrix(self) -> np.ndarray:
        if self._distances is None:
            self._distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                              len(self.collection.protein_structure_results.values())))
        return self._distances
    @staticmethod
    def metric(structure_file1: StructureFile, structure_file2: StructureFile, coord_type,
               *args) -> float:
        sequence = AminoAcids.get_rep(structure_file1.fasta)
        sequence2 = AminoAcids.get_rep(structure_file2.fasta)
        grantham_similarity = np.zeros(shape=(len(sequence)))
        for index2 in range(len(sequence)):
            row = grantham_distance_matrix_row_dict[type(sequence[index2])]
            column = grantham_distance_matrix_row_dict[type(sequence2[index2])]
            physiochemical_value = grantham_distance_matrix[row][column]
            grantham_similarity[index2] = physiochemical_value
        return float(np.sum(grantham_similarity))


class CalculateGeometricalDistance(SimilarityScore):

    @staticmethod
    def metric(structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType, *args) -> None:
        structure_coords = get_coords(structure_file1, coord_type)
        structure_inner_coords = get_coords(structure_file2, coord_type)
        R, t = find_rigid_alignment(structure_coords, structure_inner_coords)
        aligned_structure_coords = (R.dot(structure_coords.T)).T + t
        return np.sqrt((((structure_coords - aligned_structure_coords) ** 2).sum(axis=1)).mean())

    @property
    def distances_matrix(self) -> np.ndarray:
        if self._distances is None:
            self._distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                              len(self.collection.protein_structure_results.values())))
        return self._distances


        #     distances[index][index2] = self.metric(aligned_structure_coords, structure_inner_coords)
       # distance_matrix[index][index2] = self.metric(sequence1=AminoAcids.get_rep(structure.fasta),
       #                                              sequence2=AminoAcids.get_rep(
       #                                                  structure_inner.fasta))  # , coords1=aligned_structure_coords,

    #           coords2=structure_inner_coords)
    # distances[index][index2] = self.metric(structure

   # def run(self) -> None:
   #     threads = []#

#        distances = self.calculate_distances(distances)
#        [thread.wait() for thread in threads]
#        np.save(
#            file="//media/felix/ShortTerm/Research/KpsData/KpsT/pathotype_structures/grantham_only_alphafold_pathotype",
#            arr=distances)

#    def calculate_distances(self, distances):
#        for index, structures in tqdm.tqdm(enumerate(self.collection.protein_structure_results.values())):
#            for structure in structures.all_structures:
#                for index2, structures_inner in enumerate(self.collection.protein_structure_results.values()):
#                    for structure_inner in structures_inner.all_structures:
#                        structure_coords = self.get_coords(structure)
#                        structure_inner_coords = self.get_coords(structure_inner)
#                        R, t = lib.func.find_rigid_alignment(structure_coords, structure_inner_coords)
#                        aligned_structure_coords = (R.dot(structure_coords.T)).T + t
                        #  rot, score= Rotation.align_vectors(structure_coords,structure_inner_coords)
                        #  aligned_coords_inner = rot.apply(structure_inner_coords)
#                        sequence = AminoAcids.get_rep(structure.fasta)
#                        sequence2 = AminoAcids.get_rep(structure_inner.fasta)
                        #     distances[index][index2] = self.metric(aligned_structure_coords, structure_inner_coords)
#                        distances[index][index2] = self.metric(sequence1=AminoAcids.get_rep(structure.fasta),
#                                                               sequence2=AminoAcids.get_rep(
#                                                                   structure_inner.fasta))  # , coords1=aligned_structure_coords,
                    #           coords2=structure_inner_coords)
                    # distances[index][index2] = self.metric(structure.fasta.strip("X"), structure_inner.fasta.strip("X"), structure_coords, structure_inner_coords, encoding)
#        return distances


class CalculateLLMBasedDistanceMetric(SimilarityScore):

    @property
    def distances_matrix(self) -> numpy.ndarray:
        if self._distances is None:
            self._distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                              len(self.collection.protein_structure_results.values())))
        return self._distances

    #TODO Come back and fix
    @staticmethod
    def metric(structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType,
               *args) -> float:
        pass
        sequence = AminoAcids.get_rep(structure_file1.fasta)
        sequence2 = AminoAcids.get_rep(structure_file2.fasta)
        embeddings = get_encoding([sequence2], self.tokenizer, self.esm_model)
        l1_norm = torch.sum(torch.abs(encoding[0]["embedding"] - embeddings[0]["embedding"]))
        #geomerty = np.sum(np.sum(
        #    (1 / (np.cosh(0.3 * np.linalg.norm(np.array([coords1[amino_acid1], coords2[amino_acid2]]))))) for
        #    amino_acid2 in range(len(coords2)))
        #                  for amino_acid1 in range(len(coords1)))
        return  l1_norm
    #   return np.sum(1 / (0.3 * np.cosh(calculate_rmsd(coords1, coords2))) * l1_norm)
