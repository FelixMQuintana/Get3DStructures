import os
import time
from copy import deepcopy
from subprocess import Popen, PIPE

import numpy
import numpy as np
import pandas
import pandas as pd
from Bio.PDB import DSSP
from transformers import EsmModel, EsmTokenizer
import torch
import blosum as bl
from Commands.calculate.calculate import SimilarityScore
from Commands.command import FactoryBuilder, Command
from dataStructure.protein.structure.representation.atomic.atomic import Atomic
from dataStructure.protein.structure.representation.surface.computeCharges import assignChargesToNewMesh
from dataStructure.protein.structure.representation.surface.surface import MolecularSurfaceRepresentation
from dataStructure.protein.structure.structure import StructureFile
from lib.const import AminoAcids, grantham_distance_matrix_row_dict, grantham_distance_matrix, CoordinateType, \
    SimilarityDistance, Metrics, phat_matrix
#from lib.func import find_rigid_alignment, get_coords, get_encoding
from proteinAnalysis.bonds.bonds import BondProteinAnalysis
from proteinAnalysis.electrostatics.electrostatic_charge import ElectrostaticChargeProteinAnalysis
from proteinAnalysis.electrostatics.electrostatic_potential import ABPS
from proteinAnalysis.hydrophobic.hydrophobic import HydrophobicityProteinAnalysis


class CalculateGranthamDistances(SimilarityScore):

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

    def metric(self, structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType,
               *args) -> None:
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

class CalculateContactsAndSurfaceDistance(SimilarityScore):
    @property
    def distances_matrix(self) -> np.ndarray:
        if self._distances is None:
            self._distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                              len(self.collection.protein_structure_results.values())))
        return self._distances

    def get_molecular_rep(self, structure_file: StructureFile) -> MolecularSurfaceRepresentation:
        molecular_surface = MolecularSurfaceRepresentation(pdb_filename=structure_file.path)
        hydrophobicity = HydrophobicityProteinAnalysis(molecular_surface).compute()
        hbond = ElectrostaticChargeProteinAnalysis(molecular_surface).compute()
        molecular_structure_old: MolecularSurfaceRepresentation = deepcopy(molecular_surface)
        molecular_surface.repair()
        hbond_verts = assignChargesToNewMesh(molecular_surface.vertices, molecular_structure_old.vertices,
                                             hbond)
        hydrophobic_verts = assignChargesToNewMesh(molecular_surface.vertices, molecular_structure_old.vertices,
                                                   hydrophobicity)

        charges = ABPS(molecular_surface).compute()
        normalized_charges = charges / 10
        molecular_surface.mesh.add_attribute("charges")
        molecular_surface.mesh.set_attribute("charges", normalized_charges)

        molecular_surface.mesh.add_attribute("hydrophobicity")
        molecular_surface.mesh.set_attribute("hydrophobicity", hydrophobic_verts)

        molecular_surface.mesh.add_attribute("hbond")
        molecular_surface.mesh.set_attribute("hbond", hbond_verts)

        n1 = molecular_surface.vertex_normals[:, 0]
        n2 = molecular_surface.vertex_normals[:, 1]
        n3 = molecular_surface.vertex_normals[:, 2]
        molecular_surface.mesh.add_attribute("vertex_nx")
        molecular_surface.mesh.set_attribute("vertex_nx", n1)
        molecular_surface.mesh.add_attribute("vertex_ny")
        molecular_surface.mesh.set_attribute("vertex_ny", n2)
        molecular_surface.mesh.add_attribute("vertex_nz")
        molecular_surface.mesh.set_attribute("vertex_nz", n3)
        return molecular_surface

    def get_contacts(self, structure_file: StructureFile):
        atomic_rep = Atomic(structure_file.path)
        bond_analysis = BondProteinAnalysis(atomic_rep)
        bond_analysis.compute()
        tsv_file = structure_file.path.parent.joinpath("bond_contacts_" + structure_file.id+".tsv")

        contacts_df = pandas.read_csv(tsv_file, sep="\t", header=None, skiprows=2)

        #dssp = DSSP(Atomic.pdb, structure_file.path ,dssp='mkdssp')
        #dssp
        return contacts_df,

    def get_secondary_structure(self, structure_file: StructureFile) -> pd.DataFrame:
        args = [
            "mkdssp",
            str(structure_file.path.name),
            str(structure_file.path.name).split(".")[0] + ".dssp"
        ]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(structure_file.path.parent))
        stdout, stderr = p2.communicate()
        try:
            df = pd.read_csv(
                str(structure_file.path.parent) + "/" + str(structure_file.path.name).split(".")[0] + ".dssp", header=None, skiprows=28)
            hits = []
            for entry in df[0]:
                if len(entry.split()) < 10:
                    continue
                split = entry.split()
                hits.append([split[0:5]])


            #open()
    #        dssp_file = open(str(structure_file.path.parent)+"/"+str(structure_file.path.name).split(".")[0] + ".dssp", "r")
    #    # dssp_lines = dssp_file.readlines()
    #        lines = []
    #        data_start = False
    #        for line in dssp_file:
    #            if line == "_dssp_struct_summary.z_ca \n":
    #                data_start = True
     #               continue
     #           if data_start:
     #               if line == "# \n":
     #                   break
     #               else:
    #                    lines.append(line.split())
        except FileNotFoundError:
            raise FileNotFoundError("Couldnt find the structure")
     #   pd.DataFrame
     #   df = pd.read_csv(str(structure_file.path.parent)+"/"+str(structure_file.path.name).split(".")[0]
     #               + ".dssp",header=None,skiprows=28)
     #   for entry df[0]
     #   #df = pd.DataFrame(lines)
        #return df[[3,4]]
        return pd.DataFrame(hits)

    def metric(self, structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType,
               *args) -> None:

        molecular_surface1 = self.get_molecular_rep(structure_file1)
        molecular_surface2 = self.get_molecular_rep(structure_file2)
        contacts1 = self.get_contacts(structure_file1)
        contacts2 = self.get_contacts(structure_file2)
        secondary_structures1 = self.get_secondary_structure(structure_file1)
        secondary_structures2 = self.get_secondary_structure(structure_file2)
        hits_of_interest = None
        for entry in secondary_structures1.values:
            if entry[0][4] not in ['H', 'I']:
                break
            else:
                hits_of_interest = int(entry[0][0])
        contacts_of_interest = []
        for entry in contacts1[0].values:
            contact_1 = entry[2].split(':')
            contact_2 = entry[3].split(':')
            if int(contact_1[2]) <= hits_of_interest:
                contacts_of_interest.append(entry)
            elif int(contact_2[2]) <= hits_of_interest:
                contacts_of_interest.append(entry)
        contacts_of_interest = pd.DataFrame(contacts_of_interest)
        print("poop")
        #structure_file1




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

    def __init__(self):
        super().__init__()
        self._distances = None
        self.esm_model = EsmModel.from_pretrained("facebook/esm2_t6_650M_UR50D")
        self.tokenizer = EsmTokenizer.from_pretrained("facebook/esm2_t6_650M_UR50D")
        self.save_embeddings = []

    @property
    def distances_matrix(self) -> numpy.ndarray:
        if self._distances is None:
            self._distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                              len(self.collection.protein_structure_results.values())))
        return self._distances

    # TODO Come back and fix
    def metric(self, structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType,
               *args) -> float:
        sequence = AminoAcids.get_rep(structure_file1.fasta)
        sequence2 = AminoAcids.get_rep(structure_file2.fasta)
        #TODO: commented soemthing out
        #embeddings = get_encoding([sequence2], self.tokenizer, self.esm_model)
        # self.save_embeddings.extend([embeddings[0]["embedding"], embeddings[1]["embedding"]])
        #l1_norm = torch.sum(torch.abs(embeddings[1]["embedding"] - embeddings[0]["embedding"]))
        # geomerty = np.sum(np.sum(
        #    (1 / (np.cosh(0.3 * np.linalg.norm(np.array([coords1[amino_acid1], coords2[amino_acid2]]))))) for
        #    amino_acid2 in range(len(coords2)))
        #                  for amino_acid1 in range(len(coords1)))
      #  return l1_norm
    #   return np.sum(1 / (0.3 * np.cosh(calculate_rmsd(coords1, coords2))) * l1_norm)


class CalculateBlosumSimilarityScore(SimilarityScore):

    def __init__(self, structure_type = Metrics.HOMOLOGY_STRUCTURES,
                 coord_type: CoordinateType = CoordinateType.BACKBONE, out_file_name: str = None):
        super().__init__(structure_type,coord_type,out_file_name)
        self.matrix = bl.BLOSUM(80)

    @property
    def distances_matrix(self) -> numpy.ndarray:
        if self._distances is None:
            self._distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                              len(self.collection.protein_structure_results.values())))
        return self._distances

    def metric(self, structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType,
               *args):
        sequence = AminoAcids.get_rep(structure_file1.fasta)
        sequence2 = AminoAcids.get_rep(structure_file2.fasta)
        value = 0
        for index, letter in enumerate(sequence):
            try:
                value += self.matrix[letter.single_letter_rep][sequence2[index].single_letter_rep]
            except IndexError:
                value += -1
        return value
class CalculatePHATSimilarityScore(SimilarityScore):

    def __init__(self, structure_type = Metrics.HOMOLOGY_STRUCTURES,
                 coord_type: CoordinateType = CoordinateType.BACKBONE, out_file_name: str = None):
        super().__init__(structure_type,coord_type,out_file_name)
        #TODO Fix
        self.matrix = None

    @property
    def distances_matrix(self) -> numpy.ndarray:
        if self._distances is None:
            self._distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                              len(self.collection.protein_structure_results.values())))
        return self._distances

    def metric(self, structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType,
               *args):
        sequence = AminoAcids.get_rep(structure_file1.fasta)
        sequence2 = AminoAcids.get_rep(structure_file2.fasta)
        value = 0
        for index, letter in enumerate(sequence):
            try:
                value += self.matrix[letter.single_letter_rep][sequence2[index].single_letter_rep]
            except IndexError:
                value += -1
        return value

sim_commands = {
    SimilarityDistance.GRANTHAM.value: CalculateGranthamDistances,
    SimilarityDistance.BACKBONE.value: CalculateGeometricalDistance,
    SimilarityDistance.LLM.value: CalculateLLMBasedDistanceMetric,
    SimilarityDistance.BLOSUM62.value: CalculateBlosumSimilarityScore,
    SimilarityDistance.PHAT.value: CalculatePHATSimilarityScore,
    SimilarityDistance.FULL.value :CalculateContactsAndSurfaceDistance
}


class CalculateClassFactory(FactoryBuilder):
    @staticmethod
    def build(mode: SimilarityDistance, *args, **kwargs) -> Command:
        return sim_commands[mode](*args, **kwargs)
