import abc
import subprocess
import threading
from abc import ABC
import itertools
import scipy
from Bio.PDB import PDBIO
from scipy.spatial import distance_matrix
from Bio import PDB, AlignIO
import lib.func
from Commands.command import Command, FactoryBuilder
import tqdm
from dataStructure.collections import ExperimentalStructureFetcher, \
    UniProtAcessionFetcher
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure import StructureFile, HomologyStructure
from lib.const import StructureCharacteristicsMode, MotifSearchMode, MotifRefinements, AllowedExt, AminoAcids
from lib.func import get_encoding, calculate_rmsd, calculate_grantham_distance
# from tmtools import tm_align
# from tmtools.io import get_structure, get_residue_data
from scipy.spatial.transform import Rotation
import matplotlib
# from matplotlib import pyplot as plt
# matplotlib.use('TkAgg')
# import transformers
from lib.func import *
from lib.const import e_coli_k_type


class Characteristics(Command, ABC):

    def __init__(self, binding_site_database: Path):
        super().__init__()

        if not binding_site_database.exists():
            logging.info("Database doesn't appear to exist. Building it now!")
            logging.info("Building directory: %s" % binding_site_database)
            os.mkdir(binding_site_database)
        self.binding_site_database: Path = binding_site_database
        self.collection = Collection(self.working_directory, ExperimentalStructureFetcher())

    @abc.abstractmethod
    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures, *args) -> None:
        raise NotImplementedError

    def run(self) -> None:
        [(logging.info("Building database for %s." % self.binding_site_database.joinpath(structure_file.id)),
          os.mkdir(self.binding_site_database.joinpath(structure_file.id)))
         for protein_structures in self.collection.protein_structure_results.values()
         for structure_file in protein_structures.all_structures
         if not self.binding_site_database.joinpath(structure_file.id).exists()]
        [self.command(structure, protein_structures) for protein_structures in
         self.collection.protein_structure_results.values() for
         structure in protein_structures.all_structures]
        # threads: List[threading.Thread] = [threading.Thread(target=self.command(structure, protein_structures))
        #                                   for protein_structures in self.collection.protein_structure_results.values() for
        #                                   structure in protein_structures.all_structures]
        # [thread.start() for thread in threads]
        # [thread.join() for thread in threads]


class FindPockets(Characteristics):

    def __init__(self, binding_site_database: Path):
        super().__init__(binding_site_database)
        self.autosite_location = self.args[StructureCharacteristicsMode.AUTOSITE.value]

    def command(self, structure_file: StructureFile, protein_structures: ProteinStructures) -> None:
        with threading.Lock():
            logging.info("Running the command:  %s/prepare_receptor -r %s -o %s"
                         % (self.autosite_location, structure_file.path,
                            structure_file.path.name.split(".")[0] + ".pdbqt"))
            p1 = subprocess.Popen([
                self.autosite_location + "prepare_receptor -r" + structure_file.path.as_posix() + " -o " +
                structure_file.path.name.split(".")[0] + ".pdbqt"],
                shell=True)
            p1.wait()
            logging.info("Running the command: %sautosite -r %s -o %s" % (self.autosite_location,
                                                                          structure_file.path.name.split(".")[
                                                                              0] + ".pdbqt",
                                                                          self.binding_site_database.joinpath(
                                                                              structure_file.id).joinpath(
                                                                              structure_file.path.name.split(".")[0])))
            if not self.binding_site_database.joinpath(structure_file.id).joinpath(
                    structure_file.path.name.split(".")[0]).exists():
                os.mkdir(self.binding_site_database.joinpath(structure_file.id).joinpath(
                    structure_file.path.name.split(".")[0]))
            os.system("%s/autosite -r %s -o %s" % (self.autosite_location,
                                                   structure_file.path.name.split(".")[0] + ".pdbqt",
                                                   self.binding_site_database.joinpath(structure_file.id).joinpath(
                                                       structure_file.path.name.split(".")[0])))


class FixPDBFiles(Characteristics):

    def __init__(self, binding_site_database: Path):

        super().__init__(binding_site_database)
        self.collection = Collection(self.working_directory, HomologyStructureFetcher(), UniProtAcessionFetcher())

    def command(self, structure_file: HomologyStructure, protein_structures: ProteinStructures) -> None:
        pdb_parser = PDB.PDBParser()
        working_dir: Path = self.binding_site_database.joinpath(structure_file.id)

        pdb = pdb_parser.get_structure(structure_file.path.name, structure_file.path)
        io = PDBIO()
        residue_ids = []
        for index, score in enumerate(structure_file.piddt):
            if score > 70:
                break
            residue_ids.append(index)
        for index, score in enumerate(reversed(structure_file.piddt)):
            if score > 70:
                break
            residue_ids.append(len(structure_file.piddt) - index)

        residue_ids_to_remove = [id.id[1] for id in pdb.get_residues() if
                                 int(id.id[1]) in residue_ids]
        try:
            for chain in pdb[0]:
                [chain.detach_child((' ', id, ' ')) for id in residue_ids_to_remove]
        except KeyError as ex:
            for chain in pdb:
                [chain.detach_child((' ', id, ' ')) for id in residue_ids_to_remove]
        io.set_structure(pdb)
        if len([atom for atom in pdb.get_atoms()]) == 0:
            return None
        io.save(str(self.binding_site_database.joinpath(structure_file.id).joinpath(
            structure_file.path.name.removesuffix(".pdb") + "_trimmed.pdb")),
            preserve_atom_numbering=True)
        os.system(f"cp {structure_file.path.parent.joinpath(structure_file.id).with_suffix(AllowedExt.FASTA.value)} "
                  f"{working_dir}")
        os.system(f"cp {structure_file.path.parent.joinpath(structure_file.id).with_suffix(AllowedExt.JSON.value)} "
                  f"{working_dir}")


# def command(self, protein_structure: ProteinStructures) -> None:
#      with open("/mnt/ResearchLongTerm/FunSoCTrainingData/lineages.txt", "a") as e:
#        e.write(protein_structure.id)
#        [e.write(" " + item) for item in protein_structure.uniprotID.structural_data["organism"]["lineage"]]
#        e.write("\n")


class FindCustomBindingSite(Command):

    def __init__(self, align_file, coverage_region):
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())
        self.alignment = AlignIO.read(align_file, "fasta")

    @staticmethod
    def calculate_coverage(rec, start, end):
        alignment_motif = rec.seq[start:end + 1]  # 194:206]
        length_motif = len(alignment_motif)
        invalid_count = 0
        indexs = []
        for index, char in enumerate(alignment_motif):
            if char == "-":
                invalid_count += 1
                indexs.append(index)
        if invalid_count == 0:
            return alignment_motif
        if invalid_count / length_motif < .40:
            motif = ""
            previous_index = None
            print(indexs)
            for index in indexs:
                if motif == "":
                    motif += str(rec.seq[start:start + index])
                else:
                    motif += str(rec.seq[previous_index:start + index])
                previous_index = index
            motif += str(rec.seq[start + previous_index + 1:end + len(indexs) + 1])
        else:
            return None
        return motif

    def run(self) -> None:
        pdb_parser = PDB.PDBParser()
        io = PDBIO()
        accessions = strain_to_accession(Path("/home/felix/testing.json"), e_coli_k_type)
        for index, protein_structures in enumerate(self.collection.protein_structure_results.values()):
            record_of_interest = None
            for rec in self.alignment:
                if rec.id == accessions[protein_structures.id][0].replace(" ", ""):
                    record_of_interest = rec
            motif = str(self.calculate_coverage(record_of_interest, 199, 209))
            stability_motif = re.search(motif, str(protein_structures.all_structures[0].fasta))
            try:
                stability_motif = list(range(stability_motif.start() + 1, stability_motif.end() + 1))
            except AttributeError as ex:
                print(f'skipping because of {ex}')
                continue
            pdb = pdb_parser.get_structure(protein_structures.all_structures[0].path.name,
                                           protein_structures.all_structures[0].path)
            pdb_ids = [residue.id[1] for residue in pdb.get_residues()]
            for chain in pdb[0]:
                [chain.detach_child((' ', id, ' ')) for id in pdb_ids if id not in stability_motif]
            io.set_structure(pdb)
            io.save(str(self.working_directory.joinpath(
                protein_structures.all_structures[0].path.name + "_beta_sheet.pdb")))


class CalculateSubstructureDistances(Characteristics):

    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures, *args) -> None:
        pass

    def __init__(self, binding_site):
        super().__init__(binding_site)
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())
        self.metric = lib.func.granthem_distance_fixed

    def run(self) -> None:
        threads = []
        distances = np.zeros(shape=(len(self.collection.protein_structure_results.values()),
                                    len(self.collection.protein_structure_results.values())))
        distances = self.calculate_distances(distances)
        [thread.wait() for thread in threads]
        np.save(
            file="//media/felix/ShortTerm/Research/KpsData/KpsT/pathotype_structures/grantham_only_alphafold_pathotype",
            arr=distances)

    def calculate_distances(self, distances):
        for index, structures in tqdm.tqdm(enumerate(self.collection.protein_structure_results.values())):
            for structure in structures.all_structures:
                # encoding = get_encoding([structure.fasta], self.tokenizer, self.esm_model)
                for index2, structures_inner in enumerate(self.collection.protein_structure_results.values()):
                    for structure_inner in structures_inner.all_structures:
                        structure_coords = self.get_coords(structure)
                        structure_inner_coords = self.get_coords(structure_inner)
                        #  print(structure.path)
                        #  print(structure_inner.path)
                        R, t = lib.func.find_rigid_alignment(structure_coords, structure_inner_coords)
                        aligned_structure_coords = (R.dot(structure_coords.T)).T + t
                        #  rot, score= Rotation.align_vectors(structure_coords,structure_inner_coords)
                        #  aligned_coords_inner = rot.apply(structure_inner_coords)
                        sequence = AminoAcids.get_rep(structure.fasta)
                        sequence2 = AminoAcids.get_rep(structure_inner.fasta)
                        #     distances[index][index2] = self.metric(aligned_structure_coords, structure_inner_coords)
                        distances[index][index2] = self.metric(sequence1=AminoAcids.get_rep(structure.fasta),
                                                               sequence2=AminoAcids.get_rep(
                                                                   structure_inner.fasta))  # , coords1=aligned_structure_coords,
                    #           coords2=structure_inner_coords)
                    # distances[index][index2] = self.metric(structure.fasta.strip("X"), structure_inner.fasta.strip("X"), structure_coords, structure_inner_coords, encoding)
        return distances

    def get_coords(self, structure):
        pdb_parser = PDB.PDBParser()
        protein_b = pdb_parser.get_structure(structure.id, structure.path)
        coords_protein_b = []
        for x in protein_b.get_atoms():
            #    if x.name in ["N", "CA", "C", "O"]:
            if x.name in ["CA"]:
                coords_protein_b.append(list(x.coord))
        return np.array(coords_protein_b)

    @staticmethod
    def do_it(x, protein_struct, distances, metric, proteins2, encoding):
        for index, protein_structures_inner in enumerate(x):
            pdb_parser = PDB.PDBParser()
            for protein in protein_structures_inner.all_structures:
                #     if int(str(protein_struct.path).split("_")[-1].strip('.pdb')) != int(
                #             str(protein.path).split("_")[-1].strip('.pdb')):
                #    print(f"Skipping {protein.path}")
                #         continue
                #   pdb1 = pdb_parser.get_structure(protein_struct.path.name, protein_struct.path)
                #    pdb_parser2 = PDB.PDBParser()
                #    pdb2 = pdb_parser2.get_structure(protein_structures_inner.path.name, protein_structures_inner.path)
                #     s1 = get_structure(protein_struct.path)

                protein_a = pdb_parser.get_structure(protein_struct.id, protein_struct.path)
                coords_protein_a = []
                for x in protein_a.get_atoms():
                    if x.name in ["N", "CA", "C", "O"]:
                        coords_protein_a.append(list(x.coord))
                coords_protein_a = np.array(coords_protein_a)
                protein_name = protein_struct.path.name.split("_")[0:-1]
                if protein_name[0] == "WP":
                    protein_name = protein_name[0] + "_" + protein_name[1] + ".1"
                else:
                    protein_name = protein_name[0] + ".1"
                #  chain = next(s1.get_chains())
                #   coords, seq = get_residue_data(chain)
                #   s2 = get_structure(protein.path)
                protein_b = pdb_parser.get_structure(protein.id, protein.path)
                coords_protein_b = []
                for x in protein_b.get_atoms():
                    if x.name in ["N", "CA", "C", "O"]:
                        coords_protein_b.append(list(x.coord))
                coords_protein_b = np.array(coords_protein_b)
                protein_name2 = protein.path.name.split("_")[0:-1]
                if protein_name2[0] == "WP":
                    protein_name2 = protein_name2[0] + "_" + protein_name2[1] + ".1"
                else:
                    protein_name2 = protein_name2[0] + ".1"

                rot, rmsd = scipy.spatial.transform.Rotation.align_vectors(coords_protein_a, coords_protein_b)
                coords_protein_b_aligned = rot.apply(coords_protein_b)
                #   chain2 = next(s2.get_chains())
                #     coords2, seq2 = get_residue_data(chain2)
                #         alignment = tm_align(coords, coords2, seq, seq2)
                #        aligned_seq_coords1 = coords.dot(alignment.u.T) + alignment.t
                sequence = AminoAcids.get_rep(protein_struct.fasta)
                sequence2 = AminoAcids.get_rep(protein.fasta)
                distance_metric = metric(coords_protein_a, coords_protein_b_aligned)
                # distance_metric = metric(sequence,sequence2,coords_protein_a,coords_protein_b_aligned) #metric(protein_struct.fasta, protein.fasta,  coords_protein_a, coords_protein_b_aligned)#, proteins2[protein_name].crystal_fastas[0],proteins2[protein_name2].crystal_fastas[0],encoding)
                distances[index] = distance_metric

                break

    def llm_embedding_distance(self, sequence1, sequence2, coords1, coords2, encoding):
        embeddings = get_encoding([sequence2], self.tokenizer, self.esm_model)
        l1_norm = torch.sum(torch.abs(encoding[0]["embedding"] - embeddings[0]["embedding"]))
        geomerty = np.sum(np.sum(
            (1 / (np.cosh(0.3 * np.linalg.norm(np.array([coords1[amino_acid1], coords2[amino_acid2]]))))) for
            amino_acid2 in range(len(coords2)))
                          for amino_acid1 in range(len(coords1)))
        return geomerty * l1_norm
    #   return np.sum(1 / (0.3 * np.cosh(calculate_rmsd(coords1, coords2))) * l1_norm)


supported_commands = {
    StructureCharacteristicsMode.AUTOSITE: FindPockets,
    StructureCharacteristicsMode.FIND_BINDING_POCKETS: FindCustomBindingSite,
    StructureCharacteristicsMode.CALCULATE_RMSD: CalculateSubstructureDistances,
    StructureCharacteristicsMode.TRIM_PDB: FixPDBFiles,
}


class CharacteristicsFactory(FactoryBuilder):

    @staticmethod
    def build(mode: StructureCharacteristicsMode, *args, **kwargs):
        return supported_commands[mode](*args, **kwargs)

