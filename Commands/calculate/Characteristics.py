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





supported_commands = {
    StructureCharacteristicsMode.AUTOSITE: FindPockets,
    StructureCharacteristicsMode.FIND_BINDING_POCKETS: FindCustomBindingSite,
#    StructureCharacteristicsMode.CALCULATE_RMSD: CalculateSubstructureDistances,
    StructureCharacteristicsMode.TRIM_PDB: FixPDBFiles,
}





class CharacteristicsFactory(FactoryBuilder):

    @staticmethod
    def build(mode: StructureCharacteristicsMode, *args, **kwargs):
        return supported_commands[mode](*args, **kwargs)

