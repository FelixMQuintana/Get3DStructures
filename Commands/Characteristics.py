import abc
import logging
import os
import subprocess
import threading
from abc import ABC
from typing import List

from Bio.PDB import PDBIO, Selection
from matplotlib import pyplot as plt
from scipy.spatial import distance_matrix
import numpy as np
import pandas
from Bio import PDB
from Commands.command import Command, FactoryBuilder
from pathlib import Path

from dataStructure.collections import Collection, HomologyStructureFetcher, ExperimentalStructureFetcher, \
    UniProtAcessionFetcher
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure import StructureFile
from lib.const import StructureCharacteristicsMode, MotifSearchMode, MotifRefinements
import numpy as np

class Characteristics(Command, ABC):

    def __init__(self, binding_site_database: Path):
        super().__init__()

        if not binding_site_database.exists():
            logging.info("Database doesn't appear to exist. Building it now!")
            logging.info("Building directory: %s" % binding_site_database)
            os.mkdir(binding_site_database)
        self.binding_site_database: Path = binding_site_database
        self.collection = Collection(self.working_directory, ExperimentalStructureFetcher(),
                                     HomologyStructureFetcher(), )

    @abc.abstractmethod
    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures) -> None:
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


class BuildMotifStructures(Characteristics):

    def __init__(self, binding_site_database: Path):

        super().__init__(binding_site_database)
        self.collection.add_fetchers(UniProtAcessionFetcher())

    def command(self, structure_file: StructureFile, protein_structures: ProteinStructures) -> None:
        pdb_parser = PDB.PDBParser()
        print(structure_file.path)
        pdb = pdb_parser.get_structure(structure_file.path.name, structure_file.path)

        # residues = [residue.get_atoms() for residue in pdb.get_residues() if residue.id[1] in structure_file.binding_site_residues]
        io = PDBIO()
     #   try:
     #       pred_binding_sites = np.loadtxt(
     #           self.binding_site_database.joinpath(structure_file.id).joinpath("autosite_pred.txt"))
     #   except:
     #       logging.warning("Cant do this")
     #       return None
      #  print(pred_binding_sites)
        residue_ids_to_remove = [id.id[1] for id in pdb.get_residues() if
                                 int(id.id[1]) not in protein_structures.uniprotID.binding_site_residues ]
                                 #and id.id[1] not in pred_binding_sites]
        try:
            for chain in pdb[0]:
                [chain.detach_child((' ', id, ' ')) for id in residue_ids_to_remove]
        except KeyError as ex:
            for chain in pdb:
                [chain.detach_child((' ', id, ' ')) for id in residue_ids_to_remove]
        io.set_structure(pdb)
        if len([atom for atom in pdb.get_atoms()]) == 0:
            return None
        io.save(str(self.binding_site_database.joinpath(structure_file.id).joinpath(structure_file.path.name.removesuffix(".pdb") + "_motif.pdb")),
                preserve_atom_numbering=True)


class FindBindingSite(Characteristics):

    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures) -> None:
        pdb_parser = PDB.PDBParser()

        # pdb_coords = pandas.read_csv(structure_file.path, delim_whitespace=True,skiprows=2).iloc[:, 6:9].to_numpy()
        try:
            cluster_coords = pandas.read_csv(self.binding_site_database.joinpath(structure_file.id).joinpath(
                structure_file.path.name.split(".")[0]).joinpath(
                structure_file.path.name.split(".")[0] + "_cl_002.pdb"), delim_whitespace=True).iloc[:, 5:8].to_numpy()
        except FileNotFoundError:
            return
        pdb = pdb_parser.get_structure(structure_file.path.name, structure_file.path)
        # pdb_cluster_rep = pdb_parser.get_structure(id ,self.binding_site_database.joinpath(id).joinpath(
        #                                              structure_file.path.name.split(".")[0]).joinpath(structure_file.path.name.split(".")[0] + "_cl_001.pdb",
        #                                          ))
        pdb_coords = np.array([residue.center_of_mass() for residue in pdb.get_residues()])
        d_matrix = distance_matrix(pdb_coords, cluster_coords)
        np.savetxt(self.binding_site_database.joinpath(structure_file.id).joinpath("autosite_pred.txt"),
                   np.unique((d_matrix[:] < 4.5).nonzero()[0]))


class CheckBindingSiteQuality(Characteristics):
    """

    """

    def __init__(self, binding_site_database: Path):
        super().__init__(binding_site_database)
        self._correct_state = []

    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures) -> None:
        try:
            pred_binding_sites = np.loadtxt(
                self.binding_site_database.joinpath(structure_file.id).joinpath("autosite_pred.txt"))
        except:
            return None
        if protein_structure.uniprotID.binding_site_residues is None:
            return None
        percent_correctly_found = sum(
            el in pred_binding_sites for el in protein_structure.uniprotID.binding_site_residues) / sum(
            len(site) for site in protein_structure.uniprotID.binding_site_residues)
        self._correct_state.append(percent_correctly_found)

    #    all_atom_distance_matches = np.unique((d_matrix[:]<1.5).nonzero()[0])
    #    atoms = np.array([x for x in pdb.get_atoms()])
    #    matched_atoms_indexs = atoms[all_atom_distance_matches]
    #    matched_residue_ids = np.unique(np.array( [atom.parent.id[1] for atom in matched_atoms_indexs]))


class CalculateLeastRootMeanSquareDistance(Characteristics):

    def __init__(self, binding_site_database: Path):
        super().__init__(binding_site_database)
        self.per_motif_rmsd_vectors = []

    def run(self) -> None:
        [(logging.info("Building database for %s." % self.binding_site_database.joinpath(structure_file.id)),
          os.mkdir(self.binding_site_database.joinpath(structure_file.id)))
         for protein_structures in self.collection.protein_structure_results.values()
         for structure_file in protein_structures.all_structures
         if not self.binding_site_database.joinpath(structure_file.id).exists()]
        [self.command(structure, protein_structures) for protein_structures in
         self.collection.protein_structure_results.values() for
         structure in protein_structures.all_structures]
        rmsd_vector = np.array(self.per_motif_rmsd_vectors)

    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures) -> None:
        """

        Parameters
        ----------
        structure_file
        protein_structure

        Returns
        -------

        """
        pdb_parser = PDB.PDBParser()
        pdb_reference = pdb_parser.get_structure(structure_file.path.name, structure_file.path)
        pdb_coords_reference = np.array([atom.get_coord() for atom in pdb_reference.get_atoms()])
        for protein_structures in self.collection.protein_structure_results.values():
            for structure in protein_structures.all_structures:
                if structure.path == structure_file.path:
                    continue

                pdb = pdb_parser.get_structure(structure.path.name, structure.path)
                pdb_coords = np.array([atom.get_coord() for atom in pdb.get_atoms()])
                self.per_motif_rmsd_vectors.append(np.sqrt((pdb_coords_reference - pdb_coords)**2/
                                                           pdb_coords_reference.shape[0]))

supported_commands = {
    StructureCharacteristicsMode.AUTOSITE: FindPockets,
    StructureCharacteristicsMode.BUILD_MOTIFS: BuildMotifStructures,
    StructureCharacteristicsMode.FIND_BINDING_POCKETS: FindBindingSite,
    StructureCharacteristicsMode.CHECK_MOTIF_QUALITY: CheckBindingSiteQuality,
    StructureCharacteristicsMode.CALCULATE_RMSD: CalculateLeastRootMeanSquareDistance
}


class CharacteristicsFactory(FactoryBuilder):

    @staticmethod
    def build(mode: StructureCharacteristicsMode, binding_site_database: Path):
        return supported_commands[mode](binding_site_database)

    #   threads: List = [
    #        [self.thread_pool.apply_async(self.mode, [structures, structure_result.id])
    #         for structures in structure_result.all_structures]
    #        for structure_result in self._structure_results]
    #     [[thread.wait() for thread in thread_list] for thread_list in threads]
    #    plt.hist(self._correct_state )
    #    plt.title("Binding Site Residues Correctly Predicted Dist")
    #    plt.ylabel("Count")
    #    plt.xlabel("Percent Correct")
    #    plt.savefig(self.working_directory.joinpath("Correct_Binding_Sites.tif"))
    #    plt.close()
