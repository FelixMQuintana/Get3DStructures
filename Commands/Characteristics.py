import logging
import os
import subprocess
import threading
from typing import List

from Bio.PDB import PDBIO, Selection
from matplotlib import pyplot as plt
from scipy.spatial import distance_matrix
import numpy as np
import pandas
from Bio import PDB
from Commands.Structure import StructureFile
from Commands.post_processing import PostProcessing
from pathlib import Path
from lib.const import StructureCharacteristicsMode, MotifSearchMode, MotifRefinements


class Characteristics(PostProcessing):

    def __init__(self, specific_file: Path, mode: StructureCharacteristicsMode, binding_site_database: Path):
        super().__init__(specific_file)
        if mode == StructureCharacteristicsMode.AUTOSITE:
            self.mode = self.finding_ligand_binding_pockets
            self.autosite_location = self.args[StructureCharacteristicsMode.AUTOSITE.value]
        elif mode == StructureCharacteristicsMode.MOTIFS:
            self.mode = self.build_motifs
          #  self.mode = self.get_binding_site_motif
           # self.mode = self.create_binding_site_quality
            self._correct_state = []
        self.binding_site_database: Path = binding_site_database
        if not self.binding_site_database.exists():
            logging.info("Database doesn't appear to exist. Building it now!")
            logging.info("Building directory: %s" % self.binding_site_database)
            os.mkdir(self.binding_site_database)
        [(logging.info("Building database for %s." % self.binding_site_database.joinpath(structure_result.id)),
          os.mkdir(self.binding_site_database.joinpath(structure_result.id)))
         for structure_result in self._structure_results
         if not self.binding_site_database.joinpath(structure_result.id).exists()]

    def run(self) -> None:
        [self.mode(structures, r.id) for r in self._structure_results for structures in r.all_structures ]

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

    def finding_ligand_binding_pockets(self, structure_file: StructureFile, id: str):

        with threading.Lock():
            logging.info("Running the command:  %s/prepare_receptor -r %s -o %s"
                         % (self.autosite_location, structure_file.path, structure_file.path.name.split(".")[0] + ".pdbqt"))
            p1 = subprocess.Popen([self.autosite_location + "prepare_receptor -r" + structure_file.path.as_posix() + " -o " + structure_file.path.name.split(".")[0] + ".pdbqt" ],
                                  shell=True)
            p1.wait()
            logging.info("Running the command: %sautosite -r %s -o %s" % (self.autosite_location,
                                                                          structure_file.path.name.split(".")[
                                                                              0] + ".pdbqt",
                                                                          self.binding_site_database.joinpath(
                                                                              id).joinpath(
                                                                              structure_file.path.name.split(".")[0])))
            if not self.binding_site_database.joinpath(id).joinpath(structure_file.path.name.split(".")[0]).exists():
                os.mkdir(self.binding_site_database.joinpath(id).joinpath(structure_file.path.name.split(".")[0]))
            os.system("%s/autosite -r %s -o %s" % (self.autosite_location,
                                                   structure_file.path.name.split(".")[0] + ".pdbqt",
                                                   self.binding_site_database.joinpath(id).joinpath(
                                                       structure_file.path.name.split(".")[0])))

    def get_binding_site_motif(self, structure_file: StructureFile, id: str, ): #mode: MotifSearchMode, refinement: MotifRefinements) -> None:
        pdb_parser = PDB.PDBParser()

        #pdb_coords = pandas.read_csv(structure_file.path, delim_whitespace=True,skiprows=2).iloc[:, 6:9].to_numpy()
        try:
            cluster_coords = pandas.read_csv(self.binding_site_database.joinpath(id).joinpath(
                                                       structure_file.path.name.split(".")[0]).joinpath(
            structure_file.path.name.split(".")[0] + "_cl_002.pdb"), delim_whitespace=True).iloc[:, 5:8].to_numpy()
        except FileNotFoundError:
            return
        pdb = pdb_parser.get_structure(structure_file.path.name, structure_file.path)
        #pdb_cluster_rep = pdb_parser.get_structure(id ,self.binding_site_database.joinpath(id).joinpath(
         #                                              structure_file.path.name.split(".")[0]).joinpath(structure_file.path.name.split(".")[0] + "_cl_001.pdb",
         #                                          ))
        pdb_coords = np.array([residue.center_of_mass() for residue in pdb.get_residues() ])
        d_matrix = distance_matrix (pdb_coords, cluster_coords)
        np.savetxt(self.binding_site_database.joinpath(id).joinpath("autosite_pred.txt"), np.unique((d_matrix[:] < 4.5).nonzero()[0]))

    def build_motifs(self, structure_file: StructureFile, id: str):
        pdb_parser = PDB.PDBParser()
        pdb = pdb_parser.get_structure(structure_file.path.name, structure_file.path)
        #residues = [residue.get_atoms() for residue in pdb.get_residues() if residue.id[1] in structure_file.binding_site_residues]
        io = PDBIO()
        try:
            pred_binding_sites = np.loadtxt(self.binding_site_database.joinpath(id).joinpath("autosite_pred.txt"))
        except:
            logging.warning("Cant do this")
            return None
        print(pred_binding_sites)
        residue_ids_to_remove = [id.id[1] for id in pdb.get_residues() if int(id.id[1]) not in structure_file.binding_site_residues and id.id[1] not in pred_binding_sites]
        for chain in pdb[0]:
            [chain.detach_child((' ', id, ' ')) for id in residue_ids_to_remove]
        io.set_structure(pdb)
        io.save(str(structure_file.path.parent.joinpath(structure_file.path.name.removesuffix(".pdb")+"_motif.pdb" )), preserve_atom_numbering=True)

    def create_binding_site_quality(self, structure_file: StructureFile, id: str):
        try:
            pred_binding_sites = np.loadtxt(self.binding_site_database.joinpath(id).joinpath("autosite_pred.txt"))
        except:
            return None
        print(pred_binding_sites)
        print(structure_file.binding_site_residues)
        if structure_file.binding_site_residues is None:
            return None
        percent_correctly_found = sum(el in pred_binding_sites for el in structure_file.binding_site_residues) / sum(
            len(site) for site in structure_file._binding_site_residues)
        self._correct_state.append(percent_correctly_found)
    #    all_atom_distance_matches = np.unique((d_matrix[:]<1.5).nonzero()[0])
    #    atoms = np.array([x for x in pdb.get_atoms()])
    #    matched_atoms_indexs = atoms[all_atom_distance_matches]
    #    matched_residue_ids = np.unique(np.array( [atom.parent.id[1] for atom in matched_atoms_indexs]))
