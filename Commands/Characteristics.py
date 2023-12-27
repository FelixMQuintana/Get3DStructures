import abc
import logging
import os
import re
import subprocess
import threading
from abc import ABC
from typing import List
import itertools
from Bio.PDB import PDBIO, Selection
from matplotlib import pyplot as plt
from scipy.spatial import distance_matrix
import numpy as np
import pandas
from Bio import PDB, motifs, Seq
from Commands.command import Command, FactoryBuilder
from pathlib import Path
import tqdm
from Bio.Blast import NCBIWWW
from dataStructure.collections import Collection, HomologyStructureFetcher, ExperimentalStructureFetcher, \
    UniProtAcessionFetcher, UniProtFastaFetcher
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure import StructureFile, HomologyStructure
from lib.const import StructureCharacteristicsMode, MotifSearchMode, MotifRefinements, AllowedExt
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data


# import transformers

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
                                 int(id.id[1]) not in protein_structures.uniprotID.binding_site_residues]
        # and id.id[1] not in pred_binding_sites]
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
            structure_file.path.name.removesuffix(".pdb") + "_motif.pdb")),
            preserve_atom_numbering=True)


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


class GetLineage(Characteristics):

    def __init__(self, binding_site_database: Path):
        super().__init__(binding_site_database)
        self.collection = Collection(self.working_directory, UniProtAcessionFetcher(), ExperimentalStructureFetcher())

    def run(self) -> None:
        structs = [protein_structures for protein_structures in self.collection.protein_structure_results.values()]
        # print(structs)
        out = open("/media/felix/Research/protein_list2.txt", "w")
        with open("/media/felix/Research/protein_list.txt", "r") as e:
            for struct in structs:
                #   e.write(struct.id)
                print(struct.id)
                try:
                    struct.uniprotID.structural_data["organism"]["lineage"]
                except Exception:
                    continue
                if struct.uniprotID.structural_data["organism"]["lineage"][-1] == "Escherichia":
                    out.write(struct.id)
                    out.write("\n")

    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures) -> None:
        pass


# def command(self, protein_structure: ProteinStructures) -> None:
#      with open("/mnt/ResearchLongTerm/FunSoCTrainingData/lineages.txt", "a") as e:
#        e.write(protein_structure.id)
#        [e.write(" " + item) for item in protein_structure.uniprotID.structural_data["organism"]["lineage"]]
#        e.write("\n")

class ClusterGOTerms(Command):

    def __init__(self, uniprot_list: Path):
        super().__init__()
        self.collection = Collection(self.working_directory, UniProtAcessionFetcher())
        self.protein_list: Path = uniprot_list

    def run(self) -> None:
        protein_dict = {}
        for protein_structure in self.collection.protein_structure_results.values():
            for term in protein_structure.uniprotID.structural_data['uniProtKBCrossReferences']:
                if term['database'] == "GO":
                    try:
                        protein_dict[term['id']].append(protein_structure.uniprotID.id)
                    except KeyError as ex:
                        protein_dict[term['id']] = [protein_structure.uniprotID.id]

        with open(self.protein_list, "w") as output_file:
            output_file.write(str(protein_dict))


class TMScoreDatabase(Command):

    def __init__(self):
        super().__init__()
        self.collection = Collection(self.working_directory, ExperimentalStructureFetcher())

    def run(self) -> None:
        for protein_structures in self.collection.protein_structure_results.values():
            for protein_structures_inner_loop in self.collection.protein_structure_results.values():
                get_tm(protein_structures.crystal_structures[0].path,
                       protein_structures_inner_loop.crystal_structures[0].path)


class FindCustomBindingSite(Command):

    def __init__(self):
        super().__init__()
        self.collection = Collection(self.working_directory, ExperimentalStructureFetcher(), UniProtFastaFetcher())

    def run(self) -> None:
        pdb_parser = PDB.PDBParser()
        io = PDBIO()
        for protein_structures in self.collection.protein_structure_results.values():
            #  motifs_of_interest = [Seq.Seq("GXXGXGKST")]#, Seq.Seq("XXXXD")]

            try:
                print(protein_structures.all_structures[0].path)
                stability_motif = re.search(r"Y[A-Za-z]{2}D", str(protein_structures.crystal_structures[0].fasta))
                binding_motif = re.search(r"G[A-Za-z]{2}G[A-Za-z]GKST",
                                          str(protein_structures.crystal_structures[0].fasta))
                histine_pos = re.finditer(r"H", str(protein_structures.crystal_structures[0].fasta))
                #   if 0 >len(stability_motif) < 2:
                stability_motif = list(range(stability_motif.start() + 1, stability_motif.end() + 1))
                # elif len(stability_motif) > 0:
                #        stability_motif = list(range(stability_motif.start(), stability_motif.end()))

                binding_motif = list(range(binding_motif.start() + 1, binding_motif.end() + 1))
                match_motif = None
                for match in histine_pos:
                    if match.end() < 185:
                        match_motif = match

                histine_pos = list(range(match_motif.start() + 1, match_motif.end() + 1))
            except AttributeError as ex:
                print(f"Skipping because of {ex}")
                continue
            except IndexError as ex:
                print(f"Indexing error because of {ex}")
                continue
            stability_motif.extend(binding_motif)
            stability_motif.extend(histine_pos)
            working_dir: Path = self.working_directory.joinpath(protein_structures.all_structures[0].id)
            pdb = pdb_parser.get_structure(protein_structures.all_structures[0].path.name,
                                           protein_structures.all_structures[0].path)
            # pdb_cluster_rep = pdb_parser.get_structure(id ,self.binding_site_database.joinpath(id).joinpath(
            #                                              structure_file.path.name.split(".")[0]).joinpath(structure_file.path.name.split(".")[0] + "_cl_001.pdb",
            #                                          ))
            pdb_ids = [residue.id[1] for residue in pdb.get_residues()]
            for chain in pdb[0]:
                [chain.detach_child((' ', id, ' ')) for id in pdb_ids if id not in stability_motif]
            # pdb_coords = np.array([residue.center_of_mass() for residue in pdb.get_residues()])
            io.set_structure(pdb)
            io.save(str(self.working_directory.joinpath(
                protein_structures.crystal_structures[0].path.name + "_pocket.pdb")))
        # d_matrix = distance_matrix(pdb_coords, cluster_coords)

        #  m = motifs.create(motifs_of_interest,alphabet="GAVLITSMCPFYWHKRDENQ")
        #  for motif in m.search_instances(protein_structures.fasta):
        #       print(motif)
        # matches = re.finditer(r"")


class CreateCombinatorics(Characteristics):
    def __init__(self, binding_site_database: Path):
        super().__init__(binding_site_database)

    def command(self, structure_file: StructureFile, protein_structure: ProteinStructures) -> None:
        pdb_parser = PDB.PDBParser()
        #  working_dir: Path = self.binding_site_database.joinpath(structure_file.id)

        pdb = pdb_parser.get_structure(structure_file.path.name, structure_file.path)
        residues = [res for res in pdb.get_residues()]
        iter_substructures = itertools.combinations(residues, r=3)
        pdb_ids = [residue.id[1] for residue in pdb.get_residues()]

        for index, substructure in enumerate(iter_substructures):
            io = PDBIO()
            substructure_ids = [res_id.id[1] for res_id in substructure]
            for chain in pdb[0]:
                [chain.detach_child((' ', id, ' ')) for id in pdb_ids if id not in substructure_ids]
            io.set_structure(pdb)
            io.save(str(self.binding_site_database.joinpath(
                structure_file.path.name.split(".")[0] + "_" + str(index) + ".pdb")))
            pdb = pdb_parser.get_structure(structure_file.path.name, structure_file.path)

            #  for chain in pdb[0]:
    #      [chain.detach_child((' ', id, ' ')) for id in residue_ids_to_remove]


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


class CalculateSubstructureDistances(Command):

    def __init__(self, junk):
        super().__init__()
        self.collection = Collection(self.working_directory, ExperimentalStructureFetcher())

    def run(self) -> None:

        for protein_structures in self.collection.protein_structure_results.values():
            distances = np.empty((len(protein_structures.all_structures),
                                  len(protein_structures.all_structures)),dtype=np.float16)
            for index, protein_struct in enumerate(tqdm.tqdm(protein_structures.all_structures)):
                for index2, protein_structures_inner in enumerate(protein_structures.all_structures):

                    distances[index][index2] = self.calculate_vector(protein_struct,
                                                                     protein_structures_inner)
       # my_file = open("substructuct_distance.txt","w")
        np.save(file= "/home/felix/substructuct_distance.txt",arr= distances)

    def calculate_vector(self, structure1, structure2):
        s1 = get_structure(structure2.path)
        chain = next(s1.get_chains())
        coords, seq = get_residue_data(chain)
        s2 = get_structure(structure1.path)
        chain2 = next(s2.get_chains())
        coords2, seq2 = get_residue_data(chain2)
        alignment = tm_align(coords, coords2, seq, seq2)
        aligned_seq_1 = coords.dot(alignment.u.T) + alignment.t
        distances = self.calculate_distance(aligned_seq_1, coords2)
        return distances

    def calculate_distance(self, coord1, coord2):
        return np.sqrt(((coord1 - coord2) ** 2).mean())


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


class CalculateDistance(Characteristics):

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
                self.per_motif_rmsd_vectors.append(np.sqrt((pdb_coords_reference - pdb_coords) ** 2 /
                                                           pdb_coords_reference.shape[0]))


supported_commands = {
    StructureCharacteristicsMode.AUTOSITE: FindPockets,
    StructureCharacteristicsMode.BUILD_MOTIFS: BuildMotifStructures,
    StructureCharacteristicsMode.FIND_BINDING_POCKETS: FindBindingSite,
    StructureCharacteristicsMode.CHECK_MOTIF_QUALITY: CheckBindingSiteQuality,
    StructureCharacteristicsMode.CALCULATE_RMSD: CalculateSubstructureDistances,
    StructureCharacteristicsMode.TRIM_PDB: FixPDBFiles,
    StructureCharacteristicsMode.GET_LINEAGE: GetLineage,
    StructureCharacteristicsMode.CREATE_COMBINATORICS: CreateCombinatorics
}


class CharacteristicsFactory(FactoryBuilder):

    @staticmethod
    def build(mode: StructureCharacteristicsMode, *args, **kwargs):
        return supported_commands[mode](*args, **kwargs)

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
