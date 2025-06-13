import abc
import json
import logging
import os
import time
from abc import ABC
from copy import deepcopy
from pathlib import Path
from subprocess import Popen, PIPE
from typing import Optional

import Bio.PDB
import numpy
import numpy as np
import pandas
import pandas as pd
import pymesh
import scipy
from Bio.PDB import PDBIO, PDBParser
from scipy.spatial import cKDTree
from simpleicp import PointCloud, SimpleICP
from tqdm import tqdm
from pymol import cmd

from dataStructure.protein.accession import UniProtAcessionFile
from dataStructure.protein.structure.representation.atomic.atomic import Atomic
from dataStructure.protein.structure.representation.surface.surface import MolecularSurfaceRepresentation
from lib.compute_polar_coordinates import output_patch_coords
from lib.const import Metrics, CoordinateType, CalculateOptions, masif_opts
from Commands.command import Command, FactoryBuilder
from dataStructure.protein.structure.structure import StructureFile
from dataStructure.collections import ExperimentalStructureFetcher, HomologyStructureFetcher, Collection, \
    UniProtAcessionFetcher
from lib.read_data_from_surface import read_data_from_surface
from lib.read_ply import read_ply
from proteinAnalysis.bonds.bonds import BondProteinAnalysis
from proteinAnalysis.electrostatics._helper import assignChargesToNewMesh
from proteinAnalysis.electrostatics.electrostatic_charge import ElectrostaticChargeProteinAnalysis
from proteinAnalysis.electrostatics.electrostatic_potential import ABPS
from proteinAnalysis.hydrophobic.hydrophobic import HydrophobicityProteinAnalysis


class Calculate(Command, ABC):
    """

    """


fetcher_type = {
    Metrics.HOMOLOGY_STRUCTURES: HomologyStructureFetcher,
    Metrics.CRYSTAL_STRUCTURES: ExperimentalStructureFetcher

}


class SimilarityScore(Calculate, ABC):

    def __init__(self, structure_type: Metrics, coord_type: CoordinateType, out_file_name: str = None):
        super().__init__()
        if out_file_name is None:
            out_file_name = self.working_directory.joinpath("similarityArray" + str(time.time()))
        self.collection = Collection(self.working_directory, fetcher_type[structure_type]())
        self._distances: Optional[np.ndarray] = None
        self.coord_type = coord_type
        self._out_file_name: str = out_file_name

    @abc.abstractmethod
    def metric(self, structure_file1: StructureFile, structure_file2: StructureFile, coord_type: CoordinateType,
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
        ordering = open(self.working_directory.joinpath("ordering_of_distances" + str(time.time()) + ".txt"), "w")

        for index, structures in tqdm(enumerate(self.collection.protein_structure_results.values())):
            for structure in structures.all_structures:
                ordering.write(structure.id + "\n")

                for inner_index, structures_inner in enumerate(self.collection.protein_structure_results.values()):
                    for structure_inner in structures_inner.all_structures:
                        self.distances_matrix[index, inner_index] = self.metric(structure, structure_inner,
                                                                                self.coord_type)
        np.save(self._out_file_name, self.distances_matrix)
        ordering.close()


class CalculateMeshes(Calculate):

    def __init__(self,working_dir, structure):
        super().__init__(working_dir,structure)
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())

    def run(self) -> None:
        p = []
        from multiprocessing import Pool
        pool = Pool(20)
        for sturctures in tqdm(self.collection.protein_structure_results.values()):
            for structure in sturctures.all_structures:
                if structure.path.parent.joinpath(structure.path.name.split(".")[0] + ".ply").exists():
                    continue
                p.append(pool.apply_async(CalculateMeshes.metric, args=[structure]))
        [process.wait() for process in p]

    @staticmethod
    def metric(structure):
        logging.info(f"Generating molecular surface for {structure.path.name} in path {structure.path}")
        molecular_surface = MolecularSurfaceRepresentation(structure.path, resolution=1.0)  # needs to be float
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
        molecular_surface.mesh.add_attribute("hphob")
        molecular_surface.mesh.set_attribute("hphob", hydrophobic_verts)
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
        pymesh.save_mesh(
            str(structure.path.parent.joinpath(structure.path.name.split(".")[0] + ".ply")), molecular_surface.mesh,
            *molecular_surface.mesh.get_attribute_names(), use_float=True, ascii=True
        )


class CalculateMeshes_need_to_rename(Calculate):

    def __init__(self):
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())
        custom_collection = Collection(self.working_directory.parent, UniProtAcessionFetcher())
        cluster_rep = self.working_directory.parent.name

        self.residues_of_interest = set(
            custom_collection.protein_structure_results[cluster_rep].uniprotID.binding_site_residues)

    #  parser = Bio.PDB.PDBParser()
    #  structure=parser.get_structure('rep',custom_collection.protein_structure_results[cluster_rep].path)
    #  for res in structure.get_residues():
    #      print(res)

    def run(self) -> None:
        for structures in self.collection.protein_structure_results.values():
            for structure in structures.all_structures:
                count = 0
                parser = Bio.PDB.PDBParser()
                structure_file = parser.get_structure('rep', structure.path)
                for residue in structure_file.get_residues():
                    count += 1
                    print(self.residues_of_interest)
                total_residues = [i for i in range(1, count + 1)]
                total_residues = set(total_residues)
                res_to_remove = total_residues.difference(self.residues_of_interest)
                #              atom_vectors[index] = [atom.get_vector() for atom in residue.get_atoms()]

                for chain in structure_file[0]:
                    [chain.detach_child((' ', id, ' ')) for id in res_to_remove]
                io = PDBIO()
                io.set_structure(structure_file)
                io.save(str(self.working_directory.joinpath(structure.path.name + "-contacts-only.pdb")),
                        preserve_atom_numbering=True)


class CalculateSurfaceOfInterest(Calculate):

    def __init__(self, output_dir: Path, ):  # residues_of_interest: []) -> None:
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher(), UniProtAcessionFetcher())
        # custom_collection = Collection(self.working_directory.parent, UniProtAcessionFetcher())

        cluster_rep = self.working_directory.name.split('_')[-1]
        # self.residues_of_interest = self.collection.protein_structure_results[cluster_rep].uniprotID.binding_site_residues

        # self.residues_of_interest = custom_collection.protein_structure_results[cluster_rep].uniprotID.binding_site_residues
        # logging.info(f"Residues of interest: {self.residues_of_interest}")

        self.output_dir: Path = Path(output_dir)
        #  self.residues_of_interest: list = [343]
        parser = Bio.PDB.PDBParser()
        # structure=parser.get_structure('rep',self.collection .protein_structure_results[cluster_rep].path)
        # for res in structure.get_residues():
        #    print(res)
        for structures in self.collection.protein_structure_results.values():
            for structure in structures.all_structures:
                if not self.output_dir.joinpath(structure.id).exists():
                    os.mkdir(self.output_dir.joinpath(structure.id))

    # def get_idx(self, rep_res):

    @staticmethod
    def get_secondary_structure(structure_file: StructureFile, compute=False) -> pd.DataFrame:
        if compute:
            args = [
                "mkdssp",
                str(structure_file.path.name),
                str(structure_file.path.name).split(".")[0] + ".dssp"
            ]
            p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(structure_file.path.parent))
            stdout, stderr = p2.communicate()
        try:
            df = pd.read_csv(
                str(structure_file.path.parent) + "/" + str(structure_file.path.name).split(".")[0] + ".dssp",
                header=None, skiprows=28)
            hits = []
            for entry in df[0]:
                if len(entry.split()) < 10:
                    continue
                split = entry.split()
                hits.append([split[0:5]])
        except FileNotFoundError as ex:
            raise ex
        return pd.DataFrame(hits)

    def run(self):
        jobs = []
        for structures in self.collection.protein_structure_results.values():
            for structure in structures.all_structures:
                CalculateSurfaceOfInterest.compute(structures.uniprotID.binding_site_residues, structure)
                # self.get_idx(structure)
                # jobs.append(self.thread_pool.apply_async(CalculateSurfaceOfInterest.compute,
                #                                         args=[self.residues_of_interest, structure]))
        # var = [job.wait() for job in jobs]

    @staticmethod
    def compute(residues_of_interest, structure):
        # residues_of_interest = [90,91,92, 251, 256, 257]
        if len(residues_of_interest) < 1:
            return
        # print(residues_of_interest)
        if structure.path.parent.joinpath(
                structure.path.name.split(".")[0] + "-full-patches-with-vertex-pos-no-cut-off.npy").exists():
            return
        if not structure.path.parent.joinpath(structure.path.name.split(".")[0] + ".ply").exists():
            return
        mesh_vertices, mesh_faces, normals, _, _, _, _ = read_ply(str(structure.path.parent.joinpath(
            structure.path.name.split(".")[0] + ".ply")))
        # hits_of_interest =[53,54,57,58,68,69,79,80,81,84,85,118,119,129,130, 192, 282, 339,340] # [66,68,84,276,327, 328,331]

        sub_vertices, sub_faces, sub_normals, vertex_idx = CalculateSurfaceOfInterest.find_verticies(
            residues_of_interest, structure)  # hits_of_interest, structure)
        input_feat, rho, theta, mask, neigh_indices, iface_labels, verticies = read_data_from_surface(
            str(structure.path.parent.joinpath(
                structure.path.name.split(".")[0] + ".ply")), masif_opts["site"], vertex_idx)
        # CalculateSurfaceOfInterest.debug(sub_faces,np.array(neigh_indices),vertex_idx,verticies,sub_normals,theta,rho,name=structure.path.parent.name)

        numpy.save(structure.path.parent.joinpath(
            structure.path.name.split(".")[0] + "-full-patches-with-vertex-pos-no-cut-off"),
                   input_feat)

    @staticmethod
    def debug(faces, neigh_indices, vertex, vertices, normals, theta, rho, name):
        faces_of_interest = []
        for x, y, z in faces:
            if x in neigh_indices[vertex] and y in neigh_indices[vertex] and z in neigh_indices[vertex]:
                faces_of_interest.append([x, y, z])
        faces_old_vertex = np.array(faces_of_interest)
        old_vertex_to_new_dict = {}
        count = 0
        new_verticies = []
        new_normals = []
        for element in neigh_indices[vertex]:
            for ele in element:
                # old_vertex_to_new_dict[element] = count
                old_vertex_to_new_dict[ele] = count
                new_verticies.append(vertices[ele])
                new_normals.append(normals[ele])
                count += 1
        new_vertex_list = []
        for x, y, z in faces_old_vertex:
            new_vertex_list.append(
                [old_vertex_to_new_dict[x], old_vertex_to_new_dict[y], old_vertex_to_new_dict[z]])
        new_faces = np.array(new_vertex_list)
        # output_patch_coords(vertices[neigh_indices[vertex]], new_faces, normals[neigh_indices[vertex]], vertex,
        #                     neigh_indices[vertex], theta, rho)
        output_patch_coords(np.unique(np.array(new_verticies), axis=0), new_faces,
                            np.unique(np.array(new_normals), axis=0), vertex,
                            np.unique(neigh_indices[vertex].flatten()), theta, rho, optional_name=name)

    @staticmethod
    def find_verticies(hits_of_interest, structure):
        mesh_vertices, mesh_faces, normals, _, _, _, _ = read_ply(str(structure.path.parent.joinpath(
            structure.path.name.split(".")[0] + ".ply")))
        protein_structure = open(structure.path, "r").readlines()
        closest_vertices = []
        closest_faces = []
        for line in protein_structure:
            if line.startswith("ATOM"):
                if int(line.split()[5]) in hits_of_interest:
                    coords = np.array(line.split()[6:9])
                    coords = coords.astype(np.float64)
                    distance_arry = np.linalg.norm(coords - mesh_vertices, axis=1, )
                    # if np.min(distance_arry) > 5:
                    #    continue
                    closest_vertex_idx = distance_arry.argmin()
                    cond = mesh_faces[:, 0] == closest_vertex_idx
                    faces_idx = cond.nonzero()[0]
                    closest_vertices.append(closest_vertex_idx)
                    closest_faces.append(faces_idx)
        closest_vertices = np.array(closest_vertices)
        # closest_faces = np.concatenate(closest_faces)
        return mesh_vertices[closest_vertices], mesh_faces, normals, closest_vertices


from scipy.spatial.distance import cdist


class ScoreSurfacePatches(Calculate):

    def __init__(self,working_dir,structure, output_dir: Path) -> None:
        super().__init__(working_dir,structure)
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())
        self.cluster_rep = self.working_directory.parent.name.split('_')[-1]

    def run(self) -> None:
        rep_structure = self.collection.protein_structure_results[self.cluster_rep].all_structures[0]
        if not rep_structure.path.parent.joinpath(
                rep_structure.path.name.split(".")[0] + "-nearest-cavities-vertex-pos.npy").exists():
            return
        if rep_structure.path.parent.parent.joinpath("scores_dict_cavs.json").exists():
            return
        ref = np.load(str(rep_structure.path.parent.joinpath(
            rep_structure.path.name.split(".")[0] + "-nearest-cavities-vertex-pos.npy")))
        scores = {}
        for structures in tqdm(self.collection.protein_structure_results.values()):
            for structure in structures.all_structures:
                if structure.path.parent.joinpath(
                        structure.path.name.split(".")[0] + "-nearest-cavities-vertex-pos.npy").exists():
                    target = np.load(str(structure.path.parent.joinpath(
                        structure.path.name.split(".")[0] + "-nearest-cavities-vertex-pos.npy")))
                    scores[structure.id] = ScoreSurfacePatches.jaccard_similarity(ref, target)
        with open(rep_structure.path.parent.parent.joinpath("scores_dict_cavs.json"), 'w') as fileio:
            json.dump(scores, fileio)

    @staticmethod
    def jaccard_similarity(ref, target):

        comparison = ScoreSurfacePatches.build_set_comparison(ref, target)
        a_intersection_b = np.count_nonzero(
            np.linalg.norm(comparison[:, [0, 1, 2, 3]] - comparison[:, [10, 11, 12, 13]], axis=1) <= 4)
        # print(a_intersection_b)
        a_or_b = np.unique(ref, axis=0).shape[0] + np.unique(target, axis=0).shape[0] - a_intersection_b
        jaccrad_sim = a_intersection_b / np.unique(ref, axis=0).shape[0]  # a_or_b
        return jaccrad_sim

    @staticmethod
    def build_set_comparison(reference, target):
        reference = np.unique(reference, axis=0)
        target = np.unique(target, axis=0)
        # if len(reference) > len(target):
        #    reference, target = target, reference
        distance_matrix = cdist(reference[:, [4, 5, 6]], target[:, [4, 5, 6]])
        minimum_distance_index = np.argmin(distance_matrix, axis=1)
        hits, counts = np.unique(minimum_distance_index, return_counts=True)
        counts_indx = np.where(counts > 1)[0]
        neighbors = np.ones(shape=(reference.shape[0]), dtype=int, )
        neighbors = neighbors * -1
        exaustive_check_count = 0
        while (len(counts_indx) > 0 and exaustive_check_count < 5):
            dups = hits[counts_indx]
            for duplicate in dups:
                dups_indx = np.where(duplicate == minimum_distance_index)[0]
                indx_to_increase_value = np.delete(dups_indx, np.argmin(distance_matrix[dups_indx, duplicate]))
                distance_matrix[indx_to_increase_value, duplicate] += 1000
            minimum_distance_index = np.argmin(distance_matrix, axis=1)
            hits, counts = np.unique(minimum_distance_index, return_counts=True)
            counts_indx = np.where(counts > 1)[0]
            exaustive_check_count += 1

        indx = np.array([i for i in range(reference.shape[0])])
        restriction = distance_matrix[indx, minimum_distance_index] < 3
        restriction_indx = np.where(restriction == True)
        neighbors[restriction_indx] = minimum_distance_index[restriction_indx]
        no_match_indx = np.where(neighbors == -1)[0]
        matches = target[neighbors.T].reshape(neighbors.T.shape[0],
                                              1 * 10)
        matches[no_match_indx] = [-10, -10, -10, -10, -10, -10, -10, -10, -10, -10]  # reference[no_match_indx]
        feature = np.concatenate((reference[indx], matches), axis=1)
        return feature


class AlignMeshes(Calculate):
    def __init__(self, working_dir,structure, output_dir: Path) -> None:
        super().__init__(working_dir, structure)
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())

    def run(self) -> None:
        count = 0
        for structures in tqdm(self.collection.protein_structure_results.values()):
            for structure in structures.all_structures:
                structure_mesh = pymesh.load_mesh(str(structure.path.with_suffix(".ply")))
                verts = structure_mesh.vertices
                # num_patches = patches_from_structure.shape[0]
                # num_vert = 100
                # patches_from_structure =patches_from_structure.reshape(num_patches*num_vert, 10)
                # patches_from_structure = np.unique(patches_from_structure, axis=0)

                if count == 0:
                    pc_fix = PointCloud(verts,
                                        columns=["x", "y", "z"])
                    count += 1
                    mesh = pymesh.form_mesh(structure_mesh.vertices, structure_mesh.faces)
                    mesh.add_attribute("vertex_charges")
                    mesh.set_attribute("vertex_charges", structure_mesh.get_attribute("vertex_charges"))

                    mesh.add_attribute("vertex_hphob")
                    mesh.set_attribute("vertex_hphob", structure_mesh.get_attribute("vertex_hphob"))

                    mesh.add_attribute("vertex_hbond")
                    mesh.set_attribute("vertex_hbond", structure_mesh.get_attribute("vertex_hbond"))

                    # mesh.add_attribute("vertex_nx")
                    # mesh.set_attribute("vertex_nx", structure_mesh.get_attribute("vertex_nx"))
                    # mesh.add_attribute("vertex_ny")
                    # mesh.set_attribute("vertex_ny", structure_mesh.get_attribute("vertex_ny"))
                    # mesh.add_attribute("vertex_nz")
                    # mesh.set_attribute("vertex_nz", structure_mesh.get_attribute("vertex_nz"))
                    pymesh.save_mesh(str(structure.path.parent.joinpath(
                        structure.path.name.split(".")[0] + "-aligned.ply")), mesh,
                        *mesh.get_attribute_names(), use_float=True, ascii=True)
                    continue

                icp = SimpleICP()
                pc_mov = PointCloud(verts,
                                    columns=["x", "y", "z"])
                icp.add_point_clouds(pc_fix, pc_mov)
                H, X_mov_transformed, rigid_body_transformation_params, distance_residuals = icp.run()

                mesh = pymesh.form_mesh(X_mov_transformed, structure_mesh.faces)
                mesh.add_attribute("vertex_charges")
                mesh.set_attribute("vertex_charges", structure_mesh.get_attribute("vertex_charges"))

                mesh.add_attribute("vertex_hphob")
                mesh.set_attribute("vertex_hphob", structure_mesh.get_attribute("vertex_hphob"))

                mesh.add_attribute("vertex_hbond")
                mesh.set_attribute("vertex_hbond", structure_mesh.get_attribute("vertex_hbond"))

                # mesh.add_attribute("vertex_nx")
                # mesh.set_attribute("vertex_nx", structure_mesh.get_attribute("vertex_nx"))
                # mesh.add_attribute("vertex_ny")
                # mesh.set_attribute("vertex_ny", structure_mesh.get_attribute("vertex_ny"))
                # mesh.add_attribute("vertex_nz")
                # mesh.set_attribute("vertex_nz", structure_mesh.get_attribute("vertex_nz"))
                pymesh.save_mesh(
                    str(structure.path.parent.joinpath(
                        structure.path.name.split(".")[0] + "-aligned.ply")), mesh,
                    *mesh.get_attribute_names(), use_float=True, ascii=True
                )

        #     structure_mesh.vertices = X_mov_transformed

        #     pymesh.save_mesh(structure.path.parent.joinpath(
        #         structure.path.name.split(".")[0] + "-aligned.ply"), structure_mesh,
        #         *structure_mesh.mesh.get_attribute_names(), use_float=True, ascii=True)


class CalculateICP(Calculate):

    def __init__(self, output_dir: Path) -> None:
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())

    def run(self) -> None:
        count = 0
        for structures in tqdm(self.collection.protein_structure_results.values()):
            for structure in structures.all_structures:
                patches_from_structure = np.load(
                    structure.path.parent.joinpath(structure.path.name.split(".")[0] + "-patches-with-pos-n-ddc.npy"))
                num_patches = patches_from_structure.shape[0]
                num_vert = 100
                # patches_from_structure =patches_from_structure.reshape(num_patches*num_vert, 10)
                # patches_from_structure = np.unique(patches_from_structure, axis=0)

                if count == 0:
                    pc_fix = PointCloud(patches_from_structure[:, [4, 5, 6, 7, 8, 9]],
                                        columns=["x", "y", "z", "nx", "ny", "nz"])
                    count += 1
                    np.save(structure.path.parent.joinpath(
                        structure.path.name.split(".")[0] + "-patches-with-pos-n-ddc-aligned.npy"),
                        patches_from_structure)
                    continue

                icp = SimpleICP()
                pc_mov = PointCloud(patches_from_structure[:, [4, 5, 6, 7, 8, 9]],
                                    columns=["x", "y", "z", "nx", "ny", "nz"])
                icp.add_point_clouds(pc_fix, pc_mov)
                H, X_mov_transformed, rigid_body_transformation_params, distance_residuals = icp.run()
                patches_from_structure[:, [4, 5, 6]] = X_mov_transformed

                np.save(
                    structure.path.parent.joinpath(
                        structure.path.name.split(".")[0] + "-patches-with-pos-n-ddc-aligned.npy"),
                    patches_from_structure)

                # dataset.append(patches_from_structure)


class CompleteWorkFlow(Calculate):

    def __init__(self, output_dir: Path) -> None:
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())

    def run(self) -> None:
        k = 2
        aligned_tuples = []
        for structures in tqdm(self.collection.protein_structure_results.values()):
            for structure_outer in structures.all_structures:
                patches_from_structure_outer = np.load(
                    structure_outer.path.parent.joinpath(
                        structure_outer.path.name.split(".")[0] + "with-vertex-pos-aligned.npy"))
                pc_fix = PointCloud(patches_from_structure_outer[:, [4, 5, 6, 7, 8, 9]],
                                    columns=["x", "y", "z", "nx", "ny", "nz"])
                features = np.zeros((len(self.collection.protein_structure_results.values()),
                                     patches_from_structure_outer.shape[0], (k + 1) * 10))
                feature_idx = 0
                for structures_inner in self.collection.protein_structure_results.values():
                    for structure_inner in structures_inner.all_structures:
                        patches_from_structure_inner = np.load(
                            structure_inner.path.parent.joinpath(
                                structure_inner.path.name.split(".")[0] + "with-vertex-pos-aligned.npy"))
                        #    found=False
                        #    for entry in aligned_tuples:
                        #        if entry[0] == structure_inner.id and entry[1] == structure_outer.id:
                        #            patches_from_structure_inner[:,[4,5,6]] = entry[2]
                        #            found=True
                        #        elif entry[1] == structure_inner.id and entry[0] == structure_outer.id:
                        #            patches_from_structure_inner[:,[4,5,6]] = entry[3]
                        #            found=True
                        #    if not found:
                        #        icp = SimpleICP()
                        #        pc_mov = PointCloud(patches_from_structure_inner[:, [4, 5, 6, 7, 8, 9]],
                        #                        columns=["x", "y", "z", "nx", "ny", "nz"])
                        #        icp.add_point_clouds(pc_fix, pc_mov)
                        #        H, X_mov_transformed, rigid_body_transformation_params, distance_residuals = icp.run()
                        #        aligned_tuples.append ((structure_outer.id,structure_inner.id,patches_from_structure_outer[:, [4,5,6]], X_mov_transformed))
                        #        patches_from_structure_inner[:,[4,5,6]] = X_mov_transformed
                        # aligned_dict[structure_outer.id].append(structure_inner.id)
                        # aligned_dict_vals[structure_outer.id].append((structure_inner.id ,X_mov_transformed))
                        # for entry in aligned_dict[structure_outer.id]:
                        #     if entry[0] == structure_inner.id:
                        #         patches_from_structure_inner[:,[4,5,6]]=entry[1]
                        #         break
                        distance_matrix = scipy.spatial.distance.cdist(patches_from_structure_outer[:, [4, 5, 6]],
                                                                       patches_from_structure_inner[:, [4, 5, 6]])
                        neighbors = np.ones(shape=(k, patches_from_structure_outer.shape[0]), dtype=int, )
                        neighbors = neighbors * -1
                        indx = np.array([i for i in range(patches_from_structure_outer.shape[0])])
                        for neighbor in range(0, k):
                            minimum_distance_index = np.argmin(distance_matrix, axis=1)
                            # neighbors[neighbor] = minimum_distance_index
                            restriction = distance_matrix[indx, minimum_distance_index] < 1
                            restriction_indx = np.where(restriction == True)
                            neighbors[neighbor, restriction_indx] = minimum_distance_index[restriction_indx]
                            #    neighbors[neighbor] = -1
                            #    break
                            distance_matrix[indx, minimum_distance_index] = 1000000
                        # TODO Fix later here. You are cheating by overridingf the value of the last vertex to be used as a catch all for errors.
                        patches_from_structure_inner[-1] = [-10, -10, -10, -10, -10, -10, -10, -10, -10, -10]
                        matches = patches_from_structure_inner[neighbors[:, indx].T].reshape(
                            neighbors[:, indx].T.shape[0],
                            k * 10)
                        feature = np.concatenate((patches_from_structure_outer[indx], matches), axis=1)
                        features[feature_idx] = feature
                        feature_idx += 1
                np.save(
                    structure_outer.path.parent.joinpath(structure_outer.path.name.split(".")[0] + "-featurized.npy"),
                    features)


class CalculateKNN(Calculate):

    def __init__(self, output_dir: Path) -> None:
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())

    def run(self) -> None:
        k = 2
        num_proteins = len(self.collection.protein_structure_results.values())
        dist_matrix = np.zeros((num_proteins, num_proteins - 1))
        for structures in tqdm(self.collection.protein_structure_results.values()):
            for structure in structures.all_structures:
                patches_from_structure = np.load(
                    structure.path.parent.joinpath(structure.path.name.split(".")[0] + "with-vertex-pos-aligned.npy"))
                features = np.zeros((len(self.collection.protein_structure_results.values()) - 1,
                                     patches_from_structure.shape[0], (k + 1) * 10))
                feature_idx = 0
                for structures in self.collection.protein_structure_results.values():
                    for structure_inner in structures.all_structures:
                        if structure_inner.path.name == structure.path.name:
                            continue
                        patches_from_structure2 = np.load(structure_inner.path.parent.joinpath(
                            structure_inner.path.name.split(".")[0] + "with-vertex-pos-aligned.npy"))
                        distance_matrix = scipy.spatial.distance.cdist(patches_from_structure[:, [4, 5, 6]],
                                                                       patches_from_structure2[:, [4, 5, 6]])
                        neighbors = np.ones(shape=(k, patches_from_structure.shape[0]), dtype=int, )
                        neighbors = neighbors * -1
                        indx = np.array([i for i in range(patches_from_structure.shape[0])])
                        for neighbor in range(0, k):
                            minimum_distance_index = np.argmin(distance_matrix, axis=1)
                            # neighbors[neighbor] = minimum_distance_index
                            restriction = distance_matrix[indx, minimum_distance_index] < 1.5
                            restriction_indx = np.where(restriction == True)
                            neighbors[neighbor, restriction_indx] = minimum_distance_index[restriction_indx]
                            #    neighbors[neighbor] = -1
                            #    break
                            distance_matrix[indx, minimum_distance_index] = 1000000
                        # TODO Fix later here. You are cheating by overridingf the value of the last vertex to be used as a catch all for errors.
                        #    patches_from_structure2[-1] = [-10, -10, -10, -10, -10, -10, -10, -10, -10, -10]
                        # relu or linear scaling
                        matches = patches_from_structure2[neighbors[:, indx].T].reshape(neighbors[:, indx].T.shape[0],
                                                                                        k * 10)
                        feature = np.concatenate((patches_from_structure[indx], matches), axis=1)
                        features[feature_idx] = feature
                        feature_idx += 1

                np.save(structure.path.parent.joinpath(structure.path.name.split(".")[0] + "-featurized.npy"), features)

                #  for vertex in patches_from_structure:

                #      distances = np.linalg.norm(vertex-patches_from_structure2, axis=1)


class CalculateComparison(Calculate):

    def __init__(self, output_dir: Path) -> None:
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())

    def run(self) -> None:
        num_proteins = len(self.collection.protein_structure_results.values())
        dist_matrix = np.zeros((num_proteins, num_proteins))
        for index, structures in tqdm(enumerate(self.collection.protein_structure_results.values())):
            for structure in structures.all_structures:
                structure_raw_knn_matches = np.load(
                    structure.path.parent.joinpath(structure.path.name.split(".")[0] + "-featurized.npy"))
                for index_inner, protein in enumerate(structure_raw_knn_matches):
                    # vertex_norms = []
                    # for vertex in protein:
                    #   if not (vertex[[10, 11, 12, 13]] == -10).all():
                    #       vertex_norms.append(np.linalg.norm(vertex[[0,1,2,3]] - vertex[[10,11,12,13]]))
                    dist = np.linalg.norm(protein[:, [0, 1, 2, 3]] - protein[:, [10, 11, 12, 13]], axis=1)
                    dist_matrix[index][index_inner] = np.mean(dist)
                    # dist_matrix[index][index_inner] = np.mean(np.array(vertex_norms))
        np.save(str(structure.path.parent.parent) + "/dist_matrix_k_2_mean", dist_matrix)


class CalculateContacts(Calculate):

    def __init__(self, output_dir: Path) -> None:
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())
        self.output_dir: Path = Path(output_dir)
        for structures in self.collection.protein_structure_results.values():
            for structure in structures.all_structures:
                if not self.output_dir.joinpath(structure.id).exists():
                    os.mkdir(self.output_dir.joinpath(structure.id))

    def get_contacts(self, structure_file: StructureFile):
        atomic_rep = Atomic(structure_file.path)
        bond_analysis = BondProteinAnalysis(atomic_rep)
        bond_analysis.compute()
        tsv_file = structure_file.path.parent.joinpath("bond_contacts_" + structure_file.id + ".tsv")

        contacts_df = pandas.read_csv(tsv_file, sep="\t", header=None, skiprows=2)

        # dssp = DSSP(Atomic.pdb, structure_file.path ,dssp='mkdssp')
        # dssp
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
                str(structure_file.path.parent) + "/" + str(structure_file.path.name).split(".")[0] + ".dssp",
                header=None, skiprows=28)
            hits = []
            for entry in df[0]:
                if len(entry.split()) < 10:
                    continue
                split = entry.split()
                hits.append([split[0:5]])
        except FileNotFoundError as ex:
            raise ex
        return pd.DataFrame(hits)

    def run(self) -> None:
        for sturctures in self.collection.protein_structure_results.values():
            for structure in sturctures.all_structures:
                print(structure.path)
                contacts = self.get_contacts(structure)
                secondary_structure = self.get_secondary_structure(structure)
                hits_of_interest = None
                for entry in secondary_structure.values:
                    if entry[0][4] in ['S', 'T', '.', 'G', '<<']:
                        break
                    else:
                        hits_of_interest = int(entry[0][0])
                contacts_of_interest = []
                for entry in contacts[0].values:
                    contact_1 = entry[2].split(':')
                    contact_2 = entry[3].split(':')
                    if int(contact_1[2]) <= hits_of_interest:
                        contacts_of_interest.append(entry)
                    elif int(contact_2[2]) <= hits_of_interest:
                        contacts_of_interest.append(entry)
                contacts_of_interest = pd.DataFrame(contacts_of_interest)
                structure.representation = Atomic
                res_to_keep = []
                count = 0
                for residue in structure.representation.pdb.get_residues():
                    count += 1
                    for index in range(len(contacts_of_interest)):
                        if residue.id[1] == int(contacts_of_interest.loc[index, 2].split(':')[2]) and \
                                contacts_of_interest.loc[index, 1] != "vdw":
                            res_to_keep.append(residue.id[1])
                            break
                    for index in range(len(contacts_of_interest)):
                        if residue.id[1] == int(contacts_of_interest.loc[index, 3].split(':')[2]) and \
                                contacts_of_interest.loc[index, 1] != "vdw":
                            res_to_keep.append(residue.id[1])
                            break

                res_to_keep = set(res_to_keep)
                total_residues = [i for i in range(1, count + 1)]
                total_residues = set(total_residues)
                res_to_remove = total_residues.difference(res_to_keep)
                #              atom_vectors[index] = [atom.get_vector() for atom in residue.get_atoms()]

                for chain in structure.representation.pdb[0]:
                    [chain.detach_child((' ', id, ' ')) for id in res_to_remove]
                io = PDBIO()
                io.set_structure(structure.representation.pdb)
                io.save(str(self.output_dir.joinpath(structure.id).joinpath(structure.id + "-contacts-only.pdb")),
                        preserve_atom_numbering=True)
                # io.save(str(structure.path.with_name(structure.id+'-contacts-only.pdb')), preserve_atom_numbering=True)
                atom_vectors = {}
            #  for residue in structure.representation.pdb.get_residues():
            #      for index in range(len(contacts_of_interest)):
            #          if residue.id[1] == int(contacts_of_interest.loc[index, 2].split(':')[2]):
            #              atom_vectors[index] = [atom.get_vector() for atom in residue.get_atoms()]


class AlignDataset(Calculate):

    def __init__(self, working_dir ,structure,output_dir: Path) -> None:
        super().__init__(working_dir,structure)
        self.collection = Collection(self.working_directory, HomologyStructureFetcher(), UniProtAcessionFetcher())
        self.output_dir = output_dir
        if not Path(self.output_dir).exists():
            os.mkdir(self.output_dir)
        for structures in self.collection.protein_structure_results.values():
            for structure in structures.all_structures:
                if not self.output_dir.joinpath(structure.id).exists():
                    os.mkdir(self.output_dir.joinpath(structure.id))

    def run(self) -> None:
        reference_id, reference_structures = self.collection.protein_structure_results.popitem()
        cmd.load(str(reference_structures.homology_structures[-1].path), reference_structures.id)
        cmd.save(self.output_dir.joinpath(reference_structures.id).joinpath(
            str(reference_structures.homology_structures[-1].path.name.split(".pdb")[0] + ".pdb")),
            reference_structures.id)
        os.system("cp " + str(reference_structures.uniprotID.path) + " " + str(
            self.output_dir.joinpath(reference_structures.id).joinpath(reference_structures.uniprotID.path.name)))

        count = 0
        for structures in self.collection.protein_structure_results.values():
            for structure in structures.all_structures:
                print(structure.path)

                cmd.load(str(structure.path), structure.id)
                cmd.align(structure.id, reference_structures.id)
                cmd.save(self.output_dir.joinpath(structure.id).joinpath(
                    str(structure.path.name.split(".pdb")[0] + ".pdb")), structure.id)
                os.system("cp " + str(structures.uniprotID.path) + " " + str(
                    self.output_dir.joinpath(structure.id).joinpath(structures.uniprotID.path.name)))

        # for structure in self.collection.protein_structure_results.values():


class ReassignBondsViaAlignment(Calculate):

    def __init__(self, output_directory: Path) -> None:
        super().__init__()
        self.collection = Collection(self.working_directory, HomologyStructureFetcher())
        self.output_directory = Path(output_directory)

    def get_contacts(self, structure_file: StructureFile, compute):
        if compute:
            atomic_rep = Atomic(structure_file.path)
            bond_analysis = BondProteinAnalysis(atomic_rep)
            bond_analysis.compute()
        tsv_file = structure_file.path.parent.joinpath("bond_contacts_" + structure_file.id + ".tsv")

        contacts_df = pandas.read_csv(tsv_file, sep="\t", header=None, skiprows=2)

        # dssp = DSSP(Atomic.pdb, structure_file.path ,dssp='mkdssp')
        # dssp
        return contacts_df,

    def get_secondary_structure(self, structure_file: StructureFile, compute=False) -> pd.DataFrame:
        if compute:
            args = [
                "mkdssp",
                str(structure_file.path.name),
                str(structure_file.path.name).split(".")[0] + ".dssp"
            ]
            p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(structure_file.path.parent))
            stdout, stderr = p2.communicate()
        try:
            df = pd.read_csv(
                str(structure_file.path.parent) + "/" + str(structure_file.path.name).split(".")[0] + ".dssp",
                header=None, skiprows=28)
            hits = []
            for entry in df[0]:
                if len(entry.split()) < 10:
                    continue
                split = entry.split()
                hits.append([split[0:5]])
        except FileNotFoundError as ex:
            raise ex
        return pd.DataFrame(hits)

    def find_coords(self, structure, contacts_of_interest):
        structure_file_lines = open(structure.path, "r").readlines()
        features = []

        for entry in contacts_of_interest.values:
            if entry[1] == 'vdw':
                continue
            coords_for_a = None
            coords_for_b = None
            for line in structure_file_lines:
                if line.startswith("ATOM"):
                    line_split = line.split()
                    if int(line_split[5]) == int(entry[2].split(':')[2]) and line_split[2] == entry[2].split(':')[-1]:
                        coords_for_a = np.array(line_split[6:9]).astype(np.float64)
                    elif int(line_split[5]) == int(entry[3].split(':')[2]) and line_split[2] == entry[3].split(':')[-1]:
                        coords_for_b = np.array(line_split[6:9]).astype(np.float64)
            features.append([entry[1], entry[2].split(':')[1], coords_for_a[0], coords_for_a[1], coords_for_a[2],
                             entry[3].split(':')[1],
                             coords_for_b[0], coords_for_b[1], coords_for_b[2]])

        return features

    def run(self) -> None:
        structs = []
        for structures in self.collection.protein_structure_results.values():
            for structure in structures.all_structures:
                contacts = self.get_contacts(structure, False)
                secondary_structure = self.get_secondary_structure(structure, False)
                hits_of_interest = None
                for entry in secondary_structure.values:
                    if entry[0][4] in ['S', 'T', '.', 'G', '<<']:
                        break
                    else:
                        hits_of_interest = int(entry[0][0])
                contacts_of_interest = []
                for entry in contacts[0].values:
                    contact_1 = entry[2].split(':')
                    contact_2 = entry[3].split(':')
                    if int(contact_1[2]) <= hits_of_interest:
                        contacts_of_interest.append(entry)
                    elif int(contact_2[2]) <= hits_of_interest:
                        contacts_of_interest.append(entry)
                contacts_of_interest = pd.DataFrame(contacts_of_interest)
                features = self.find_coords(structure, contacts_of_interest)
                features_numpy = np.array(features)
                np.save(structure.path.parent.joinpath(structure.id + "-bond-features"), features_numpy)


class GetClosestVerticiesToMesh(Calculate):

    def __init__(self,working_dir,structure, output_directory_path):
        super().__init__(working_dir,structure)
        self.output_dir = output_directory_path
        self.collection = Collection(self.working_directory, HomologyStructureFetcher(), UniProtAcessionFetcher())

    def run(self) -> None:
        for structures in self.collection.protein_structure_results.values():
            for protein_structure in structures.all_structures:
                #load structure
                protein_structure_coords = PDBParser(QUIET=True).get_structure(protein_structure.path.name, protein_structure.path)
                residues_coords = []
                vertices_of_interest = []
                for residue in structures.uniprotID.binding_site_residues:
                    protein_residue = protein_structure_coords[0]['A'][(' ',residue,' ')]
                    residue_coords = []
                    for atom in protein_residue:
                        residue_coords.append(atom.coord)
                    residue_coords = np.array(residue_coords)
                    #residues_coords.append(residue_coords)

                # Step 1: get protein structure and find the closets set of cavities, if no cavities exist, try to compute them
                    if not Path("/home/felix/cavities/").joinpath(protein_structure.path.parent.parent.parent.name+"_cavity").exists():
                        os.mkdir(Path("/home/felix/cavities/").joinpath(protein_structure.path.parent.parent.parent.name+"_cavity"))
                    logging.info("Generating cavities")
                    path=Path("/home/felix/cavities/").joinpath(protein_structure.path.parent.parent.parent.name+"_cavity").joinpath(protein_structure.id)
                    if not path.exists():
                        os.mkdir(path)
                    path_to_save = path.joinpath(protein_structure.id)
                    residue_info = protein_residue.resname+str(protein_residue.id[1])+"A-"
                    Popen(args=["/home/felix/Software/bin/Get_Cleft", "-p" ,protein_structure.path, "-a" ,residue_info, "-o", str(path_to_save)]).communicate(timeout=300)
                    meshes_collection = Collection( Path("/home/felix/cavities/").joinpath(protein_structure.path.parent.parent.parent.name+"_cavity").joinpath(protein_structure.id), ExperimentalStructureFetcher())
                    meshes = meshes_collection.protein_structure_results.popitem()
                    min_distance = 100000
                    closest_mesh=None
                    for mesh in meshes[1].all_structures:
                        mesh_structure_coords = PDBParser(QUIET=True).get_structure(mesh.path.name, mesh.path)
                        mesh_atoms = [atom for atom in mesh_structure_coords.get_atoms()]
                        coords = np.array([atom.coord for atom in mesh_atoms])
                        distances = cdist(coords,residue_coords,  metric='euclidean')
                        #distances = np.linalg.norm(coords-residue_coords,axis=1)
                        mesh_min_distance = distances.mean() #distances[np.argmin(distances)]
                        if mesh_min_distance < min_distance:
                            min_distance = mesh_min_distance
                            closest_mesh = mesh_structure_coords
                    # Step 2 find closest vertices to mesh
                #mesh_vertices, mesh_faces, normals, _, _, _, _ = read_ply(str(protein_structure.path.parent.joinpath(
                #    protein_structure.path.name.split(".")[0] + ".ply")))
                    mesh = pymesh.load_mesh(str(protein_structure.path.parent.joinpath(protein_structure.path.name.split(".")[0] + ".ply")))
                    vertices = mesh.vertices
                    tree= cKDTree(vertices)
                    nearest_points=[]
                    for atom in closest_mesh.get_atoms():
                        target_coords =atom.coord
                        indices = tree.query_ball_point(target_coords, r=3.0)
                        nearby_points = indices
                        nearest_points.extend(nearby_points)
                    nearest_points =list(np.unique(np.array(nearest_points)))
                    vertices_of_interest.extend(nearest_points)
                if len(vertices_of_interest) < 1:
                    continue
                vertices_of_interest = np.unique(np.array(vertices_of_interest))
                input_feat, rho, theta, mask, neigh_indices, iface_labels, verticies = read_data_from_surface(
                    str(protein_structure.path.parent.joinpath(
                        protein_structure.path.name.split(".")[0] + ".ply")), masif_opts["site"], vertices_of_interest)
                mesh_vertices, mesh_faces, normals, _, _, _, _ = read_ply(str(protein_structure.path.parent.joinpath(
                    protein_structure.path.name.split(".")[0] + ".ply")))
                CalculateSurfaceOfInterest.debug(mesh_faces,np.array(neigh_indices),vertices_of_interest,verticies,normals,theta,rho,name=protein_structure.path.parent.name)
                numpy.save(protein_structure.path.parent.joinpath(
                    protein_structure.path.name.split(".")[0] + "-nearest-cavities-vertex-pos"),
                    input_feat)

                #protein_structure.path.parent.joinpath(protein_structure.path.name.split(".")[0] + ".ply")




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
        self.clustering_algorithm(self.data, self.ordering_file)
