"""
Get representative structure and binding site data if available. If binding site data is unavailable,
use ML-based tool to indentify some potential sites. Afterward, align representative to the members,
check if there is a suitable alignment to the binding site from the representative, if so, do site comparison.
Report back the cluster's performance.
"""
import json
import logging
import csv
import os
import subprocess
from pathlib import Path
from typing import Dict, Optional, List

import numpy as np
import scipy.spatial.distance
from Bio import PDB
from Bio.PDB import PDBIO, Selection, Superimposer
from pymol import cmd

from query import AlphaFoldQuery

logging.basicConfig(filename="get_scores_run",
                    filemode='a',
                    format='%(asctime)s,%(msecs)03d %(name)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.DEBUG)

logging.info("Running script")

logger = logging.getLogger('binding')
#WKDIR_PATH = "/mnt/research/cytotoxic/"
WKDIR_PATH = "/mnt/research/uniref50_ec_mismatches_for_patches/"
CSV_FILE_NAME = "_cluster_scores.csv"


def flow():
    cluster_data = get_clusters()
    for representative, members in cluster_data.items():
        iterate_through_cluster(representative, members)


def identify_member_binding_site(representative_binding_site_structure_data, aligned_member_structure_data,
                                 aligned_member_path):
    # Parse the structure from a PDB file

    # Select the C-alpha atoms from both structures
    selection = Selection.unfold_entities(representative_binding_site_structure_data,
                                          'R')  # Use only ATOM records (discard HETATM and hetatm)
    ca1 = [res for res in selection]
    if len(ca1) < 1:
        return None
    selection = Selection.unfold_entities(aligned_member_structure_data, 'R')
    ca2 = [res for res in selection]

    #ca2 = [atom for atom in selection if atom.get_name() == "CA"]

    # Calculate RMSD using Superimposer

    #  coord1 = np.array([a.get_coord() for a in ca1])  # Coordinates of C-alpha atoms from structure1
    #  coord2 = np.array([a.get_coord() for a in ca2])  # Coordinates of C-alpha atoms from structure2
    best_matches = []
    for reference_amino_acid in ca1:
        max_dist = 10000
        best_match_residue = None
        for res_num, query_amino_acid in enumerate(ca2):
            reference_amino_acid_pos = np.array([atom.get_vector().get_array() for atom in reference_amino_acid.get_atoms() if atom.id in ["CA","CB","N","O"]])
            query_amino_acid_pos = np.array([atom.get_vector().get_array() for atom in query_amino_acid.get_atoms() if atom.id in ["CA","CB","N","O"]])
            try:
                dist = np.linalg.norm(query_amino_acid_pos-reference_amino_acid_pos,ord='fro')
            except ValueError as ex:
                logger.warning(f"Could compare {ex}")
                continue
            #reference_amino_acid_pos =reference_amino_acid.get_vector().get_array()
            #query_amino_acid_pos =query_amino_acid.get_vector().get_array()
            #dist = np.linalg.norm(query_amino_acid_pos-reference_amino_acid_pos)
            if dist < max_dist:
                max_dist = dist
                best_match_residue = res_num+1
        if best_match_residue not in best_matches and best_match_residue is not None:
            best_matches.append(best_match_residue)


    #best_rmsd = 10000
    #best_ca_indx = None
    #for iteration in range(0, len(ca2), len(ca1)):
    #    sup = Superimposer()
    #    start = iteration
    #    ca_of_interest = ca2[start:start + len(ca1)]
    #    if len(ca_of_interest) < len(ca1):
    #        extra_res = len(ca1) - len(ca_of_interest)
    #        start -= extra_res
    #    ca_of_interest = ca2[start:iteration + len(ca1)]
    #    sup.set_atoms(ca1, ca_of_interest)  # Superimpose on the basis of C-alphas
    #    if sup.rms < best_rmsd:
    #        best_rmsd = sup.rms
    #        best_ca_indx = range(start, iteration + len(ca1))
    residues_to_remove = []
    io = PDBIO()
    for residue in aligned_member_structure_data.get_residues():
        if int(residue.id[1]) not in best_matches:
            residues_to_remove.append(int(residue.id[1]))
    try:
        for chain in aligned_member_structure_data[0]:
            [chain.detach_child((' ', id, ' ')) for id in residues_to_remove]
    except KeyError as ex:
        for chain in aligned_member_structure_data:
            [chain.detach_child((' ', id, ' ')) for id in residues_to_remove]
    io.set_structure(aligned_member_structure_data)
    io.save(str(aligned_member_path.joinpath(aligned_member_path.name + "_binding_site_data.pdb")),
            preserve_atom_numbering=False)
    return aligned_member_path.joinpath(aligned_member_path.name + "_binding_site_data.pdb")

    # for residue in representative_binding_site_structure_data.get_residues():
    #      atom_coodinates = []
    #      for atom in residue.get_atoms():
    #          atom_coords = atom.get_vector().get_array()
    #          atom_coodinates.append(atom_coords)
    #      atom_coodinates = np.array(atom_coodinates)
    #      for aligned_residue in aligned_member.get_residues():
    #          aligned_coords = []
    #          for aligned_atom in aligned_residue.get_atoms():
    #              atom_coords = aligned_atom.get_vector().get_array()
    #              aligned_coords.append(atom_coords)
    #          aligned_coords = np.array(aligned_coords)
    #      pairwise_distance_matrix = scipy.spatial.distance.cdist(atom_coodinates, aligned_coords)


def align_protein(representative_structure, query_protein_structure):
    cmd.load(str(representative_structure), representative_structure.parent.name)
    cmd.load(str(query_protein_structure), query_protein_structure.parent.name)
    cmd.align(query_protein_structure.parent.name,representative_structure.parent.name)
    cmd.save(query_protein_structure.parent.joinpath(query_protein_structure.name.replace(".pdb", "")
                                                     + "-" + representative_structure.name),query_protein_structure.parent.name)
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(query_protein_structure.name + "_aligned",
                                         query_protein_structure.parent.joinpath(
                                             query_protein_structure.name.replace(".pdb", "")
                                             + "-" + representative_structure.name))
    structure_path = query_protein_structure.parent.joinpath(
        query_protein_structure.name.replace(".pdb", "")
        + "-" + representative_structure.name)
    cmd.remove(representative_structure.parent.name)
    cmd.remove(query_protein_structure.parent.name)
    #cmd.remove(query_protein_structure.parent.joinpath(query_protein_structure.name.replace(".pdb", "")
    #                                                 + "-" + representative_structure.name))
    return structure, structure_path


def compare_binding_sites(representative_binding_site, query_binding_site):
    """

    Parameters
    ----------
    representative_binding_site
    query_binding_site

    Returns
    -------

    """
    with open(query_binding_site, 'r') as file1:
        lines= file1.readlines()
    lines=lines[:-2]
    with open(query_binding_site, 'w') as file1:
        for line in lines:
            file1.write(line)
    #with open(query_binding_site, 'r') as file1:
    #    file1.write("TER\n")
    p1 = subprocess.Popen(args=['java', 'AssignChemicalFeatures', representative_binding_site])
    p1.wait()
    p2 = subprocess.Popen(args=['java', 'AssignChemicalFeatures', query_binding_site])
    p2.wait()
    p3 = subprocess.Popen(args=['glosa', '-s1', representative_binding_site, '-s1cf',
                                representative_binding_site.parent.joinpath(
                                    representative_binding_site.name.split('.pdb')[0] + "-cf.pdb"), '-s2',
                                query_binding_site, '-s2cf',
                                query_binding_site.parent.joinpath(
                                    query_binding_site.name.split('.pdb')[0] + "-cf.pdb"),
                                "-itercf", "2",
                                "-o", "0"])
    p3.wait()
    with open(representative_binding_site.parent.joinpath(representative_binding_site.name+"-"+query_binding_site.name+"ga-score.txt"), "r") as file1:
        data =file1.read().splitlines()
        target=data[0].split(" ")[-1]
    return float(target)


def iterate_through_cluster(representative, members):
    try:
        representative_structure_path, representative_binding_site_structure_data, representative_binding_site_path = get_representative_data(
        representative)
    except RuntimeError as ex:
        return
    member_score_list = {}
    with open(representative_binding_site_path, 'r') as file1:
        lines = file1.readlines()
    lines = lines[:-2]
    with open(representative_binding_site_path, 'w') as file1:
        for line in lines:
            file1.write(line)
    for member in members:
        if member == representative:
            continue
        try:
            protein_data_path, binding_site_data_query = get_protein_structure(member)
        except RuntimeError as ex:
            continue
        aligned_member_structure_data, aligned_member_path = align_protein(representative_structure_path,
                                                                           protein_data_path)
        member_structureal_binding_site = identify_member_binding_site(representative_binding_site_structure_data,
                                                                       aligned_member_structure_data,
                                                                       member)
        if member_structureal_binding_site is not None:

            score = compare_binding_sites(representative_binding_site_path, member_structureal_binding_site)
        else:
            score = 0
            logger.info(f"score is 0 for {member}")
        member_score_list[member] = score
    scores_path = WKDIR_PATH + representative.name + CSV_FILE_NAME
    with open(scores_path, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in member_score_list.items():
            writer.writerow([key, value])
    logger.info(f"CSV file '{scores_path}' created successfully.")


def get_clusters() -> Dict:
    #dirs = os.listdir(Path(WKDIR_PATH).joinpath("cluster_data"))
    dirs = os.listdir(Path(WKDIR_PATH).joinpath("clusters"))
    clusters = {}
    for cluster in dirs:
        clusters[Path(WKDIR_PATH).joinpath("clusters").joinpath(cluster).joinpath(cluster.split("_")[-1])] = []
        for root, dirs, files in os.walk(Path(WKDIR_PATH).joinpath("clusters").joinpath(cluster)):
            for file in files:
                if file.endswith(".json"):
                    clusters[Path(WKDIR_PATH).joinpath("clusters").joinpath(cluster).joinpath(
                        cluster.split("_")[-1])].append(Path(root).joinpath(file).parent)
    return clusters


def get_representative_data(representative: Path):
    """
    This method is to get the representative data from uniprot, identify the binding site from metadata,
    capture protein structure and keep note of where the binding site infromation for the structure is.
    Returns
    -------

    """
    protein_stucture_path, protein_structure = get_protein_structure(representative)
    binding_site_structure_data, binding_site_path = get_binding_site(representative, protein_structure)
    return protein_stucture_path, binding_site_structure_data, binding_site_path


def get_protein_structure(representative):
    uniprot_id = representative.name
    potential_structure=representative.joinpath("AF-" + uniprot_id + "-F1-model_v4.pdb")
#    folder_of_interest = representative.parent.parent.parent.joinpath("tmp").joinpath(uniprot_id)
#    potential_structure = folder_of_interest.joinpath("AF-" + uniprot_id + "-F1-model_v4.pdb")
    if potential_structure.exists():
        return potential_structure, PDB.PDBParser(QUIET=True).get_structure(uniprot_id, potential_structure)
    else:
        raise RuntimeError(f"{uniprot_id}")
    result = _get_structure_from_af_website(representative)
    if result is None:
        result = _predict_protein_structure(uniprot_id)
    return None, result


def _get_structure_from_af_website(representative):
    """
    Grabs protein structure from AF's DB website: https://alphafold.ebi.ac.uk/
    Parameters
    ----------
    uniprot_id:

    Returns
    -------

    """
    pass
    # AlphaFoldQuery(representative).query(representative.name)


def _predict_protein_structure(uniprot_id):
    """
    Predicts protein structure using ColabFold.
    Parameters
    ----------
    uniprot_id

    Returns
    -------

    """
    pass


def get_binding_site(representative, protein_structure):
    binding_site_residue_data = get_binding_site_from_uniprot(representative)
    binding_site_path = None
    if len(binding_site_residue_data) < 1:
        binding_site_residue_data = predict_binding_site(protein_structure)
    binding_site_data, binding_site_path = fetch_binding_site_structure_data(binding_site_residue_data,
                                                                             protein_structure, representative)

    return binding_site_data, binding_site_path


def fetch_binding_site_structure_data(binding_site_residue_data, protein_structure, representative):
    residues_to_remove = []
    io = PDBIO()
    for residue in protein_structure.get_residues():
        if int(residue.id[1]) not in binding_site_residue_data:
            residues_to_remove.append(int(residue.id[1]))
    try:
        for chain in protein_structure[0]:
            [chain.detach_child((' ', id, ' ')) for id in residues_to_remove]
    except KeyError as ex:
        for chain in protein_structure:
            [chain.detach_child((' ', id, ' ')) for id in residues_to_remove]
    io.set_structure(protein_structure)
    io.save(str(representative.joinpath(representative.name + "_binding_site_data.pdb")), preserve_atom_numbering=False)
    return protein_structure, representative.joinpath(representative.name + "_binding_site_data.pdb")


def get_binding_site_from_uniprot(representative):
    uniprot_data = representative.joinpath(representative.name + ".json")
    if uniprot_data.exists():
        binding_site_residues_from_uniprot = read_uniprot_entry_data(uniprot_data)
        return binding_site_residues_from_uniprot


def read_uniprot_entry_data(uniprot_data):
    uniprot_data = json.load(open(uniprot_data, 'r'))
    acceptable_sites = ["ACT_SITE", "BINDING", "Binding site", "Active site", "Site"]
    if uniprot_data["entryType"] == "UniProtKB unreviewed (TrEMBL)":
        return []
    try:
        features: List = uniprot_data["features"]
    except KeyError:
        return []
    binding_site_residues = []
    for feature in features:
        if feature["type"] in acceptable_sites:
            try:
                binding_site_residues.extend(range(int(feature["begin"]), int(feature["end"]) + 1))
            except KeyError as ex:
                logger.warning("Getting features failed with first method. Attempting second method.")
                binding_site_residues.extend(range(int(feature["location"]["start"]["value"]),
                                                   int(feature["location"]["end"]["value"]) + 1))

    return binding_site_residues


def predict_binding_site(protein_structure_path):
    """

    Parameters
    ----------
    protein_structure_path

    Returns
    -------

    """
    return []


if __name__ == '__main__':
    flow()
