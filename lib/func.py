import json
import logging
import os
import re
from pathlib import Path
from typing import Dict, List

import Bio.SeqIO
import pandas
from sklearn.mixture import GaussianMixture
import numpy as np
import transformers
from transformers import EsmTokenizer, EsmModel
from Bio import SeqIO
from sklearn.decomposition import PCA

from dataStructure.collections import Collection, HomologyStructureFetcher
from lib.const import grantham_distance_matrix_row_dict, grantham_distance_matrix
from scipy.cluster.hierarchy import dendrogram, linkage
from lib.const import e_coli_k_type

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
def change_directory(directory: Path):
    """

    :param directory:
    :param skip:
    :return:
    """
    if not directory.exists():
        os.mkdir(directory)
    logging.debug("Changing directory to path %s" % directory)
    os.chdir(directory)
    return True


def load_json(filename: Path) -> Dict:
    """

    Parameters
    ----------
    filename

    Returns
    -------

    """
    try:
        with open(filename, "r") as f:
            data = json.load(f)
    except FileNotFoundError as ex:
        raise ex
    return data


def split_fasta_file(file_name: str, output_dir):
    for record in SeqIO.parse(file_name, "fasta"):
        if not os.path.exists(Path(output_dir).joinpath(record.name)):
            os.mkdir(Path(output_dir).joinpath(record.name))
        with open(Path(output_dir).joinpath(record.name).joinpath(record.name + ".fasta"), "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")


def get_files(file_path: Path, regular_expression: str) -> List[Path]:
    return list(file_path.rglob(regular_expression))


#   file_path.glob('**/*')
#   files: List[Path] = []
#   for file in os.listdir(file_path):

#    return [Path(str(file)) for file in os.listdir(file_path)]

def create_accession_list_for_genbank_IDs(fasta_file: Path, output_file_path: Path):
    """
    Creates accession_list  from fasta file genbank entries at output file path

    Parameters
    ----------
    output_file_path: Path where to output accession_list
    fasta_file: fasta file containing entries of interest

    -------

    """
    write_file = open(output_file_path, "w")
    with open(fasta_file, "r") as read_handle:
        for record in SeqIO.parse(read_handle, "fasta"):
            logging.debug("Writing %s" % record.id)
            write_file.write(record.id + "\n")
    write_file.close()


def parse_fasta_to_gp():
    with open("/home/felix/kspT_200_230_parsed.fasta", "w") as output_handle:
        for record in SeqIO.parse("/home/felix/kspT_200_230_parsed.gp", "gb"):
            SeqIO.write(record, output_handle, format="fasta")


def build_db(path: Path, accession_db):
    accession = open(accession_db, "r")
    for line in accession:
        accession = line.split("\n")[0]
        accession = accession.split(".1")[0]
        accession = ''.join(accession.split('_'))
        try:
            os.mkdir(path.joinpath(accession))
        except FileExistsError as ex:
            print(f"Caution: {ex}")
    x = os.scandir(path)  # "/media/felix/Research/KpsData/KpsT/substructures")
    for entry in x:
        p = path.joinpath(''.join(entry.name.split(".")[0].replace("_", "")))
        print(p)
        os.system("mv %s %s" % (entry.path, str(p.joinpath(entry.name))))


def find_rigid_alignment(A, B):
    """
    See: https://en.wikipedia.org/wiki/Kabsch_algorithm
    2-D or 3-D registration with known correspondences.
    Registration occurs in the zero centered coordinate system, and then
    must be transported back.
        Args:
        -    A: Numpy array of shape (N,D) -- Point Cloud to Align (source)
        -    B: Numpy array of shape (N,D) -- Reference Point Cloud (target)
        Returns:
        -    R: optimal rotation
        -    t: optimal translation
    """
    a_mean = A.mean(axis=0)
    b_mean = B.mean(axis=0)
    A_c = A - a_mean
    B_c = B - b_mean
    # Covariance matrix
    H = A_c.T.dot(B_c)
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    # Rotation matrix
    R = V.dot(U.T)
    # Translation vector
    t = b_mean - R.dot(a_mean)
    return R, t


def calculate_rmsd(A, B):
    return np.sqrt((((A - B) ** 2).sum(axis=1)).mean())


# def calculate_rmsd(coord1: np.ndarray, coord2: np.ndarray):
#    """#

#    Parameters
#    ----------
#    coord1
#    coord2#

#    Returns
#    -------

#    """
#    return np.sqrt(((coord1 - coord2) ** 2).mean())

def calculate_grantham_distance(sequence1, sequence2, coords1, coords2):
    """

    Parameters
    ----------
    sequence1
    sequence2

    Returns
    -------

    """
    grantham_similarity = np.zeros(shape=(len(sequence1), len(sequence2)))
    for index2 in range(len(sequence1)):
        for index in range(len(sequence2)):
            row = grantham_distance_matrix_row_dict[type(sequence1[index2])]
            column = grantham_distance_matrix_row_dict[type(sequence2[index])]
            physiochemical_value = grantham_distance_matrix[row][column]
            grantham_similarity[index2][index] = physiochemical_value
    return np.sum(np.sum(1 / (np.cosh(0.3 * np.linalg.norm(np.array([coords1[amino_acid1], coords2[amino_acid2]])))) * (
            1 - grantham_similarity[amino_acid1][amino_acid2] / 215) for amino_acid2 in range(len(sequence2))) for
                  amino_acid1 in range(len(sequence1)))


def granthem_distance(sequence1, sequence2, ):
    grantham_similarity = np.zeros(shape=(len(sequence1), len(sequence2)))
    for index2 in range(len(sequence1)):
        for index in range(len(sequence2)):
            row = grantham_distance_matrix_row_dict[type(sequence1[index2])]
            column = grantham_distance_matrix_row_dict[type(sequence2[index])]
            physiochemical_value = grantham_distance_matrix[row][column]
            grantham_similarity[index2][index] = physiochemical_value
    return np.sum(grantham_similarity)
def granthem_distance_fixed(sequence1, sequence2, ):
    grantham_similarity = np.zeros(shape=(len(sequence1)))
    for index2 in range(len(sequence1)):
        row = grantham_distance_matrix_row_dict[type(sequence1[index2])]
        column = grantham_distance_matrix_row_dict[type(sequence2[index2])]
        physiochemical_value = grantham_distance_matrix[row][column]
        grantham_similarity[index2] = physiochemical_value
    return np.sum(grantham_similarity)

def grantham_distance_backbone(sequence1, sequence2, coords1, coords2):
    """

    Parameters
    ----------
    sequence1
    sequence2

    Returns
    -------

    """
    grantham_similarity = np.zeros(shape=(len(sequence1), len(sequence2)))
    for index2 in range(len(sequence1)):
        for index in range(len(sequence2)):
            row = grantham_distance_matrix_row_dict[type(sequence1[index2])]
            column = grantham_distance_matrix_row_dict[type(sequence2[index])]
            physiochemical_value = grantham_distance_matrix[row][column]
            grantham_similarity[index2][index] = physiochemical_value

    return np.sum(np.sum(1 / (np.cosh(0.3 * np.linalg.norm(np.array([coords1[amino_acid1], coords2[amino_acid2]])))) * (
            1 - grantham_similarity[int(amino_acid1 / 4)][int(amino_acid2 / 4)] / 215) for amino_acid2 in
                         range(len(coords2))) for amino_acid1 in range(len(coords1)))


def prep_seq(seq: str):
    """
    Credit:

    Adding spaces between AAs and replace rare AA [UZOB] to X.
    ref: https://huggingface.co/Rostlab/prot_bert.

    :param seq: a string of AA sequence.

    Returns:
        String representing the input sequence where U,Z,O and B has been replaced by X.
    """
    seq_spaced = " ".join(seq)
    seq_input = re.sub(r"[UZOB]", "X", seq_spaced)
    return seq_input


def get_encoding(sequences: [str], tokenizer, model):
    """

    Returns
    -------

    """
    record_dict = {}
    sequences = [prep_seq(sequence) for sequence in sequences]
    encodings = tokenizer(sequences, return_tensors="pt", padding=True, truncation=True).to("cuda")
    for index, seq in enumerate(sequences):
        record_dict[index] = {
            "input_ids": encodings["input_ids"][index],
            "attention_mask": encodings["attention_mask"][index]
        }
    for index, record in enumerate(record_dict.values()):
        pooling_output = model(record["input_ids"].reshape(1, -1),
                               record["attention_mask"].reshape(1, -1)).pooler_output
        record_dict[index]["embedding"] = pooling_output
    return record_dict


def calc_pca(matrix: np.ndarray, n_components: int = None):
    pca = PCA(n_components="mle")
    return pca.fit_transform(matrix)


def unsupervised_gmm_clustering_and_majority_voting(matrix):
    labels_matrix = np.empty(shape=matrix.shape)
    for index, feature_in_matrix in enumerate(matrix):
        stop_criteron = False
        count = 1
        aic_lowest_score = 100000
        bic_lowest_score = 100000
        # while not stop_criteron:
        gmm = GaussianMixture(4)
        gmm.fit(feature_in_matrix.reshape(-1, 1))
        #  aic = gmm.aic(feature_in_matrix.reshape(-1,1))
        #   bic = gmm.bic(feature_in_matrix.reshape(-1,1))
        #   if aic < aic_lowest_score and bic < bic_lowest_score:
        #        aic_lowest_score = aic
        #       bic_lowest_score = bic
        #       count+=1
        #    else:
        #        break
        labels_matrix[index] = gmm.predict(feature_in_matrix.reshape(-1, 1))
    return labels_matrix


def cluster_gmm(matrix, index):
    # matrix = calc_pca(matrix)
    gmm = GaussianMixture(index)
    gmm.fit(matrix)
    return gmm.predict(matrix)


def plot_gmm(matrix, labels):
    plt.scatter(range(len(matrix.transpose())), np.sum(matrix, axis=1).reshape(-1, 1), c=labels, s=40,
                cmap='viridis')


# def cluster_purity(labels):
#    for label in labels:

def strain_to_accession(mapping_json_file: Path, strains: dict) -> dict:
    json_object = json.load(open(mapping_json_file))
    accessions = {}
    for strain in strains.keys():
        for accession, strain_id in json_object.items():
            if strain_id == strain:
                accessions[accession] = (strain,strains[strain].value)
                # accessions.append(accession)
                break
    return accessions

#for index, item in enumerate(self.collection.protein_structure_results.values()):
#    for structure in item.all_structures:
#         structure_cluster_results[structure.id]=dendo["color_list"][index]


def accession_representative_fetcher(accession_ids, mapping_file):
    csv = pandas.read_csv(mapping_file,sep="\t",header=None)
    accessions_mapping = {}
    for accession in accession_ids.keys():
        representative = csv[csv[1] == accession][0]
        if len(representative) > 1:
            representative = representative.array[0]
        else:
            representative = representative.item()
        if accessions_mapping.get(representative) is None:
            accessions_mapping[representative] = [(accession, accession_ids[accession])]
        else:
            accessions_mapping[representative].append((accession, accession_ids[accession]))
    return accessions_mapping


def map_clustering(accessions_mapping, collection, dendogram):
    for index, structures in enumerate(collection.protein_structure_results.values()):
        for structure in structures.all_structures:
#        #print(structure.id.replace("WP","WP_") +".1")
            if accessions_mapping.get(structure.id.replace("WP","WP_") +".1") is not None:
                accessions_mapping[structure.id.replace("WP", "WP_") + ".1"].append(dendogram["color_list"][index])
    return accessions_mapping
def temp_fix_load(data_path,working_path):
    data = np.load(data_path)
    collection = Collection(working_path, HomologyStructureFetcher())
    labels= cluster_gmm(data,4)
    labels_fixed = [[],[],[],[]]
    for index, item in enumerate(collection.protein_structure_results.values()):
        for structure in item.all_structures:
            labels_fixed[labels[index]].append(structure.id)
    return labels_fixed

def build_dendrogram(data):
    return dendrogram(linkage(data,method="ward"))

def get_dendrogram(collection: Collection):
    accessions = strain_to_accession(Path("/home/felix/testing.json"),e_coli_k_type)
    accessions = accession_representative_fetcher(accessions, "/home/felix/kpsT_cluster_cluster.tsv", )
    data = np.load("/media/felix/ShortTerm/Research/KpsData/KpsT/pathotype_structures/geometry_only_alphafold_pathotype.npy")
    labels = []
    for structs in collection.protein_structure_results.values():
        for struct in structs.all_structures:
            labels.append(struct.id +" "+ e_coli_k_type[struct.id].value[0])
    dendo = dendrogram(linkage(data),labels=labels, color_threshold=6, above_threshold_color='grey', orientation="left")
    plt.axhline(y=6, c='grey', lw=1, linestyle='dashed')
    return map_clustering(accessions, collection,dendo)


def create_strain_fastas(accessions, working_dir: Path):
    for record in Bio.SeqIO.parse("/home/felix/Research/tmpKps/kpsT_all_seqs.fasta","fasta"):
        if record.id in accessions.keys():
            try:
                os.mkdir(working_dir.joinpath(record.id))
            except:
                continue
            Bio.SeqIO.write(record, working_dir.joinpath(record.id).joinpath(accessions[record.id][0]+".fasta"),"fasta")
            #working_dir.joinpath(record.id).joinpath(accessions[recor)

def create_pathotype_fasta(accessions, working_dir: Path):
    records_to_save = []
    for record in Bio.SeqIO.parse("/home/felix/Research/tmpKps/kpsT_all_seqs.fasta","fasta"):
        if record.id in accessions.keys():
            record.description = ""
            #record.description = accessions[record.id][0] + " " + accessions[record.id][1][0]
            record.id = accessions[record.id][0] + " " + accessions[record.id][1][0]
            records_to_save.append(record)

    Bio.SeqIO.write(records_to_save, working_dir.joinpath("kpsT_pathotype_sequences.fasta"),"fasta")