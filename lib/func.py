import json
import logging
import os
from pathlib import Path
from typing import Dict, List
import transformers

from Bio import SeqIO


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
            write_file.write(record.id+"\n")
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
        os.mkdir(path.joinpath(accession))
    x = os.scandir(path)#"/media/felix/Research/KpsData/KpsT/substructures")
    for entry in x:
        p = path.joinpath(''.join(entry.name.split("_")[0:-1]))
        os.system("mv %s %s"% (entry.path, str(p.joinpath(entry.name))))

