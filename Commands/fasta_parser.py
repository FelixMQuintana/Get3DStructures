import json
import os
from pathlib import Path
from typing import Optional
import ast
import Bio.Align
import pandas as pd

from Commands.command import Command
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
from modeller import *


def find_strain_from_accession_id(database_of_strains: Path):#, database_of_accession_ids: Path):
    """

    Parameters
    ----------
    database_of_strains
    database_of_accession_ids

    Returns
    -------

    """
    #accession_db = open(database_of_accession_ids, "r")
    mapping = {}
    #for accession_id in accession_db:
    for record in SeqIO.parse(database_of_strains, "gb"):
        record = RecordWrapper(record)
            #if accession_id.strip("\n") == record.record.id:
        mapping[str(record.record.id)] = record.strain
             #   continue
    json.dump(mapping, open("/home/felix/testing.json", "w"))


class FastaData(Command):

    def __init__(self, fasta_file: str, gene_name: str, lower_bound: int, upper_bound: int, out_file_path: str):
        super().__init__()
        self._fasta_file_path = fasta_file
        self._gene_name = gene_name
        self._upper_bound = upper_bound
        self._lower_bound = lower_bound
        self._outfile_name = Path(out_file_path)

    def run(self) -> None:
        records: [RecordWrapper] = []
        for record in SeqIO.parse(self.working_directory.joinpath(self._fasta_file_path), "gb"):
            wrapper = RecordWrapper(record)
            if self._gene_name not in wrapper.region_name:
                continue
            elif self._lower_bound > len(wrapper.protein_length) or len(wrapper.protein_length) > self._upper_bound:
                continue

            #   print(len(wrapper.protein_length))
            records.append(wrapper)

        for index_1, record in enumerate(records):
            for index_2, record_repeat in enumerate(records):
                if index_1 == index_2:
                    continue
                if record.strain == record_repeat.strain and record.strain != None:
                    x = Bio.Align.PairwiseAligner()
                    x.open_gap_score = -12
                    x.extend_gap_score = -4
                    matrix = substitution_matrices.load("BLOSUM62")
                    x.substitution_matrix = matrix
                    alignment_ref = x.align(record.record.seq, records[0].record.seq)
                    alignment_query = x.align(record_repeat.record.seq, records[0].record.seq)
                    # alignment = x.align(record.record.seq, record_repeat.record.seq)
                    if alignment_query.score > alignment_ref.score:
                        try:
                            records.remove(record)
                        except ValueError:
                            print("Already removed")
                    else:
                        records.remove(record_repeat)
        #
        #   if alignment[0].counts()[2] == 0:
        #       print("Match!")
        #    print(alignment[0])
        #    else:
        #        print("mismatch")
        #        print(alignment[0])
        #        x.align(record.record.seq, records[0].record.seq)
        print(len(records))

        with open(self._outfile_name, "w") as output_handle:
            for record in records:
                SeqIO.write(record.record, output_handle, format="gb")


class ParseMapping(Command):

    def __init__(self, json_file: str, out_file_path: str):
        super().__init__()
        self._json_file: Path = Path(json_file)
        self._out_file_path: Path = Path(out_file_path)

    def run(self) -> None:
        json_object = json.load(open(self._json_file, "r"))["results"]
        mapping_dict = {}
        for entry in json_object:
            mapping_dict[entry["from"]] = entry["to"]["primaryAccession"]
        with open(self._out_file_path, "w") as out:
            for value in mapping_dict.values():
                out.write(value + "\n")
        with open(
                self._out_file_path.parent.joinpath("mapping_file_" + self._out_file_path.name.split(".")[0] + ".txt"),
                "w") as out:
            out.write(str(mapping_dict))


class MappPDBToGenBankHits(Command):

    def __init__(self, json_file: str, cluster_representative_file: str):
        super().__init__()

        self._mapping_dict = json_file
        self._cluster_file = Path(cluster_representative_file)

    def run(self) -> None:
        with open(self._mapping_dict) as f:
            data = f.read()
        self._mapping_dict = ast.literal_eval(data)
        cluster_file = pd.read_csv(self._cluster_file, header=None, sep="\t")
        cluster_represntatives = cluster_file[0].unique()
        mapping_file_found_representatives = []
        for key in self._mapping_dict:
            cluster_rep_hits = cluster_file.index[cluster_file[0] == key].tolist()
            if len(cluster_rep_hits) == 0:
                cluster_member = cluster_file.index[cluster_file[1] == key]
                associated_represent = cluster_file[0][cluster_member].unique().tolist()
                mapping_file_found_representatives.append((associated_represent[0], self._mapping_dict[key]))
            else:
                cluster_rep_list = (
                    cluster_file[0][cluster_file.index[cluster_file[0] == key]].unique().tolist()[0],
                    self._mapping_dict[key])

                mapping_file_found_representatives.append(cluster_rep_list)
        for key, value in mapping_file_found_representatives:
            if not os.path.exists("/media/felix/Research/KpsData/KpsS/" + key):
                os.mkdir("/media/felix/Research/KpsData/KpsS/" + key)
                os.system("cp -r /media/felix/Research/KpsData/KpsS/accession_hits/" + value +
                          " /media/felix/Research/KpsData/KpsS/" + key)

        # cluster_rep_hits = cluster_file[0][mapping_file_found_representatives].unique()

        # print(mapping_file_found_representatives)


class RecordWrapper:
    """
    Difference between strain and isolate:  Isolates are viruses/bacteria isolated from a host, strain is a genome
    variant that expresses a different phenotype.

    """

    def __init__(self, record: SeqIO.SeqRecord):
        self.record = record

    @property
    def strain(self) -> Optional[str]:
        """

        Returns
        -------

        """
        try:
            return self.record.features[0].qualifiers["strain"][0]
        except KeyError:
            return None

    @property
    def isolate(self) -> Optional[str]:
        """

        Returns
        -------

        """
        try:
            return self.record.features[0].qualifiers["isolate"][0]
        except KeyError:
            return None

    @property
    def region_name(self):
        try:
            return self.record.features[2].qualifiers["region_name"][0]
        except KeyError:
            return []
        except IndexError:
            return []

    @property
    def protein_length(self):
        try:
            return self.record.features[1].location
        except KeyError:
            return None
