import abc
import logging
import multiprocessing
import os
import threading
from pathlib import Path
from typing import List

#import modeller
import typer
import textwrap
from Commands.command import Command, FactoryBuilder
from dataStructure.collections import Collection, HomologyStructureFetcher, ExperimentalStructureFetcher, \
    UniProtAcessionFetcher, UniProtFastaFetcher
from dataStructure.protein.accession import UniProtIDFastaFile
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure import HomologyStructure
from lib.const import COLABFOLD_WORKING_DIRECTORY, ColabFoldOptions, \
    ACCESSIONS_LOOKUP_TABLE, COLABFOLDResponses, ALPHA_FOLD_STRUCTURE_EXT, ALPHA_FOLD_PAE_EXT, APP_DIRECTORY, \
    AllowedExt, StructureBuildMode
from query import AlphaFoldQuery, UNIPROT_RESPONSE, PDBQuery, UniProtIDQuery, FastaQuery
import dask.dataframe as df
#from modeller import *
#from modeller.automodel import *


def generate_alpha_fold_structures(fasta_path: Path):
    os.system(COLABFOLD_WORKING_DIRECTORY + ColabFoldOptions.USE_GPU.value + ColabFoldOptions.MSA_MODE.value %
              COLABFOLDResponses.MMSEQS2_UNIREF_ENV.value +  # ColabFoldOptions.AMBER.value +
              ColabFoldOptions.NUM_ENSEMBLE.value % "3" + ColabFoldOptions.PAIR_MODE.value %
              COLABFOLDResponses.UNPAIRED_PAIRED.value + fasta_path.as_posix() + " " + fasta_path.parent.as_posix())


class GetStructureData(Command, abc.ABC):

    def __init__(self, uniprot_list: Path):
        super().__init__()
        self._uniprot_list = self.get_list(uniprot_list)
        logging.debug("Building database directory")
        for uniprot_id in self._uniprot_list:
            if not os.path.exists(self.working_directory.joinpath(uniprot_id)):
                logging.debug("Built directory %s" % self.working_directory.joinpath(uniprot_id))
                os.mkdir(self.working_directory.joinpath(uniprot_id))

    @staticmethod
    def get_list(list_file: Path) -> List[str]:
        if not list_file.exists():
            raise FileNotFoundError(f"File doesn't seem to exist: {list_file}")
        with open(list_file, "r") as file:
            possible_uniprot_ids: List[str] = file.read().splitlines()
        # return [UniProtID(uni_id.strip("\n"), working_directory) for uni_id in possible_uniprot_ids]
        return possible_uniprot_ids

    @abc.abstractmethod
    def command(self) -> None:
        raise NotImplementedError

    def run(self) -> None:
        self.command()


class GetAccessionData(GetStructureData):

    def command(self) -> None:
        threads = [threading.Thread(target=UniProtIDQuery(accession_id + AllowedExt.JSON.value,
                                                          self.working_directory.joinpath(accession_id)).query,
                                    args=[accession_id]) for accession_id in self._uniprot_list]
        [thread.start() for thread in threads]
        [thread.join() for thread in threads]


class GetFastaData(GetStructureData):

    def command(self) -> None:
        threads = [threading.Thread(target=FastaQuery(self.working_directory.joinpath(accession_id)).query,
                                    args=[accession_id + AllowedExt.FASTA.value]) for accession_id in
                   self._uniprot_list]
        [thread.start() for thread in threads]
        [thread.join() for thread in threads]


class FindAlphaFoldStructures(GetStructureData):

    def command(self) -> None:
        self.look_up_alpha_fold_structures()

    def look_up_alpha_fold_structures(self):
        """

        :return:
        """
        complete_dataframe: df.DataFrame = df.read_csv(
            Path(APP_DIRECTORY) / "data" / ACCESSIONS_LOOKUP_TABLE,
            header=None)
        # uniprot_id_strs = [uniprot.id for uniprot in self._uniprot_id_query_list]
        queried_dataframe: df.DataFrame = complete_dataframe[complete_dataframe[0].apply(
            lambda x: x in self._uniprot_list)]
        queried_dataframe = queried_dataframe.compute()
        threads = [(threading.Thread(target=AlphaFoldQuery(self.working_directory.joinpath(accession)).query,
                                     args=[alpha_fold_model_name + ALPHA_FOLD_STRUCTURE_EXT % version + "." +
                                           self.structure_type])
                    , threading.Thread(target=AlphaFoldQuery(self.working_directory.joinpath(accession)).query,
                                       args=[alpha_fold_model_name + ALPHA_FOLD_PAE_EXT % version]))
                   for accession, alpha_fold_model_name, version in
                   zip(queried_dataframe[0], queried_dataframe[3], queried_dataframe[4])]
        [(thread[0].start(), thread[1].start()) for thread in threads]
        [(thread[0].join(), thread[1].join()) for thread in threads]
        return None


class ComputeAlphaFoldStructures(GetStructureData):

    def __init__(self, uniprot_list: Path):
        super().__init__(uniprot_list)
        self.collection = Collection(self.working_directory, HomologyStructureFetcher(),
                                     UniProtFastaFetcher())

    def command(self) -> None:
        [self.compute_structures(protein_structures) for protein_structures in
         self.collection.protein_structure_results.values() if
         not self.working_directory.joinpath(protein_structures.id.split('.')[0]).joinpath("log.txt").exists()]

    @staticmethod
    def compute_structures(protein_structures: ProteinStructures):
        if len(protein_structures.all_structures) == 0:  # and protein_structures.fasta:
            typer.secho(f"Generating AlphaFold Structures for UniProtID {protein_structures.id}", fg=typer.colors.BLUE)
            generate_alpha_fold_structures(protein_structures.fasta_file.path)


class FindExperimentalStructures(GetStructureData):

    def command(self) -> None:
        [self.get_pdb_structures(protein_structures) for protein_structures in
         self.collection.protein_structure_results.values()]

    def __init__(self, uniprot_list: Path):
        super().__init__(uniprot_list)
        self.collection = Collection(self.working_directory, UniProtAcessionFetcher())

    def get_pdb_structures(self, protein_structures: ProteinStructures) -> None:
        """

        Parameters
        ----------
        protein_structures
        """
        if protein_structures.uniprotID.structural_data.get(UNIPROT_RESPONSE.STRUCTURE.value) is None:
            return None

        threads = [
            threading.Thread(target=PDBQuery(self.working_directory.joinpath(protein_structures.id)).query,
                             args=[pdb_code + "." + self.structure_type]) for
            pdb_code in
            protein_structures.uniprotID.structural_data.get(UNIPROT_RESPONSE.STRUCTURE.value)]
        [t.start() for t in threads]
        [t.join() for t in threads]


def create_pir_file(fasta_file: UniProtIDFastaFile):
    logging.info("Creating PIR %s file for fasta file %s" % (fasta_file.path.with_suffix(".pir"), fasta_file.path))
    with open(fasta_file.path.with_suffix(".pir"), "w") as out_handle:
        out_handle.write(">P1;" + fasta_file.id.split(".")[0] + "\n")
        out_handle.write("sequence::1::" + str(len(fasta_file.fasta)) + ":::::\n")

        out_handle.write(textwrap.fill(fasta_file.fasta + "*", width=60))





supported_commands = {
    StructureBuildMode.ACCESSION_DATA: GetAccessionData,
    StructureBuildMode.FASTA_DATA: GetFastaData,
    StructureBuildMode.PDB_STRUCTURES: FindExperimentalStructures,
    StructureBuildMode.ALPHAFOLD_STRUCTURES: FindAlphaFoldStructures,
    StructureBuildMode.COMPUTE_ALPHA_STRUCTURES: ComputeAlphaFoldStructures,
}


class StructureFactory(FactoryBuilder):

    @staticmethod
    def build(mode: StructureBuildMode, *args, **kwargs) -> Command:
        return supported_commands[mode](*args, **kwargs)
