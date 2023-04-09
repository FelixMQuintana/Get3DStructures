import logging
import os
import threading
from pathlib import Path
from typing import Optional, List
import typer
from fasta_reader import read_fasta

from Commands.command import Command, UniProtID
from lib.const import COLABFOLD_WORKING_DIRECTORY, COLABFOLD_OPTIONS, \
    ACCESSIONS_LOOKUP_TABLE, COLABFOLDResponses, ALPHA_FOLD_STRUCTURE_EXT, ALPHA_FOLD_PAE_EXT, APP_DIRECTORY, \
    ALLOWED_EXT
from query import AlphaFoldQuery, UNIPROT_RESPONSE, PDBQuery
import dask.dataframe as df
from rich.progress import track


class StructureFile:

    def __init__(self, path: Path):
        self._path: Path = path

    @property
    def path(self) -> Path:
        return self._path


class HomologyStructure(StructureFile):
    """

    """

    def __init__(self, path: Path):
        super().__init__(path)
        self._piddt: Optional[List[float]] = None

    @property
    def piddt(self) -> List:
        if self._piddt is None:
            read_alpha_fold_file_obj = open(self.path.with_suffix(ALLOWED_EXT.CIF.value), "r")
            current_res = 0
            plddt = []
            for line in read_alpha_fold_file_obj:
                if line.startswith("ATOM") and line.split()[2] == "C":
                    if int(line.split()[8]) > current_res:
                        plddt.append(float(line.split()[14]))
                        current_res = int(line.split()[8])
            self._piddt = plddt
        return self._piddt


class CrystalStructure(StructureFile):
    """

    """


def generate_alpha_fold_structures(uniprotid: UniProtID):
    os.system(COLABFOLD_WORKING_DIRECTORY + COLABFOLD_OPTIONS.USE_GPU.value + COLABFOLD_OPTIONS.MSA_MODE.value %
              COLABFOLDResponses.MMSEQS2_UNIREF_ENV.value + COLABFOLD_OPTIONS.AMBER.value +
              COLABFOLD_OPTIONS.NUM_ENSEMBLE.value % "3" + COLABFOLD_OPTIONS.PAIR_MODE.value %
              COLABFOLDResponses.UNPAIRED_PAIRED.value + uniprotid.path.as_posix() + " " + uniprotid.path.parent.as_posix())


class Structure(Command):
    def __init__(self, uniprot_id: Optional[str] = None,
                 uniprot_ids_list: Optional[Path] = None):
        super().__init__()
        if not self.working_directory.exists():
            logging.warning("Database working directory doesn't exist. Building now")
            os.mkdir(self.working_directory)
        if uniprot_id is not None:
            self._uniprot_id_query_list: List[UniProtID] = [UniProtID(uniprot_id, self.working_directory)]
        else:
            self._uniprot_id_query_list: List[UniProtID] = self.get_list(uniprot_ids_list, self.working_directory)
            typer.secho("Building database directory")
            for uniprot_id in self._uniprot_id_query_list:
                if not os.path.exists(self.working_directory.joinpath(uniprot_id.id)):
                    logging.debug("Built directory %s" % self.working_directory.joinpath(uniprot_id.id))
                    os.mkdir(self.working_directory.joinpath(uniprot_id.id))
        self.active_threads = []

    @staticmethod
    def get_list(list_file: Path, working_directory: Path) -> List[UniProtID]:
        if not list_file.exists():
            raise FileNotFoundError(f"File doesn't seem to exist: {list_file}")
        with open(list_file, "r") as file:
            possible_uniprot_ids: List[str] = file.readlines()
        return [UniProtID(uni_id.strip("\n"), working_directory) for uni_id in possible_uniprot_ids]

    def run(self) -> None:
    #    alpha_fold_threads = self.look_up_alpha_fold_structures()
        fasta_threads: List[threading.Thread] = [threading.Thread(target=uniprot.query_fasta, args=[]) for uniprot in
                                                 self._uniprot_id_query_list]
        accession_data_threads: List[threading.Thread] = [threading.Thread(target=uniprot.query_accession_data, args=[])
                                                          for uniprot in self._uniprot_id_query_list]
        [(fasta_thread.start(), accession_data_thread.start()) for fasta_thread, accession_data_thread in
         track(zip(fasta_threads, accession_data_threads))]
        [(fasta_thread.join(), accession_data_thread.join()) for fasta_thread, accession_data_thread in
         track(zip(fasta_threads, accession_data_threads))]
    #    [(threads[0].join(), threads[1].join()) for threads in track(alpha_fold_threads)]
        [self.get_structures(uniprot) for uniprot in self._uniprot_id_query_list]
        [process.wait() for process in self.active_threads]

    def look_up_alpha_fold_structures(self):
        """

        :param list_of_uniprot_ids:
        :param database:
        :return:
        """
        complete_dataframe: df.DataFrame = df.read_csv(
            Path(APP_DIRECTORY) / "data" / ACCESSIONS_LOOKUP_TABLE,
            header=None)
        uniprot_id_strs = [uniprot.id for uniprot in self._uniprot_id_query_list]
        queried_dataframe: df.DataFrame = complete_dataframe[complete_dataframe[0].apply(
            lambda x: x in uniprot_id_strs)]
        queried_dataframe = queried_dataframe.compute()
        threads = [(threading.Thread(target=AlphaFoldQuery(self.working_directory.joinpath(accession)).query,
                                     args=[alpha_fold_model_name + ALPHA_FOLD_STRUCTURE_EXT % version + "." +
                                           self.structure_type])
                    , threading.Thread(target=AlphaFoldQuery(self.working_directory.joinpath(accession)).query,
                                       args=[alpha_fold_model_name + ALPHA_FOLD_PAE_EXT % version]))
                   for accession, alpha_fold_model_name, version in
                   zip(queried_dataframe[0], queried_dataframe[3], queried_dataframe[4])]
        [(thread[0].start(), thread[1].start()) for thread in threads]
        return threads

    def get_structures(self, uniprotid: UniProtID):
        self.get_pdb_structures(uniprotid)
        pdbs = [file for file in os.listdir(self.working_directory.joinpath(uniprotid.id)) if
                file.endswith("." + self.structure_type)]
        fatsa_len = [[len(item.sequence) for item in read_fasta(self.working_directory.joinpath(uniprotid.id).joinpath(str(file)))][0] for file in os.listdir(self.working_directory.joinpath(uniprotid.id)) if file.endswith(".fasta")][0]
        if len(pdbs) == 0 and fatsa_len <= 1000:
            print(threading.active_count())
            typer.secho(f"Generating AlphaFold Structures for UniProtID {uniprotid.id}", fg=typer.colors.BLUE)
            self.active_threads.append(self.thread_pool.apply_async(generate_alpha_fold_structures, [uniprotid]))

    def get_pdb_structures(self, uniprot_id: UniProtID) -> None:
        """

        :param uniprot_id:
        :return:
        """
        threads = [
            threading.Thread(target=PDBQuery(self.working_directory.joinpath(uniprot_id.id)).query,
                             args=[pdb_code + "." + self.structure_type]) for
            pdb_code in
            uniprot_id.structural_data.get(UNIPROT_RESPONSE.STRUCTURE.value)]
        [t.start() for t in threads]
        [t.join() for t in threads]
