import os
from pathlib import Path
from typing import Optional, List, Dict

from Commands.command import Command, UniProtID
from lib.const import COLABFOLD_WORKING_DIRECTORY, ALLOWED_EXT, COLABFOLD_OPTIONS, \
    ACCESSIONS_LOOKUP_TABLE, COLABFOLDResponses, ALPHA_FOLD_STRUCTURE_EXT, ALPHA_FOLD_PAE_EXT
from lib.func import change_directory
from query import AlphaFoldQuery, UNIPROT_RESPONSE, PDBQuery
import dask.dataframe as df


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
    def piddt(self) -> List[float]:
        return self._piddt

    @piddt.setter
    def piddt(self, piddt_list: List[float]) -> None:
        self._piddt = piddt_list


class CrystalStructure(StructureFile):
    """

    """


def generate_alpha_fold_structures(uniprotid: UniProtID):
    os.system(COLABFOLD_WORKING_DIRECTORY + COLABFOLD_OPTIONS.USE_GPU.value + COLABFOLD_OPTIONS.MSA_MODE.value %
              COLABFOLDResponses.MMSEQS2_UNIREF_ENV.value + COLABFOLD_OPTIONS.AMBER.value +
              COLABFOLD_OPTIONS.NUM_ENSEMBLE.value % "5" + COLABFOLD_OPTIONS.PAIR_MODE.value %
              COLABFOLDResponses.UNPAIRED_PAIRED.value + uniprotid.id + ALLOWED_EXT.FASTA.value + " ./")


def get_pdb_structures(uniprot_pdb_structure_names: Dict) -> None:
    """

    :param uniprot_pdb_structure_names:
    :param database_location:
    :return:
    """
    [PDBQuery().query(pdb_code + ALLOWED_EXT.CIF.value) for pdb_code in
     uniprot_pdb_structure_names.get(UNIPROT_RESPONSE.STRUCTURE.value)]


class Structure(Command):
    def __init__(self, working_directory: str, uniprot_id: Optional[str] = None,
                 uniprot_ids_list: Optional[str] = None, skip: bool = False):
        super().__init__(working_directory)
        if uniprot_id is not None:
            self._uniprot_id_query_list: List[UniProtID] = [UniProtID(uniprot_id)]
        else:
            self._uniprot_id_query_list: List[UniProtID] = self.get_list(uniprot_ids_list)
        self._skip: bool = skip

    @staticmethod
    def get_list(list_file: str) -> List[UniProtID]:
        if not os.path.isfile(list_file):
            raise FileNotFoundError(f"File doesn't seem to exist: {list_file}")
        with open(list_file, "r") as file:
            possible_uniprot_ids: List[str] = file.readlines()
        return [UniProtID(uni_id.strip("\n")) for uni_id in possible_uniprot_ids]

    def run(self) -> None:
        self.look_up_alpha_fold_structures()
        [self.get_structures(uniprot) for uniprot in self._uniprot_id_query_list if
         change_directory(self.working_directory.joinpath(uniprot.id), skip=self._skip)]

    def look_up_alpha_fold_structures(self):
        """

        :param list_of_uniprot_ids:
        :param database:
        :return:
        """
        complete_dataframe: df.DataFrame = df.read_csv(
            self.working_directory.joinpath(ACCESSIONS_LOOKUP_TABLE).absolute(),
            header=None)
        uniprot_id_strs = [uniprot.id for uniprot in self._uniprot_id_query_list]
        queried_dataframe: df.DataFrame = complete_dataframe[complete_dataframe[0].apply(
            lambda x: x in uniprot_id_strs)]
        [(AlphaFoldQuery().query(alpha_fold_model_name + ALPHA_FOLD_STRUCTURE_EXT % version),
          AlphaFoldQuery().query(alpha_fold_model_name + ALPHA_FOLD_PAE_EXT % version))
         for accession, alpha_fold_model_name, version in
         zip(queried_dataframe[0], queried_dataframe[3], queried_dataframe[4]) if
         change_directory(directory=self.working_directory.joinpath(accession), skip=self._skip)]

    def get_structures(self, uniprotid: UniProtID):
        get_pdb_structures(uniprotid.query())
        pdbs = [file for file in os.listdir(self.working_directory.joinpath(uniprotid.id)) if
                file.endswith(ALLOWED_EXT.CIF.value)]
        if len(pdbs) == 0:
            generate_alpha_fold_structures(uniprotid)
