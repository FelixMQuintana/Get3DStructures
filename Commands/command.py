import logging
from multiprocessing import pool
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Dict

import urllib3.exceptions
from fasta_reader import read_fasta

from lib.const import ALLOWED_EXT, CONFIG_PATH, CONFIG_OPTIONS, SUPPORTED_STRUCTURE_TYPES
from lib.func import load_json
from query import FastaQuery, UniProtIDQuery


class UniProtID:
    def __init__(self, uniprot_id: str, database_location: Path):
        self._id: str = self.verify(uniprot_id)
        self._base: Path = database_location.joinpath(uniprot_id)
        self._path: Optional[Path] = None
        self._uniprot_structural_data: Optional[Dict] = None
        self._path = self._base.joinpath(self._id + ALLOWED_EXT.FASTA.value)

    @property
    def id(self) -> str:
        return self._id

    @property
    def fasta(self) -> str:
        iterator = read_fasta(self._path)
        fasta_item = iterator.read_item()
        return fasta_item.sequence

    @staticmethod
    def verify(id: str):
        return id

    @property
    def path(self) -> Path:
        return self._base.joinpath(self._id + ALLOWED_EXT.FASTA.value)

    @property
    def structural_data(self) -> Dict:
        if self._uniprot_structural_data is None:
            return {}
        return self._uniprot_structural_data

    def query_accession_data(self) -> None:
        uni_query = UniProtIDQuery(self._id, self._base)
        try:
            uni_query.query(self._id)
        except urllib3.exceptions.MaxRetryError:
            raise FileNotFoundError (f"Maximum retries to query {self.id}", self.id)
        except urllib3.exceptions.SSLError:
            raise FileNotFoundError (f"Maximum retries to query {self.id}", self.id)
        self._uniprot_structural_data = uni_query.parse_response()

    def query_fasta(self) -> None:
        logging.debug("Querying fasta for uniprotID: %s" % self.id)
        FastaQuery(self._base).query(self._id + ALLOWED_EXT.FASTA.value)


class Command(ABC):
    """

    """

    def __init__(self) -> None:
        config_json = load_json(Path(CONFIG_PATH))
        self.thread_pool = pool.ThreadPool(config_json[CONFIG_OPTIONS.THREAD_COUNT.value])
        self.working_directory: Path = Path(config_json[CONFIG_OPTIONS.DATABASE_LOCATION.value])
        self.structure_type: str = config_json[CONFIG_OPTIONS.STRUCTURE_TYPE.value]
        self.args = config_json

    @abstractmethod
    def run(self) -> None:
        """

        :return:
        """
        raise NotImplemented
