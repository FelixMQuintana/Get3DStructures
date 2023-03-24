import logging
from multiprocessing import pool
from abc import ABC, abstractmethod
from pathlib import Path
import threading
from typing import Optional, Dict

import typer

from lib.const import ALLOWED_EXT, CONFIG_PATH, THREAD_COUNT
from lib.func import load_json
from query import FastaQuery, UniProtIDQuery


class UniProtID:
    def __init__(self, uniprot_id: str, database_location: Path):
        self._id: str = self.verify(uniprot_id)
        self._base: Path = database_location.joinpath(uniprot_id)
        self._path: Optional[Path] = None
        self._uniprot_structural_data: Optional[Dict] = None
    @property
    def id(self) -> str:
        return self._id

    @staticmethod
    def verify(id: str):
        return id

    @property
    def path(self) -> Path:
        return self._base.joinpath(self._id + ALLOWED_EXT.FASTA.value)

    @property
    def structural_data(self) -> Dict:
        return self._uniprot_structural_data

    def query_accession_data(self) -> None:
        uni_query = UniProtIDQuery(self._id, self._base)
        uni_query.query(self._id)
        self._uniprot_structural_data = uni_query.parse_response()

    def query_fasta(self) -> None:
        logging.debug("Querying fasta for uniprotID: %s" % self.id)
        FastaQuery(self._base).query(self._id + ALLOWED_EXT.FASTA.value)
        self._path = self._base.joinpath(self._id + ALLOWED_EXT.FASTA.value)


class Command(ABC):
    """

    """

    def __init__(self, working_directory: Path) -> None:
        self.working_directory: Path = Path(working_directory)
        config_json = load_json(Path(CONFIG_PATH))
        self.thread_pool = pool.ThreadPool(config_json[THREAD_COUNT])

    @abstractmethod
    def run(self) -> None:
        """

        :return:
        """
        raise NotImplemented
