import json
import logging
from multiprocessing import pool
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Dict

import urllib3.exceptions
from fasta_reader import read_fasta

from lib.const import ALLOWED_EXT, CONFIG_PATH, CONFIG_OPTIONS
from lib.func import load_json
from query import FastaQuery, UniProtIDQuery


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

    def exception_handler(self, args):
        logging.warning(f'Thread failed: {args.exc_value}')
        # if isinstance(args.exc_value, FileNotFoundError):
        #    self._uniprot_id_query_list.remove(args.exc_value.args[1])


class FactoryBuilder(ABC):

    @staticmethod
    @abstractmethod
    def build(*args, **kwargs) -> Command:
        raise NotImplemented
