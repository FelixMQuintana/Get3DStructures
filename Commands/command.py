"""

"""
import logging
from multiprocessing import Pool
from abc import ABC, abstractmethod
from pathlib import Path

import urllib3.exceptions

from lib.const import CONFIG_PATH, ConfigOptions
from lib.func import load_json


class Command(ABC):
    """

    """

    def __init__(self,working_dir, structure_type) -> None:
        #config_json = load_json(Path(CONFIG_PATH))
        self.thread_pool =Pool(28) # Pool(processes=config_json[ConfigOptions.THREAD_COUNT.value])
        self.working_directory: Path = Path(working_dir)# Path(config_json[ConfigOptions.DATABASE_LOCATION.value])
        self.structure_type: str = structure_type #config_json[ConfigOptions.STRUCTURE_TYPE.value]
        #self.args = config_json

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
