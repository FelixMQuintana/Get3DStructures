from abc import ABC, abstractmethod
from pathlib import Path

from lib.const import ALLOWED_EXT
from uniprot import FastaQuery, UniProtIDQuery


class UniProtID:
    def __init__(self, uniprot_id: str):
        self._id: str = self.verify(uniprot_id)

    @property
    def id(self) -> str:
        return self._id

    @staticmethod
    def verify(id: str):
        return id

    def query(self) -> dict:
        FastaQuery().query(self._id + ALLOWED_EXT.FASTA.value)
        uni_query = UniProtIDQuery(self._id + ALLOWED_EXT.JSON.value)
        return uni_query.parse_response(uni_query.query(self._id))


class Command(ABC):
    """

    """
    def __init__(self, working_directory: str) -> None:
        self.working_directory: Path = Path(working_directory)


    @abstractmethod
    def run(self) -> None:
        """

        :return:
        """
        raise NotImplemented
