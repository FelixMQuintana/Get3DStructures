import abc
import json
from enum import Enum
from typing import Optional, List, Dict
import sys

import requests
from requests.sessions import HTTPAdapter
from requests.adapters import Retry


class UNIPROT_RESPONSE(Enum):
    ACCESSION = "accession"
    STRUCTURE = "structure"
    DB_REFERENCES = "dbReferences"


class HTMLQuery:
    """

    """

    def __init__(self):
        """
        """
        self._response_history: List[Dict] = []
        self._session = self.create_http_session()

    @property
    def response_history(self):
        if self._response_history is None:
            raise Exception("History is None")
        return self._response_history

    @property
    @abc.abstractmethod
    def html_base(self) -> str:
        """

        :return:
        """
    @staticmethod
    def create_http_session():
        """
        A requests session configured with retries.
        """

        http_ = requests.Session()

        # Retry has been set for all server related errors
        retry_ = Retry(total=3, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
        adaptor = HTTPAdapter(max_retries=retry_)
        http_.mount('https://', adaptor)

        return http_

    def query(self, query: str) -> Optional[Dict]:
        """

        :param query:
        :return:
        """

        print("Querying:" + str(self.html_base) + query)

        response = self._session.get(str(self.html_base) + query, headers={"Accept": "application/json"})
        return self._query_cleanup(response)

    def _query_cleanup(self, response) -> Optional[Dict]:
        """

        :param response:
        :return:
        """
        if response.status_code == 404:
            UserWarning(f"File was not found: {response.status_code}", )
        elif not response.ok:
            response.raise_for_status()
            sys.exit(1)
        self._response_history = response
        open(response.request.url.split(self.html_base)[1], 'wb').write(response.content)
        return None


class UniProtIDQuery(HTMLQuery):
    """

    """
    def __init__(self, meta_data_file_name: str):
        super().__init__()
        self._meta_data_file_name: str = meta_data_file_name

    def _query_cleanup(self, response) -> Dict:
        """

        :param response:
        :return:
        """

        if not response.ok:
            response.raise_for_status()
            sys.exit(1)
        self._response_history = response.json()[0]
        return response.json()[0]

    @property
    def html_base(self) -> str:
        return "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession="

    def _parse_response(self, response: Dict) -> Dict:
        with open(self._meta_data_file_name, 'w') as out:
            json.dump(response, out)
        accession = response.get(UNIPROT_RESPONSE.ACCESSION.value)
        data: List[Dict] = response.get(UNIPROT_RESPONSE.DB_REFERENCES.value)
        if data is None:
            raise UserWarning(f"No known model! {data}")
        pdb_results: List[str] = [entry.get('id') for entry in data if entry.get('type') == "PDB"]
        return {UNIPROT_RESPONSE.ACCESSION.value: accession,
                UNIPROT_RESPONSE.STRUCTURE.value: pdb_results}

    def parse_response(self, response: Optional[Dict] = None) -> Dict:
        """

        :return:
        """
        return self._parse_response(response)

class FastaQuery(HTMLQuery):
    """

    """
    @property
    def html_base(self) -> str:
        return "https://rest.uniprot.org/uniprotkb/"

class PDBQuery(HTMLQuery):
    """

    """
    @property
    def html_base(self) -> str:
        return "https://files.rcsb.org/download/"


class AlphaFoldQuery(HTMLQuery):
    """

    """

    @property
    def html_base(self) -> str:
        return "https://alphafold.ebi.ac.uk/files/"
