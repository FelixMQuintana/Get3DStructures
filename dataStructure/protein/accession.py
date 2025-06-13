import json
import logging
from pathlib import Path
from typing import Optional, Dict, List

from fasta_reader import read_fasta

from dataStructure.protein.structure.structure import StructureFile
from lib.const import AllowedExt
from query import FastaQuery

acceptable_sites = ["ACT_SITE", "BINDING", "Motif", "Binding site", "Active site", "Site"]


def query_fasta(self) -> None:
    logging.debug("Querying fasta for uniprotID: %s" % self.id)
    FastaQuery(self._base).query(self._id + AllowedExt.FASTA.value)


class UniProtIDFastaFile(StructureFile):
    def __init__(self, path: Path):
        super().__init__(path)
        self._uniprot_structural_data: Optional[Dict] = None
        self._sequence: Optional[str] = None

    @property
    def id(self) -> str:
        return self.path.parent.name

    @property
    def fasta(self) -> str:
        try:
            if self._sequence is None:
                iterator = read_fasta(self._path)
                fasta_item = iterator.read_item()
                self._sequence = fasta_item.sequence
        except StopIteration as ex:
            logging.warning("File is empty returning: \"\"")
            self._sequence = ""
        return self._sequence


class UniProtAcessionFile(StructureFile):

    def __init__(self, path: Path):
        super().__init__(path)
        self._uniprot_structural_data = None
        self._binding_site_residues = None
        self._go_terms = None
    @property
    def structural_data(self) -> Dict:
        if self._uniprot_structural_data is None:
            try:
                with open(self.path) as f:
                    self._uniprot_structural_data = json.load(f)
            except FileNotFoundError:
                raise FileNotFoundError(f"Couldn't open {self.path}")
        return self._uniprot_structural_data

    @property
    def go_terms(self):
        if self._go_terms is None:
            go_terms = []
            try:
                for entry in self.structural_data['uniProtKBCrossReferences']:
                    if entry["id"].startswith("GO:"):
                        go_terms.append(entry["id"])
                self._go_terms = go_terms
            except KeyError as ex:
                return None
        return self._go_terms

    # def query_accession_data(self) -> None:
    #    uni_query = UniProtIDQuery(self._id, self._base)
    #    try:
    #        uni_query.query(self._id + ".json")
    #    except urllib3.exceptions.MaxRetryError:
    #        raise FileNotFoundError(f"Maximum retries to query {self.id}", self.id)
    #    except urllib3.exceptions.SSLError:
    #        raise FileNotFoundError(f"Maximum retries to query {self.id}", self.id)
    #    self._uniprot_structural_data = uni_query.parse_response()

    @property
    def binding_site_residues(self):
        if self._binding_site_residues is None:

            #     try:
            raw_uniprot_data = json.load(open(self._path, "r"))  # [0]
            #      except:
            #          logging.warning(f"First path didn't work{self._path}. Trying "
            #                          f"{str(self._path).replace('Repaired_structs', 'raw_structs')}")
            #          raw_uniprot_data = json.load(
            #              open(str(self._path).replace('Repaired_structs', 'raw_structs'), "r"))[0]
            try:
                features: List = raw_uniprot_data["features"]
            except KeyError:
                return []
            self._binding_site_residues = []
            for feature in features:
                if feature["type"] in acceptable_sites:
                    try:
                        self._binding_site_residues.extend(range(int(feature["begin"]), int(feature["end"]) + 1))
                    except KeyError as ex:
                        logging.warning("Getting features failed with first method. Attempting second method.")
                        self._binding_site_residues.extend(range(int(feature["location"]["start"]["value"]),
                                                                 int(feature["location"]["end"]["value"]) + 1))
        return self._binding_site_residues
