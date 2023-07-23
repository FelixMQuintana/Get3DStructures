"""

"""
import json
import logging
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO

from lib.const import ALLOWED_EXT


class SupportedFileType:
    """

    """

    def __init__(self, path: Path):
        self._path = path

    @property
    def path(self) -> Path:
        return self._path



class GOLabelsFile(SupportedFileType):
    """

    """


class StructureFile(SupportedFileType):

    def __init__(self, path: Path):
        super().__init__(path)
        self._path: Path = path
        self._sequence = ""
        self._binding_site_residues = None

    @property
    def id(self):
        return self._path.parent.name

    @property
    def fasta(self) -> str:
        if len(self._sequence) == 0:
            pdb_parser = SeqIO.parse(self._path, "pdb-atom")
            for record in pdb_parser:
                self._sequence += record.seq
        return self._sequence


class HomologyStructure(StructureFile):
    """

    """

    def __init__(self, path: Path):
        super().__init__(path)
        self._piddt: Optional[List[float]] = None

    @property
    def piddt(self, ) -> List:
        if self._piddt is None:
            read_alpha_fold_file_obj = open(self.path.with_suffix(ALLOWED_EXT.PDB.value), "r")
            #   current_res = 0
            plddt = []
            for line in read_alpha_fold_file_obj:
                if line.startswith("ATOM"):  # and line.split()[2] == "C":
                    #         if int(line.split()[5]) > current_res or int(line.split()[5])==0:
                    plddt.append(float(line.split()[-2]))
            #             current_res = int(line.split()[5])
            self._piddt = plddt
        return self._piddt


class CrystalStructure(StructureFile):
    """

    """
