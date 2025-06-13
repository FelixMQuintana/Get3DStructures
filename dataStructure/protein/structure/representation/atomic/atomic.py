from pathlib import Path

from Bio.PDB import PDBParser, Structure

from dataStructure.protein.structure.representation.representation import Representation


class Atomic(Representation):

    def __init__(self, pdb_filename: Path):
        super().__init__()
        self._file_path = pdb_filename
        self._pdb = PDBParser(QUIET=True).get_structure(pdb_filename.name, pdb_filename)

    @property
    def pdb(self) -> Structure:
        return self._pdb

    @property
    def file_path(self):
        return self._file_path