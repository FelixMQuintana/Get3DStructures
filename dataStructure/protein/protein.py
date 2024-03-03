from typing import List

from dataStructure.protein.accession import UniProtAcessionFile, UniProtIDFastaFile
from dataStructure.protein.structure import CrystalStructure, HomologyStructure, StructureFile


class ProteinStructures:

    def __init__(self, structure_files: List[StructureFile]):
        self._crystal_structures: List[CrystalStructure] = []
        self._homology_structures: List[HomologyStructure] = []
        self._id = None
        for structure_file in structure_files:
            if isinstance(structure_file, CrystalStructure):
                self._crystal_structures.append(structure_file)
            elif isinstance(structure_file, HomologyStructure):
                self._homology_structures.append(structure_file)
            elif isinstance(structure_file, UniProtAcessionFile):
                self._accession: UniProtAcessionFile = structure_file
            elif isinstance(structure_file, UniProtIDFastaFile):
                self._uniprot_fasta: UniProtIDFastaFile = structure_file
            if self._id is None:
                self._id = structure_file.id

    @property
    def id(self) -> str:
        return self._id

    @property
    def fasta(self) -> str:
        return self._accession.fasta

    @property
    def crystal_structures(self) -> List[CrystalStructure]:
        return self._crystal_structures

    @property
    def homology_structures(self) -> List[HomologyStructure]:
        return self._homology_structures

    @property
    def all_structures(self) -> List[StructureFile]:
        return [*self._crystal_structures, *self._homology_structures]

    @property
    def crystal_structure_count(self) -> int:
        return len(self._crystal_structures)

    @property
    def homology_structure_count(self) -> int:
        return len(self._homology_structures)

    @property
    def homology_fastas(self) -> List[str]:
        return [struct.fasta for struct in self.homology_structures]

    @property
    def crystal_fastas(self) -> List[str]:
        return [struct.fasta for struct in self.crystal_structures]

    @property
    def uniprotID(self, ) -> UniProtAcessionFile:
        """

        Returns
        -------

        """
        return self._accession

    @property
    def fasta_file(self) -> UniProtIDFastaFile:
        """

        Returns
        -------

        """
        return self._uniprot_fasta
