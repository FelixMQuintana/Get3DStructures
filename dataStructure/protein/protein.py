from typing import List

from dataStructure.protein.accession import UniProtAcessionFile, UniProtIDFastaFile
from dataStructure.protein.structure.structure import CrystalStructure, HomologyStructure, StructureFile, MeshFile
from lib.func import Metrics


class ProteinStructures:

    def __init__(self, structure_files: List[StructureFile]):
        self._crystal_structures: List[CrystalStructure] = []
        self._homology_structures: List[HomologyStructure] = []
        self._id = None
        self.mesh = None
        for structure_file in structure_files:
            if isinstance(structure_file, CrystalStructure):
                self._crystal_structures.append(structure_file)
            elif isinstance(structure_file, HomologyStructure):
                self._homology_structures.append(structure_file)
            elif isinstance(structure_file, UniProtAcessionFile):
                self._accession: UniProtAcessionFile = structure_file
            elif isinstance(structure_file, UniProtIDFastaFile):
                self._uniprot_fasta: UniProtIDFastaFile = structure_file
            elif isinstance(structure_file, MeshFile):
                self.mesh: MeshFile = structure_file
            if self._id is None:
                self._id = structure_file.id
        self._annotation = None

    def add_metric(self, metric: Metrics):
        self._annotation = metric.add_metric(self.id)

    @property
    def all_files(self) -> List[StructureFile]:
        structs = [self.fasta_file, self.mesh, self._accession]
        structs.extend(self.homology_structures)
        structs.extend(self.crystal_structures)
        cleaned_data = [x for x in structs if x is not None]

        return cleaned_data

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

    @property
    def annotations(self) -> List:
        return self._annotation