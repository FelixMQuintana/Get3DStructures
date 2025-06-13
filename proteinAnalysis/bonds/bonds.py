from subprocess import Popen, PIPE

from dataStructure.protein.structure.representation.atomic.atomic import Atomic
from proteinAnalysis.proteinAnalysis import ProteinAnalysisOperations


class BondProteinAnalysis(ProteinAnalysisOperations):
    @staticmethod
    def registered_methods() -> dict:
        return {Atomic.__name__: "atomic", }

    def atomic(self):
        args = [
            "python",
            "/home/felix/Downloads/getcontacts/get_static_contacts.py",
            "--structure",
            self._structure_representation.file_path,
            "--itypes",
            "all",
            "--output",
            "bond_contacts_"+str(self._structure_representation.file_path.parent.name)+".tsv",
        ]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(self._structure_representation.file_path.parent))
        stdout, stderr = p2.communicate()
        assert len(stderr) == 0
