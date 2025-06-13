import numpy as np
from Bio.PDB import PDBParser, Selection

from dataStructure.protein.structure.representation.surface.surface import MolecularSurfaceRepresentation
from proteinAnalysis.electrostatics._helper import compute_satisfied_CO_HN, compute_charge_helper
from proteinAnalysis.proteinAnalysis import ProteinAnalysisOperations


class ElectrostaticChargeProteinAnalysis(ProteinAnalysisOperations):

    @classmethod
    def registered_methods(cls) -> dict:
        return {MolecularSurfaceRepresentation.__name__: "surface"}

    def surface(self):
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure(self._structure_representation.file_name,
                                      self._structure_representation.file_name)
        residues = {}
        for res in struct.get_residues():
            chain_id = res.get_parent().get_id()
            if chain_id == "":
                chain_id = " "
            residues[(chain_id, res.get_id())] = res
        atoms = Selection.unfold_entities(struct, "A")
        satisfied_CO, satisfied_HN = compute_satisfied_CO_HN(atoms)

        charge = np.array([0.0] * len(self._structure_representation.vertices))
        # Go over every vertex
        for ix, name in enumerate(self._structure_representation.names):
            fields = name.split("_")
            chain_id = fields[0]
            if chain_id == "":
                chain_id = " "
            if fields[2] == "x":
                fields[2] = " "
            res_id = (" ", int(fields[1]), fields[2])
            aa = fields[3]
            atom_name = fields[4]
            # Ignore atom if it is BB and it is already satisfied.
            if atom_name == "H" and res_id in satisfied_HN:
                continue
            if atom_name == "O" and res_id in satisfied_CO:
                continue
            # Compute the charge of the vertex
            charge[ix] = compute_charge_helper(
                atom_name, residues[(chain_id, res_id)], self._structure_representation.vertices[ix]
            )

        return charge
