import numpy as np

from dataStructure.protein.structure.representation.surface.surface import MolecularSurfaceRepresentation
from lib.const import kd_scale
from proteinAnalysis.proteinAnalysis import ProteinAnalysisOperations


class HydrophobicityProteinAnalysis(ProteinAnalysisOperations):
    @classmethod
    def registered_methods(cls) -> dict:
        return {MolecularSurfaceRepresentation.__name__: "surface"}

    def surface(self):
        hp = np.zeros(len(self._structure_representation.names))
        for ix, name in enumerate(self._structure_representation.names):
            aa = name.split("_")[3]
            hp[ix] = kd_scale[aa]
        return hp
