import abc

from dataStructure.protein.structure.representation.representation import Representation


class ProteinAnalysisOperations(abc.ABC):
    """
    Registration of all operations that can be performed on a give StructureRepresentation object.
    """

    def __init__(self, structure_rep: Representation):
        self._structure_representation: Representation = structure_rep

    @staticmethod
    @abc.abstractmethod
    def registered_methods() -> dict:
        raise NotImplementedError

    def compute(self):
        for subclass in Representation.__subclasses__():
            if subclass.__name__ == type(self._structure_representation).__name__:
                return getattr(self, self.registered_methods()[subclass.__name__])()
