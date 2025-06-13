import unittest
import os
from typing import List

from dataStructure.collections import Collection
from pathlib import Path


import unittest
from pathlib import Path

from dataStructure.protein.structure.structure import SupportedFileType, StructureFile, HomologyStructure, MeshFile
import unittest
from pathlib import Path
from dataStructure.collections import Collection
# ... (other imports)

class TestCollection(unittest.TestCase):

    def setUp(self):
        self.working_directory = Path(__file__).parent / "data"

    # Happy path tests

    def test_add_fetchers_with_valid_input(self):
        # Arrange
        class MockDataFetcher:
            def get_data(self, file_path: Path) -> List[StructureFile]:
                return [HomologyStructure(Path(str(file_path) +"/alphafold_structure.pdb")), HomologyStructure(Path(str(file_path)+"/replicated_alphafold_structure.pdb"))]

        mock_data_fetcher = MockDataFetcher()
        collection = Collection(self.working_directory, mock_data_fetcher)

        # Act
        collection.add_fetchers(mock_data_fetcher)

        # Assert
        print(collection.protein_structure_results.values())
        self.assertEqual(len(collection.protein_structure_results.values()), 2)

    def test_remove_redundant_proteins_with_valid_input(self):
        # Arrange
        class MockDataFetcher:
            def get_data(self, file_path: Path) -> List[StructureFile]:
                return [HomologyStructure(Path("test_homology.pdb")), HomologyStructure(Path("test_homology2.pdb"))]

        mock_data_fetcher = MockDataFetcher()
        collection = Collection(self.working_directory, mock_data_fetcher)

        # Act
        collection.add_fetchers(mock_data_fetcher)
        collection.remove_redundant_proteins()

        # Assert
        self.assertEqual(len(collection.protein_structure_results), 1)

    # Edge case tests

    def test_add_fetchers_with_empty_working_directory(self):
        # Arrange
        class MockDataFetcher:
            def get_data(self, file_path: Path) -> List[StructureFile]:
                return []

        mock_data_fetcher = MockDataFetcher()
        working_directory = Path("/dev/null")

        # Act
        collection = Collection(working_directory, mock_data_fetcher)
        collection.add_fetchers(mock_data_fetcher)

        # Assert
        self.assertEqual(collection.protein_structure_results, {})

    def test_remove_redundant_proteins_with_empty_working_directory(self):
        # Arrange
        class MockDataFetcher:
            def get_data(self, file_path: Path) -> List[StructureFile]:
                return []

        mock_data_fetcher = MockDataFetcher()
        working_directory = Path("/dev/null")

        # Act
        collection = Collection(working_directory, mock_data_fetcher)
        collection.remove_redundant_proteins()

        # Assert
        self.assertEqual(collection.protein_structure_results, {})

if __name__ == "__main__":
    unittest.main()



#def test_collections():
#    print("Testing collections")
#    project_root_dir = os.path.dirname(os.path.abspath(__file__))
#    my_collection = collections.Collection(Path(project_root_dir+ "/data/"), collections.HomologyStructureFetcher())
##    print(my_collection.protein_structure_results)
 #   assert len(my_collection.protein_structure_results.values()) == 2
#    my_collection.remove_redundant_proteins()
#    assert len(my_collection.protein_structure_results.values()) == 1
