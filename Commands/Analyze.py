import os
from pathlib import Path
from typing import List, Optional
import statistics
import matplotlib.pyplot as plt

from Commands.Structure import CrystalStructure, HomologyStructure, get_structure_files
from Commands.command import UniProtID, Command


class StructureResults:

    def __init__(self, uniprot: UniProtID, structures: tuple[List[CrystalStructure], List[HomologyStructure]]) -> None:
        self._accession: UniProtID = uniprot
        self._crystal_structures: List[CrystalStructure] = structures[0]
        self._homology_structures: List[HomologyStructure] = structures[1]

    @property
    def id(self) -> str:
        return self._accession.id

    @property
    def crystal_structures(self) -> List[CrystalStructure]:
        return self._crystal_structures

    @property
    def homology_structures(self) -> List[HomologyStructure]:
        return self._homology_structures

    @property
    def crystal_structure_count(self) -> int:
        return len(self._crystal_structures)

    @property
    def homology_structure_count(self) -> int:
        return len(self._homology_structures)


def read_af_piddt(alpha_fold_structure_path: HomologyStructure) -> List[float]:
    read_alpha_fold_file_obj = open(alpha_fold_structure_path.path, "r")
    return [float(line.split()[10]) for line in read_alpha_fold_file_obj
            if line.startswith("ATOM") and line.split()[2] == "C"]


def make_piddt_plot(piddt_list: List[float]):
    plt.plot(piddt_list)
    plt.show()

class Analyze(Command):

    def __init__(self, working_directory: str, all: bool, piddt_plots:bool) -> None:
        super().__init__(working_directory)
        self._structure_results: List[StructureResults] = \
            [StructureResults(UniProtID(self.working_directory.joinpath(str(directories)).as_uri()),
                              get_structure_files(self.working_directory.joinpath(str(directories))))
             for directories in os.listdir(self.working_directory)]
        if all:
            self._mode = self.get_structures_metrics
        elif piddt_plots:
            self._mode = self.plot_piddts

    def run(self) -> None:
        self._mode()

    def plot_piddts(self):
        [[make_piddt_plot(alpha_fold_structure.piddt) for alpha_fold_structure in structure.homology_structures] for structure in self._structure_results]

    def get_structures_metrics(self) -> None:
        number_of_uniprot_ids = len(self._structure_results)
        number_of_crystal = sum([structure.crystal_structure_count for structure in self._structure_results])
        number_of_homology = sum([structure.homology_structure_count for structure in self._structure_results])
        print(f"Number of uniprotIDs {number_of_uniprot_ids}")
        print(f"Number of crystals {number_of_crystal} with mean "
              f"{statistics.mean([structure.crystal_structure_count for structure in self._structure_results])} and stdev:"
              f"{statistics.stdev([structure.crystal_structure_count for structure in self._structure_results])}"
              f" Mode: {statistics.mode([structure.crystal_structure_count for structure in self._structure_results])} "
              f"and median {statistics.median([structure.crystal_structure_count for structure in self._structure_results])}"
              f"Max {max([structure.crystal_structure_count for structure in self._structure_results])}")
        plt.hist([structure.crystal_structure_count for structure in self._structure_results], bins=10)
        plt.show()
        print(f"Number of homology {number_of_homology} and mean "
              f"{statistics.mean([structure.homology_structure_count for structure in self._structure_results])}"
              f" and stdev: {statistics.stdev([structure.homology_structure_count for structure in self._structure_results])} "
              f"and median: {statistics.median([structure.homology_structure_count for structure in self._structure_results])}")
        plt.hist([structure.homology_structure_count for structure in self._structure_results], bins=10)
        plt.show()
        print([[statistics.mean(structure.piddt) for structure in structure.homology_structures] for structure in self._structure_results])
        # explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

        # fig, ax = plt.subplots()
        # ax.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        #   shadow=True, startangle=90)
        # plt.show()
