from typing import List
import statistics

import matplotlib
import matplotlib.pyplot as plt

from Commands.Structure import HomologyStructure
from Commands.post_processing import PostProcessing
from lib.func import change_directory

matplotlib.use('TkAgg')


def make_piddt_plot(alpha_fold_structure: HomologyStructure):
    plt.plot(range(len(alpha_fold_structure.piddt)), alpha_fold_structure.piddt)
    plt.xlabel("Residue ID")
    plt.ylabel("PLDDT")
    plt.savefig(alpha_fold_structure.path.with_suffix(".TIF"))
    plt.close()


class Analyze(PostProcessing):

    def __init__(self, working_directory: str, all_files: bool, specific_file: str, piddt_plots: bool) -> None:
        super().__init__(working_directory=working_directory, all_files=all_files, specific_file=specific_file)
        if piddt_plots:
            self._mode = self.plot_piddts

    def run(self) -> None:
        self._mode()

    def plot_piddts(self):
        [[make_piddt_plot(alpha_fold_structure) for alpha_fold_structure in structure.homology_structures if
          change_directory(self.working_directory.joinpath(structure.id), skip=False)] for
         structure in self._structure_results]

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
        print([[statistics.mean(structure.piddt) for structure in structure.homology_structures] for structure in
               self._structure_results])
