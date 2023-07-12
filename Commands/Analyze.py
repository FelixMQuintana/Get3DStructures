import json
import logging
from typing import List, Optional, Dict
import statistics

import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

import numpy

from Commands.Structure import HomologyStructure, CrystalStructure, StructureFile
from Commands.post_processing import PostProcessing, StructureResults
from lib.const import AnalysisMode, CountStatistics, SequenceLengthStatistics, Metrics, HomologyStructureStatistics
from lib.func import change_directory
#import pandas as pd
matplotlib.use('TkAgg')
log = logging.getLogger()
def make_piddt_plot(alpha_fold_structure: HomologyStructure):
    plt.plot(range(len(alpha_fold_structure.piddt)), alpha_fold_structure.piddt)
    plt.xlabel("Residue ID")
    plt.ylabel("PLDDT")
    plt.savefig(alpha_fold_structure.path.with_suffix(".TIF"))
    plt.close()


class Analyze(PostProcessing):

    def __init__(self, specific_file: Optional[Path], mode: AnalysisMode) -> None:
        super().__init__(specific_file=specific_file)
        self._mode = self.get_structures_metrics
        #if mode.value == AnalysisMode.PLDDT:
        #    self._mode = self.plot_piddts
        #else:
        #    self._mode = self.get_structures_metrics

    def run(self) -> None:
       #self.get_go_terms()
        self._mode()
    #   job = [self._mode()  for structures in self._structure_results]
      # print((sum(job))/len(job))

    def plot_piddts(self):
        [[make_piddt_plot(alpha_fold_structure) for alpha_fold_structure in structure.homology_structures if
          change_directory(self.working_directory.joinpath(structure.id))] for
         structure in self._structure_results]

    #    def get_ligand_binding_sites(self, ):


    def get_structures_metrics(self) -> None:

    #    self.remove_small_structures()
        number_of_uniprot_ids = len(self._structure_results)
      #  number_of_crystal = sum([structure.crystal_structure_count for structure in self._structure_results])
        number_of_homology = sum([structure.homology_structure_count for structure in self._structure_results])
        plddt_results = [(uniprotid.id, statistics.mean(
            [statistics.mean(homology_structure.piddt) for homology_structure in uniprotid.homology_structures]))
                         for uniprotid in self._structure_results if len(uniprotid.homology_structures) > 0]
        metrics_dict = {
            Metrics.UNIPROTID_COUNT: number_of_uniprot_ids,
      #      Metrics.CRYSTAL_STRUCTURES: {
      #          CountStatistics.COUNT.value: number_of_crystal,
      #          CountStatistics.MEAN.value: statistics.mean([structure.crystal_structure_count for structure in self._structure_results]),
      #          CountStatistics.ST_DEV.value: statistics.stdev([structure.crystal_structure_count for structure in self._structure_results]),
      #          CountStatistics.MAX.value: max([structure.crystal_structure_count for structure in self._structure_results]),
      #          CountStatistics.MIN.value: min([structure.crystal_structure_count for structure in self._structure_results]),
      #          CountStatistics.MODE.value: statistics.mode([structure.crystal_structure_count for structure in self._structure_results]),
      #          SequenceLengthStatistics.MEAN.value: statistics.mean([len(fasta) for structures in self._structure_results if structures.crystal_structure_count>0 for fasta in structures.crystal_fastas]),
      #          SequenceLengthStatistics.ST_DEV.value: statistics.stdev([len(fasta) for structures in self._structure_results if structures.crystal_structure_count>0 for fasta in structures.crystal_fastas]),
      #          SequenceLengthStatistics.MAX.value: max([max([len(fasta) for fasta in structures.crystal_fastas]) for structures in self._structure_results if structures.crystal_structure_count > 0]),
      #          SequenceLengthStatistics.MIN.value: min([min([len(fasta) for fasta in structures.crystal_fastas]) for structures in self._structure_results if structures.crystal_structure_count > 0])
      #      },
            Metrics.HOMOLOGY_STRUCTURES: {
                CountStatistics.COUNT.value: number_of_homology,
                CountStatistics.MEAN.value: statistics.mean([structure.homology_structure_count for structure in self._structure_results]),
                CountStatistics.ST_DEV.value: statistics.stdev([structure.homology_structure_count for structure in self._structure_results]),
                CountStatistics.MAX.value: max([structure.homology_structure_count for structure in self._structure_results]),
                CountStatistics.MIN.value: min([structure.homology_structure_count for structure in self._structure_results]),
                CountStatistics.MODE.value: statistics.mode([structure.homology_structure_count for structure in self._structure_results]),
                SequenceLengthStatistics.MEAN.value: statistics.mean([statistics.mean([len(fasta) for fasta in structures.homology_fastas]) for structures in self._structure_results if structures.homology_structure_count>0]),
                SequenceLengthStatistics.ST_DEV.value: statistics.stdev([len(fasta) for structures in self._structure_results if structures.homology_structure_count>0 for fasta in structures.homology_fastas]),
                SequenceLengthStatistics.MAX.value: max([max([len(fasta) for fasta in structures.homology_fastas]) for structures in self._structure_results if structures.homology_structure_count > 0]),
                SequenceLengthStatistics.MIN.value: min([min([len(fasta) for fasta in structures.homology_fastas]) for structures in self._structure_results if structures.homology_structure_count > 0]),
                HomologyStructureStatistics.PLDDT_MEAN.value: statistics.mean([plddt[1] for plddt in plddt_results])
            },
            "raw_data":{
                "plddt":plddt_results,
                "homolgy_sequence_len": [(stucts.path.name, len(stucts.fasta))
                                        for structures in self._structure_results if structures.homology_structure_count > 0 for stucts in structures.homology_structures],
        #        "crystal_sequence_len": [(stucts.path.name, len(stucts.fasta))
        #                                 for structures in self._structure_results if structures.crystal_structure_count > 0 for stucts in structures.crystal_structures],

            }

            # fasta_len = [len(Bio.SeqIO.read(self.working_directory.joinpath(uniprotid.id).joinpath(str(file)), "fasta").seq) for file in os.listdir(self.working_directory.joinpath(uniprotid.id)) if file.endswith(".fasta")][0]
        }
        json_object = json.dumps(metrics_dict)
        with open(self.working_directory.joinpath("metrics.json"),"w" ) as outfile:
            outfile.write(json_object)
        plt.hist([plddt[1] for plddt in plddt_results], bins=[60,70,75,80,85,90,100])
        plt.title("PLDDT Distribution")
        plt.ylabel("Count")
        plt.xlabel("PLDDT Score")
        plt.savefig(self.working_directory.joinpath("PLDDT_distribution.tif"))
        plt.close()
       # plt.hist([len(fatsa) for structures in self._structure_results if structures.crystal_structure_count > 0 for fatsa in structures.crystal_fastas], bins=[25,50,100,150,200,300,400,500])
       # plt.title("Crystal Fasta Distribution")
       # plt.ylabel("Count")
       # plt.xlabel("Fasta Length")
       # plt.savefig(self.working_directory.joinpath("crystal_fasta_distribution.tif"))
       # plt.close()
        plt.hist([len(fatsa) for structures in self._structure_results if structures.homology_structure_count > 0 for fatsa in structures.homology_fastas], bins=[25,50,100,150,200,300,400,500])
        plt.title("Homology Fasta Distribution")
        plt.ylabel("Count")
        plt.xlabel("Fasta Length")
        plt.savefig(self.working_directory.joinpath("homology_fasta_distribution.tif"))
        plt.close()
    #    plt.hist([structure.crystal_structure_count for structure in self._structure_results], bins=100)
    #    plt.title("Crystal Structure Distribution")
    #    plt.ylabel("Count")
    #    plt.xlabel("Number of Structures")
    #    plt.savefig(self.working_directory.joinpath("crystal_structure_distribution.tif"))
      #  plt.close()
        plt.hist([structure.homology_structure_count for structure in self._structure_results], bins=100)
        plt.title("Homology Structure Distribution")
        plt.ylabel("Count")
        plt.xlabel("Number of Structures")
        plt.savefig(self.working_directory.joinpath("homology_structure_distribution.tif"))
        plt.close()
        #print(f"Number of uniprotIDs {number_of_uniprot_ids}")
        #print(f"Number of crystals {number_of_crystal} with mean "
        #      f"{statistics.mean([structure.crystal_structure_count for structure in self._structure_results])} and stdev:"
        #      f"{statistics.stdev([structure.crystal_structure_count for structure in self._structure_results])}"
        #      f" Mode: {statistics.mode([structure.crystal_structure_count for structure in self._structure_results])} "
        #      f"and median {statistics.median([structure.crystal_structure_count for structure in self._structure_results])}"
        #      f"Max {max([structure.crystal_structure_count for structure in self._structure_results])}")
        #plt.hist([structure.crystal_structure_count for structure in self._structure_results], bins=10)
        #plt.show()
        #print(f"Number of homology {number_of_homology} and mean "
        #      f"{statistics.mean([structure.homology_structure_count for structure in self._structure_results])}"
        #      f" and stdev: {statistics.stdev([structure.homology_structure_count for structure in self._structure_results])} "
        #      f"and median: {statistics.median([structure.homology_structure_count for structure in self._structure_results])}")
        #plt.hist([structure.homology_structure_count for structure in self._structure_results], bins=10)
        #plt.show()
        #print([[statistics.mean(structure.piddt) for structure in structure.homology_structures] for structure in
        #       self._structure_results])
      #  plddt_results = [(uniprotid.id, statistics.mean(
      #      [statistics.mean(homology_structure.piddt) for homology_structure in uniprotid.homology_structures]))
      #                   for uniprotid in self._structure_results if len(uniprotid.homology_structures) > 0]
      #  plt.close()
    #    plt.plot(range(len(plddt_results)), [plddt[1] for plddt in plddt_results])
    #    plt.show()
    #    plt.savefig(self.working_directory.name + "PerStructurePLDDT.tif")
        #print(statistics.mean([plddt[1] for plddt in plddt_results]))
        #with open(self.working_directory.name + "PLDDT_RESULTS.TXT", "w") as f:
        #    [f.write(plddt_result[0] + " " + str(plddt_result[1]) + "\n") for plddt_result in plddt_results]

    def get_go_terms(self) -> Dict:
        my_dict= {}
        for structures in self._structure_results:
            terms = []
            try:
                references = structures.uniprotID.structural_data["dbReferences"]
                go="type"
            except FileNotFoundError as ex:
                continue
            except KeyError as ex:
                try:
                    references = structures.uniprotID.structural_data["uniProtKBCrossReferences"]
                    go = "database"
                except:
                    print(structures.uniprotID.path)
            except TypeError:
                references = structures.uniprotID.structural_data[0]["dbReferences"]
               # raise RuntimeError
            for term_dict in references:
                if term_dict[go] == "GO":
                    terms.append(term_dict["id"])
            my_go_terms_encoding = numpy.isin(self.go_terms, terms)
            my_go_terms_encoding = my_go_terms_encoding.astype(numpy.float32)
            for structure in structures.all_structures:
                my_dict[structure.path.name] = my_go_terms_encoding
        return my_dict

    def check_published_binding_sites(self, uniprot_id_structures: StructureFile):
        if uniprot_id_structures.binding_site_residues is None:
            return 1
        return 0
        #print(uniprot_id_structures.binding_site_residues)

    def remove_small_structures(self) -> None:
        for uniprot_id in self._structure_results:
            print(uniprot_id.id)
            if len(uniprot_id.fasta) < 50:
                self._structure_results.remove(uniprot_id)
                continue
            for struct in uniprot_id.all_structures:
                print(struct.path.name)
                if len(struct.fasta) < 50:
                    log.info(f"Removing {struct.path} from dataset due to small size with AA length of "
                             f"{len(struct.fasta)}")
                #    self._structure_results.remove(uniprot_id)
                    if isinstance(struct, HomologyStructure):
                        uniprot_id.homology_structures.remove(struct)
                    elif isinstance(struct, CrystalStructure):
                        uniprot_id.crystal_structures.remove(struct)
                    self._structure_results.append(uniprot_id)