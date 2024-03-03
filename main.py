"""
Database: used for my asshole
"""
import json
import os.path
from typing import Optional

from pathlib import Path

# from Commands.Analyze import Analyze
from Commands.Characteristics import CharacteristicsFactory, ClusterGOTerms, FindCustomBindingSite
from Commands.Structure.Structure import StructureFactory
from Commands.Structure.model import HomologyModel
from Commands.repair import RepairStructures
from Commands.fasta_parser import FastaData, ParseMapping, MappPDBToGenBankHits
from Commands.similarity_metrics import SimilarityMetricBuilder
from lib.const import AnalysisMode, APP_DIRECTORY, StructureCharacteristicsMode, ConfigOptions, \
    SupportedStructureFileTypes, StructureBuildMode, SimilarityDistance
# from Commands.repair import RepairPDB
import typer
import logging
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
app = typer.Typer(pretty_exceptions_show_locals=False)


@app.command()
def find_characteristics(mode: StructureCharacteristicsMode =
                         typer.Argument(..., help="Mode for characteristics to be calculated."),
                         binding_site_database: Path = typer.Argument(..., help="Please specify where to store the "
                                                                                "database of binding sites")):
    """
    This command is for finding characteristics of proteins.
    """
    if mode == StructureCharacteristicsMode.AUTOSITE:
        auto_site_directory = typer.prompt(text="Autosite binaries directory?", default="~/ADFRsuite-1.0/bin/")
        with open(APP_DIRECTORY + "/config.json", "r") as json_config:
            my_json = json.load(json_config)
            my_json[StructureCharacteristicsMode.AUTOSITE.value] = auto_site_directory
        with open(APP_DIRECTORY + "/config.json", "w") as json_config:
            json.dump(my_json, json_config)
    command = CharacteristicsFactory().build(mode, Path(os.path.abspath(binding_site_database)))
    command.run()


@app.command()
def protein_similarity(mode: SimilarityDistance = typer.Argument(..., help="Mode of how to calculate "
                                                                           "pairwise-distances of structures"),
                       distances_file_name: str = typer.Argument(...,
                                                                 help="Distances file name to be saved as numpy file")):
    """
    This command is for calculating pairwise distances of structures from a dataset of structures.
    """
    command = SimilarityMetricBuilder().build(mode, distances_file_name)
    command.run()


@app.command()
def cluster_go_terms(uniprot_id_list: Path = typer.Option(..., )):
    command = ClusterGOTerms(uniprot_id_list)
    command.run()


@app.command()
def tm_score():
    command = FindCustomBindingSite("nothingyet")
    command.run()


@app.command()
def get_structures(mode: StructureBuildMode =
                   typer.Argument(..., help="Mode for getting structures")
                   , uniprot_id_list: Path = typer.Option(...,
                                                          help="Path to UniProtID list to find structures for."),
                   sequence_db: Optional[Path] = None):
    """
    This command is to populate structures based on provided uniprotID(s). Structures are pulled from the AF2 database,
    PDB databank. If no structure is available, ColabFold will run locally to generate the structure.
    """

    command = StructureFactory().build(mode, uniprot_id_list)
    command.run()


# @app.command()
# def homology_model(uniprot_id_list: Path = typer.Option(...,
#                                                          help="Path to UniProtID list to find structures for."),
#                   target_sequence: Optional[Path] = typer.Option(None, help="For homology modeling")):
#    HomologyModel(uniprot_id_list, target_sequence).run()
@app.command()
def repair(new_dataset_location: Path = typer.Argument(..., help="Where you'd like to save the new dataset?")):
    """This command is to repair structures from a specified database or structure file. The repaired database is
    stored in new_dataset_location"""

    command = RepairStructures(Path(os.path.abspath(new_dataset_location)))
    command.run()


@app.command()
def analyze(mode: AnalysisMode = typer.Argument(..., help="Mode for analysis.", rich_help_panel="Mode")):
    """
    This command is to analyze a given database of PDB. Supported analysis modes: plddt - calculates plddt values of
    homology structures made from AF2; STATS - provides stats on database provided such as number of crystal and
    homology structures and average values per uniProtID..
    """
    command = Analyze(mode)
    command.run()


@app.command()
def read_sequence_data(file_name: str = typer.Argument(..., help="Name of fasta file of interest."),
                       gene_name: str = typer.Argument(...,
                                                       help="Gene name to search for in GP file under region name"),
                       lower_bound: int = typer.Argument(..., help="lower bound for AA sequence size )temporary "
                                                                   "implementation needs tp be more robust"),
                       upper_bound: int = typer.Argument(...,
                                                         help="Upper bound for AA sequence size (temporary "
                                                              "implementation"),
                       out_file_path: str = typer.Argument(..., help="Output fasta file path and name for output")):
    FastaData(file_name, gene_name, lower_bound, upper_bound, out_file_path).run()


@app.command()
def parse_mapping(file_name: str = typer.Argument(..., help="Name of fasta file of interest."),
                  out_file_path: str = typer.Argument(..., help="Output fasta file path and name for output")):
    ParseMapping(file_name, out_file_path).run()


@app.command()
def map_represntatives_pdb(mapping_json_file: str = typer.Argument(..., help="Name of fasta file of interest."),
                           cluster_representatives_tsv_file: str
                           = typer.Argument(..., help="Output fasta file path and name for output")):
    MappPDBToGenBankHits(mapping_json_file, cluster_representatives_tsv_file).run()


@app.callback()
def main(database: Path = typer.Argument(..., help="Path of database of interest for different modes."),
         number_of_threads: int = typer.Option(default=os.cpu_count(),
                                               #  prompt="Specify the number of threads for this process?",
                                               help="Option specifying number of threads to use for this process. "
                                                    "Default is number of logical cores."),
         mode: SupportedStructureFileTypes = typer.Argument(SupportedStructureFileTypes.CIF,
                                                            help="Structure type to perform analysis on.")):
    database = os.path.abspath(Path(database))
    config = {
        ConfigOptions.DATABASE_LOCATION.value: database,
        ConfigOptions.THREAD_COUNT.value: number_of_threads,
        ConfigOptions.STRUCTURE_TYPE.value: mode.value
    }
    with open(APP_DIRECTORY + "/config.json", "w") as json_config:
        json.dump(config, json_config)


if __name__ == '__main__':
    app()
