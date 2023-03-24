"""
Database: used for my asshole
"""
import json
import os.path
from typing import Optional

from pathlib import Path

from Commands.Analyze import Analyze
from Commands.Characteristics import Characteristics
from Commands.Structure import Structure
from lib.const import ALLOWED_EXT, AnalysisMode, APP_DIRECTORY, StructureCharacteristicsMode, THREAD_COUNT
from Commands.repair import RepairPDB
import typer
import logging
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
app = typer.Typer(pretty_exceptions_show_locals=False)


@app.command()
def find_characteristics(database: Path = typer.Argument(..., help="Path to database."),
                         mode: StructureCharacteristicsMode =
                         typer.Argument(..., help="Mode for characteristics to be calculated."),
                         structure_file: Optional[Path] =
                         typer.Option(None,
                                      help="Specific structure file to find characteristics based on mode selected")):
    """
    This command is for finding characteristics of proteins.
    """
    command = Characteristics(Database, structure_file, mode)
    command.run()


@app.command()
def get_structures(#database: Path = typer.Argument(..., help="Path of where to create a database."),
                   uniprot_id: Optional[str] = typer.Option(None, help="Specific UnitProtID to get structures for."),
                   uniprot_id_list: Path = typer.Option(APP_DIRECTORY + "/data/UniProtIDs.txt",
                                                        help="Path to UniProtID list to find structures for.")):
    """
    This command is to populate structures based on provided uniprotID(s). Structures are pulled from the AF2 database,
    PDB databank. If no structure is available, ColabFold will run locally to generate the structure.
    """
    command = Structure(Database, uniprot_id, Path(os.path.abspath(uniprot_id_list)))
    command.run()


@app.command()
def repair(database: Path = typer.Argument(..., help="Path of database to perform repair on."),
           new_dataset_location: Path = typer.Option(Path(APP_DIRECTORY + "/data/database/"),
                                                     prompt=f"Is this where you'd like to save the new dataset?"),
           specific_file: Optional[Path] = typer.Option(None, help="Specific file to perform repair on.")):
    """This command is to repair structures from a specified database or structure file. The repaired database is
    stored in new_dataset_location"""
    if specific_file is not None and specific_file.suffix not in {ALLOWED_EXT.CIF.value, ALLOWED_EXT.PDB.value}:
        raise FileExistsError(f"File type provided {specific_file.suffix} isn't supported. Supported types are: "
                              f"{ALLOWED_EXT.PDB.value, ALLOWED_EXT.CIF.value}")
    command = RepairPDB(Database, specific_file, new_dataset_location)
    command.run()


@app.command()
def analyze(mode: AnalysisMode = typer.Argument(..., help="Mode for analysis.", rich_help_panel="Mode"),
            database: Path = typer.Argument(..., help="Path of database to perform analysis on."),
            file: Optional[Path] = typer.Option(..., help=f"Provide path to file of interest for analysis. "
                                                          f"The supported file types "
                                                          f"{ALLOWED_EXT.PDB.value, ALLOWED_EXT.CIF.value}")):
    """
    This command is to analyze a given database of PDB. Supported analysis modes: plddt - calculates plddt values of
    homology structures made from AF2; STATS - provides stats on database provided such as number of crystal and
    homology structures and average values per uniProtID..
    """
    if file is not None and file.suffix not in {ALLOWED_EXT.CIF.value, ALLOWED_EXT.PDB.value}:
        raise FileExistsError(f"File type provided {file.suffix} isn't supported. Supported types are: "
                              f"{ALLOWED_EXT.PDB.value, ALLOWED_EXT.CIF.value}")
    command = Analyze(Database, file, mode)
    command.run()




@app.callback()
def main(database: Path = typer.Argument(..., help="Path of database of interest for different modes."),
         number_of_threads: int = typer.Option(os.cpu_count(),
                                               prompt="Specify the number of threads for this process?",
                                               help="Option specifying number of threads to use for this process. "
                                                    "Default is number of logical cores.")):
    global Database
    Database = os.path.abspath(Path(database))

    config = {
        "database_location": Database,
        THREAD_COUNT: number_of_threads
    }
    with open(APP_DIRECTORY+"/config.json", "w") as json_config:
        json.dump(config, json_config)

if __name__ == '__main__':
    app()
