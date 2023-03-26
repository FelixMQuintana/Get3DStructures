from enum import Enum

import typer

ACCESSIONS_LOOKUP_TABLE = "accession_ids.csv"
APP_NAME = "StructCompare"
APP_DIRECTORY = typer.get_app_dir(APP_NAME)
CONFIG_PATH = APP_DIRECTORY + "/config.json"
COLABFOLD_WORKING_DIRECTORY = "/home/felix/Software/colabfold_batch/colabfold-conda/bin/colabfold_batch"
LOGFILE = "log.txt"


class AnalysisMode(Enum):
    PLDDT = "plddt"
    STATS = "stats"

class CONFIG_OPTIONS(Enum):
    DATABASE_LOCATION = "database_location"
    THREAD_COUNT = "thread_count"
    STRUCTURE_TYPE = "structure_type"

class StructureCharacteristicsMode(Enum):
    AUTOSITE = "Autosite"
    # LIGANDBINDING


class SUPPORTED_MODES(Enum):
    STRUCTURE = "structure"
    REPAIR_STRUCTURES = "repair-structure"
    ANALYZE_DATA = "analyze-dataset"
    FIND_BINDING_SITES = "find-binding-sites"


class SUPPORTED_STRUCTURE_TYPES(Enum):
    PDB = "pdb"
    CIF = "cif"


class COLABFOLD_OPTIONS(Enum):
    MSA_MODE = " --msa-mode %s "
    USE_GPU = " --use-gpu-relax "
    NUM_ENSEMBLE = " --num-ensemble %s"
    AMBER = " --amber "
    PAIR_MODE = " --pair-mode %s "


class COLABFOLDResponses(Enum):
    UNPAIRED_PAIRED = "unpaired_paired"
    MMSEQS2_UNIREF_ENV = "mmseqs2_uniref_env"


class ALLOWED_EXT(Enum):
    CIF = ".cif"
    PDB = ".pdb"
    JSON = ".json"
    FASTA = ".fasta"


class FUNSOCS(Enum):
    ANTIBIOTIC_RESISTANCE = "AntibioticResistance"
    BACTERIAL_COUNTER_SIGNALING = "BacterialCounterSignaling"
    COUNTER_IMMUNOGLOBULIN = "CounterImmunoglobulin"
    CYTOTOXICITY = "Cytotoxicity"
    DEGRADE_ECM = "DegradeECM"
    DEVELOPMENT_IN_HOST = "DevelopmentInHost"
    DISABLE_ORGAN = "DisableOrgan"


ALPHA_FOLD_STRUCTURE_EXT = "-model_v%s"
ALPHA_FOLD_PAE_EXT = "-predicted_aligned_error_v%s" + ALLOWED_EXT.JSON.value
