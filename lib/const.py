from enum import Enum

import typer

ACCESSIONS_LOOKUP_TABLE = "accession_ids.csv"
APP_NAME = "StructCompare"
APP_DIRECTORY = typer.get_app_dir(APP_NAME)
CONFIG_PATH = APP_DIRECTORY + "/config.json"
COLABFOLD_WORKING_DIRECTORY = "/home/felix/Data/colabfold_batch/colabfold-conda/bin/colabfold_batch"
LOGFILE = "log.txt"


class AnalysisMode(Enum):
    PLDDT = "plddt"
    STATS = "stats"


class ConfigOptions(Enum):
    DATABASE_LOCATION = "database_location"
    THREAD_COUNT = "thread_count"
    STRUCTURE_TYPE = "structure_type"


class StructureCharacteristicsMode(Enum):
    AUTOSITE = "Autosite"
    BUILD_MOTIFS = "build_motif"
    FIND_BINDING_POCKETS = "find_binding_pocket"
    CHECK_MOTIF_QUALITY = "check_motif_quality"
    CALCULATE_RMSD = "calculate_rmsd"
    TRIM_PDB = "trim_pdb"
    GET_LINEAGE = "get_lineage"
    CREATE_COMBINATORICS = "create_combinatorics"
    # LIGANDBINDING


class StructureBuildMode(Enum):
    ACCESSION_DATA = "accession_data"
    FASTA_DATA = "fasta_data"
    PDB_STRUCTURES = "pdb_structures"
    ALPHAFOLD_STRUCTURES = "alphafold_structures"
    COMPUTE_ALPHA_STRUCTURES = "compute_alphafold_structures"
    HomologyModel = "homology_model"


class MotifSearchMode(Enum):
    CLUSTER = "cluster"
    FEATURE_POINTS = "feature_points"


class MotifRefinements(Enum):
    ALL_ATOM = "all_atom"
    RESIDUE = "residue"


class SupportedModes(Enum):
    STRUCTURE = "structure"
    REPAIR_STRUCTURES = "repair-structure"
    ANALYZE_DATA = "analyze-dataset"
    FIND_BINDING_SITES = "find-binding-sites"


class SupportedStructureFileTypes(Enum):
    PDB = "pdb"
    CIF = "cif"


class SupportedFileTypeRegex(Enum):
    EXPERIMENTAL_STRUCTURE = "[!AF,!sp]*.pdb"
    HOMOLOGY_STRUCTURE = "[AF,sp]*.pdb"
    CSV_FILE = "*.csv"
    FASTA_FILE = '*.fasta'
    JSON_FILE = '*.json'


class ColabFoldOptions(Enum):
    MSA_MODE = " --msa-mode %s "
    USE_GPU = " --use-gpu-relax "
    NUM_ENSEMBLE = " --num-ensemble %s"
    AMBER = " --amber "
    PAIR_MODE = " --pair-mode %s "


class COLABFOLDResponses(Enum):
    UNPAIRED_PAIRED = "unpaired_paired"
    MMSEQS2_UNIREF_ENV = "mmseqs2_uniref_env"


class AllowedExt(Enum):
    CIF = ".cif"
    PDB = ".pdb"
    JSON = ".json"
    FASTA = ".fasta"


class FunSoCs(Enum):
    ANTIBIOTIC_RESISTANCE = "AntibioticResistance"
    BACTERIAL_COUNTER_SIGNALING = "BacterialCounterSignaling"
    COUNTER_IMMUNOGLOBULIN = "CounterImmunoglobulin"
    CYTOTOXICITY = "Cytotoxicity"
    DEGRADE_ECM = "DegradeECM"
    DEVELOPMENT_IN_HOST = "DevelopmentInHost"
    DISABLE_ORGAN = "DisableOrgan"


class Metrics:
    UNIPROTID_COUNT = "uniprotID_count"
    CRYSTAL_STRUCTURES = "crystal_structures"
    HOMOLOGY_STRUCTURES = "homology_structures"


class Statistics(Enum):
    """"""


class CountStatistics(Statistics):
    COUNT = "count"
    MEAN = "count_mean"
    ST_DEV = "count_stdev"
    MAX = "count_max"
    MIN = "count_min"
    MODE = "count_mode"


class HomologyStructureStatistics(Statistics):
    PLDDT_MEAN = "plddt_mean"


class SequenceLengthStatistics(Statistics):
    MEAN = "sequence_length_mean"
    ST_DEV = "sequence_length_stdev"
    MAX = "sequence_length_max"
    MIN = "sequence_length_min"


ALPHA_FOLD_STRUCTURE_EXT = "-model_v%s"
ALPHA_FOLD_PAE_EXT = "-predicted_aligned_error_v%s" + AllowedExt.JSON.value
