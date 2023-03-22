from enum import Enum

ACCESSIONS_LOOKUP_TABLE = "accession_ids.csv"

COLABFOLD_WORKING_DIRECTORY = "/home/felix/Software/colabfold_batch/colabfold-conda/bin/colabfold_batch"


class AnalysisMode(Enum):
    PLDDT = "plddt"
    STATS = "stats"


class SUPPORTED_MODES(Enum):
    STRUCTURE = "structure"
    REPAIR_STRUCTURES = "repair-structure"
    ANALYZE_DATA = "analyze-dataset"
    FIND_BINDING_SITES = "find-binding-sites"


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


ALPHA_FOLD_STRUCTURE_EXT = "-model_v%s" + ALLOWED_EXT.CIF.value
ALPHA_FOLD_PAE_EXT = "-predicted_aligned_error_v%s" + ALLOWED_EXT.JSON.value
