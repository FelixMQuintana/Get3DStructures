from enum import Enum

ACCESSIONS_LOOKUP_TABLE = "accession_ids.csv"

ALPHA_FOLD_EXT = "-model_v%s.pdb"

COLABFOLD_WORKING_DIRECTORY = "/home/felix/Software/colabfold_batch/colabfold-conda/bin/colabfold_batch"


class SUPPORTED_MODES(Enum):
    STRUCTURE = "structure"
    ANALYZE_DATA = "analyze-dataset"

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
