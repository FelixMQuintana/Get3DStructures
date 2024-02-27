import abc
import dataclasses
import enum
import logging
from enum import Enum

import numpy
import numpy as np
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


class SimilarityDistance(Enum):
    GRANTHAM = "grantham"
    BACKBONE = "backbone"
    LLM = "llm"

class ClusterAnalysisType(Enum):
    RAND_INDEX = "rand-index"
    DUNN_INDEX = "dunn-index"

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
    EXPERIMENTAL_STRUCTURE = "^[!AF,!sp]*.pdb"
    HOMOLOGY_STRUCTURE = "*ptm_model*.pdb"
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


class AminoAcidRep:
    """

    """


@dataclasses.dataclass
class ThreeLetterRep(AminoAcidRep):
    value: str


@dataclasses.dataclass
class SingleLetterRep(AminoAcidRep):
    value: str


class AminoAcid(abc.ABC):
    """

    """
    single_letter_rep: str
    three_letter_rep: str

    def __init__(self, rep_value: AminoAcidRep):
        self._rep = rep_value

    @property
    def rep(self) -> AminoAcidRep:
        return self._rep

    @classmethod
    def check(cls, potential_code: str):
        if potential_code == cls.single_letter_rep:
            return cls(SingleLetterRep(cls.single_letter_rep))
        elif potential_code == cls.three_letter_rep:
            return cls(ThreeLetterRep(cls.three_letter_rep))
        else:
            return None


class ALANINE(AminoAcid):
    single_letter_rep = "A"
    three_letter_rep = "ALA"


class ARGININE(AminoAcid):
    single_letter_rep = "R"
    three_letter_rep = "ARG"


class ASPARAGINE(AminoAcid):
    single_letter_rep = "N"
    three_letter_rep = "ASN"


class ASPARTATE(AminoAcid):
    single_letter_rep = "D"
    three_letter_rep = "ASP"


class CYSTEINE(AminoAcid):
    single_letter_rep = "C"
    three_letter_rep = "CYS"


class GLUTAMATE(AminoAcid):
    single_letter_rep = "E"
    three_letter_rep = "GLU"


class GLUTAMINE(AminoAcid):
    single_letter_rep = "Q"
    three_letter_rep = "GLN"


class GLYCINE(AminoAcid):
    single_letter_rep = "G"
    three_letter_rep = "GLY"


class HISTIDINE(AminoAcid):
    single_letter_rep = "H"
    three_letter_rep = "HIS"


class ISOLEUCINE(AminoAcid):
    single_letter_rep = "I"
    three_letter_rep = "ILE"


class LEUCINE(AminoAcid):
    single_letter_rep = "L"
    three_letter_rep = "LEU"


class LYSINE(AminoAcid):
    single_letter_rep = "K"
    three_letter_rep = "LYS"


class METHIONINE(AminoAcid):
    single_letter_rep = "M"
    three_letter_rep = "MET"


class PHENYLALANINE(AminoAcid):
    single_letter_rep = "F"
    three_letter_rep = "PHE"


class PROLINE(AminoAcid):
    single_letter_rep = "P"
    three_letter_rep = "PRO"


class SERINE(AminoAcid):
    single_letter_rep = "S"
    three_letter_rep = "SER"


class THREONINE(AminoAcid):
    single_letter_rep = "T"
    three_letter_rep = "THR"


class TRYPTOPHAN(AminoAcid):
    single_letter_rep = "W"
    three_letter_rep = "TRP"


class TYROSINE(AminoAcid):
    single_letter_rep = "Y"
    three_letter_rep = "TYR"


class VALINE(AminoAcid):
    single_letter_rep = "V"
    three_letter_rep = "VAL"


class AminoAcids(Enum):
    """

    """
    ALANINE = ALANINE
    ARGININE = ARGININE
    ASPARAGINE = ASPARAGINE
    ASPARTATE = ASPARTATE
    CYSTEINE = CYSTEINE
    GLUTAMATE = GLUTAMATE
    GLUTAMINE = GLUTAMINE
    GLYCINE = GLYCINE
    HISTIDINE = HISTIDINE
    ISOLEUCINE = ISOLEUCINE
    LEUCINE = LEUCINE
    LYSINE = LYSINE
    METHIONINE = METHIONINE
    PHENYLALANINE = PHENYLALANINE
    PROLINE = PROLINE
    SERINE = SERINE
    THREONINE = THREONINE
    TRYPTOPHAN = TRYPTOPHAN
    TYROSINE = TYROSINE
    VALINE = VALINE

    @classmethod
    def get_rep(cls, amino_acid_string: str):
        # logging.debug("Finding type of string")
        rep = []
        for letter in amino_acid_string:
            for amino_acid in cls.__members__.values():
                result = amino_acid.value.check(letter)
                if result:
                    rep.append(result)
                    break
        return rep

    #    amino_acid_string[0]


grantham_distance_matrix_row_dict = {
    SERINE: 0,
    ARGININE: 1,
    LEUCINE: 2,
    PROLINE: 3,
    THREONINE: 4,
    ALANINE: 5,
    VALINE: 6,
    GLYCINE: 7,
    ISOLEUCINE: 8,
    PHENYLALANINE: 9,
    TYROSINE: 10,
    CYSTEINE: 11,
    HISTIDINE: 12,
    GLUTAMINE: 13,
    ASPARAGINE: 14,
    LYSINE: 15,
    ASPARTATE: 16,
    GLUTAMATE: 17,
    METHIONINE: 18,
    TRYPTOPHAN: 19
}
grantham_distance_matrix_column_dict = {
    SERINE: 0,
    ARGININE: 1,
    LEUCINE: 2,
    "PRO": 3,
    "THR": 4,
    "ALA": 5,
    "VAL": 6,
    "GLY": 7,
    "ILE": 8,
    PHENYLALANINE: 9,
    TYROSINE: 10,
    "CYS": 11,
    "HIS": 12,
    "GLN": 13,
    "ASN": 14,
    "LYS": 15,
    "ASP": 16,
    "GLU": 17,
    "MET": 18,
    "TRP": 19
}

grantham_distance_matrix = np.array([
    # SER # ARG # LEU # PRO # THR # ALA # VAL # GLY # ILE # PHE # TYR # CYS # HIS # GLN # ASN # LYS # ASP # GLU # MET # TRP
    np.array(
        [1.00, 0.49, 0.33, 0.66, 0.73, 0.54, 0.42, 0.74, 0.34, 0.28, 0.33, 0.48, 0.59, 0.68, 0.79, 0.44, 0.70, 0.63,
         0.37, 0.18]),  # SER
    np.array(
        [0.49, 1.00, 0.53, 0.52, 0.67, 0.48, 0.55, 0.42, 0.55, 0.55, 0.64, 0.16, 0.87, 0.80, 0.60, 0.88, 0.55, 0.75,
         0.58, 0.53]),  # ARG
    np.array(
        [0.33, 0.53, 1.00, 0.54, 0.57, 0.55, 0.85, 0.36, 0.98, 0.90, 0.83, 0.08, 0.54, 0.47, 0.29, 0.50, 0.20, 0.36,
         0.93, 0.72]),  # LEU
    np.array(
        [0.66, 0.52, 0.54, 1.00, 0.82, 0.87, 0.68, 0.80, 0.56, 0.47, 0.49, 0.21, 0.64, 0.65, 0.58, 0.52, 0.50, 0.57,
         0.60, 0.32]),  # PRO
    np.array(
        [0.73, 0.67, 0.57, 0.82, 1.00, 0.73, 0.68, 0.73, 0.59, 0.52, 0.57, 0.31, 0.78, 0.80, 0.70, 0.64, 0.60, 0.70,
         0.62, 0.40]),  # THR
    np.array(
        [0.54, 0.48, 0.55, 0.87, 0.73, 1.00, 0.70, 0.72, 0.56, 0.47, 0.48, 0.09, 0.60, 0.58, 0.48, 0.51, 0.41, 0.50,
         0.61, 0.31]),  # ALA
    np.array(
        [0.42, 0.55, 0.85, 0.68, 0.68, 0.70, 1.00, 0.49, 0.87, 0.77, 0.74, 0.11, 0.61, 0.55, 0.38, 0.55, 0.29, 0.44,
         0.90, 0.59]),  # VAL
    np.array(
        [0.74, 0.42, 0.36, 0.80, 0.73, 0.72, 0.49, 1.00, 0.37, 0.29, 0.32, 0.26, 0.54, 0.60, 0.63, 0.41, 0.56, 0.54,
         0.41, 0.14]),  # GLY
    np.array(
        [0.34, 0.55, 0.98, 0.56, 0.59, 0.56, 0.87, 0.37, 1.00, 0.90, 0.85, 0.08, 0.56, 0.49, 0.31, 0.53, 0.22, 0.38,
         0.95, 0.72]),  # ILE
    np.array(
        [0.28, 0.55, 0.90, 0.47, 0.52, 0.47, 0.77, 0.29, 0.90, 1.00, 0.90, 0.05, 0.53, 0.46, 0.27, 0.53, 0.18, 0.35,
         0.87, 0.81]),  # PHE
    np.array(
        [0.33, 0.64, 0.83, 0.49, 0.57, 0.48, 0.74, 0.32, 0.85, 0.90, 1.00, 0.10, 0.61, 0.54, 0.33, 0.60, 0.26, 0.43,
         0.83, 0.83]),  # TYR
    np.array(
        [0.48, 0.16, 0.08, 0.21, 0.31, 0.09, 0.11, 0.26, 0.08, 0.05, 0.10, 1.00, 0.19, 0.28, 0.35, 0.06, 0.28, 0.21,
         0.09, 0.00]),  # CYS
    np.array(
        [0.59, 0.87, 0.54, 0.64, 0.78, 0.60, 0.61, 0.54, 0.56, 0.53, 0.61, 0.19, 1.00, 0.89, 0.68, 0.85, 0.62, 0.81,
         0.60, 0.47]),  # HIS
    np.array(
        [0.68, 0.80, 0.47, 0.65, 0.80, 0.58, 0.55, 0.60, 0.49, 0.46, 0.54, 0.28, 0.89, 1.00, 0.79, 0.75, 0.72, 0.87,
         0.53, 0.40]),  # GLN
    np.array(
        [0.79, 0.60, 0.29, 0.58, 0.70, 0.48, 0.38, 0.63, 0.31, 0.27, 0.33, 0.35, 0.68, 0.79, 1.00, 0.56, 0.89, 0.80,
         0.34, 0.19]),  # ASN
    np.array(
        [0.44, 0.88, 0.50, 0.52, 0.64, 0.51, 0.55, 0.41, 0.53, 0.53, 0.60, 0.06, 0.85, 0.75, 0.56, 1.00, 0.53, 0.74,
         0.56, 0.49]),  # LYS
    np.array(
        [0.70, 0.55, 0.20, 0.50, 0.60, 0.41, 0.29, 0.56, 0.22, 0.18, 0.26, 0.28, 0.62, 0.72, 0.89, 0.53, 1.00, 0.79,
         0.26, 0.16]),  # ASP
    np.array(
        [0.63, 0.75, 0.36, 0.57, 0.70, 0.50, 0.44, 0.54, 0.38, 0.35, 0.43, 0.21, 0.81, 0.87, 0.80, 0.74, 0.79, 1.00,
         0.41, 0.29]),  # GLU
    np.array(
        [0.37, 0.58, 0.93, 0.60, 0.62, 0.61, 0.90, 0.41, 0.95, 0.87, 0.83, 0.09, 0.60, 0.53, 0.34, 0.56, 0.26, 0.41,
         1.00, 0.69]),  # MET
    np.array(
        [0.18, 0.53, 0.72, 0.32, 0.40, 0.31, 0.59, 0.14, 0.72, 0.81, 0.83, 0.00, 0.47, 0.40, 0.19, 0.49, 0.16, 0.29,
         0.69, 1.00])  # TRP
])

class ktypes(Enum):
    k_plus = "k+",
    k1 = "k1",
    k2 = "k2",
    k100 = "k100",
    k52 = "k52",
    k5 = "k5",
    k14="k14"
    k15="k15"

e_coli_k_type = {
    "MNCRE44": ktypes.k_plus,
    "APEC O1": ktypes.k1,
    "CE10": ktypes.k_plus,
    "G749": ktypes.k1,
    "IAI39": ktypes.k1,
    "IHE3034": ktypes.k_plus,
    "APEC IMT5155": ktypes.k1,
    "MCJCHV-1":ktypes.k_plus,
    "NC101": ktypes.k1,
    "NMEC O18": ktypes.k1,
    "O18": ktypes.k1,
    "RS218": ktypes.k1,
    "S88": ktypes.k1,
    "SF-088": ktypes.k1,
    "SF-166": ktypes.k1,
    "SF-173": ktypes.k1,
    "SF-468": ktypes.k1,
    "UM146": ktypes.k_plus,
    "UTI89": ktypes.k1,
    "VR50": ktypes.k1,
    "ST648": ktypes.k_plus,
    "Ecol_743":ktypes.k_plus,
    "536": ktypes.k15,
    "Ecol_448":ktypes.k_plus,
    "EC958":ktypes.k100,
    "Ecol_745":ktypes.k_plus,
    "MVAST0167":ktypes.k100,
    "NRG 857C":ktypes.k_plus,
    "LF82":ktypes.k_plus,
    "SE15":ktypes.k_plus,
    "JJ1887":ktypes.k14,
    "ED1a": ktypes.k_plus,
    "REL606":ktypes.k_plus,
    "K-15KW01":ktypes.k_plus,
    "042":ktypes.k_plus,
    "NA114":ktypes.k_plus,
    "ZH193":ktypes.k5,
    "ZH063":ktypes.k5,
    "SaT040": ktypes.k5,
    "CD306":ktypes.k5,
    "ABU 83972":ktypes.k5,
    "JJ2434": ktypes.k5
}