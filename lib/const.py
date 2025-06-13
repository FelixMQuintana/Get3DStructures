import abc
import dataclasses

from enum import Enum

import numpy as np
import typer


ACCESSIONS_LOOKUP_TABLE = "accession_ids.csv"
APP_NAME = "StructCompare"
APP_DIRECTORY = typer.get_app_dir(APP_NAME)
CONFIG_PATH = APP_DIRECTORY + "/config.json"
COLABFOLD_WORKING_DIRECTORY = "/home/felix/Data/colabfold_batch/colabfold-conda/bin/colabfold_batch"
LOGFILE = "log.txt"


class METRIC(Enum):
    PDB_ID = "PDB ID"
    ACC_ID = "Accession ID"
    PLDDT = "plddt"
    Annotations = "Annotations"


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
    HOMOLOGY_MODEL = "homology_model"


class SimilarityDistance(Enum):
    GRANTHAM = "grantham"
    BACKBONE = "geometric_backbone"
    LLM = "llm"
    BLOSUM62 = "blosum62"
    PHAT = "phat"
    FULL = "full"


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
  #  EXPERIMENTAL_STRUCTURE = "*.pdb"
    EXPERIMENTAL_STRUCTURE = "[!AF,!sp]*.pdb"
    HOMOLOGY_STRUCTURE = "*.pdb"
    CSV_FILE = "*.csv"
    FASTA_FILE = '*.fasta'
    JSON_FILE = '*.json'
    PLY_FILE = "*.ply"
    ATOM_STRUCTURE = "*.pdb"


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


class CalculateOptions(Enum):
    SIMILARITY_DISTANCE = "sim-dist"
    CLUSTER = "cluster"
    GET_REP = "get-rep"
    GET_MESH = "get-mesh"
    ALIGN_DATA = "align-data"
    REASSIGN_CONTACTS = "reassign-contacts"
    MOLECULAR_PARTITION = "molecular-partition"
    GET_ICP = "get-icp"
    GET_KNN = "get-knn"
    COMPARE = "compare"
    COMPARE_PATCHES="compare-patches"
    CLOSEST_MESH = "closest-mesh"

class Metrics:
    UNIPROT_ID_COUNT = "uniprotID_count"
    CRYSTAL_STRUCTURES = "crystal_structures"
    HOMOLOGY_STRUCTURES = "homology_structures"


class CoordinateType(Enum):
    BACKBONE = "backbone"
    C_ALPHA = "c-alpha"


coordinate_type_conversion_dict = {
    CoordinateType.BACKBONE: ["N", "CA", "C", "O"],
    CoordinateType.C_ALPHA: ["CA"]

}


class Statistics(Enum):
    """"""


class FoldingEnvironment(Enum):
    LOCAL = "local"
    DOCKER = "docker"


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
    k14 = "k14"
    k15 = "k15"


e_coli_k_type = {
    "MNCRE44": ktypes.k_plus,
    "APEC O1": ktypes.k1,
    "CE10": ktypes.k_plus,
    "G749": ktypes.k1,
    "IAI39": ktypes.k1,
    "IHE3034": ktypes.k_plus,
    "APEC IMT5155": ktypes.k1,
    "MCJCHV-1": ktypes.k_plus,
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
    "Ecol_743": ktypes.k_plus,
    "536": ktypes.k15,
    "Ecol_448": ktypes.k_plus,
    "EC958": ktypes.k100,
    "Ecol_745": ktypes.k_plus,
    "MVAST0167": ktypes.k100,
    "NRG 857C": ktypes.k_plus,
    "LF82": ktypes.k_plus,
    "SE15": ktypes.k_plus,
    "JJ1887": ktypes.k14,
    "ED1a": ktypes.k_plus,
    "REL606": ktypes.k_plus,
    "K-15KW01": ktypes.k_plus,
    "042": ktypes.k_plus,
    "NA114": ktypes.k_plus,
    "ZH193": ktypes.k5,
    "ZH063": ktypes.k5,
    "SaT040": ktypes.k5,
    "CD306": ktypes.k5,
    "ABU 83972": ktypes.k5,
    "JJ2434": ktypes.k5
}
strain_names = [
    "MG1655",
    "MNCRE44",
    "APEC O1",
    "CE10",
    "G749",
    "IAI39",
    "IHE3034",
    "APEC IMT5155",
    "MCJCHV-1",
    "NC101",
    "NMEC O18",
    "O18",
    "RS218",
    "S88",
    "SF-088",
    "SF-166",
    "SF-173",
    "SF-468",
    "UM146",
    "UTI89",
    "VR50",
    "ST648",
    "Ecol_743",
    "536",
    "Ecol_448",
    "EC958",
    "Ecol_745",
    "MVAST0167",
    "NRG 857C",
    "LF82",
    "SE15",
    "JJ1887",
    "ED1a",
    "606",
    "K-15KW01",
    "042",
    "NA114",
    "ZH193",
    "ZH063",
    "SaT040",
    "CD306",
    "ABU 83972",
    "JJ2434",
    "UMN026",
    "MS6198",
    "ST648",
    "NMECO18",
    "789",
    "11128",
    "11368",
    "CFSAN027343",
    "E2865",
    "FORC_028",
    "RM8426",
    "2011C-3911",
    "150",
    "180-PT54",
    "1130",
    "28RC1",
    "ATCC43889",
    "EC4115",
    "EDL933",
    "FRIK944",
    "FRIK2069",
    "FRIK2455",
    "JEONG-1266",
    "Sakai",
    "SRCC1675",
    "SS52",
    "TW14359",
    "Xuzhou21",
    "RM13514",
    "2013C-4465",
    "RM13516",
    "042",
    "90-9272",
    "H10407",
    "214-4",
    "90-9269",
    "ATCC 43886",
    "90-9281",
    "103605",
    "2014EL-1345-2",
    "E2348",
    "CB9615",
    "RM12579"
]

phat_matrix = {
    "A": {"A": 5, "R": -6, "N": -2, "D": -5, "C": 1, "Q": -3, "E": -5, "G": 1, "H": -3, "I": 0,
          "L": -1, "K": -7, "M": -1, "F": -1, "P": -3, "S": 2, "T": 0, "W": -4, "Y": -3, "V": 1},
    "R": {"A": -6, "R": 9, "N": -3, "D": -7, "C": -8, "Q": -2, "E": -6, "G": -5, "H": -4, "I": -6,
          "L": -6, "K": -1, "M": -6, "F": -7, "P": -7, "S": -6, "T": -6, "W": -7, "Y": -6, "V": -7},
    "N": {"A": -2, "R": -3, "N": 11, "D": 2, "C": -2, "Q": 2, "E": 0, "G": -1, "H": 4, "I": -3,
          "L": -3, "K": -2, "M": -2, "F": -1, "P": -4, "S": 1, "T": -1, "W": -5, "Y": 2, "V": -3},
    "D": {"A": -5, "R": -7, "N": 2, "D": 12, "C": -7, "Q": 0, "E": 6, "G": -2, "H": -1, "I": -5,
          "L": -5, "K": -5, "M": -5, "F": -5, "P": -5, "S": -4, "T": -5, "W": -7, "Y": -4, "V": -5},
    "C": {"A": 1, "R": -8, "N": -2, "D": -7, "C": 7, "Q": -5, "E": -7, "G": -2, "H": -7, "I": -3,
          "L": -2, "K": -10, "M": -2, "F": 0, "P": -8, "S": 1, "T": -1, "W": -4, "Y": -1, "V": -2},
    "Q": {"A": -3, "R": -2, "N": 2, "D": 0, "C": -5, "Q": 9, "E": 1, "G": -2, "H": 2, "I": -3,
          "L": -3, "K": -1, "M": -1, "F": -2, "P": -3, "S": -1, "T": -3, "W": 1, "Y": 0, "V": -3},
    "E": {"A": -5, "R": -6, "N": 0, "D": 6, "C": -7, "Q": 9, "E": 12, "G": -3, "H": -1, "I": -5,
          "L": -5, "K": -4, "M": -5, "F": -5, "P": -5, "S": -3, "T": -5, "W": -7, "Y": -2, "V": -5},
    "G": {"A": 1, "R": -5, "N": -1, "D": -2, "C": -2, "Q": -2, "E": -3, "G": 9, "H": -4, "I": -2,
          "L": -2, "K": -5, "M": -1, "F": -2, "P": -3, "S": 1, "T": -1, "W": -5, "Y": -3, "V": -2},
    "H": {"A": -3, "R": -4, "N": 4, "D": -1, "C": -7, "Q": 2, "E": -1, "G": -4, "H": 11, "I": -5,
          "L": -4, "K": -5, "M": -4, "F": -2, "P": -6, "S": -2, "T": -4, "W": -3, "Y": 3, "V": -5},
    "I": {"A": 0, "R": -6, "N": -3, "D": -5, "C": -3, "Q": -3, "E": -5, "G": -2, "H": -5, "I": 5,
          "L": 2, "K": -7, "M": 3, "F": 0, "P": -4, "S": -2, "T": -1, "W": -4, "Y": -3, "V": 3},
    "L": {"A": -1, "R": -6, "N": -3, "D": -5, "C": -2, "Q": -3, "E": -5, "G": -2, "H": -4, "I": 2,
          "L": 4, "K": -7, "M": 2, "F": 1, "P": -5, "S": -2, "T": -1, "W": -3, "Y": -2, "V": 1},
    "K": {"A": -7, "R": -1, "N": -2, "D": -5, "C": -10, "Q": -1, "E": -4, "G": -5, "H": -5, "I": -7,
          "L": -7, "K": 5, "M": -6, "F": -7, "P": -4, "S": -5, "T": -6, "W": -8, "Y": -4, "V": -8},
    "M": {"A": -1, "R": -6, "N": -2, "D": -5, "C": -2, "Q": -1, "E": -5, "G": -1, "H": -4, "I": 3,
          "L": 2, "K": -6, "M": 6, "F": 0, "P": -5, "S": -2, "T": 0, "W": -4, "Y": -2, "V": 1},
    "F": {"A": -1, "R": -7, "N": -1, "D": -5, "C": 0, "Q": -2, "E": -5, "G": -2, "H": -2, "I": 0,
          "L": 1, "K": -7, "M": 0, "F": 6, "P": -5, "S": -2, "T": -2, "W": 0, "Y": 4, "V": -1},
    "P": {"A": -3, "R": -7, "N": -4, "D": -5, "C": -8, "Q": -3, "E": -5, "G": -3, "H": -6, "I": -4,
          "L": -5, "K": -4, "M": -5, "F": -5, "P": 13, "S": -3, "T": 4, "W": -6, "Y": -5, "V": -4},
    "S": {"A": 2, "R": -6, "N": 1, "D": -4, "C": 1, "Q": -1, "E": -3, "G": 1, "H": -2, "I": -2,
          "L": -2, "K": -5, "M": -2, "F": -2, "P": -3, "S": 6, "T": 1, "W": -5, "Y": -2, "V": -2},
    "T": {"A": 0, "R": -6, "N": -1, "D": -5, "C": -1, "Q": -3, "E": -5, "G": -1, "H": -4, "I": -1,
          "L": -1, "K": -6, "M": 0, "F": -2, "P": 4, "S": 1, "T": 3, "W": -7, "Y": -3, "V": 0},
    "W": {"A": -4, "R": -7, "N": -5, "D": -7, "C": 4, "Q": 1, "E": -7, "G": -5, "H": -3, "I": -4,
          "L": -3, "K": -8, "M": -4, "F": 0, "P": -6, "S": -5, "T": -7, "W": 11, "Y": 1, "V": -4},
    "Y": {"A": -3, "R": -6, "N": 2, "D": -4, "C": -1, "Q": 0, "E": -2, "G": -3, "H": 3, "I": -3,
          "L": -2, "K": -4, "M": -2, "F": 4, "P": -5, "S": -2, "T": -3, "W": 1, "Y": 11, "V": -3},
    "V": {"A": 1, "R": -7, "N": -3, "D": -5, "C": -2, "Q": -3, "E": -5, "G": -2, "H": -5, "I": 3,
          "L": 1, "K": -8, "M": 1, "F": -1, "P": -4, "S": -2, "T": 0, "W": -4, "Y": -3, "V": 4}
}

# Pablo Gainza - LPDI STI EPFL 2018-2019
# Released under an Apache License 2.0


# radii for atoms in explicit case.
radii = {}
radii["N"] = "1.540000"
radii["N"] = "1.540000"
radii["O"] = "1.400000"
radii["C"] = "1.740000"
radii["H"] = "1.200000"
radii["S"] = "1.800000"
radii["P"] = "1.800000"
radii["Z"] = "1.39"
radii["X"] = "0.770000"  ## Radii of CB or CA in disembodied case.
# This  polar hydrogen's names correspond to that of the program Reduce.
polarHydrogens = {}
polarHydrogens["ALA"] = ["H"]
polarHydrogens["GLY"] = ["H"]
polarHydrogens["SER"] = ["H", "HG"]
polarHydrogens["THR"] = ["H", "HG1"]
polarHydrogens["LEU"] = ["H"]
polarHydrogens["ILE"] = ["H"]
polarHydrogens["VAL"] = ["H"]
polarHydrogens["ASN"] = ["H", "HD21", "HD22"]
polarHydrogens["GLN"] = ["H", "HE21", "HE22"]
polarHydrogens["ARG"] = ["H", "HH11", "HH12", "HH21", "HH22", "HE"]
polarHydrogens["HIS"] = ["H", "HD1", "HE2"]
polarHydrogens["TRP"] = ["H", "HE1"]
polarHydrogens["PHE"] = ["H"]
polarHydrogens["TYR"] = ["H", "HH"]
polarHydrogens["GLU"] = ["H"]
polarHydrogens["ASP"] = ["H"]
polarHydrogens["LYS"] = ["H", "HZ1", "HZ2", "HZ3"]
polarHydrogens["PRO"] = []
polarHydrogens["CYS"] = ["H"]
polarHydrogens["MET"] = ["H"]

hbond_std_dev = np.pi / 3

# Dictionary from an acceptor atom to its directly bonded atom on which to
# compute the angle.
acceptorAngleAtom = {}
acceptorAngleAtom["O"] = "C"
acceptorAngleAtom["O1"] = "C"
acceptorAngleAtom["O2"] = "C"
acceptorAngleAtom["OXT"] = "C"
acceptorAngleAtom["OT1"] = "C"
acceptorAngleAtom["OT2"] = "C"
# Dictionary from acceptor atom to a third atom on which to compute the plane.
acceptorPlaneAtom = {}
acceptorPlaneAtom["O"] = "CA"
# Dictionary from an H atom to its donor atom.
donorAtom = {}
donorAtom["H"] = "N"
# Hydrogen bond information.
# ARG
# ARG NHX
# Angle: NH1, HH1X, point and NH2, HH2X, point 180 degrees.
# radii from HH: radii[H]
# ARG NE
# Angle: ~ 120 NE, HE, point, 180 degrees
donorAtom["HH11"] = "NH1"
donorAtom["HH12"] = "NH1"
donorAtom["HH21"] = "NH2"
donorAtom["HH22"] = "NH2"
donorAtom["HE"] = "NE"

# ASN
# Angle ND2,HD2X: 180
# Plane: CG,ND2,OD1
# Angle CG-OD1-X: 120
donorAtom["HD21"] = "ND2"
donorAtom["HD22"] = "ND2"
# ASN Acceptor
acceptorAngleAtom["OD1"] = "CG"
acceptorPlaneAtom["OD1"] = "CB"

# ASP
# Plane: CB-CG-OD1
# Angle CG-ODX-point: 120
acceptorAngleAtom["OD2"] = "CG"
acceptorPlaneAtom["OD2"] = "CB"

# GLU
# PLANE: CD-OE1-OE2
# ANGLE: CD-OEX: 120
# GLN
# PLANE: CD-OE1-NE2
# Angle NE2,HE2X: 180
# ANGLE: CD-OE1: 120
donorAtom["HE21"] = "NE2"
donorAtom["HE22"] = "NE2"
acceptorAngleAtom["OE1"] = "CD"
acceptorAngleAtom["OE2"] = "CD"
acceptorPlaneAtom["OE1"] = "CG"
acceptorPlaneAtom["OE2"] = "CG"

# HIS Acceptors: ND1, NE2
# Plane ND1-CE1-NE2
# Angle: ND1-CE1 : 125.5
# Angle: NE2-CE1 : 125.5
acceptorAngleAtom["ND1"] = "CE1"
acceptorAngleAtom["NE2"] = "CE1"
acceptorPlaneAtom["ND1"] = "NE2"
acceptorPlaneAtom["NE2"] = "ND1"

# HIS Donors: ND1, NE2
# Angle ND1-HD1 : 180
# Angle NE2-HE2 : 180
donorAtom["HD1"] = "ND1"
donorAtom["HE2"] = "NE2"

# TRP Donor: NE1-HE1
# Angle NE1-HE1 : 180
donorAtom["HE1"] = "NE1"

# LYS Donor NZ-HZX
# Angle NZ-HZX : 180
donorAtom["HZ1"] = "NZ"
donorAtom["HZ2"] = "NZ"
donorAtom["HZ3"] = "NZ"

# TYR acceptor OH
# Plane: CE1-CZ-OH
# Angle: CZ-OH 120
acceptorAngleAtom["OH"] = "CZ"
acceptorPlaneAtom["OH"] = "CE1"

# TYR donor: OH-HH
# Angle: OH-HH 180
donorAtom["HH"] = "OH"
acceptorPlaneAtom["OH"] = "CE1"

# SER acceptor:
# Angle CB-OG-X: 120
acceptorAngleAtom["OG"] = "CB"

# SER donor:
# Angle: OG-HG-X: 180
donorAtom["HG"] = "OG"

# THR acceptor:
# Angle: CB-OG1-X: 120
acceptorAngleAtom["OG1"] = "CB"

# THR donor:
# Angle: OG1-HG1-X: 180
donorAtom["HG1"] = "OG1"

import tempfile

masif_opts = {}
# Default directories
masif_opts["raw_pdb_dir"] = "data_preparation/00-raw_pdbs/"
masif_opts["pdb_chain_dir"] = "data_preparation/01-benchmark_pdbs/"
masif_opts["ply_chain_dir"] = "data_preparation/01-benchmark_surfaces/"
masif_opts["tmp_dir"] = tempfile.gettempdir()
masif_opts["ply_file_template"] = masif_opts["ply_chain_dir"] + "/{}_{}.ply"

# Surface features
masif_opts["use_hbond"] = True
masif_opts["use_hphob"] = True
masif_opts["use_apbs"] = True
masif_opts["compute_iface"] = True
# Mesh resolution. Everything gets very slow if it is lower than 1.0
masif_opts["mesh_res"] = 1.0
masif_opts["feature_interpolation"] = True

# Coords params
masif_opts["radius"] = 12.0

# Neural network patch application specific parameters.
masif_opts["ppi_search"] = {}
masif_opts["ppi_search"]["training_list"] = "lists/training.txt"
masif_opts["ppi_search"]["testing_list"] = "lists/testing.txt"
masif_opts["ppi_search"]["max_shape_size"] = 200
masif_opts["ppi_search"]["max_distance"] = 12.0  # Radius for the neural network.
masif_opts["ppi_search"][
    "masif_precomputation_dir"
] = "data_preparation/04b-precomputation_12A/precomputation/"
masif_opts["ppi_search"]["feat_mask"] = [1.0] * 5
masif_opts["ppi_search"]["max_sc_filt"] = 1.0
masif_opts["ppi_search"]["min_sc_filt"] = 0.5
masif_opts["ppi_search"]["pos_surf_accept_probability"] = 1.0
masif_opts["ppi_search"]["pos_interface_cutoff"] = 1.0
masif_opts["ppi_search"]["range_val_samples"] = 0.9  # 0.9 to 1.0
masif_opts["ppi_search"]["cache_dir"] = "nn_models/sc05/cache/"
masif_opts["ppi_search"]["model_dir"] = "nn_models/sc05/all_feat/model_data/"
masif_opts["ppi_search"]["desc_dir"] = "descriptors/sc05/all_feat/"
masif_opts["ppi_search"]["gif_descriptors_out"] = "gif_descriptors/"
# Parameters for shape complementarity calculations.
masif_opts["ppi_search"]["sc_radius"] = 12.0
masif_opts["ppi_search"]["sc_interaction_cutoff"] = 1.5
masif_opts["ppi_search"]["sc_w"] = 0.25

# Neural network patch application specific parameters.
masif_opts["site"] = {}
masif_opts["site"]["training_list"] = "lists/training.txt"
masif_opts["site"]["testing_list"] = "lists/testing.txt"
masif_opts["site"]["max_shape_size"] = 100
masif_opts["site"]["n_conv_layers"] = 3
masif_opts["site"]["max_distance"] = 9.0  # Radius for the neural network.
masif_opts["site"][
    "masif_precomputation_dir"
] = "data_preparation/04a-precomputation_9A/precomputation/"
masif_opts["site"]["range_val_samples"] = 0.9  # 0.9 to 1.0
masif_opts["site"]["model_dir"] = "nn_models/all_feat_3l/model_data/"
masif_opts["site"]["out_pred_dir"] = "output/all_feat_3l/pred_data/"
masif_opts["site"]["out_surf_dir"] = "output/all_feat_3l/pred_surfaces/"
masif_opts["site"]["feat_mask"] = [1.0] * 5

# Neural network ligand application specific parameters.
masif_opts["ligand"] = {}
masif_opts["ligand"]["assembly_dir"] = "data_preparation/00b-pdbs_assembly"
masif_opts["ligand"]["ligand_coords_dir"] = "data_preparation/00c-ligand_coords"
masif_opts["ligand"][
    "masif_precomputation_dir"
] = "data_preparation/04a-precomputation_12A/precomputation/"
masif_opts["ligand"]["max_shape_size"] = 200
masif_opts["ligand"]["feat_mask"] = [1.0] * 5
masif_opts["ligand"]["train_fract"] = 0.9 * 0.8
masif_opts["ligand"]["val_fract"] = 0.1 * 0.8
masif_opts["ligand"]["test_fract"] = 0.2
masif_opts["ligand"]["tfrecords_dir"] = "data_preparation/tfrecords"
masif_opts["ligand"]["max_distance"] = 12.0
masif_opts["ligand"]["n_classes"] = 7
masif_opts["ligand"]["feat_mask"] = [1.0, 1.0, 1.0, 1.0, 1.0]
masif_opts["ligand"]["costfun"] = "dprime"
masif_opts["ligand"]["model_dir"] = "nn_models/all_feat/"
masif_opts["ligand"]["test_set_out_dir"] = "test_set_predictions/"

kd_scale = {}
kd_scale["ILE"] = 4.5
kd_scale["VAL"] = 4.2
kd_scale["LEU"] = 3.8
kd_scale["PHE"] = 2.8
kd_scale["CYS"] = 2.5
kd_scale["MET"] = 1.9
kd_scale["ALA"] = 1.8
kd_scale["GLY"] = -0.4
kd_scale["THR"] = -0.7
kd_scale["SER"] = -0.8
kd_scale["TRP"] = -0.9
kd_scale["TYR"] = -1.3
kd_scale["PRO"] = -1.6
kd_scale["HIS"] = -3.2
kd_scale["GLU"] = -3.5
kd_scale["GLN"] = -3.5
kd_scale["ASP"] = -3.5
kd_scale["ASN"] = -3.5
kd_scale["LYS"] = -3.9
kd_scale["ARG"] = -4.5
