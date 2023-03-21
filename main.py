import argparse
from argparse import Namespace
import os.path
import re
from typing import List, Dict

from Commands.Analyze import Analyze
from Commands.Structure import Structure
from Commands.command import Command
from Commands.repair import RepairPDB
from lib.const import SUPPORTED_MODES, ALLOWED_EXT
from query import UniProtIDQuery
import json

# DEFAULT_STRUCTURE_PATH = "/Users/felixquintana/Downloads/FunSoCTrainingData/"
ALLOWED_REGEX = "[a-zA-Z0-9]{6}"


def verify_uniprot_ids_format(possible_uniprot_id: str) -> None:
    """

    :param possible_uniprot_id:
    :return:
    """
    if not re.search(ALLOWED_REGEX, possible_uniprot_id):
        raise Exception(f"Input {possible_uniprot_id}, does not follow expected regex.")


def digest_uniprot_ids(input_file: str) -> List[str]:
    """

    :param input_file:
    :return:
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"File doesn't seem to exist: {input_file}")
    with open(input_file, "r") as file:
        possible_uniprot_ids: List[str] = file.readlines()
        possible_uniprot_ids = [id.strip("\n") for id in possible_uniprot_ids]
        for possible_uniprot_id in possible_uniprot_ids:
            verify_uniprot_ids_format(possible_uniprot_id)
    return possible_uniprot_ids


def save_responses(responses: List[Dict]):
    """

    :param responses:
    :return:
    """
    with open('responses.txt', 'w') as out:
        json.dump(responses, out)


def query_uniprot_id(uniprot_id: str) -> Dict:
    """

    :type uniprot_id: str
    :param uniprot_id:
    :return:
    """

    uniprot = UniProtIDQuery(uniprot_id + ALLOWED_EXT.JSON.value)
    return uniprot.parse_response(uniprot.query(uniprot_id))


def get_args() -> tuple[Namespace, list[str]]:
    parser = argparse.ArgumentParser(prog="GenerateStructures",
                                     description="Program is designed to capture structures from a specified set of "
                                                 "UniProtIDs ")
    # TODO add exlusive groups
    parser.add_argument("--wd", metavar="WORKING-DIRECTORY", required=True, type=str,
                        help="WORKING-DIRECTORY for any operation of the program")
    parser.add_argument('--skip', metavar="SKIP_BOOL", default=False, type=bool)
    subparser = parser.add_subparsers(help="Sub-Commands", dest="mode")
    structures = subparser.add_parser(SUPPORTED_MODES.STRUCTURE.value, help="Called for getting structures.")
    structures_to_check = structures.add_mutually_exclusive_group(required=True)
    structures_to_check.add_argument('--id', metavar="UniProtID",
                                     help="Specify UniProtID to look up protein structures.")
    structures_to_check.add_argument('--file', metavar="UniProtID-List", help="UnitProtID List ")
    structures.set_defaults(mode=SUPPORTED_MODES.STRUCTURE.value)
    analyze_dataset = subparser.add_parser(SUPPORTED_MODES.ANALYZE_DATA.value,
                                           help="This is called when you want an "
                                                "introspective of the dataset generated")
    analyze_dataset.add_argument("--all", type=bool, default=True)
    analyze_dataset.add_argument("--file", metavar="PDB_FILE", type=str, help="Path of specific PDB file you'd like "
                                                                              "to investigate")
    analyze_dataset.add_argument("--mode", type=str, choices=["stats", "plddt"], help="Can be set as ")
    analyze_dataset.add_argument("--stats", type=bool, default=False, help="Returns statistics on dataset specified "
                                                                           "by WORKING_DIRECTORY")
    analyze_dataset.add_argument("--plddt", type=bool, help="Called to analyze PLDDT plots.")
    analyze_dataset.set_defaults(mode=SUPPORTED_MODES.ANALYZE_DATA.value)
    fixer = subparser.add_parser(SUPPORTED_MODES.REPAIR_STRUCTURES.value, help="Repair pdb files")
    fixer.add_argument("--all", type=bool, default=True, help="Set as True by default. Set True if you want all "
                                                              "PDB Files in WORKING_DIRECTORY to be repaired.")
    fixer.add_argument("--file", metavar="PDB_FILE", type=str, help="Called when you want fix a specific PDB_FILE",
                       default=False)
    fixer.add_argument("--ds", metavar="DATASET-LOCATION", type=str, help="Where to save new CIF files.", required=True)
    lbs_finder = subparser.add_parser(SUPPORTED_MODES.FIND_BINDING_SITES.value, help="This is called when you want to "
                                                                        "locate binding site")
    # parser.add_argument('--file', metavar='F', required=False, type=str, help="File of UniProtIDs")
    # parser.add_argument('--database', type=str, required=True, help="Location where PDB structure "
    #                                                                "database should be held. IF folders do "
    #                                                                "not exist for a given accession, this "
    #                                                                "code will make it.")
    # parser.add_argument('--id', metavar='S', type=str, help="UniProtID", )
    return parser.parse_known_args()


class CommandDigest:

    def __init__(self, namespace: Namespace):
        self._args: Namespace = namespace

    def digest(self) -> Command:
        return MODE_OPTIONS[self._args.mode](self._args)

    @staticmethod
    def structure(args: Namespace) -> Command:
        return Structure(args.wd, args.id, args.file, args.skip)

    @staticmethod
    def analyze_data(args: Namespace) -> Command:
        return Analyze(args.wd, args.all, args.file, args.plddt)

    @staticmethod
    def repair(args: Namespace) -> Command:
        return RepairPDB(args.wd, args.all, args.file, args.ds)
    @staticmethod
    def find_characteristics(args: Namespace) -> Command:
        return
MODE_OPTIONS = {SUPPORTED_MODES.STRUCTURE.value: CommandDigest.structure,
                SUPPORTED_MODES.STRUCTURE.ANALYZE_DATA.value: CommandDigest.analyze_data,
                SUPPORTED_MODES.REPAIR_STRUCTURES.value: CommandDigest.repair}

if __name__ == '__main__':
    namespace, extra = get_args()
    digester = CommandDigest(namespace)
    command: Command = digester.digest()
    command.run()
