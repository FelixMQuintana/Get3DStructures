"""

"""
import abc
import os
from abc import ABC
from pathlib import Path
from typing import List, Optional, Type

from dataStructure.protein.structure.representation.representation import Representation
from Bio import SeqIO
from lib.const import AllowedExt


class SupportedFileType:
    """

    """

    def __init__(self, path: Path):
        self._path = path

    @property
    def path(self) -> Path:
        return self._path

    @property
    def id(self):
        return self._path.parent.name


class GOLabelsFile(SupportedFileType):
    """

    """


class MeshFile(SupportedFileType):
    """

    """


class StructureFile(SupportedFileType, abc.ABC):

    def __init__(self, path: Path, ):  # representation_type:Representation):
        super().__init__(path)
        self._path: Path = path
        #   if self.path.stat().st_size == 0:
        #        os.system(f"rm {self.path}")
        self._sequence = ""
        self._binding_site_residues = None
        self._representation = Representation

    @property
    def fasta(self) -> str:
        if len(self._sequence) == 0:
            pdb_parser = SeqIO.parse(self._path, "pdb-atom")
            for record in pdb_parser:
                self._sequence += record.seq
        return self._sequence

    @property
    def representation(self) -> Type[Representation]:
        return self._representation


# @property
# def representation(self) -> StructureRepresentation:
#     if self._representation is None:
#         raise RuntimeError("Representation not set")
#     else:
#         return self._representation

# @representation.setter
# def representation(self, representation: StructureRepresentation):
#    if self._representation is None:
#        self._representation = representation(self._path)


class HomologyStructure(StructureFile, ABC):
    """
        Homology Structure File Type
    """

    #    @property
    #    def representation(self) -> StructureRepresentation:
    #        return Atomic(self.path)

    def __init__(self, path: Path):
        super().__init__(path)
        self._piddt: Optional[List[float]] = None
        self.preprocess()

    @property
    def piddt(self) -> List:
        if self._piddt is None:
            read_alpha_fold_file_obj = open(self.path.with_suffix(AllowedExt.PDB.value), "r")
            plddt = []
            current_res_num = 0
            for line in read_alpha_fold_file_obj:
                if line.startswith("ATOM"):
                    if current_res_num < 1000 and line.split()[5].isnumeric():
                        if int(line.split()[5]) > current_res_num:
                            plddt.append(float(line.split()[-2]))
                            current_res_num = int(line.split()[5])
                    else:
                        #     print(line)
                        if int(line.split()[4][1:]) > current_res_num:
                            plddt.append(float(line.split()[-2]))
                            current_res_num = int(line.split()[4][1:])
            self._piddt = plddt
        return self._piddt

    def preprocess(self):
        read_alpha_fold_file_obj = open(self.path.with_suffix(AllowedExt.PDB.value), "r")
        read_lines = read_alpha_fold_file_obj.readlines()
        read_alpha_fold_file_obj.close()
        write_alpha_fold_file_obj = open(self.path.with_suffix(AllowedExt.PDB.value), "w")

        for line in read_lines:
            if line.startswith("DBREF"):
                continue
            else:
                write_alpha_fold_file_obj.write(line)
        write_alpha_fold_file_obj.close()

    def trim(self, cutoff=70):
        read_alpha_fold_file_obj = open(self.path.with_suffix(AllowedExt.PDB.value), "r")
        read_lines = read_alpha_fold_file_obj.readlines()
        read_alpha_fold_file_obj.close()
        write_alpha_fold_file_obj = open(self.path.with_suffix(AllowedExt.PDB.value), "w")
        # read_alpha_fold_file_obj = open(self.path.with_suffix(AllowedExt.PDB.value), "r")
        plddt = []
        not_middle_of_structure = True
        current_res_num = 0
        canidates = []
        for line in read_lines:
            if line.startswith("ATOM"):
                if current_res_num < 1000 and line.split()[5].isnumeric():
                    if int(line.split()[5]) > current_res_num:  # and not_middle_of_structure:
                        if float(line.split()[-2]) < cutoff:
                            canidates.append(False)
                            continue
                else:
                    if int(line.split()[4][1:]) > current_res_num:
                        if float(line.split()[-2]) < cutoff:  # and not_middle_of_structure:
                            canidates.append(False)
                            continue
                canidates.append(True)
        print(canidates)
        forward_copy = iter(canidates)
        backward_copy = iter(reversed(canidates))
        index = 0
        while not next(forward_copy):
            index += 1
        forward_copy_index = index
        index = 0
        while not next(backward_copy):
            index += 1
        backward_copy_index = index
        for index, entry in enumerate(canidates):
            if forward_copy_index < index < len(canidates)-backward_copy_index:
                canidates[index] = True

        print(canidates)
        index=-1
        for line in read_lines:
            if line.startswith("ATOM"):
                index+=1
                if canidates[index]:
                    write_alpha_fold_file_obj.write(line)
            else:
                write_alpha_fold_file_obj.write(line)
                #if current_res_num < 1000 and line.split()[5].isnumeric():
                #    if int(line.split()[5]) > current_res_num:  # and not_middle_of_structure:
                #        if float(line.split()[-2]) < cutoff:
                #            canidates.append(False)
                #            continue
                #else:
                #    if int(line.split()[4][1:]) > current_res_num:
                #        if float(line.split()[-2]) < cutoff:  # and not_middle_of_structure:
                #            canidates.append(False)
                #            continue
                #canidates.append(True)


        write_alpha_fold_file_obj.close()
            # write_alpha_fold_file_obj.write(line)


class CrystalStructure(StructureFile, ABC):
    """

    """
