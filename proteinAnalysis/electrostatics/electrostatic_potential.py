from pathlib import Path
from subprocess import PIPE, Popen

import numpy as np

from dataStructure.protein.structure.representation.surface.surface import MolecularSurfaceRepresentation
from proteinAnalysis.proteinAnalysis import ProteinAnalysisOperations


class ABPS(ProteinAnalysisOperations):

    @staticmethod
    def registered_methods() -> dict:
        return {MolecularSurfaceRepresentation.__name__: "surface"}

    def surface(self):

        """
        computeAPBS.py: Wrapper function to compute the Poisson Boltzmann electrostatics for a surface using APBS.
        Pablo Gainza - LPDI STI EPFL 2019
        This file is part of MaSIF.
        Released under an Apache License 2.0
        """
        vertices = self._structure_representation.vertices
        pdb_file = Path(self._structure_representation.file_name)
        # def computeAPBS(vertices, pdb_file, tmp_file_base):

        #   fields = tmp_file_base.split("/")[0:-1]
        #  directory = "/".join(fields) + "/"
        # filename_base = tmp_file_base.split("/")[-1]
        pdbname = pdb_file.name.strip(".pdb")
        args = [
            "pdb2pqr",
            "--ff=PARSE",
            "--whitespace",
            "--noopt",
            "--apbs-input",
            pdbname,
            pdb_file.name,
            pdb_file.name.strip(".pdb") + ".pqr"
            #  filename_base,
        ]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(pdb_file.parent))
        stdout, stderr = p2.communicate()

        args = ["apbs", pdbname]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(pdb_file.parent))
        stdout, stderr = p2.communicate()

        vertfile = open(str(pdb_file.parent) + "/" + pdbname + ".csv", "w")
        for vert in vertices:
            vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))
        vertfile.close()

        args = [
            "multivalue",
            pdbname + ".csv",
            pdbname + ".pqr.dx",
            pdbname + "_out.csv",
        ]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(pdb_file.parent))
        stdout, stderr = p2.communicate()

        # Read the charge file
        chargefile = open(str(pdb_file.parent) + "/" + pdbname + "_out.csv")
        charges = np.array([0.0] * len(vertices))
        for ix, line in enumerate(chargefile.readlines()):
            charges[ix] = float(line.split(",")[3])

        # remove_fn = os.path.join(directory, filename_base)
        # os.remove(remove_fn)
        # os.remove(remove_fn + '.csv')
        # os.remove(remove_fn + '.dx')
        # os.remove(remove_fn + '.in')
        # os.remove(remove_fn + '-input.p')
        # os.remove(remove_fn + '_out.csv')

        return charges