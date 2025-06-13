import os
import random
from pathlib import Path
from subprocess import Popen, PIPE

import numpy as np
import pymesh
from numpy.linalg import norm
from numpy.matlib import repmat

from dataStructure.protein.structure.representation.representation import Representation
from lib.func import output_pdb_as_xyzrn, read_msms


def compute_ms(pdb_file, protonate=True):
    randnum = random.randint(1, 10000000)
    file_base = "./msms_" + str(randnum)
    out_xyzrn = file_base + ".xyzrn"
    if protonate:
        output_pdb_as_xyzrn(pdb_file, out_xyzrn)
    else:
        raise RuntimeError("Error - pdb2xyzrn is deprecated.")
        # Now run MSMS on xyzrn file
    FNULL = open(os.devnull, 'w')
    args = ["msms", "-density", "3.0", "-hdensity", "3.0", "-probe", \
            "1.5", "-if", out_xyzrn, "-of", file_base, "-af", file_base]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()
    vertices, faces, normals, names = read_msms(file_base)
    areas = {}
    ses_file = open(file_base + ".area")
    next(ses_file)  # ignore header line
    for line in ses_file:
        fields = line.split()
        areas[fields[3]] = fields[1]

        # Remove temporary files.
    os.remove(file_base + '.area')
    os.remove(file_base + '.xyzrn')
    os.remove(file_base + '.vert')
    os.remove(file_base + '.face')
    return vertices, faces, normals, areas, names


"""
fixmesh.py: Regularize a protein surface mesh. 
- based on code from the PyMESH documentation. 
"""


def fix_mesh(mesh, resolution, detail="normal"):
    bbox_min, bbox_max = mesh.bbox
    diag_len = norm(bbox_max - bbox_min)
    if detail == "normal":
        target_len = diag_len * 5e-3
    elif detail == "high":
        target_len = diag_len * 2.5e-3
    elif detail == "low":
        target_len = diag_len * 1e-2
    target_len = resolution
    # print("Target resolution: {} mm".format(target_len));
    # PGC 2017: Remove duplicated vertices first
    mesh, _ = pymesh.remove_duplicated_vertices(mesh, 0.001)

    count = 0
    print("Removing degenerated triangles")
    mesh, __ = pymesh.remove_degenerated_triangles(mesh, 100)
    mesh, __ = pymesh.split_long_edges(mesh, target_len)
    num_vertices = mesh.num_vertices
    while True:
        mesh, __ = pymesh.collapse_short_edges(mesh, 1e-6)
        mesh, __ = pymesh.collapse_short_edges(mesh, target_len,
                                               preserve_feature=True)
        mesh, __ = pymesh.remove_obtuse_triangles(mesh, 150.0, 100)
        if mesh.num_vertices == num_vertices:
            break
        num_vertices = mesh.num_vertices
        count += 1
        if count > 10: break

    mesh = pymesh.resolve_self_intersection(mesh)
    mesh, __ = pymesh.remove_duplicated_faces(mesh)
    mesh = pymesh.compute_outer_hull(mesh)
    mesh, __ = pymesh.remove_duplicated_faces(mesh)
    mesh, __ = pymesh.remove_obtuse_triangles(mesh, 179.0, 5)
    mesh, __ = pymesh.remove_isolated_vertices(mesh)
    mesh, _ = pymesh.remove_duplicated_vertices(mesh, 0.001)

    return mesh


"""
compute_normal.py: Compute the normals of a closed shape.
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF, based on previous matlab code by Gabriel Peyre, converted to Python by Pablo Gainza
"""

###
eps = 1.0e-6


def compute_normal(vertex, face):
    """
    compute_normal - compute the normal of a triangulation
    vertex: 3xn matrix of vertices
    face: 3xm matrix of face indices.

      normal,normalf = compute_normal(vertex,face)

      normal(i,:) is the normal at vertex i.
      normalf(j,:) is the normal at face j.

    Copyright (c) 2004 Gabriel Peyr
    Converted to Python by Pablo Gainza LPDI EPFL 2017
    """

    vertex = vertex.T
    face = face.T
    nface = np.size(face, 1)
    nvert = np.size(vertex, 1)
    normal = np.zeros((3, nvert))
    # unit normals to the faces
    normalf = crossp(
        vertex[:, face[1, :]] - vertex[:, face[0, :]],
        vertex[:, face[2, :]] - vertex[:, face[0, :]],
    )
    sum_squares = np.sum(normalf ** 2, 0)
    d = np.sqrt(sum_squares)
    d[d < eps] = 1
    normalf = normalf / repmat(d, 3, 1)
    # unit normal to the vertex
    normal = np.zeros((3, nvert))
    for i in np.arange(0, nface):
        f = face[:, i]
        for j in np.arange(3):
            normal[:, f[j]] = normal[:, f[j]] + normalf[:, i]

    # normalize
    d = np.sqrt(np.sum(normal ** 2, 0))
    d[d < eps] = 1
    normal = normal / repmat(d, 3, 1)
    # enforce that the normal are outward
    vertex_means = np.mean(vertex, 0)
    v = vertex - repmat(vertex_means, 3, 1)
    s = np.sum(np.multiply(v, normal), 1)
    if np.sum(s > 0) < np.sum(s < 0):
        # flip
        normal = -normal
        normalf = -normalf
    return normal.T


def crossp(x, y):
    # x and y are (m,3) dimensional
    z = np.zeros((x.shape))
    z[0, :] = np.multiply(x[1, :], y[2, :]) - np.multiply(x[2, :], y[1, :])
    z[1, :] = np.multiply(x[2, :], y[0, :]) - np.multiply(x[0, :], y[2, :])
    z[2, :] = np.multiply(x[0, :], y[1, :]) - np.multiply(x[1, :], y[0, :])
    return z


class MolecularSurfaceRepresentation(Representation):
    """

    """

    def __init__(self, pdb_filename: Path, resolution: float = 0.5, ):
        """
        ----------
        resolution: resolution in angstroms of mesh
        pdb_filename
        """
        super().__init__()
        self._file_name = pdb_filename
        self._resolution = resolution
        self._mesh = None
        self._vertices, self._faces, self._normals, self._areas, self._names = compute_ms(self.file_name)

    @property
    def vertices(self):
        return self._vertices

    @property
    def faces(self):
        return self._faces

    @property
    def file_name(self):
        return self._file_name

    @property
    def resolution(self):
        return self._resolution

    @property
    def names(self):
        return self._names

    @property
    def mesh(self) -> pymesh.Mesh:
        if self._mesh is None:
            self.repair()
        else:
            return self._mesh

    @property
    def vertex_normals(self):
        if self._mesh is None:
            self.repair()
        else:
            return self._vertex_normals

    def repair(self):
        mesh = pymesh.form_mesh(self._vertices, self._faces)
        fixed_mesh = fix_mesh(mesh, self._resolution)
        self._vertex_normals = compute_normal(fixed_mesh.vertices, fixed_mesh.faces)
        self._mesh = fixed_mesh
        self._vertices = fixed_mesh.vertices
        self._faces = fixed_mesh.faces
