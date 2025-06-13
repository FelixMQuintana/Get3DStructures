import numpy as np
from Bio.PDB import calc_angle, Vector, calc_dihedral, NeighborSearch
from scipy.spatial import KDTree

from lib.const import acceptorPlaneAtom, acceptorAngleAtom, donorAtom, hbond_std_dev, polarHydrogens


def compute_angle_deviation(a, b, c, theta):
    return abs(calc_angle(Vector(a), Vector(b), Vector(c)) - theta)


# Compute the angle deviation from a plane
def compute_plane_deviation(a, b, c, d):
    dih = calc_dihedral(Vector(a), Vector(b), Vector(c), Vector(d))
    dev1 = abs(dih)
    dev2 = np.pi - abs(dih)
    return min(dev1, dev2)


# angle_deviation from ideal value. TODO: do a more data-based solution
def compute_angle_penalty(angle_deviation):
    # Standard deviation: hbond_std_dev
    return max(0.0, 1.0 - (angle_deviation / (hbond_std_dev)) ** 2)


def is_polar_hydrogen(atom_name, res):
    if atom_name in polarHydrogens[res.get_resname()]:
        return True
    else:
        return False


def is_acceptor_atom(atom_name, res):
    if atom_name.startswith("O"):
        return True
    else:
        if res.get_resname() == "HIS":
            if atom_name == "ND1" and "HD1" not in res:
                return True
            if atom_name == "NE2" and "HE2" not in res:
                return True
    return False


# Compute the list of backbone C=O:H-N that are satisfied. These will be ignored.
def compute_satisfied_CO_HN(atoms):
    ns = NeighborSearch(atoms)
    satisfied_CO = set()
    satisfied_HN = set()
    for atom1 in atoms:
        res1 = atom1.get_parent()
        if atom1.get_id() == "O":
            neigh_atoms = ns.search(atom1.get_coord(), 2.5, level="A")
            for atom2 in neigh_atoms:
                if atom2.get_id() == "H":
                    res2 = atom2.get_parent()
                    # Ensure they belong to different residues.
                    if res2.get_id() != res1.get_id():
                        # Compute the angle N-H:O, ideal value is 180 (but in
                        # helices it is typically 160) 180 +-30 = pi
                        angle_N_H_O_dev = compute_angle_deviation(
                            res2["N"].get_coord(),
                            atom2.get_coord(),
                            atom1.get_coord(),
                            np.pi,
                        )
                        # Compute angle H:O=C, ideal value is ~160 +- 20 = 8*pi/9
                        angle_H_O_C_dev = compute_angle_deviation(
                            atom2.get_coord(),
                            atom1.get_coord(),
                            res1["C"].get_coord(),
                            8 * np.pi / 9,
                        )
                        ## Allowed deviations: 30 degrees (pi/6) and 20 degrees
                        #       (pi/9)
                        if (
                                angle_N_H_O_dev - np.pi / 6 < 0
                                and angle_H_O_C_dev - np.pi / 9 < 0.0
                        ):
                            satisfied_CO.add(res1.get_id())
                            satisfied_HN.add(res2.get_id())
    return satisfied_CO, satisfied_HN


def compute_charge_helper(atom_name, res, v):
    res_type = res.get_resname()
    # Check if it is a polar hydrogen.
    if is_polar_hydrogen(atom_name, res):
        donor_atom_name = donorAtom[atom_name]
        a = res[donor_atom_name].get_coord()  # N/O
        b = res[atom_name].get_coord()  # H
        # Donor-H is always 180.0 degrees, = pi
        angle_deviation = compute_angle_deviation(a, b, v, np.pi)
        angle_penalty = compute_angle_penalty(angle_deviation)
        return 1.0 * angle_penalty
    # Check if it is an acceptor oxygen or nitrogen
    elif is_acceptor_atom(atom_name, res):
        acceptor_atom = res[atom_name]
        b = acceptor_atom.get_coord()
        try:
            a = res[acceptorAngleAtom[atom_name]].get_coord()
        except:
            return 0.0
        # 120 degress for acceptor
        angle_deviation = compute_angle_deviation(a, b, v, 2 * np.pi / 3)
        # TODO: This should not be 120 for all atoms, i.e. for HIS it should be
        #       ~125.0
        angle_penalty = compute_angle_penalty(angle_deviation)
        plane_penalty = 1.0
        if atom_name in acceptorPlaneAtom:
            try:
                d = res[acceptorPlaneAtom[atom_name]].get_coord()
            except:
                return 0.0
            plane_deviation = compute_plane_deviation(d, a, b, v)
            plane_penalty = compute_angle_penalty(plane_deviation)
        return -1.0 * angle_penalty * plane_penalty
        # Compute the
    return 0.0
def assignChargesToNewMesh(new_vertices, old_vertices, old_charges):
    dataset = old_vertices
    testset = new_vertices
    new_charges = np.zeros(len(new_vertices))

    if True:#seeder_opts["feature_interpolation"]:
        num_inter = 4  # Number of interpolation features
        # Assign k old vertices to each new vertex.
        kdt = KDTree(dataset)
        dists, result = kdt.query(testset, k=num_inter)
        # Square the distances (as in the original pyflann)
        dists = np.square(dists)
        # The size of result is the same as new_vertices
        for vi_new in range(len(result)):
            vi_old = result[vi_new]
            dist_old = dists[vi_new]
            # If one vertex is right on top, ignore the rest.
            if dist_old[0] == 0.0:
                new_charges[vi_new] = old_charges[vi_old[0]]
                continue

            total_dist = np.sum(1 / dist_old)
            for i in range(num_inter):
                new_charges[vi_new] += (
                    old_charges[vi_old[i]] * (1 / dist_old[i]) / total_dist
                )
 #   else:
        # Assign k old vertices to each new vertex.
 #       kdt = KDTree(dataset)
 #       dists, result = kdt.query(testset)
 #       new_charges = old_charges[result]
    return new_charges