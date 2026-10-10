import numpy as np
from Bio.PDB.vectors import Vector, calc_dihedral
from Bio.SVDSuperimposer import SVDSuperimposer

##### CHIRALITY ######

def triple_product(coords):
    
    if not all(atom in coords for atom in ("N", "CA", "C", "CB")):
        return None

    N = coords["N"]
    CA = coords["CA"]
    C = coords["C"]
    CB = coords["CB"]

    vec1 = N - CA
    vec2 = C - CA
    vec3 = CB - CA
    return float(np.dot(vec1, np.cross(vec2, vec3))) # (N - CA) . ((C - CA) x (CB - CA)); Positive = L, Negative = D


def chirality_label(coords, planar_cutoff):
    # Label of 'L', 'D', 'planar' (too flat to determine) or None (glycine or missing atoms)
    value = triple_product(coords)
    if value is None:
        return None
    if abs(value) < planar_cutoff:
        return "planar"
    elif value > 0:
        return "L"
    else:
        return "D"


#### SECONDARY STRUCTURE ######

def phi_psi(residues, i):
    if i == 0 or i == len(residues) - 1:
        return None, None

    prev_res = residues[i - 1]
    this_res = residues[i]
    next_res = residues[i + 1]

    # If atoms missing flanking or in the residue, or there's a gap between residues, return nothing
    if "C" not in prev_res or "N" not in next_res or not all(atom in this_res for atom in ("N", "CA", "C")):
        return None, None        
    if np.linalg.norm(this_res["N"] - prev_res["C"]) > 2.0 or np.linalg.norm(next_res["N"] - this_res["C"]) > 2.0:
        return None, None

    phi = np.degrees(calc_dihedral(
        Vector(prev_res["C"]), Vector(this_res["N"]),
        Vector(this_res["CA"]), Vector(this_res["C"]),
    ))
    psi = np.degrees(calc_dihedral(
        Vector(this_res["N"]), Vector(this_res["CA"]),
        Vector(this_res["C"]), Vector(next_res["N"]),
    ))
    return phi, psi


def is_within_range(angle_1, angle_2, range_1, range_2):
    return range_1[0] <= angle_1 <= range_1[1] and range_2[0] <= angle_2 <= range_2[1]


def determine_secondary_structure(phi, psi):
    if phi is None or psi is None:
        return "terminus"
    # Check if the angles fit an alpha-helix, for either an L-residue or a D-residue
    if is_within_range(phi, psi, (-100, -30), (-77, -5)) or is_within_range(-phi, -psi, (-100, -30), (-77, -5)):
        return "a-helix"
    # Check if the angles fit a beta-sheet, for either an L-residue or a D-residue
    if is_within_range(phi, psi, (-170, -90), (90, 170)) or is_within_range(-phi, -psi, (-170, -90), (90, 170)):
        return "b-sheet"
    return "loop"


##### RMSD ######

def paired_ca_coordinates(crystal_residues, predicted_residues):
    n_residues = min(len(crystal_residues), len(predicted_residues))
    crystal_atoms = []
    predicted_atoms = []
    for i in range(n_residues):
        crystal_res = crystal_residues[i]
        predicted_res = predicted_residues[i]
        if "CA" in crystal_res and "CA" in predicted_res:
            crystal_atoms.append(crystal_res["CA"])
            predicted_atoms.append(predicted_res["CA"])
    return np.array(crystal_atoms), np.array(predicted_atoms)


def shape_rmsd(crystal_xyz, predicted_xyz):
    # Superimpose AF3 predicted D-peptide onto crystal D-peptide to get rmsd of the D-peptide binder
    superimposer = SVDSuperimposer()
    superimposer.set(crystal_xyz, predicted_xyz)
    superimposer.run()
    return superimposer.get_rms()


def pose_rmsd(crystal_protein_xyz, predicted_protein_xyz, crystal_peptide_xyz, predicted_peptide_xyz):
    # Superimpose AF3 predicted protein-peptide complex with crystal protein-peptide complex to find the rmsd of binding pose
    # Superimpose AF3 predicted protein with crystal protein to find rmsd of the protein target
    superimposer = SVDSuperimposer()
    superimposer.set(crystal_protein_xyz, predicted_protein_xyz)
    superimposer.run()
    rotation, translation = superimposer.get_rotran()
    moved_peptide = np.dot(predicted_peptide_xyz, rotation) + translation

    distances = np.linalg.norm(moved_peptide - crystal_peptide_xyz, axis=1)   # distance in Angstroms per atom
    rmsd = np.sqrt((distances ** 2).mean())
    return rmsd, superimposer.get_rms()  # returns (pose_rmsd, protein_rmsd)
