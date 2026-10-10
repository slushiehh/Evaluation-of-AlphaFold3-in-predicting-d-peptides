import json

import numpy as np
import pandas as pd

from common import (
    HYDROPHOBICITY,
    get_protein_residues,
    get_residues,
    get_structure,
    map_af3_ligand,
    residue_coords,
    to_l_name,
)
from config import (
    COMPLEXES,
    PLANAR_CUTOFF,
    RESULTS_DIR,
    TERMINAL_RESIDUES,
    af3_results_folder,
    crystal_path,
)
from structure_calculations import (
    chirality_label,
    determine_secondary_structure,
    paired_ca_coordinates,
    phi_psi,
    pose_rmsd,
    shape_rmsd,
)

#### READING AF3 FILES #####

def load_crystal(complex_):
    structure = get_structure(crystal_path(complex_))

    peptide_residues = get_residues(structure, complex_["peptide_chain"], complex_["peptide_resnames"])
    peptide_coords = [residue_coords(residue) for residue in peptide_residues]

    protein_residues = get_protein_residues(structure, complex_["protein_chain"])
    protein_ca = np.array([residue["CA"].get_coord() for residue in protein_residues])

    # {"peptide_residues": peptide Bio.PDB residues, "peptide_coords": dict of each atom's coords, "protein_ca": array of each alpha carbon coords of protein target}
    return {"peptide_residues": peptide_residues, "peptide_coords": peptide_coords, "protein_ca": protein_ca}


def read_input_smiles(results_folder):
    # Read the AF3 job json file to get the SMILES string again
    data_file = next(results_folder.glob("*_data.json"))
    with open(data_file) as file:
        data = json.load(file)

    for sequence_entry in data["sequences"]:
        if "ligand" in sequence_entry:
            return sequence_entry["ligand"]["smiles"]


def parse_seed_sample(sample_seed_folders):
    # AF3 results folders "seed-x_sample-y" parsed into (x, y)
    seed_part, sample_part = sample_seed_folders.split("_")  
    seed = int(seed_part.split("-")[1])               
    sample = int(sample_part.split("-")[1])           
    return seed, sample


def read_confidences(sample_dir):
    # AF3 confidence score of one sample
    confidence_file = next(sample_dir.glob("*_summary_confidences.json"))

    with open(confidence_file) as file:
        confidence_data = json.load(file)

    return {"ranking_score": confidence_data["ranking_score"], "iptm": confidence_data["iptm"], 
            "ptm": confidence_data["ptm"], "has_clash": confidence_data["has_clash"]}




#### EVALUATE EACH RESIDUE'S FEATURES #####

def build_feature_rows(crystal_residues, crystal_coords, predicted_coords):
    # Add a row with values for feature columns for each residue (each row are values for a residue)

    rows = []
    n_residues = min(len(crystal_coords), len(predicted_coords)) 

    for i in range(n_residues):
        crystal_label = chirality_label(crystal_coords[i], PLANAR_CUTOFF)
        predicted_label = chirality_label(predicted_coords[i], PLANAR_CUTOFF)

        # glycine or missing atoms
        if crystal_label is None or predicted_label is None:
            continue  
        
        # if scalar triple product < planar_cutoff
        if predicted_label == "planar":
            inverted = None  
        else:
            if crystal_label != predicted_label:
                inverted = 1
            else:
                inverted = 0

        phi, psi = phi_psi(crystal_coords, i)
        secondary_structure = determine_secondary_structure(phi, psi)

        if i < TERMINAL_RESIDUES or i >= n_residues - TERMINAL_RESIDUES:
            terminal_residue = True
        else:
            terminal_residue = False

        residue_name = crystal_residues[i].get_resname()
        rows.append({
            "position": i,
            "residue_name": residue_name,
            "crystal_label": crystal_label,
            "predicted_label": predicted_label,
            "inverted": inverted,
            "hydrophobicity": HYDROPHOBICITY[to_l_name(residue_name)],
            "terminal_residue": terminal_residue,
            "secondary_structure": secondary_structure,
        })

    return rows   # [{position, residue_name, crystal_label, predicted_label, inverted, hydrophobicity, terminal_residue, secondary_structure}]



#### EVALUATE EACH SAMPLE ######

def evaluate_sample(sample_dir, crystal, input_smiles):
    af3_model = get_structure(next(sample_dir.glob("*_model.cif")))
    predicted_coords = map_af3_ligand(af3_model, crystal["peptide_residues"], input_smiles)

    rows = build_feature_rows(crystal["peptide_residues"], crystal["peptide_coords"], predicted_coords)

    crystal_ca, predicted_ca = paired_ca_coordinates(crystal["peptide_coords"], predicted_coords)
    predicted_protein_ca = np.array([residue["CA"].get_coord() for residue in get_protein_residues(af3_model, "A")]) # chain A is the protein in the af3 job in 1.make_af3_job.py
    binding_pose, protein_fit = pose_rmsd(crystal["protein_ca"], predicted_protein_ca, crystal_ca, predicted_ca)

    n_inverted = sum(row["inverted"] == 1 for row in rows)
    n_planar = sum(row["predicted_label"] == "planar" for row in rows)

    sample_row = {
        "n_residues": len(rows),
        "n_inverted": n_inverted,
        "n_planar": n_planar,
        "chirality_violation_rate": n_inverted / (len(rows) - n_planar),  # planars shouldn't count
        "shape_rmsd": shape_rmsd(crystal_ca, predicted_ca),
        "binding_pose_rmsd": binding_pose,
        "protein_rmsd": protein_fit
    }
    sample_row.update(read_confidences(sample_dir)) # adds {"ranking_score":, "iptm":, "ptm":, "has_clash":} of the sample to the dictionary
    return rows, sample_row    
    # sample row = {"n_residues":, "n_inverted":, "n_planar":, "chirality_violation_rate":, "shape_rmsd":,  "binding_pose_rmsd":, "protein_rmsd":, "ranking_score":, "iptm":, "ptm":, "has_clash":}


#### MAIN SCRIPT #####

def main():
    residue_rows = []
    sample_rows = []

    for complex_ in COMPLEXES:
        results_folder = af3_results_folder(complex_)
        if not results_folder.exists():
            print("\nNo AF3 results for ", complex_["pdb_id"], "-  ", results_folder, "doesn't exist")
            continue

        crystal = load_crystal(complex_)
        input_smiles = read_input_smiles(results_folder)

        for sample_dir in sorted(results_folder.glob("seed-*_sample-*")):
            seed, sample = parse_seed_sample(sample_dir.name)
            rows, sample_row = evaluate_sample(sample_dir, crystal, input_smiles)

            for row in rows:
                row.update({"pdb_id": complex_["pdb_id"], "seed": seed, "sample": sample})
            sample_row.update({"pdb_id": complex_["pdb_id"], "seed": seed, "sample": sample})

            residue_rows.extend(rows)
            sample_rows.append(sample_row)
 
    pd.DataFrame(residue_rows).to_csv(RESULTS_DIR / "each_residue_results.csv", index=False)
    pd.DataFrame(sample_rows).to_csv(RESULTS_DIR / "each_sample_results.csv", index=False)


if __name__ == "__main__":
    main()
