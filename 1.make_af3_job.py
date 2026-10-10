import json

from common import (
    download_crystal_structure,
    get_one_letter_sequence,
    get_protein_residues,
    get_residues,
    get_structure,
    peptide_to_smiles,
)
from config import (
    AF3_INPUT_DIR,
    COMPLEXES,
    MODEL_SEEDS,
    crystal_path,
    job_name,
)


def make_alphafold3_job(protein_sequence, smiles, job_name, save_path):
    job = {
        "name": job_name,
        "modelSeeds": MODEL_SEEDS,  
        "sequences": [
            {"protein": {"id": "A", "sequence": protein_sequence}},  # L-protein target
            {"ligand": {"id": "B", "smiles": smiles}},  # D-peptide binder
        ],
        "dialect": "alphafold3",
        "version": 4,
    }
    with open(save_path, "w") as file:
        json.dump(job, file, indent=1)
    print(" Made", save_path, "for AlphaFold3")


#### MAIN SCRIPT #####

def setup(complex_):
    pdb_id = complex_["pdb_id"]
    chain = complex_["peptide_chain"]
    crystal_file = crystal_path(complex_)
    only_resnames = complex_["peptide_resnames"]

    if not crystal_file.exists():
        download_crystal_structure(pdb_id, crystal_file)

    crystal_structure = get_structure(crystal_file)
    crystal_residues = get_residues(crystal_structure, chain, only_resnames)
    sequence = get_one_letter_sequence(crystal_residues)

    print("\n D-peptide is " + pdb_id + " chain " + chain + " (" + str(len(crystal_residues)) + " residues): ")

    smiles = peptide_to_smiles(crystal_residues) 
    print("  SMILES:", smiles)

    print("  Peptide sequence: ", sequence)
    protein_residues = get_protein_residues(crystal_structure, complex_["protein_chain"])
    protein_sequence = get_one_letter_sequence(protein_residues)
    print("  Protein sequence: ", protein_sequence)

    name = job_name(complex_)
    make_alphafold3_job(protein_sequence, smiles, name, AF3_INPUT_DIR / (name + ".json"))


def main():
    for complex_ in COMPLEXES:
        setup(complex_)


if __name__ == "__main__":
    main()
