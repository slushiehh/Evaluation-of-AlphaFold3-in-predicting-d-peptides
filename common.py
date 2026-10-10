from io import StringIO

import requests
from Bio.PDB import PDBIO, MMCIFParser, PDBParser
from Bio.PDB.Chain import Chain
from rdkit import Chem
 
# Kyte-Doolittle Hydrophobicity scale
HYDROPHOBICITY = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
    "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
    "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

# Needed since RDKit only reads standard residue names
D_TO_L_THREE = {
    "DAL": "ALA", "DAR": "ARG", "DSG": "ASN", "DAS": "ASP", "DCY": "CYS",
    "DGN": "GLN", "DGL": "GLU", "DHI": "HIS", "DIL": "ILE", "DLE": "LEU",
    "DLY": "LYS", "DPN": "PHE", "DPR": "PRO", "DSN": "SER", "DTH": "THR",
    "DTR": "TRP", "DTY": "TYR", "DVA": "VAL",
}

# Needed since AF3 reads one letter residues
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",

    "DAL": "A", "DAR": "R", "DSG": "N", "DAS": "D", "DCY": "C",
    "DGN": "Q", "DGL": "E", "DHI": "H", "DIL": "I", "DLE": "L",
    "DLY": "K", "DPN": "F", "DPR": "P", "DSN": "S", "DTH": "T",
    "DTR": "W", "DTY": "Y", "DVA": "V",
}


def to_l_name(resname):
    return D_TO_L_THREE.get(resname, resname)


####### GET STRUCTURES AND RESIDUES ##########

def download_crystal_structure(pdb_id, save_path):
    url = "https://files.rcsb.org/download/" + pdb_id + ".cif"
    response = requests.get(url, timeout=30)
    response.raise_for_status()  # check if PDB ID url exists
    with open(save_path, "w") as file:
        file.write(response.text)
    print(" Downloaded", pdb_id, "to", save_path)


def get_structure(filepath):
    # Choose the right parser for either .cif or pdb text files
    if str(filepath).lower().endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("structure", filepath)


def get_residues(structure, chain_id=None, only_resnames=None):
    model = structure[0]  # first crystal structure

    if chain_id is None:
        chain = next(model.get_chains())
    else:
        chain = model[chain_id]
    residues = []

    for residue in chain:
        residue_name = residue.get_resname()
        if residue_name == "HOH":
            continue  # skip water molecules
        if only_resnames is not None and residue_name not in only_resnames:
            continue  # skip residues that are not in our list
        if all(atom in residue for atom in ("N", "CA", "C")):
            residues.append(residue)
    return residues


def get_one_letter_sequence(residues):
    letters = []
    for residue in residues:
        resname = residue.get_resname()
        letter = THREE_TO_ONE.get(resname)
        if letter is None:
            print(" Unknown residue '" + resname + "' - using 'X' ")
            letter = "X"
        letters.append(letter)
    return "".join(letters)


def get_protein_residues(structure, chain_id):
    kept = []
    previous_number = None
    for residue in structure[0][chain_id]:
        if residue.id[0] != " ":
            continue  # HETATM: water, ligands, the D-residues of 5N8T
        if residue.get_resname() not in THREE_TO_ONE:
            print("\nWarning: Skipping unknown residue ", residue)
            continue
        if "CA" not in residue:
            print("\nWarning: Skipping residue without alpha carbon ", residue)
            continue
        if previous_number != None and residue.id[1] != previous_number + 1:
            print("\nWarning: There is a gap in the chain between residue ", residue.id[1], " and ", previous_number)
        kept.append(residue)
        previous_number = residue.id[1]
    return kept


def residue_coords(residue):
    # Each Bio.PDB residue has atom coordinates in a dictionary {'N': array([x, y, z]), 'CA': ..., 'C': ..., 'CB': ...}.
    return {atom.get_name(): atom.get_coord() for atom in residue}


############# CRYSTAL PEPTIDE TO SMILES ##################
def peptide_pdb_block(residues):
    chain = Chain("A")

    for position, residue in enumerate(residues):
        residue_copy = residue.copy()
        residue_copy.resname = to_l_name(residue_copy.resname) 
        residue_copy.id = (" ", position + 1, " ")  # blank hetero flag = ATOM record
        chain.add(residue_copy)

    pdb_writer = PDBIO()
    pdb_writer.set_structure(chain)
    text = StringIO()
    pdb_writer.save(text)
    return text.getvalue()


def peptide_to_mol(residues):
    molecule = Chem.MolFromPDBBlock(peptide_pdb_block(residues), removeHs=False)
    if molecule is None:
        raise ValueError("RDKit could not read the peptide")

    Chem.AssignStereochemistryFrom3D(molecule)
    return molecule


def peptide_to_smiles(residues):
    molecule = peptide_to_mol(residues)   
    return Chem.MolToSmiles(molecule)


###### LIGAND ATOMS TO RESIDUES ########

def map_af3_ligand(af3_structure, crystal_residues, input_smiles):
    reference = Chem.MolFromPDBBlock(peptide_pdb_block(crystal_residues)) 
    ligand = Chem.MolFromSmiles(input_smiles)  

    af3_atoms = [atom for chain in af3_structure[0] for residue in chain if residue.get_resname().startswith("LIG") for atom in residue]

    af3_elements = [a.element.capitalize() for a in af3_atoms]
    smiles_elements = [a.GetSymbol() for a in ligand.GetAtoms()]
    if af3_elements != smiles_elements:
        raise ValueError(str(len(af3_elements)) + " AF3 atoms vs " + str(len(smiles_elements)) + " SMILES atoms don't line up")

    parameters = Chem.AdjustQueryParameters.NoAdjustments() 
    parameters.makeBondsGeneric = True  
    match = ligand.GetSubstructMatch(Chem.AdjustQueryProperties(reference, parameters))

    residues = [{} for _ in crystal_residues]

    for reference_index, ligand_index in enumerate(match):
        residue_info = reference.GetAtomWithIdx(reference_index).GetPDBResidueInfo()
        position = residue_info.GetResidueNumber() - 1
        atom_name = residue_info.GetName().strip()
        residues[position][atom_name] = af3_atoms[ligand_index].get_coord()

    return residues

