from pathlib import Path

### PATHS ###
CRYSTAL_DIR = Path("data/crystal")
AF3_INPUT_DIR = Path("data/af3_inputs")
AF3_RESULTS_DIR = Path("data/af3_results")
RESULTS_DIR = Path("data/results")
FIGURES_DIR = Path("data/figures")

### SETTINGS ###
MODEL_SEEDS = [1, 2, 3, 4, 5]  # each seed creates 5 AF3 samples

TERMINAL_RESIDUES = 2  # How many residues at the start and end of the peptide count as "terminal"

# Threshold in cubic angstroms where alpha-carbon is too flat to be labelled either L or D 
PLANAR_CUTOFF = 1.0

COMPLEX_COLOURS = {"7YH8": "red", "5N8T": "blue", "3IWY":"green"}    # Colour for data points of each L-protein D-peptide complex pdb id

# L-Protein D-peptide complexes to predict
COMPLEXES = [
    {
        "pdb_id": "7YH8",
        "peptide_chain": "B",
        "protein_chain": "A",
        "peptide_resnames": None  # Specifies which of the chain is the D-peptide, here all of chain B is the peptide
    },
    {
        "pdb_id": "5N8T",
        "peptide_chain": "A",  
        "protein_chain": "A",
        "peptide_resnames": ["DLE", "DTR", "DGN", "DHI", "DGL", "DAL", "DTH", "DLY"] # only one chain and only 9/13 peptide residues are in chain A, so need to specify manually
    },
    {
        "pdb_id": "3IWY",
        "peptide_chain": "B",     # Use author chain bc Biopython looks at auth chain
        "protein_chain": "A",
        "peptide_resnames": None
    }
    
]

def job_name(complex_):
    return complex_["pdb_id"] + "_" + complex_["peptide_chain"] + "_dpeptide_complex"


def crystal_path(complex_):
    return CRYSTAL_DIR / (complex_["pdb_id"] + "_crystal.cif")


def af3_results_folder(complex_):
    return AF3_RESULTS_DIR / job_name(complex_)