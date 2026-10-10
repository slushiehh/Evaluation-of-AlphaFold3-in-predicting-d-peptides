import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from config import COMPLEX_COLOURS, FIGURES_DIR, RESULTS_DIR


def save_figure(name):
    plt.savefig(FIGURES_DIR / name, dpi=200, bbox_inches="tight")
    plt.close()


#### CHIRALITY VIOLATION RATE + RMSDs #####

def summary_table(_residues, samples):
    summary = samples.groupby("pdb_id")[["chirality_violation_rate", "shape_rmsd", "binding_pose_rmsd"]].agg(["mean", "std"])

    print("\n--------- Summary per complex ----------")
    print(summary)


#### AF3 CONFIDENCE vs CHIRALITY VIOLATION RATE #####

def confidence_vs_violation(_residues, samples):
    fig, ax = plt.subplots(figsize=(8, 7))
    for pdb_id, group in samples.groupby("pdb_id"):
        ax.scatter(group["ranking_score"], group["chirality_violation_rate"], color=COMPLEX_COLOURS[pdb_id], edgecolor="white", linewidth=1.5, label=pdb_id, zorder=3)

    ax.set_xlabel("AF3 ranking score (higher = more confident)")
    ax.set_ylabel("Fraction of residues inverted")
    ax.set_ylim(0, 1.05)                           
    ax.legend(title="PDB ID", title_fontsize=8, fontsize=6)
    save_figure("confidence_vs_violation_rate.png")
    


#### CHIRALITY VIOLATION RATE BY SECONDARY STRUCTURE #####

def feature_tables(residues, _samples):
    print("\n------ Chirality violation rate by secondary structure --------")
    print(residues.groupby(["pdb_id", "secondary_structure"]).agg(chirality_violation_rate=("inverted", "mean"), n_residue_samples=("inverted", "count")))



#### LOGISTIC REGRESSION #####

def logistic_regression(residues, _samples):
    data = residues.dropna(subset=["inverted"]).copy()
    data["inverted"] = data["inverted"].astype(int)
    data["terminal_residue"] = data["terminal_residue"].astype(int)
    data["residue_id"] = data["pdb_id"] + "_" + data["position"].astype(str)     # eg. 5N8T_1 for first residue of 5N8T

    formula = "inverted ~ hydrophobicity + terminal_residue + secondary_structure + pdb_id"  # inverted relationship to hydrophobicity + terminal_residue + secondary_structure and pdb_id
    model = smf.logit(formula, data=data)
    fitted = model.fit(cov_type="cluster", cov_kwds={"groups": data["residue_id"]})  # cluster by each residue
    table = pd.DataFrame({"odds_ratio": np.exp(fitted.params), "p_value": fitted.pvalues})

    print("\n-------- Logistic regression -----------")
    print(table.round(2))



#### MAIN SCRIPT #####

def main():
    residues = pd.read_csv(RESULTS_DIR / "each_residue_results.csv")
    samples = pd.read_csv(RESULTS_DIR / "each_sample_results.csv")
    analyses = [summary_table, confidence_vs_violation, feature_tables, logistic_regression]

    for analysis in analyses:
        analysis(residues, samples)


if __name__ == "__main__":
    main()
