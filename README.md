# NSCLC Metabolomics & Survival Analysis Pipeline

This repository contains the data preparation, quality audit, and statistical exploration pipeline for studying plasma metabolomics, clinical covariates, and survival outcomes in NSCLC patients and healthy controls.

---

## Repository Contents

The repository currently includes two core Python scripts:

| Script | Description |
| :--- | :--- |
| [`data_corrections.py`](data_corrections.py) | Corrects 6 known cell-level defects in the primary patient spreadsheet (duplicate patient IDs, malformed dates, decimal formatting, misclassified cohort, and transposed dates) before any downstream parsing occurs. Outputs audit logs and hard assertions. |
| [`survival_data_exploration_oct15.py`](survival_data_exploration_oct15.py) | Handles raw multi-level Excel data ingestion, clinical standardization, SCLC exclusion, survival endpoint derivation, Probabilistic Quotient Normalization (PQN), differential abundance testing (FDR-corrected rank-sum), paired pre/post treatment analysis, Kaplan-Meier curves, and Cox proportional hazards regression. |

---

## Environment Setup & Installation

We recommend using **Conda** to manage dependencies.

### 1. Create and Activate Conda Environment

```bash
# Create a dedicated conda environment with Python 3.11/3.12
conda create -n nsclc_metabolomics python=3.11 -y

# Activate the environment
conda activate nsclc_metabolomics
```

### 2. Install Required Packages

Install all necessary statistical, bioinformatic, and plotting libraries via `conda` (using `conda-forge`) and `pip`:

```bash
# Core numerical, scientific, and visualization libraries
conda install -c conda-forge \
    numpy \
    pandas \
    scipy \
    statsmodels \
    matplotlib \
    seaborn \
    openpyxl \
    tabulate \
    -y

# Survival analysis and plot text placement
pip install lifelines adjustText
```

#### Package Summary

- **`pandas` & `openpyxl`**: Multi-index Excel parsing, table manipulation, and Excel workbook reporting.
- **`numpy` & `scipy`**: Numerical matrix computations, median quotient normalization, and non-parametric hypothesis tests (`ranksums`, `wilcoxon`, `spearmanr`, `fisher_exact`).
- **`statsmodels`**: Benjamini-Hochberg False Discovery Rate (FDR) multitest corrections.
- **`lifelines`**: Survival analysis (`KaplanMeierFitter`, `CoxPHFitter`, log-rank tests).
- **`matplotlib` & `seaborn`**: Publication-ready styling and plot rendering (volcano plots, boxplots, heatmaps, Kaplan-Meier curves).
- **`adjustText`**: Dynamic collision avoidance for text labels in volcano plots.
- **`tabulate`**: Terminal and markdown table formatting for provenance/audit reports.

---

## Data Preparation

The pipeline expects the study Excel workbook in the working directory:
```
Combien GeoC& Metabo Epigenetic 08_26_2026 AM an AH without identifiers Responders vs Non-Responders .xlsx
```
*(Alternatively, update `FILE_PATH` near the top of `survival_data_exploration_oct15.py` to point to your local copy).*


---

## How to Run

### 1. Full Exploration Pipeline (Default: PQN Normalized)
Runs the entire analytical battery across 254 metabolites, generates all diagnostic figures, and saves master statistical workbooks:

```bash
python survival_data_exploration_oct15.py
```
*(or explicitly with `--pqn`)*:
```bash
python survival_data_exploration_oct15.py --pqn
```

### 2. Sensitivity Analysis (Raw / Non-Normalized Metabolites)
To run all analyses on raw, un-normalized ion intensities for comparison:

```bash
python survival_data_exploration_oct15.py --no-pqn
# or
python survival_data_exploration_oct15.py --raw
```

### 3. Fast Execution: Table only
To regenerate baseline clinical demographics and cohort characteristics (Table 1) without running the full plot generation suite:

```bash
# With PQN normalization
python survival_data_exploration_oct15.py --table-one

# With raw intensities
python survival_data_exploration_oct15.py --table-one --no-pqn
```

---

## Output Artifacts

Execution automatically generates two organized output directories:

- **`figs_oct15/`**: 
  - `dilution_factor_distribution.png`: PQN diagnostic dilution factor distribution comparing cancer vs. control.
  - `volcano_cancer_vs_control.png`: FDR-adjusted volcano plot of differential metabolites.
  - `paired_pre_post_treatment.png`: Trajectory plots for chemotherapy/immunotherapy response.
  - `km_curves_*.png`: Kaplan-Meier survival curves stratified by high/low metabolite expression.
  - `all_statistical_tests.xlsx`: Multi-sheet Excel workbook compiling all statistical outputs.
- **`outputs_oct15/`**: Tabular CSV summaries, Table 1 outputs (`table_one_baseline.csv`), and data correction audit logs (`data_corrections_oct15.csv`).
