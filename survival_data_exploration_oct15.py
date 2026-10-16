# Use the /Users/tnandi/Codes/ai_env conda env
from bdb import set_trace
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from pathlib import Path
from scipy.stats import ranksums, spearmanr, wilcoxon, fisher_exact
from statsmodels.stats.multitest import fdrcorrection
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from pdb import set_trace
import warnings
# Import data corrections module that corrects for known source spreadsheet defects
# (duplicate patient IDs, invalid dates, decimal commas, misclassified cohorts, and transposed dates)
import data_corrections as _dc

warnings.filterwarnings("ignore")

PUB_PLOT_PARAMS = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "Computer Modern Roman"],
    "font.size": 12,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 12,
}
plt.rcParams.update(PUB_PLOT_PARAMS)

BASE_DIR = Path(__file__).resolve().parent
FILE_PATH = str(BASE_DIR / "Combien GeoC& Metabo Epigenetic 08_26_2026 AM an AH without identifiers Responders vs Non-Responders .xlsx")

FIG_DIR = BASE_DIR / "figs_oct15"
OUTPUT_DIR = BASE_DIR / "outputs_oct15"
FIG_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

import re as _re
_KEGG_PREFIX = _re.compile(r'^(C\d{5}(?:\s*[/\-]\s*C\d{5})*)_')

survival_time_threshold = 30000 # Use 3000 to remove the high survival samples that interfere with survival modeling task

# --- Zero / non-detection handling -------------------------------------------
# A zero is a non-detection (left-censored), not a measured concentration.
#
# ZERO_MAX  routes a feature between the two differential tests: below this
#           fraction of zeros in BOTH arms it is tested on concentration, above
#           it in either arm it is tested on detection rate instead.
# pct_zeros excludes features from WITHIN-CANCER analyses only, computed on the
#           cancer cohort. Cancer-vs-control is deliberately NOT filtered on
#           sparsity, so that differential detection is captured rather than
#           discarded; only never-detected features are dropped there, because
#           no statistic exists for them.
ZERO_MAX = 0.5
pct_zeros = 0.5   # max fraction of zeros among CANCER patients for within-cancer work
SOCIAL_MODE_MAX = 0.5  # skip a social indicator whose modal value exceeds this share

# --- Metabolite Normalization Configuration ----------------------------------
# Set USE_PQN = True to apply Probabilistic Quotient Normalization (PQN)
# Set USE_PQN = False to analyze raw peak intensities without normalization
# Can also be overridden via command-line flags: --pqn or --no-pqn / --raw
USE_PQN = True

LOG2FC_EPS_PQN = 4.28  # half of global min non-zero PQN-normalized metabolite value
LOG2FC_EPS_RAW = 2.96  # half of global min non-zero raw metabolite value
LOG2FC_EPS = LOG2FC_EPS_PQN if USE_PQN else LOG2FC_EPS_RAW

# Volcano significance gates, applied to every volcano figure.
# The fold change is specified on the LINEAR scale and converted, since 1.25-fold
# is the interpretable quantity; the plot axis is log2.
VOLCANO_FDR      = 0.1
VOLCANO_FC       = 1.25
VOLCANO_LOG2FC   = float(np.log2(VOLCANO_FC))   # 0.3219

# The methylation assay writes non-detections as ~1e-22 rather than 0, and every
# marker shows a clean ~18-order-of-magnitude gap between that floor and its real
# 1e-3..5 range. Anything at or below this counts as not detected.
METH_DETECT_MIN = 1e-6
# Flatten-mangled follow-up draw timing columns (bogus Baseline_C00669_ prefix, double spaces)
FU_DRAW_DATE_COL = 'Baseline_C00669_Follow up Date blood was drawn'
FU_AFTER_IO_COL    = 'Baseline_C00669_Follow Up Blood Draw before  (No) or After (yes)  Immunotherapy'
FU_AFTER_CHEMO_COL = 'Baseline_C00669_Follow Up Blood Draw before  (No) or After (yes)  Chemotherapy'
FU_AFTER_SBRT_COL  = 'Baseline_C00669_Follow Up Blood Draw before  (No) or After (yes)  SBRT'

# Full names for the social/geospatial indicator codes
SOCIAL_LABELS = {
    'ADI':    'Area Deprivation Index',
    'EKI':    'Environmental burden index',
    'PMC':    'Particulate matter, PM2.5 (µg/m³)',
    'RHI':    'Respiratory hazard index',
    'TBC':    'Cigarette smoking rate (% of adults)',
    'CCG':    'Lung cancer diagnosis rate (per 100,000)',
    'CRV':    'Violent crime (per 100,000 residents)',
    'EDB':    'High school graduation rate (%)',
    'EKG':    'Green space proximity (% of area)',
    'EKR':    'Proximity to roads, railways, airports (index)',
    'FAI':    'Food insecurity (% of residents)',
    'HBS':    'Severe housing cost burden (% of units)',
    'HDX':    'Hardship Index (score)',
    'HTA':    'Single-parent households (%)',
    'ICR':    'Lifetime inhalation cancer risk (per million)',
    'INC':    'Median household income ($)',
    'LNG':    'COPD prevalence (% of adults)',
    'MCD':    'Medicaid coverage (% of residents)',
    'NVC':    'No vehicle available (% of households)',
    'POV':    'Poverty rate (% of residents)',
    'UMP':    'Unemployment rate (%)',
    'UNS':    'Uninsured rate (% of residents)',
    'PCT-A':  'Asian (% of residents)',
    'PCT-B':  'Non-Hispanic Black (% of residents)',
    'PCT-W':  'Non-Hispanic White (% of residents)',
    'PCT-H':  'Hispanic or Latino (% of residents)',
    'PCT-NA': 'Native American (% of residents)',
    'PCT-P':  'Pacific Islander / Native Hawaiian (%)',
    'PCT-F':  'Female (% of residents)',
    'PCT-M':  'Male (% of residents)',
}

# DNAm promoter methylation markers 
METH_COLS_ALL = ['CDO1. P', 'TAC1. P', 'SOX17.P', 'HOXA7.', 'HOXA9.P', 'ZFP42.P']
METH_LABELS = {
    'CDO1. P': 'CDO1 promoter', 'TAC1. P': 'TAC1 promoter', 'SOX17.P': 'SOX17 promoter',
    'HOXA7.': 'HOXA7 promoter', 'HOXA9.P': 'HOXA9 promoter', 'ZFP42.P': 'ZFP42 promoter',
}

# Tumor-infiltrate panel: populated for only 9/80 patients
IMMUNE_COLS = ['Macrophages', 'B Cells', 'NK Cells', 'T Cells (Total)', 'CD4', 'CD8']


def met_label(col, maxlen=None):
    """Strip 'Baseline_' prefix and leading KEGG ID(s) from a metabolite column name."""
    s = str(col)
    if s.startswith('Baseline_'):
        s = s[len('Baseline_'):]
    m = _KEGG_PREFIX.match(s)
    if m:
        s = s[m.end():]
    if maxlen and len(s) > maxlen:
        s = s[:maxlen - 1] + '…'
    return s
FIG_DIR.mkdir(exist_ok=True)


class DataLoader:
    def __init__(self, file_path):
        self.file_path = file_path
        
    def load_data(self):
        print(f"Loading data from {self.file_path}...")
        xls = pd.ExcelFile(self.file_path)
        df_raw = pd.read_excel(
            xls, 
            'Cancer, controls and follow ups',
            header=[0, 1, 2],
            keep_default_na=False,
            na_values=['N/A', '#N/A', 'NULL']
        )
        
        new_columns = []
        for col in df_raw.columns:
            level0, level1, level2 = col
            if ("Unnamed" in str(level0) or pd.isna(level0)) and pd.notna(level2) and "Unnamed" not in str(level2):
                clean_name = str(level2).strip()
            else:
                parts = [str(p).strip() for p in col if "Unnamed" not in str(p) and pd.notna(p)]
                clean_name = "_".join(parts)
            new_columns.append(clean_name)
        df_raw.columns = new_columns
        return df_raw

    # Clean the data for downstream survival and metabolomic analyses:
    # standardizes column names and clinical covariates (Sex, Race, Smoking, Stage),
    # excludes SCLC histology, parses dates, calculates Time to Event, and filters cohorts.
    def clean_data(self, df):
        print("Cleaning data...")
        df = df.copy()

        def _report(label, df_before, df_after, cohort_col=None):
            lost = len(df_before) - len(df_after)
            if cohort_col and cohort_col in df_before.columns and cohort_col in df_after.columns:
                lost_c = (df_before[cohort_col] == 'CANCER').sum() - (df_after[cohort_col] == 'CANCER').sum()
                lost_ctrl = (df_before[cohort_col] == 'No').sum() - (df_after[cohort_col] == 'No').sum()
                print(f"  [{label}] removed {lost} rows  (cancer={lost_c}, control={lost_ctrl})  → {len(df_after)} remaining")
            else:
                print(f"  [{label}] removed {lost} rows  → {len(df_after)} remaining")

        n_start = len(df)
        print(f"  [Start] {n_start} total rows loaded")

        rename_map = {
            'Vital Status\n0=Alive\n1=Dead': 'Event',
            'Race \nB= Black\nW= White\nA= Asian': 'Race',
            'Double Deidentifier': 'Patient ID',
            'Stage (AJCC CLINICAL STAGE)': 'Stage'
        }
        df.rename(columns=rename_map, inplace=True)

        if 'Histological Classification' in df.columns:
            df_before = df.copy()
            df = df[~df['Histological Classification'].astype(str).str.upper().str.contains('SCLC', na=False)]
            _report('SCLC exclusion', df_before, df)
            df['Histological Classification'] = df['Histological Classification'].astype(str).str.strip().str.title().replace({'Nan': np.nan, 'Na': np.nan, 'None': np.nan})

        if 'Sex' in df.columns:
            df['Sex'] = df['Sex'].astype(str).str.strip().str.upper().map({'MALE': 'Male', 'M': 'Male', 'FEMALE': 'Female', 'F': 'Female'})

        if 'Smoking Status' in df.columns:
            df['Smoking Status'] = df['Smoking Status'].astype(str).str.strip().str.upper().map({'CURRENT': 'Current', 'FORMER': 'Former', 'NEVER': 'Never'})

        if 'Race' in df.columns:
            df['Race'] = df['Race'].astype(str).str.strip().str.upper().map({'W': 'White', 'B': 'Black', 'A': 'Asian'})

        if 'Stage' in df.columns:
            # Preserve original stage label before collapsing
            stage_raw = df['Stage'].astype(str).str.strip().str.upper()
            df['Stage_orig'] = stage_raw.map(
                {'STAGE I': 'Stage I', 
                'STAGE II': 'Stage II',
                 'STAGE III': 'Stage III', 
                 'STAGE IV': 'Stage IV'}
            )
            df['Stage'] = stage_raw.replace(
                {'STAGE IV': 'Metastatic', 
                'STAGE III': 'Non-Metastatic', 
                'STAGE II': 'Non-Metastatic', 
                'STAGE I': 'Non-Metastatic', 
                'NAN': np.nan, 'NA': np.nan, 'NONE': np.nan}
            )

        date_cols = ['Date of Initial Diagnosis', 'Date of death', 'Last Day of Follow-up', 'Date of baseline blood draw (visit)']
        for col in date_cols:
            if col in df.columns:
                df[col] = pd.to_datetime(df[col], errors='coerce')

        df['Event'] = pd.to_numeric(df['Event'], errors='coerce').fillna(0).astype(int)

        # Survival clock ORIGIN: baseline blood draw (visit) when metabolites were measured
        # (fallback to Date of Initial Diagnosis if missing)
        fu_draw_cols = [c for c in df.columns if 'Follow up Date blood was drawn' in str(c)]
        if fu_draw_cols:
            df[fu_draw_cols[0]] = pd.to_datetime(df[fu_draw_cols[0]], errors='coerce')
            last_contact = df[['Last Day of Follow-up', fu_draw_cols[0]]].max(axis=1)
        else:
            last_contact = df['Last Day of Follow-up']

        end_date = df["Date of death"].where(df["Event"] == 1, last_contact)
        origin_date = df['Date of baseline blood draw (visit)'].fillna(df['Date of Initial Diagnosis'])
        df['Time to Event'] = (end_date - origin_date).dt.days

        # Determine cohort
        df['Cohort'] = 'No'
        if 'Cancer YES or NO' in df.columns:
            is_cancer = df['Cancer YES or NO'].astype(str).str.upper().str.contains('CANCER', na=False)
            df.loc[is_cancer, 'Cohort'] = 'CANCER'

        n_cancer_before = (df['Cohort'] == 'CANCER').sum()
        n_ctrl_before   = (df['Cohort'] == 'No').sum()
        print(f"  [After cohort labelling] cancer={n_cancer_before}, control={n_ctrl_before}")

        # Filter cancer patients: valid Time to Event and survival threshold
        cancer_mask   = df['Cohort'] == 'CANCER'
        neg_time_mask    = cancer_mask & (df['Time to Event'] <= 0)
        over_thresh_mask = cancer_mask & (df['Time to Event'] > survival_time_threshold)
        missing_tte_mask = cancer_mask & df['Time to Event'].isna()
        n_neg_time    = neg_time_mask.sum()
        n_over_thresh = over_thresh_mask.sum()
        n_missing_tte = missing_tte_mask.sum()
        print(f"  [Cancer filters] missing TTE={n_missing_tte}, TTE≤0={n_neg_time}, TTE>{survival_time_threshold}d={n_over_thresh}")
        # why were these rejected
        print("rejected:", df[(df['Cohort'] == 'CANCER') & (neg_time_mask | over_thresh_mask | missing_tte_mask)].index)
        # set_trace()

        # valid_cancer = cancer_mask & ~(neg_time_mask | over_thresh_mask | missing_tte_mask)
        valid_cancer = cancer_mask & ~(over_thresh_mask | missing_tte_mask)

        df_before = df.copy()
        df = df[(df['Cohort'] == 'No') | valid_cancer].copy()
        _report('Survival filter (cancer only)', df_before, df, cohort_col='Cohort')
        # set_trace()

        n_cancer_final = (df['Cohort'] == 'CANCER').sum()
        n_ctrl_final   = (df['Cohort'] == 'No').sum()
        print(f"  [Final] cancer={n_cancer_final}, control={n_ctrl_final}, total={len(df)}")

        return df
        
    def get_feature_groups(self, df):
        all_cols = df.columns
        exclude_keywords = ['Date', 'Follow Up', 'Draw', 'Immunotherapy', 'Chemotherapy', 'SBRT']
        metabolite_cols = [
            c for c in all_cols 
            if c.startswith('Baseline_C') 
            and not any(keyword in c for keyword in exclude_keywords)
        ]
        
        # Ensure numeric for metabolites to calculate sparsity properly
        for c in metabolite_cols:
             df[c] = pd.to_numeric(df[c], errors='coerce')

        # Two scoped feature lists, because the right exclusion depends on the
        # comparison being made:
        #
        #   metabolite_cols     WITHIN-CANCER work. Zero fraction computed on the
        #                       cancer cohort only -- controls are the wrong
        #                       denominator for an early-vs-late or survival split,
        #                       and a feature that is mostly zeros among cases makes
        #                       such a test degenerate (nearly all ties at zero).
        #                       Computed on the cohort rather than per arm so that
        #                       small arms do not drive exclusions.
        #
        #   metabolite_cols_cc  CANCER-vs-CONTROL. Deliberately NOT filtered on
        #                       sparsity: a difference in detection rate between the
        #                       groups is itself a finding, and a threshold would
        #                       discard exactly the features showing it. Only
        #                       never-detected features are removed, since no
        #                       statistic is defined for them.
        is_cancer = df['Cohort'] == 'CANCER'
        zero_cancer = (df.loc[is_cancer, metabolite_cols] == 0).mean()
        detected_any = (df[metabolite_cols] > 0).sum() > 0

        metabolite_cols_cc = [c for c in metabolite_cols if detected_any[c]]
        metabolite_cols = [c for c in metabolite_cols
                           if detected_any[c] and zero_cancer[c] < pct_zeros]

        print(f"  Metabolite features: {len(all_cols[[str(c).startswith('Baseline_C') for c in all_cols]])} raw "
              f"-> {len(metabolite_cols_cc)} for cancer-vs-control (detected at least once), "
              f"{len(metabolite_cols)} for within-cancer (<{pct_zeros:.0%} zeros among cases)")

        methylation_cols = [c for c in all_cols if c in ['CDO1. P', 'TAC1. P', 'SOX17.P', 'HOXA7.', 'HOXA9.P', 'ZFP42.P']]

        social_indicators = ['ADI', 'EKI', 'PMC', 'RHI', 'TBC', 'CCG', 'CRV', 'EDB', 'EKG', 'EKR', 'FAI', 'HBS', 'HDX', 'HTA', 'ICR', 'INC', 'LNG', 'MCD', 'NVC', 'POV', 'UMP', 'UNS', 'PCT-A', 'PCT-B', 'PCT-W', 'PCT-H', 'PCT-NA', 'PCT-P', 'PCT-F', 'PCT-M']
        social_cols = [c for c in social_indicators if c in all_cols]

        return metabolite_cols, methylation_cols, social_cols, metabolite_cols_cc

def plot_clinical_distribution(df):
    """4x4 clinical plots for CANCER patients."""
    print("Plotting clinical distributions for Cancer patients...")
    # Filter cancer cohort that meets survival criteria
    df_cancer = df[df['Cohort'] == 'CANCER'].copy()
    
    numeric_cols = ['Age at the time of diagnosis', 'BMI', 'Pack Year History', 'Size 1 of tumor (cms)', 'Time to Event']
    for c in numeric_cols:
        if c in df_cancer.columns:
            df_cancer[c] = pd.to_numeric(df_cancer[c], errors='coerce')
            
    PANEL_FS  = 25   # axis title
    TICK_FS   = 20   # tick labels
    LABEL_FS  = 25   # axis labels

    fig, axes = plt.subplots(4, 4, figsize=(32, 26))
    axes = axes.flatten()

    plot_idx = 0

    def _style(ax, title):
        ax.set_title(title, fontsize=PANEL_FS, fontweight='bold', pad=6)
        ax.tick_params(axis='both', labelsize=TICK_FS)
        ax.set_xlabel(ax.get_xlabel(), fontsize=LABEL_FS)
        ax.set_ylabel(ax.get_ylabel(), fontsize=LABEL_FS)

    def plot_dist(col, title):
        nonlocal plot_idx
        if col in df_cancer.columns:
            sns.histplot(df_cancer[col].dropna(), kde=True, ax=axes[plot_idx], color='skyblue')
            _style(axes[plot_idx], title)
            plot_idx += 1

    def plot_count(col, title):
        nonlocal plot_idx
        if col in df_cancer.columns:
            sns.countplot(y=df_cancer[col].rename(title).dropna(), ax=axes[plot_idx], color='grey')
            _style(axes[plot_idx], title)
            axes[plot_idx].set_ylabel('')
            plot_idx += 1

    # Distributions
    plot_dist('Age at the time of diagnosis', 'Distribution of Age at diagnosis')
    plot_dist('Time to Event', 'Distribution of Time to Event')
    plot_dist('BMI', 'Distribution of BMI')
    plot_dist('Pack Year History', 'Distribution of Pack Year History')
    plot_dist('Size 1 of tumor (cms)', 'Distribution of Tumor Size (cms)')
    
    # Counts
    plot_count('Event', 'Count of Event')
    plot_count('Race', 'Count of Race')
    plot_count('Sex', 'Count of Sex')
    plot_count('Histological Classification', 'Count of Histological Classification')
    plot_count('Stage', 'Count of Stage')
    plot_count('Smoking Status', 'Count of Smoking Status')
    plot_count('Primary Tumor Location', 'Count of Primary Tumor Location')
    
    # Comorbidity counts
    comorb_cols = ['HTN', 'CHF', 'CAD', 'DM', 'Thyroid', 'Emphysema', 'Asthma', 'COPD']
    comorb_counts = {}
    for c in comorb_cols:
        if c in df_cancer.columns:
            # Count "Yes" or "Hypo" 
            series = df_cancer[c].astype(str).str.strip().str.upper()
            comorb_counts[c] = series[series.isin(['YES', 'HYPO'])].shape[0]
            
    if comorb_counts:
        comorb_df = pd.DataFrame(list(comorb_counts.items()), columns=['Comorbidity', 'Count']).sort_values('Count', ascending=False)
        sns.barplot(x='Count', y='Comorbidity', data=comorb_df, ax=axes[plot_idx], color='grey')
        _style(axes[plot_idx], 'Comorbidity Counts')
        axes[plot_idx].set_ylabel('')
        plot_idx += 1

    # Hide unused axes
    for i in range(plot_idx, len(axes)):
        axes[i].axis('off')
        
    plt.tight_layout(pad=2.0)
    plt.savefig(FIG_DIR / 'overview_clinical_exploration.png', dpi=400, bbox_inches='tight')
    plt.close()


def rank_sum_test(df, features, is_log=False, group_col='Cohort', group1='CANCER', group2='No',
                  zero_max=ZERO_MAX):
    """Differential test between two groups, routed by each feature's zero structure.

    A zero in this assay is a non-detection (left-censored), not a measured
    concentration of zero. Testing a feature that is mostly zeros in one arm on
    the magnitude scale produces a fold change driven entirely by the detection
    floor -- that is where the +/-15 log2FC values in the volcano came from. So:

      * < `zero_max` zeros in BOTH arms -> Wilcoxon rank-sum on concentration.
      * otherwise, if detected at least once overall -> Fisher exact on the
        detection rate. log2FC is left NaN; the effect is an odds ratio.
      * never detected in either arm -> no statistic exists, feature is skipped
        (this is not a sparsity filter: there is nothing to test).

    One BH correction is applied across both tiers, since they answer the same
    question by different means.
    """
    results = []

    df1 = df[df[group_col] == group1]
    df2 = df[df[group_col] == group2]

    for f in features:
        raw1 = pd.to_numeric(df1[f], errors='coerce').dropna()
        raw2 = pd.to_numeric(df2[f], errors='coerce').dropna()
        if len(raw1) < 5 or len(raw2) < 5:
            continue

        det1, det2 = int((raw1 > 0).sum()), int((raw2 > 0).sum())
        z1, z2 = 1.0 - det1 / len(raw1), 1.0 - det2 / len(raw2)

        common = {
            'Feature': f,
            'Zero_Pct_1': z1, 'Zero_Pct_2': z2,
            'Detect_1': det1, 'Detect_2': det2,
            'Detect_Pct_1': det1 / len(raw1), 'Detect_Pct_2': det2 / len(raw2),
        }

        if z1 < zero_max and z2 < zero_max:
            vals1 = np.log1p(raw1) if is_log else raw1
            vals2 = np.log1p(raw2) if is_log else raw2
            _, p_val = ranksums(vals1, vals2)
            med1, med2 = vals1.median(), vals2.median()
            results.append({
                **common,
                'Test': 'concentration',
                'p_value': p_val,
                'Direction': "Increased" if med1 > med2 else "Decreased",
                'Median_Diff': med1 - med2,
                'log2FC': np.log2((raw1.median() + LOG2FC_EPS) / (raw2.median() + LOG2FC_EPS)),
                'Odds_Ratio': np.nan,
            })
        elif det1 + det2 > 0:
            table = [[det1, len(raw1) - det1], [det2, len(raw2) - det2]]
            odds, p_val = fisher_exact(table)
            results.append({
                **common,
                'Test': 'detection',
                'p_value': p_val,
                'Direction': "Increased" if common['Detect_Pct_1'] > common['Detect_Pct_2'] else "Decreased",
                'Median_Diff': np.nan,
                'log2FC': np.nan,          # undefined against a detection floor
                'Odds_Ratio': odds,
            })

    res_df = pd.DataFrame(results)
    if not res_df.empty:
        res_df['FDR'] = fdrcorrection(res_df['p_value'])[1]
        n_c = int((res_df['Test'] == 'concentration').sum())
        n_d = int((res_df['Test'] == 'detection').sum())
        print(f"    tested {len(res_df)} features: {n_c} concentration, {n_d} detection "
              f"(skipped {len(features) - len(res_df)} untestable)")

    return res_df

def plot_boxplots(df, feat_df, title, filename, is_log=False, standardize=False, hue_col='Cohort', hue_order=None):
    """Generates parameterized box plots mapping raw feature metadata."""
    if feat_df is None or feat_df.empty:
        print(f"No significant features to plot for {filename}.")
        return
        
    print(f"Plotting {filename}...")
    
    # Preprocess dataframe specifically for plotting
    plot_data = []
    
    for _, row in feat_df.iterrows():
        f = row['Feature']
        fdr    = row.get('FDR', np.nan)
        pval   = row.get('p_value', np.nan)
        log2fc = row.get('log2FC', np.nan)
        stat_parts = []
        if pd.notna(pval):   stat_parts.append(f"p={pval:.1e}")
        if pd.notna(fdr):    stat_parts.append(f"q={fdr:.1e}")
        if pd.notna(log2fc): stat_parts.append(f"log2FC={log2fc:.2f} ({2 ** log2fc:.2f}x)")
        stat_str = "\n".join(stat_parts)

        # Clean up labels to prevent them from dominating the plot space
        clean_f = met_label(f)
        if len(clean_f) > 20:
            clean_f = clean_f.replace('_', '\n')

        label = f"{clean_f}\n{stat_str}" if stat_str else clean_f
        
        vals = pd.to_numeric(df[f], errors='coerce')
        if is_log:
            vals = np.log1p(vals)
        if standardize and vals.std() > 0:
            vals = (vals - vals.mean()) / vals.std()
            
        temp_df = pd.DataFrame({
            'Value': vals,
            'Cohort_Split': df[hue_col],
            'Metabolite': label
        }).dropna()
        plot_data.append(temp_df)
        
    if not plot_data: return
    plot_df = pd.concat(plot_data)
    
    plt.figure(figsize=(16, 8))
    
    if hue_order is None:
        hue_order = plot_df['Cohort_Split'].unique()
        
    sns.boxplot(
        data=plot_df, x='Metabolite', y='Value', hue='Cohort_Split', 
        hue_order=hue_order,
        palette='Set1' if hue_col != 'Cohort' else {'CANCER': 'tab:blue', 'No': 'tab:orange'},
        showfliers=False, width=0.6, linewidth=1.0
    )
    
    plt.title(title, fontsize=12)
    plt.xticks(rotation=45, ha='right')
    ylab = "Value"
    if is_log: ylab = "Value (log transformed)"
    if standardize: ylab = "Value (Z-scored)"
    plt.ylabel(ylab)
    plt.xlabel("Metabolite" if "DNAm" not in title else "Site")
    plt.legend(title=hue_col, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()

def plot_volcano(res_df, title, filename, p_col='FDR', diff_col='log2FC',
                 p_thresh=VOLCANO_FDR, diff_thresh=VOLCANO_LOG2FC,
                 group_a='Group 1', group_b='Group 2'):
    """Volcano plot: x = log2FC (with pseudocount), y = -log10(FDR). Standard metabolomics convention."""
    from adjustText import adjust_text as _adjust_text
    if res_df is None or res_df.empty:
        return

    df_v = res_df.copy()
    if diff_col not in df_v.columns:
        diff_col = 'Median_Diff'
    if p_col not in df_v.columns:
        p_col = 'p_value'

    # Detection-tier features have no defined fold change (one arm sits on the
    # detection floor), so they cannot go on a magnitude axis. Drop them here and
    # report the count -- plot_detection_tier draws them properly.
    n_detect = 0
    if 'Test' in df_v.columns:
        n_detect = int((df_v['Test'] == 'detection').sum())
        df_v = df_v[df_v['Test'] != 'detection']
    df_v = df_v[df_v[diff_col].notna() & df_v[p_col].notna()]
    if df_v.empty:
        print(f"  {filename}: no concentration-tier features to plot")
        return

    # No display cap on the y-axis. A floor on the q-value (e.g. clip at 1e-10)
    # collapses every hit below it onto one line -- here that was 50 features
    # spanning q = 2e-20 to 4e-11 all drawn at y = 10, which destroys the ranking
    # among exactly the features worth naming.
    df_v['_logy'] = -np.log10(df_v[p_col].replace(0, np.nextafter(0, 1)))

    sig_inc = (df_v[p_col] < p_thresh) & (df_v[diff_col] > diff_thresh)
    sig_dec = (df_v[p_col] < p_thresh) & (df_v[diff_col] < -diff_thresh)
    not_sig = ~(sig_inc | sig_dec)

    fig, ax = plt.subplots(figsize=(11, 8))
    ax.scatter(df_v.loc[not_sig,  diff_col], df_v.loc[not_sig,  '_logy'],
               color='#aaaaaa', alpha=0.45, s=22, label='Not significant', zorder=2)
    ax.scatter(df_v.loc[sig_inc, diff_col], df_v.loc[sig_inc, '_logy'],
               color='#d73027', alpha=0.9, s=50, label=f'Higher in {group_a}', zorder=3)
    ax.scatter(df_v.loc[sig_dec, diff_col], df_v.loc[sig_dec, '_logy'],
               color='#4575b4', alpha=0.9, s=50, label=f'Higher in {group_b}', zorder=3)

    # Label the 20 most significant overall, then force in the significant-increased
    # side. Without this the increased hits are never labelled whenever the decreased
    # side supplies all 20 top-ranked features, which is the usual case here.
    top = pd.concat([df_v[sig_inc | sig_dec].sort_values(p_col).head(20),
                     df_v[sig_inc].sort_values(p_col).head(10)]).drop_duplicates(subset='Feature')
    texts = []
    for _, row in top.iterrows():
        lbl = met_label(row['Feature'], maxlen=38)
        texts.append(ax.text(row[diff_col], row['_logy'], lbl, fontsize=8, zorder=4))
    if texts:
        _adjust_text(texts, ax=ax,
                     arrowprops=dict(arrowstyle='->', color='#555555', lw=0.6,
                                     shrinkA=2, shrinkB=3))

    ax.axhline(-np.log10(p_thresh), color='black', linestyle='--', lw=0.9, alpha=0.6,
               label=f'FDR = {p_thresh}')
    # Label the gate with the linear fold change as well -- "0.32" is not a quantity
    # anyone reads off a volcano, "1.25-fold" is.
    if diff_col == 'log2FC':
        diff_lbl = (rf"|$\log_2$(FC)| = {diff_thresh:.2f}"
                    rf"  ({2 ** diff_thresh:.2f}-fold)")
    else:
        diff_lbl = f'|median diff| = {diff_thresh:.2f}'
    ax.axvline( diff_thresh, color='gray', linestyle='--', lw=0.8, alpha=0.4,
               label=diff_lbl)
    ax.axvline(-diff_thresh, color='gray', linestyle='--', lw=0.8, alpha=0.4)

    # y = -log10(FDR) cannot go below 0 (that would be FDR > 1), so anchor the axis
    # exactly there rather than letting the default margin open a region no point
    # can occupy.
    y_max = float(df_v['_logy'].max())
    ax.set_ylim(bottom=0.0, top=y_max * 1.08)

    x_lbl = r"$\log_2$(FC)" if diff_col == 'log2FC' else 'Median Difference (log1p)'
    ax.set_xlabel(x_lbl, fontsize=13)
    ax.set_ylabel(r'-$\log_{10}$(FDR)', fontsize=13)
    ax.set_title(title, fontsize=13, fontweight='bold')
    if n_detect:
        fig.text(0.01, 0.005,
                 f"{n_detect} further feature(s) are mostly non-detected in one group; a fold "
                 f"change is undefined against a detection floor, so they are tested on detection "
                 f"rate and shown separately rather than plotted here.",
                 fontsize=8.5, color='#444444', va='bottom')
    ax.legend(fontsize=10, framealpha=0.85)
    ax.tick_params(labelsize=11)
    plt.tight_layout()
    plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    n_fdr_only = int((df_v[p_col] < p_thresh).sum())
    print(f"  Volcano saved: {filename}  |  FDR<{p_thresh} & |log2FC|>{diff_thresh:.3f} "
          f"({2 ** diff_thresh:.2f}-fold): ↑{sig_inc.sum()} ↓{sig_dec.sum()}  "
          f"(FDR<{p_thresh} alone: {n_fdr_only})")

def plot_detection_tier(res_df, title, filename, group_a='Group 1', group_b='Group 2',
                        fdr_alpha=0.05, n_max=30):
    """Detection-rate figure for features that cannot carry a fold change.

    These are the features sitting on the detection floor in one arm. The effect
    is the difference in the proportion of samples where the metabolite was seen
    at all, so that is what is drawn -- a paired dot per feature with a connector,
    rather than a magnitude on a log2FC axis where the value would be an artefact
    of the pseudocount.
    """
    if res_df is None or res_df.empty or 'Test' not in res_df.columns:
        return
    d = res_df[(res_df['Test'] == 'detection') & (res_df['FDR'] < fdr_alpha)].copy()
    if d.empty:
        print(f"  {filename}: no significant detection-tier features")
        return

    d['_gap'] = (d['Detect_Pct_1'] - d['Detect_Pct_2']).abs()
    d = d.nlargest(min(n_max, len(d)), '_gap').sort_values('Detect_Pct_1')
    y = np.arange(len(d))

    fig, ax = plt.subplots(figsize=(10, max(4.0, 0.34 * len(d) + 2.0)))
    ax.hlines(y, d['Detect_Pct_2'] * 100, d['Detect_Pct_1'] * 100,
              color='#bbbbbb', lw=1.8, zorder=1)
    ax.scatter(d['Detect_Pct_2'] * 100, y, s=62, color='#4575b4', zorder=3, label=group_b)
    ax.scatter(d['Detect_Pct_1'] * 100, y, s=62, color='#d73027', zorder=3, label=group_a)

    ax.set_yticks(y)
    ax.set_yticklabels([met_label(f, maxlen=34) for f in d['Feature']], fontsize=9)
    ax.set_xlabel('Samples with the metabolite detected (%)', fontsize=12)
    ax.set_xlim(-2, 102)
    ax.set_title(title, fontsize=13, fontweight='bold', loc='left')
    ax.legend(fontsize=10, frameon=False, loc='lower right')
    ax.grid(axis='x', color='#eeeeee', lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    fig.text(0.01, 0.005,
             f"Features mostly non-detected in one group, so a fold change is undefined. "
             f"Fisher exact on detection rate, FDR < {fdr_alpha:g} "
             f"(corrected jointly with the concentration-tier tests). "
             f"Showing {len(d)} of {int(((res_df['Test'] == 'detection') & (res_df['FDR'] < fdr_alpha)).sum())} "
             f"significant, ranked by detection-rate gap.",
             fontsize=8.5, color='#444444', va='bottom')

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  {filename}: {len(d)} features plotted")


def generate_km_plots(df, features, is_categorical=False):
    """Generates KM survival plots for each feature and returns log-rank statistics."""
    km_dir = FIG_DIR / 'km_plots'
    km_dir.mkdir(exist_ok=True)
    
    log_rank_results = []
    
    # Filter for cancer cohort 
    df_cancer = df[df['Cohort'] == 'CANCER'].copy()
    
    for f in features:
        if f not in df_cancer.columns:
            continue
            
        if not is_categorical:
            vals = pd.to_numeric(df_cancer[f], errors='coerce')
            if vals.notna().sum() > 5: # Require at least a few instances
                median_val = vals.median()
                high_mask = vals > median_val
                low_mask = vals <= median_val
                
                valid_mask = vals.notna() & df_cancer['Time to Event'].notna() & df_cancer['Event'].notna()
                if high_mask[valid_mask].sum() == 0 or low_mask[valid_mask].sum() == 0:
                    continue
                    
                T_high = df_cancer.loc[high_mask & valid_mask, 'Time to Event']
                E_high = df_cancer.loc[high_mask & valid_mask, 'Event']
                T_low = df_cancer.loc[low_mask & valid_mask, 'Time to Event']
                E_low = df_cancer.loc[low_mask & valid_mask, 'Event']
                
                try:
                    results = logrank_test(T_high, T_low, event_observed_A=E_high, event_observed_B=E_low)
                    p_val = results.p_value
                except Exception:
                    continue
                    
                log_rank_results.append({'Feature': f, 'P_Value': p_val, 'Type': 'Continuous (> Median)'})
                    
                plt.figure(figsize=(8, 6))
                kmf_high = KaplanMeierFitter()
                kmf_low = KaplanMeierFitter()
                
                kmf_high.fit(T_high, event_observed=E_high, label=f'> Median ({median_val:.2f})')
                kmf_low.fit(T_low, event_observed=E_low, label=f'<= Median ({median_val:.2f})')
                
                ax = kmf_high.plot()
                kmf_low.plot(ax=ax)
                
                plt.title(f"{f}\nLog-Rank p-value: {p_val:.4f}")
                plt.ylabel("Survival Probability")
                plt.xlabel("Time (Days)")
                plt.ylim(0, 1)
                plt.xlim(xmin=0)
                plt.tight_layout()
                
                clean_filename = str(f).replace('/', '_').replace(' ', '_').replace('\n', '')
                plt.savefig(km_dir / f"{clean_filename}.png", dpi=300)
                plt.close()
                
        else:
            temp_df = df_cancer[[f, 'Time to Event', 'Event']].copy()
            # String parsing for categorical to prevent Nan crashes
            temp_df[f] = temp_df[f].astype(str).str.strip().str.title()
            temp_df = temp_df.replace({'Nan': np.nan, 'Na': np.nan, 'None': np.nan}).dropna()
            
            if temp_df.empty:
               continue
               
            categories = temp_df[f].unique()
            if len(categories) < 2:
                continue
                
            try:
                results = multivariate_logrank_test(temp_df['Time to Event'], temp_df[f], temp_df['Event'])
                p_val = results.p_value
            except Exception:
                continue
                
            log_rank_results.append({'Feature': f, 'P_Value': p_val, 'Type': 'Categorical'})
                
            plt.figure(figsize=(8, 6))
            ax = plt.gca()
            for cat in categories:
                mask = temp_df[f] == cat
                if mask.sum() == 0:
                    continue
                kmf = KaplanMeierFitter()
                kmf.fit(temp_df.loc[mask, 'Time to Event'], event_observed=temp_df.loc[mask, 'Event'], label=str(cat))
                kmf.plot(ax=ax)
                
            plt.title(f"{f}\nLog-Rank p-value: {p_val:.4f}")
            plt.ylabel("Survival Probability")
            plt.xlabel("Time (Days)")
            plt.ylim(0, 1)
            plt.xlim(xmin=0)
            plt.tight_layout()
            
            clean_filename = str(f).replace('/', '_').replace(' ', '_').replace('\n', '')
            plt.savefig(km_dir / f"{clean_filename}.png", dpi=300)
            plt.close()
            
    return pd.DataFrame(log_rank_results)

def _order_profiles(M):
    """Hierarchical-clustering order for the rows of a rectangular matrix, using
    each row's correlation profile as its feature vector."""
    if M.shape[0] < 3:
        return list(M.index)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        V = M.fillna(0.0).values
        return [M.index[i] for i in leaves_list(linkage(V, method='average', metric='euclidean'))]
    except Exception as exc:
        print(f"  (clustering failed, keeping input order: {exc})")
        return list(M.index)


def cross_corr_matrices(df_sub, social_cols, met_cols, min_n=10):
    """Spearman rho and BH-FDR q for every social x metabolite pair (rows=social)."""
    soc = {s: pd.to_numeric(df_sub[s], errors='coerce')
           for s in social_cols if s in df_sub.columns}
    met = {m: np.log1p(pd.to_numeric(df_sub[m], errors='coerce'))
           for m in met_cols if m in df_sub.columns}

    rr, pp = {}, {}
    for s, sv in soc.items():
        row_r, row_p = {}, {}
        for m, mv in met.items():
            mask = sv.notna() & mv.notna()
            if mask.sum() < min_n or sv[mask].nunique() < 3 or mv[mask].nunique() < 3:
                row_r[m], row_p[m] = np.nan, np.nan
            else:
                r, p = spearmanr(sv[mask], mv[mask])
                row_r[m], row_p[m] = r, p
        rr[s], pp[s] = row_r, row_p

    rho = pd.DataFrame(rr).T
    pval = pd.DataFrame(pp).T
    flat = pval.values.ravel().astype(float)
    ok = ~np.isnan(flat)
    qf = np.full_like(flat, np.nan)
    if ok.sum():
        qf[ok] = fdrcorrection(flat[ok])[1]
    q = pd.DataFrame(qf.reshape(pval.shape), index=pval.index, columns=pval.columns)
    return rho, q


def plot_cross_corr_heatmap(rho, q, title, filename, fdr_alpha=0.05, footnote=""):
    """Social x metabolite heatmap.

    This matrix is rectangular and asymmetric, so unlike the within-block
    correlation figures there is no redundant half to mask -- every cell is a
    distinct pair and all of them are drawn.
    """
    if rho.empty or rho.shape[1] == 0:
        print(f"  {filename}: nothing to plot")
        return

    rho = rho.loc[_order_profiles(rho)]
    rho = rho[_order_profiles(rho.T)]
    q = q.loc[rho.index, rho.columns]

    n_r, n_c = rho.shape
    sig = (q.values < fdr_alpha) & rho.notna().values
    annotate = n_c <= 32
    annot = (np.where(sig, np.vectorize(lambda v: f"{v:.2f}")(rho.fillna(0).values), "")
             if annotate else False)

    fig, ax = plt.subplots(figsize=(max(9.0, 0.46 * n_c + 5.0), max(7.0, 0.42 * n_r + 3.0)))
    sns.heatmap(rho, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                annot=annot, fmt="" if annotate else None,
                annot_kws={"size": 7.5} if annotate else None,
                xticklabels=[met_label(c, maxlen=28) for c in rho.columns],
                yticklabels=rho.index,
                linewidths=0.4 if n_c < 60 else 0, linecolor='white',
                cbar_kws={'shrink': 0.5, 'label': r'Spearman $\rho$',
                          'ticks': [-1, -0.5, 0, 0.5, 1]}, ax=ax)

    # Ring the FDR-significant cells so they are findable in the unannotated view.
    for i, j in zip(*np.where(sig)):
        ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='#111111', lw=1.1))

    ax.set_title(title, fontsize=13, pad=14, loc='left')
    ax.tick_params(left=False, bottom=False, labelsize=8.5 if n_c > 40 else 10)
    plt.setp(ax.get_xticklabels(), rotation=90 if n_c > 40 else 45, ha='right' if n_c <= 40 else 'center')
    plt.setp(ax.get_yticklabels(), rotation=0)
    for s in ax.spines.values():
        s.set_visible(False)

    note = (f"Rectangular matrix (social x metabolite): asymmetric, so no half is redundant "
            f"and every cell is shown. {int(sig.sum())}/{int(rho.notna().values.sum())} pairs "
            f"pass FDR < {fdr_alpha:g} (outlined). Rows and columns ordered independently by "
            f"hierarchical clustering. ")
    fig.text(0.01, 0.008, note + footnote, fontsize=8.5, color='#444444', va='bottom')

    plt.tight_layout(rect=[0, 0.045, 1, 1])
    plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  {filename}: {n_r}x{n_c}, {int(sig.sum())} pairs FDR<{fdr_alpha:g}")


def plot_cross_corr_lines(df_sub, rho, q, label, suffix, n_top=20):
    """Line grids for the most positive and most negative social x metabolite pairs.

    Both series are z-scored onto a single shared axis rather than plotted on twin
    y-axes: with two independent y-scales any pair can be made to look coupled or
    decoupled by the choice of limits, so a shared standardised axis is the only
    version where the visual agrees with the coefficient. Samples are ordered by
    the social indicator so a monotone relationship reads as a trend.
    """
    flat = (rho.stack(dropna=True).rename('rho').reset_index()
               .rename(columns={'level_0': 'soc', 'level_1': 'met'}))
    if flat.empty:
        return
    flat['q'] = [q.loc[r.soc, r.met] for r in flat.itertuples()]

    def grid(pairs, title, filename):
        if pairs.empty:
            return
        rows = int(np.ceil(len(pairs) / 5))
        fig, axes = plt.subplots(rows, 5, figsize=(25, 3.9 * rows), squeeze=False)
        axes = axes.flatten()
        for idx, row in pairs.reset_index(drop=True).iterrows():
            ax = axes[idx]
            x = pd.to_numeric(df_sub[row.soc], errors='coerce')
            y = np.log1p(pd.to_numeric(df_sub[row.met], errors='coerce'))
            d = pd.DataFrame({'x': x, 'y': y}).dropna().sort_values('x')
            if len(d) < 3:
                ax.axis('off')
                continue
            z = lambda s: (s - s.mean()) / (s.std(ddof=0) or 1.0)
            k = range(len(d))
            ax.plot(k, z(d.x), color='#1f4e79', lw=1.6, label=row.soc)
            ax.plot(k, z(d.y), color='#c0504d', lw=1.6, alpha=0.85,
                    label=met_label(row.met, maxlen=22))
            ax.axhline(0, color='#cccccc', lw=0.7, zorder=0)
            ax.set_title(f"rho = {row.rho:+.2f}   q = {row.q:.2g}   (n = {len(d)})", fontsize=10)
            ax.set_xticks([])
            ax.set_ylabel("z-score", fontsize=8)
            ax.legend(fontsize=7.5, loc='upper left', frameon=False)
            for s in ('top', 'right'):
                ax.spines[s].set_visible(False)
        for i in range(len(pairs), len(axes)):
            axes[i].axis('off')
        plt.suptitle(title, fontsize=16)
        fig.text(0.01, 0.005,
                 "Both series z-scored on one shared axis (not twin y-axes). Samples ordered "
                 "along x by the social indicator, so a monotone association reads as a trend.",
                 fontsize=9, color='#444444', va='bottom')
        plt.tight_layout(rect=[0, 0.02, 1, 0.965])
        plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  {filename}: {len(pairs)} pairs")

    grid(flat.nlargest(n_top, 'rho'),
         f"Top {n_top} positive social-metabolite associations ({label})",
         f'line_social_metabolite_corr_positive_{suffix}.png')
    grid(flat.nsmallest(n_top, 'rho'),
         f"Top {n_top} negative social-metabolite associations ({label})",
         f'line_social_metabolite_corr_negative_{suffix}.png')


def analyze_social_metabolite_correlations(df, social_cols, met_cols, cohort='No'):
    """Spearman correlation between social indicators and metabolites.

    cohort: 'No' = controls, 'CANCER' = cases.
    """
    label = 'Cases' if cohort == 'CANCER' else 'Controls'
    suffix = 'cases' if cohort == 'CANCER' else 'controls'
    print(f"Testing social vs metabolite correlations in {label}...")
    df_control = df[df['Cohort'] == cohort].copy()

    # FDR is computed once over ALL social x metabolite tests, then the panel is
    # narrowed for legibility -- so the q-values shown are not inflated by the
    # fact that only the strongest metabolites are displayed.
    rho_mat, q_mat = cross_corr_matrices(df_control, social_cols, met_cols)
    if not rho_mat.empty:
        n_top_met = 30
        strength = rho_mat.abs().max(axis=0).sort_values(ascending=False)
        keep = list(strength.head(n_top_met).index)
        plot_cross_corr_heatmap(
            rho_mat[keep], q_mat[keep],
            f"Social indicators x top {len(keep)} metabolites by max |rho| "
            f"({label}, n={len(df_control)})",
            f'heatmap_social_metabolite_corr_{suffix}.png',
            footnote=f"Showing the {len(keep)} of {rho_mat.shape[1]} metabolites with the "
                     f"largest |rho| against any indicator (range "
                     f"{strength.iloc[len(keep)-1]:.2f}-{strength.iloc[0]:.2f}); the rest are "
                     f"weaker still. Metabolites selected on the same coefficients displayed, "
                     f"so these are upper bounds, not typical values. FDR computed over all "
                     f"{int(rho_mat.notna().values.sum())} tests, not just those shown.")

        plot_cross_corr_lines(df_control, rho_mat, q_mat, label, suffix)

    results = []
    for soc in social_cols:
        if soc not in df_control.columns: continue
        soc_vals = pd.to_numeric(df_control[soc], errors='coerce')
        if soc_vals.notna().sum() < 5 or soc_vals.dropna().nunique() < 2: 
            continue
        
        for met in met_cols:
            if met not in df_control.columns: continue
            met_vals = np.log1p(pd.to_numeric(df_control[met], errors='coerce'))
            
            mask = soc_vals.notna() & met_vals.notna()
            if mask.sum() < 5: continue
            
            rho, pval = spearmanr(soc_vals[mask], met_vals[mask])
            results.append({
                'Social_Indicator': soc,
                'Metabolite': met,
                'Spearman_Rho': rho,
                'P_Value': pval
            })
            
    res_df = pd.DataFrame(results)
    if res_df.empty:
        print("No valid correlations found.")
        return
        
    res_df = res_df.dropna(subset=['P_Value'])
    res_df['FDR'] = fdrcorrection(res_df['P_Value'])[1]
    
    res_df['Abs_Rho'] = res_df['Spearman_Rho'].abs()
    res_df = res_df.sort_values(by=['P_Value', 'FDR'], ascending=[True, True])
    
    print("\n" + "="*80)
    print(f"TOP 15 SOCIAL-METABOLITE CORRELATIONS ({label.upper()} ONLY)")
    print("="*80)
    print(res_df.head(15).drop(columns=['Abs_Rho']).to_markdown(index=False))
    print("="*80 + "\n")
    
    # Visualization
    top_9 = res_df.head(9)
    if not top_9.empty:
        fig, axes = plt.subplots(3, 3, figsize=(15, 15))
        axes = axes.flatten()
        
        for idx, row in top_9.reset_index(drop=True).iterrows():
            soc = row['Social_Indicator']
            met = row['Metabolite']
            rho = row['Spearman_Rho']
            p_val = row['P_Value']
            fdr = row['FDR']
            
            x = pd.to_numeric(df_control[soc], errors='coerce')
            y = np.log1p(pd.to_numeric(df_control[met], errors='coerce'))
            plot_df = pd.DataFrame({soc: x, 'Metabolite (log1p)': y}).dropna()
            
            sns.regplot(data=plot_df, x=soc, y='Metabolite (log1p)', ax=axes[idx], color='tab:orange', scatter_kws={'alpha': 0.7})
            
            clean_met = met_label(met, maxlen=30)
            
            axes[idx].set_title(f"Rho = {rho:.2f} | p = {p_val:.3e}", fontsize=11)
            axes[idx].set_ylabel(clean_met)
            
        for i in range(len(top_9), len(axes)):
            axes[i].axis('off')
            
        plt.tight_layout()
        plt.savefig(FIG_DIR / f'scatter_social_metabolite_correlations_{suffix}.png', dpi=300)
        plt.close()

def effect_size_stats(a_raw, b_raw):
    """Distribution-free effect sizes matching a rank-sum test of A against B.

    Cliff's delta (identical to the rank-biserial correlation) is the natural effect
    size for Wilcoxon rank-sum: it is exactly the statistic the test ranks on, bounded
    in [-1, 1], and unlike a fold change it cannot blow up when a median lands on a
    detection floor. AUC is the same quantity rescaled to [0, 1] -- the probability
    that a randomly drawn A exceeds a randomly drawn B. Hodges-Lehmann is the median
    of all pairwise A-B differences: the location shift the rank-sum test is actually
    estimating, in native units, and far more stable at small n than a median gap.
    """
    a = pd.to_numeric(pd.Series(a_raw), errors='coerce').dropna().values
    b = pd.to_numeric(pd.Series(b_raw), errors='coerce').dropna().values
    out = {'Cliffs_delta': np.nan, 'AUC': np.nan, 'Hodges_Lehmann': np.nan}
    if len(a) == 0 or len(b) == 0:
        return out
    diff = a[:, None] - b[None, :]
    out['Cliffs_delta'] = float((np.sum(diff > 0) - np.sum(diff < 0)) / diff.size)
    out['AUC'] = float((np.sum(diff > 0) + 0.5 * np.sum(diff == 0)) / diff.size)
    out['Hodges_Lehmann'] = float(np.median(diff))
    return out


def detection_test(a_raw, b_raw, detect_min=0.0):
    """Fisher exact test on detection rate (value > detect_min) between two groups.

    For features where a sizeable share of samples sit on the assay floor, part of any
    rank-sum signal is 'detected vs not detected' rather than 'more vs less'. This
    separates the two so they can be reported apart.
    """
    a = pd.to_numeric(pd.Series(a_raw), errors='coerce').dropna().values
    b = pd.to_numeric(pd.Series(b_raw), errors='coerce').dropna().values
    out = {'P_detect': np.nan, 'OR_detect': np.nan}
    if len(a) == 0 or len(b) == 0:
        return out
    da, db = int((a > detect_min).sum()), int((b > detect_min).sum())
    if (da + db) == 0 or (da == len(a) and db == len(b)):
        return out                      # no variation in detection: nothing to test
    odds, p = fisher_exact([[da, len(a) - da], [db, len(b) - db]])
    out['P_detect'], out['OR_detect'] = float(p), float(odds)
    return out


def fold_change_stats(a_raw, b_raw, eps=None, zero_max=ZERO_MAX, n_boot=0, seed=0,
                      as_ratio=True, detect_min=0.0):
    """Median fold change of group A over group B, with detection-aware suppression."""
    if eps is None:
        eps = LOG2FC_EPS
    """
    A zero is a non-detection, not a measured zero. Once a group is mostly
    non-detected its median sits on the detection floor and the ratio measures the
    floor, not biology -- the same reasoning that makes rank_sum_test route those
    features to a Fisher test. So when either group is >= zero_max zeros the fold
    change is returned as NaN and the caller shows detection counts instead.

    n_boot > 0 adds a percentile bootstrap CI on the fold change (resampling both
    groups independently and recomputing the median ratio).

    as_ratio=False is for features where a ratio is not the right effect measure --
    methylation (values 0-5, where the eps=2.96 pseudocount would dominate and force
    every FC to ~1.0) and clinical variables like age or BMI (a ratio of ages is not
    an interpretable quantity). Those report an absolute median difference instead.
    """
    a = pd.to_numeric(pd.Series(a_raw), errors='coerce').dropna().values
    b = pd.to_numeric(pd.Series(b_raw), errors='coerce').dropna().values
    out = {'N_A': len(a), 'N_B': len(b),
           'Detect_A': int((a > detect_min).sum()), 'Detect_B': int((b > detect_min).sum()),
           'Median_A': np.nan, 'Median_B': np.nan, 'Median_Diff': np.nan,
           'log2FC': np.nan, 'FC': np.nan, 'FC_CI_low': np.nan, 'FC_CI_high': np.nan,
           'FC_label': 'FC n/a'}
    if len(a) == 0 or len(b) == 0:
        return out

    out['Median_A'], out['Median_B'] = float(np.median(a)), float(np.median(b))
    out['Median_Diff'] = out['Median_A'] - out['Median_B']

    if not as_ratio:
        # Ratio is not meaningful for this feature type; report the absolute shift.
        out['FC_label'] = f"Δ median = {out['Median_Diff']:+.3g}"
        return out

    z_a = 1.0 - out['Detect_A'] / len(a)
    z_b = 1.0 - out['Detect_B'] / len(b)

    if z_a >= zero_max or z_b >= zero_max:
        out['FC_label'] = (f"FC n/a (detected {out['Detect_A']}/{len(a)} "
                           f"vs {out['Detect_B']}/{len(b)})")
        return out

    out['log2FC'] = float(np.log2((out['Median_A'] + eps) / (out['Median_B'] + eps)))
    out['FC'] = float(2 ** out['log2FC'])

    label = f"FC {out['FC']:.2f}x"
    if n_boot:
        rng = np.random.default_rng(seed)
        bs = np.empty(n_boot)
        for i in range(n_boot):
            ma = np.median(rng.choice(a, len(a)))
            mb = np.median(rng.choice(b, len(b)))
            bs[i] = (ma + eps) / (mb + eps)
        out['FC_CI_low'], out['FC_CI_high'] = (float(x) for x in np.percentile(bs, [2.5, 97.5]))
        label += f" [{out['FC_CI_low']:.2f}-{out['FC_CI_high']:.2f}]"
    out['FC_label'] = label

    # A non-trivial minority of non-detections still distorts a median ratio, so
    # say so rather than letting the number stand unqualified.
    if max(z_a, z_b) > 0.2:
        out['FC_label'] += f" (nondetect {out['Detect_A']}/{len(a)} vs {out['Detect_B']}/{len(b)})"
    return out


def plot_combined_boxplot_panel(df_src, feat_p, group_col, group_order, palette,
                                title, filename, n_cols=4, value_label=None,
                                footnote=None, log_transform=True, log_axis=False,
                                nd_floor=1e-6):
    """One figure, one boxplot panel per feature, ordered by ascending p-value.

    feat_p is a list of (feature, p_value) or (feature, p_value, fc_label) already
    filtered and sorted by the caller; df_src must carry group_col plus every
    feature column. When present, fc_label is shown under the p-value.

    log_transform=False plots the values as they are. Use it for feature types whose
    native range is already readable -- methylation runs 0-6, where log1p is the
    identity near zero and only squashes the tail. It also keeps the axis on the same
    scale as the Δ median annotation, which fold_change_stats computes on raw values.
    The rank-sum p-values are unaffected either way: log1p is monotone, so the ranks,
    and therefore the p-values, are identical.

    log_axis=True keeps the plotted values raw but draws the y-axis on a log scale, so
    the tick labels stay in native units while a right-skewed bulk near zero is still
    legible. Non-positive values cannot be drawn on a log axis; if any are present the
    axis falls back to linear rather than silently dropping points.
    """
    if value_label is None:
        if log_transform:
            value_label = 'Value (log1p)'
        else:
            value_label = 'Value (log scale)' if log_axis else 'Value'
    if not feat_p:
        print(f"  {filename}: no features meet the threshold — not plotted.")
        return

    n_feats = len(feat_p)
    n_cols  = min(n_cols, n_feats)
    n_rows  = int(np.ceil(n_feats / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 3.6, n_rows * 3.8), squeeze=False)
    axes_flat = axes.flatten()

    for idx, rec in enumerate(feat_p):
        feat, p_val = rec[0], rec[1]
        fc_label = rec[2] if len(rec) > 2 else None
        ax = axes_flat[idx]
        vals = pd.to_numeric(df_src[feat], errors='coerce')
        pdf = pd.DataFrame({
            'Value': np.log1p(vals) if log_transform else vals,
            'Group': df_src[group_col],
        }).dropna()
        if pdf.empty:
            ax.axis('off')
            continue
        sns.boxplot(data=pdf, x='Group', y='Value', order=group_order,
                    palette=palette, showfliers=False, ax=ax)
        sns.stripplot(data=pdf, x='Group', y='Value', order=group_order,
                      color='black', alpha=0.5, size=3, ax=ax)
        # NB: not named `title` -- that is the figure-level parameter used by suptitle below.
        panel_title = f"{met_label(feat, maxlen=30)}\np = {p_val:.3g}"
        if fc_label:
            panel_title += f"\n{fc_label}"
        ax.set_title(panel_title, fontsize=8.5, fontweight='bold')
        if log_axis:
            # The methylation assay reports non-detections as ~1e-22 rather than 0, so a
            # bare log axis spans ~22 decades and squashes the real 1e-3..5 range into a
            # sliver. Anchor the axis just below the smallest genuinely detected value and
            # state how many points sit on the floor below it.
            real = pdf.loc[pdf['Value'] > nd_floor, 'Value']
            n_nd = int((pdf['Value'] <= nd_floor).sum())
            if len(real):
                ax.set_yscale('log')
                ax.set_ylim(bottom=real.min() * 0.5)
                if n_nd:
                    ax.text(0.5, 0.02, f'{n_nd}/{len(pdf)} below detection (not shown)',
                            transform=ax.transAxes, ha='center', va='bottom',
                            fontsize=7, color='#a03030')
            else:
                print(f"    [warn] {filename}: {met_label(feat)} has no value above the "
                      f"detection floor; left on a linear axis.")
        ax.set_xlabel('')
        ax.set_ylabel(value_label if idx % n_cols == 0 else '', fontsize=8)
        ax.tick_params(labelsize=7.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for i in range(n_feats, len(axes_flat)):
        axes_flat[i].axis('off')

    fig.suptitle(title, fontsize=12, fontweight='bold', y=1.005)
    plt.tight_layout()
    if footnote:
        fig.text(0.01, -0.012, footnote, fontsize=8.5, color='#444444', va='top')
    plt.savefig(FIG_DIR / filename, dpi=250, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}  ({n_feats} features)")


def analyze_survival_boxplots(df, features, met_cols=None, meth_cols=None,
                              panel_p_max=0.1, n_boot=2000):
    """Per-feature boxplots for three splits, plus a combined panel per split.

    Splits: dead low- vs high-survival (median TTE among the deceased), dead vs
    alive, and early (Stage I-II) vs advanced (Stage III-IV) stage. Each writes a
    folder of single-feature figures; the combined panels are restricted to
    met_cols (metabolites) at raw p <= panel_p_max.
    """
    print(f"Generating median-split survival boxplots for {len(features)} features...")
    
    dir_dd = FIG_DIR / 'survival_boxplots_dead_vs_dead'
    dir_da = FIG_DIR / 'survival_boxplots_dead_vs_alive'
    dir_ea = FIG_DIR / 'survival_boxplots_early_vs_advanced'
    dir_dd.mkdir(exist_ok=True)
    dir_da.mkdir(exist_ok=True)
    dir_ea.mkdir(exist_ok=True)

    df_cancer = df[df['Cohort'] == 'CANCER'].copy()

    df_dead = df_cancer[df_cancer['Event'] == 1].copy()
    df_alive = df_cancer[df_cancer['Event'] == 0].copy()

    med_dead = df_dead['Time to Event'].median()

    df_dead['Survival_Group'] = np.where(df_dead['Time to Event'] <= med_dead, f'Dead (<= {med_dead:.0f}d)', f'Dead (> {med_dead:.0f}d)')
    df_alive['Survival_Group'] = 'Alive'

    df_plot_base = pd.concat([df_dead, df_alive])
    df_plot_base['Status'] = df_plot_base['Event'].map({1: 'Dead', 0: 'Alive'})

    # Early (Stage I-II) vs advanced (Stage III-IV); unknown stage is dropped by the helper.
    df_stage = _stage_early_late(df)
    df_stage['StageLabel'] = df_stage['StageGroup'].map({'Early': 'Early (I-II)',
                                                        'Late':  'Advanced (III-IV)'})
    n_early = int((df_stage['StageGroup'] == 'Early').sum())
    n_adv   = int((df_stage['StageGroup'] == 'Late').sum())
    print(f"  Stage groups: Early (I-II) n={n_early} | Advanced (III-IV) n={n_adv}")

    _met_set  = set(met_cols)  if met_cols  is not None else set()
    _meth_set = set(meth_cols) if meth_cols is not None else set()

    group_order_dd = [f'Dead (<= {med_dead:.0f}d)', f'Dead (> {med_dead:.0f}d)']
    group_order_da = ['Dead', 'Alive']
    group_order_ea = ['Early (I-II)', 'Advanced (III-IV)']

    results = []
    
    for f in features:
        if f not in df_cancer.columns: continue
        
        vals_dead_low = pd.to_numeric(df_dead.loc[df_dead['Time to Event'] <= med_dead, f], errors='coerce').dropna()
        vals_dead_high = pd.to_numeric(df_dead.loc[df_dead['Time to Event'] > med_dead, f], errors='coerce').dropna()
        vals_alive = pd.to_numeric(df_alive[f], errors='coerce').dropna()
        vals_dead_all = pd.to_numeric(df_dead[f], errors='coerce').dropna()
        
        vals_dead_low_log = np.log1p(vals_dead_low)
        vals_dead_high_log = np.log1p(vals_dead_high)
        vals_alive_log = np.log1p(vals_alive)
        vals_dead_all_log = np.log1p(vals_dead_all)
        
        p_dead = np.nan
        if len(vals_dead_low) > 3 and len(vals_dead_high) > 3:
            _, p_dead = ranksums(vals_dead_low_log, vals_dead_high_log)
            
        p_dead_vs_alive = np.nan
        if len(vals_dead_all) > 3 and len(vals_alive) > 3:
            _, p_dead_vs_alive = ranksums(vals_dead_all_log, vals_alive_log)

        vals_early = pd.to_numeric(df_stage.loc[df_stage['StageGroup'] == 'Early', f], errors='coerce').dropna()
        vals_adv   = pd.to_numeric(df_stage.loc[df_stage['StageGroup'] == 'Late', f], errors='coerce').dropna()
        p_early_vs_adv = np.nan
        if len(vals_early) > 3 and len(vals_adv) > 3:
            _, p_early_vs_adv = ranksums(np.log1p(vals_early), np.log1p(vals_adv))

        # Fold change per split. A is always the first-named group, so the sign
        # convention matches the boxplot's left-to-right order.
        # Only metabolites get a ratio; methylation and clinical variables report an
        # absolute median shift (see fold_change_stats). Methylation also needs a
        # non-zero detection floor: its non-detections are recorded as ~1e-22, not 0,
        # so a `> 0` test would call every sample detected.
        as_ratio = (met_cols is None) or (f in _met_set)
        dmin = METH_DETECT_MIN if f in _meth_set else 0.0

        arms = {
            'Dead_vs_Dead':      (df_dead.loc[df_dead['Time to Event'] <= med_dead, f],
                                  df_dead.loc[df_dead['Time to Event'] >  med_dead, f]),
            'Dead_vs_Alive':     (df_dead[f], df_alive[f]),
            'Early_vs_Advanced': (df_stage.loc[df_stage['StageGroup'] == 'Early', f],
                                  df_stage.loc[df_stage['StageGroup'] == 'Late', f]),
        }
        row = {
            'Feature': f,
            'P_Dead_vs_Dead': p_dead,
            'P_Dead_vs_Alive': p_dead_vs_alive,
            'P_Early_vs_Advanced': p_early_vs_adv,
        }
        fc_by_tag = {}
        for tag, (va, vb) in arms.items():
            fc = fold_change_stats(va, vb, n_boot=n_boot, as_ratio=as_ratio, detect_min=dmin)
            fc.update(effect_size_stats(va, vb))
            fc.update(detection_test(va, vb, detect_min=dmin))
            fc_by_tag[tag] = fc
            for k in ('FC', 'log2FC', 'FC_CI_low', 'FC_CI_high',
                      'Median_A', 'Median_B', 'Median_Diff',
                      'Cliffs_delta', 'AUC', 'Hodges_Lehmann',
                      'P_detect', 'OR_detect',
                      'Detect_A', 'Detect_B', 'N_A', 'N_B'):
                row[f'{k}_{tag}'] = fc[k]
        fc_dd, fc_da, fc_ea = (fc_by_tag['Dead_vs_Dead'], fc_by_tag['Dead_vs_Alive'],
                               fc_by_tag['Early_vs_Advanced'])
        results.append(row)
        
        clean_filename = str(f).replace('/', '_').replace(' ', '_').replace('\n', '')
        
        # PLOT 1: Dead vs Dead
        plot_dd = pd.DataFrame({
            'Value': np.log1p(pd.to_numeric(df_dead[f], errors='coerce')),
            'Survival_Group': df_dead['Survival_Group']
        }).dropna()
        
        if not plot_dd.empty:
            plt.figure(figsize=(8, 6))
            group_order_dd = [f'Dead (<= {med_dead:.0f}d)', f'Dead (> {med_dead:.0f}d)']
            sns.boxplot(data=plot_dd, x='Survival_Group', y='Value', 
                        palette=['#d95f02', '#1b9e77'], showfliers=False, order=group_order_dd)
            sns.stripplot(data=plot_dd, x='Survival_Group', y='Value', 
                          dodge=True, color='black', alpha=0.5, size=4, order=group_order_dd)
            plt.title(f"{f}\nDead Low vs High Survival (P-value: {p_dead:.3e})\n{fc_dd['FC_label']}",
                      fontsize=11)
            plt.ylabel("Value (log1p)")
            plt.xlabel("")
            plt.tight_layout()
            plt.savefig(dir_dd / f"{clean_filename}.png", dpi=300)
            plt.close()
            
        # PLOT 2: Dead vs Alive
        plot_da = pd.DataFrame({
            'Value': np.log1p(pd.to_numeric(df_plot_base[f], errors='coerce')),
            'Status': df_plot_base['Status']
        }).dropna()
        
        if not plot_da.empty:
            plt.figure(figsize=(8, 6))
            sns.boxplot(data=plot_da, x='Status', y='Value',
                        palette=['#e7298a', '#7570b3'], showfliers=False, order=group_order_da)
            sns.stripplot(data=plot_da, x='Status', y='Value',
                          dodge=True, color='black', alpha=0.5, size=4, order=group_order_da)
            plt.title(f"{f}\nDead vs Alive (P-value: {p_dead_vs_alive:.3e})\n{fc_da['FC_label']}",
                      fontsize=11)
            plt.ylabel("Value (log1p)")
            plt.xlabel("")
            plt.tight_layout()
            plt.savefig(dir_da / f"{clean_filename}.png", dpi=300)
            plt.close()

        # PLOT 3: Early (Stage I-II) vs Advanced (Stage III-IV)
        plot_ea = pd.DataFrame({
            'Value': np.log1p(pd.to_numeric(df_stage[f], errors='coerce')),
            'StageLabel': df_stage['StageLabel']
        }).dropna()

        if not plot_ea.empty:
            plt.figure(figsize=(8, 6))
            sns.boxplot(data=plot_ea, x='StageLabel', y='Value',
                        palette=['#66a61e', '#a6761d'], showfliers=False, order=group_order_ea)
            sns.stripplot(data=plot_ea, x='StageLabel', y='Value',
                          dodge=True, color='black', alpha=0.5, size=4, order=group_order_ea)
            plt.title(f"{f}\nEarly (I-II) vs Advanced (III-IV) (P-value: {p_early_vs_adv:.3e})"
                      f"\n{fc_ea['FC_label']}", fontsize=11)
            plt.ylabel("Value (log1p)")
            plt.xlabel("")
            plt.tight_layout()
            plt.savefig(dir_ea / f"{clean_filename}.png", dpi=300)
            plt.close()

    res_df = pd.DataFrame(results)
    if not res_df.empty:
        # Rank by Dead survival Time
        res_df_valid = res_df.dropna(subset=['P_Dead_vs_Dead']).copy()
        if not res_df_valid.empty:
            res_df_valid['FDR_Dead_vs_Dead'] = fdrcorrection(res_df_valid['P_Dead_vs_Dead'])[1]
            top_20 = res_df_valid.sort_values('P_Dead_vs_Dead').head(20)
            
            print("\n" + "="*80)
            print("TOP 20 FEATURES: DEAD (LOW SURVIVAL) vs DEAD (HIGH SURVIVAL)")
            print("="*80)
            print(top_20.to_markdown(index=False))
            print("="*80 + "\n")
            
        res_df_da = res_df.dropna(subset=['P_Dead_vs_Alive']).copy()
        if not res_df_da.empty:
            res_df_da['FDR_Dead_vs_Alive'] = fdrcorrection(res_df_da['P_Dead_vs_Alive'])[1]
            top_20_da = res_df_da.sort_values('P_Dead_vs_Alive').head(20)
            
            print("\n" + "="*80)
            print("TOP 20 FEATURES: DEAD vs ALIVE")
            print("="*80)
            print(top_20_da.to_markdown(index=False))
            print("="*80 + "\n")

        res_df_ea = res_df.dropna(subset=['P_Early_vs_Advanced']).copy()
        if not res_df_ea.empty:
            res_df_ea['FDR_Early_vs_Advanced'] = fdrcorrection(res_df_ea['P_Early_vs_Advanced'])[1]
            top_20_ea = res_df_ea.sort_values('P_Early_vs_Advanced').head(20)

            print("\n" + "="*80)
            print("TOP 20 FEATURES: EARLY (STAGE I-II) vs ADVANCED (STAGE III-IV)")
            print("="*80)
            print(top_20_ea.to_markdown(index=False))
            print("="*80 + "\n")

        # Carry the FDR columns back onto the returned frame so the Excel sheet has them.
        for sub, col in ((res_df_valid, 'FDR_Dead_vs_Dead'),
                         (res_df_da,    'FDR_Dead_vs_Alive'),
                         (res_df_ea,    'FDR_Early_vs_Advanced')):
            if not sub.empty and col in sub.columns:
                res_df = res_df.merge(sub[['Feature', col]], on='Feature', how='left')

    # ── Combined panels, raw p <= panel_p_max ────────────────────────────────
    # Metabolites carry a fold change; methylation carries an absolute median shift
    # (see fold_change_stats), so the effect-size label differs between the two.
    def _emit_panels(feature_set, kind_label, specs, log_transform=True, log_axis=False):
        n_feat   = int(res_df['Feature'].isin(feature_set).sum())
        expected = panel_p_max * n_feat
        note = (f"Raw rank-sum p-values, not FDR-corrected. {n_feat} {kind_label} screened; "
                f"~{expected:.1f} would be expected below p = {panel_p_max} by chance alone.")
        for pcol, tag, src, gcol, gorder, pal, base_title, fname in specs:
            sel = res_df[res_df['Feature'].isin(feature_set) & (res_df[pcol] <= panel_p_max)]
            sel = sel.sort_values(pcol)
            feat_p = []
            for _, rw in sel.iterrows():
                fc, lo, hi = rw[f'FC_{tag}'], rw[f'FC_CI_low_{tag}'], rw[f'FC_CI_high_{tag}']
                if pd.isna(fc):
                    md = rw.get(f'Median_Diff_{tag}', np.nan)
                    lbl = (f"Δ median = {md:+.4g}" if pd.notna(md) else "effect n/a")
                else:
                    lbl = f"FC {fc:.2f}x"
                    if pd.notna(lo):
                        lbl += f" [{lo:.2f}-{hi:.2f}]"
                cd = rw.get(f'Cliffs_delta_{tag}', np.nan)
                if pd.notna(cd):
                    lbl += f"   δ = {cd:+.2f}"
                # Detection counts and the Fisher p stay in Survival_Boxplot_Stats
                # (Detect_A/B, P_detect, OR_detect) rather than on the panel.
                feat_p.append((rw['Feature'], rw[pcol], lbl))
            eff = ('FC = ' + gorder[0] + ' / ' + gorder[1]) if feature_set is met_set \
                  else ('Δ = ' + gorder[0] + ' − ' + gorder[1])
            plot_combined_boxplot_panel(
                src, feat_p, gcol, gorder, pal,
                f'{base_title}\n{len(feat_p)} {kind_label} at raw p ≤ {panel_p_max}   ({eff})',
                fname, n_cols=4, footnote=note,
                log_transform=log_transform, log_axis=log_axis)

    met_set = _met_set

    if met_cols is not None and not res_df.empty:
        _emit_panels(met_set, 'metabolites', [
            ('P_Dead_vs_Dead', 'Dead_vs_Dead', df_dead, 'Survival_Group', group_order_dd,
             ['#d95f02', '#1b9e77'],
             f'Metabolites: deceased patients, low vs high survival (split at median {med_dead:.0f} d)',
             'boxplot_combined_dead_low_vs_high_survival_p01.png'),
            ('P_Dead_vs_Alive', 'Dead_vs_Alive', df_plot_base, 'Status', group_order_da,
             ['#e7298a', '#7570b3'],
             'Metabolites: deceased vs alive cancer patients',
             'boxplot_combined_dead_vs_alive_p01.png'),
            ('P_Early_vs_Advanced', 'Early_vs_Advanced', df_stage, 'StageLabel', group_order_ea,
             ['#66a61e', '#a6761d'],
             f'Metabolites: early (Stage I-II, n={n_early}) vs advanced (Stage III-IV, n={n_adv})',
             'boxplot_combined_early_vs_advanced_p01.png'),
        ])

    if meth_cols is not None and not res_df.empty:
        _emit_panels(set(meth_cols), 'DNAm markers', [
            ('P_Dead_vs_Dead', 'Dead_vs_Dead', df_dead, 'Survival_Group', group_order_dd,
             ['#d95f02', '#1b9e77'],
             f'DNA methylation: deceased patients, low vs high survival '
             f'(split at median {med_dead:.0f} d)',
             'boxplot_combined_dnam_dead_low_vs_high_survival_p01.png'),
            ('P_Dead_vs_Alive', 'Dead_vs_Alive', df_plot_base, 'Status', group_order_da,
             ['#e7298a', '#7570b3'],
             'DNA methylation: deceased vs alive cancer patients',
             'boxplot_combined_dnam_dead_vs_alive_p01.png'),
            ('P_Early_vs_Advanced', 'Early_vs_Advanced', df_stage, 'StageLabel', group_order_ea,
             ['#66a61e', '#a6761d'],
             f'DNA methylation: early (Stage I-II, n={n_early}) vs advanced (Stage III-IV, n={n_adv})',
             'boxplot_combined_dnam_early_vs_advanced_p01.png'),
        ], log_transform=False, log_axis=True)

    return res_df

def analyze_group_metabolites(d, met_cols, group_col, group_order, title, file_stem,
                              folder_name, n_boot=2000, p_max=0.1, panel_p_max=0.05,
                              palette=('#7570b3', '#d95f02'), per_feature_plots=True,
                              as_ratio=True, detect_min=0.0, log_transform=True,
                              log_axis=False, make_volcano=True, kind_label='metabolites',
                              total_signal_note=True):
    """Two-group feature comparison: volcano, per-feature boxplots, combined panel.

    `d` must already be subset to the patients of interest and carry `group_col` with
    exactly the two values in `group_order` (group A first -- the direction convention
    for every fold change and effect size reported).

    The defaults describe metabolites. Methylation needs as_ratio=False (values 0-6, so
    the eps=2.96 pseudocount would pin every ratio to 1.0), detect_min=METH_DETECT_MIN
    (non-detections are stored as ~1e-22, not 0), log_axis=True with log_transform=False,
    and make_volcano=False (a log2FC axis is meaningless without a ratio).
    """
    g1, g2 = group_order
    n1 = int((d[group_col] == g1).sum())
    n2 = int((d[group_col] == g2).sum())
    print("\n" + "="*70)
    print(title.upper())
    print("="*70)
    print(f"  {g1} n={n1} | {g2} n={n2}")
    if n1 < 5 or n2 < 5:
        print("  Insufficient samples -- skipping.")
        return None

    # Total measured signal often differs between groups in this dataset; when it does,
    # every fold change below carries that global offset on top of any real difference.
    ratio, p_tot = np.nan, np.nan
    if total_signal_note:
        num = d[met_cols].apply(pd.to_numeric, errors='coerce')
        tot = num.sum(axis=1)
        ta, tb = tot[d[group_col] == g1], tot[d[group_col] == g2]
        ratio = float(tb.median() / ta.median()) if ta.median() else np.nan
        _, p_tot = ranksums(ta, tb)
        print(f"  [note] total ion signal {g2}/{g1} = {ratio:.2f} (rank-sum p={p_tot:.3g}); "
              f"fold changes inherit this offset.")

    results = rank_sum_test(d, met_cols, is_log=True, group_col=group_col,
                            group1=g1, group2=g2)
    if results is None or results.empty:
        print("  No results.")
        return None

    print(f"  FDR<{p_max}: {int((results['FDR'] < p_max).sum())} | "
          f"FDR<0.05: {int((results['FDR'] < 0.05).sum())} | "
          f"raw p<={p_max}: {int((results['p_value'] <= p_max).sum())} "
          f"(~{p_max*len(results):.0f} expected by chance) | "
          f"min FDR={results['FDR'].min():.3f}")

    if make_volcano:
        plot_volcano(results, f"{kind_label.capitalize()} ({g1} n={n1} vs {g2} n={n2})",
                     f"volcano_metabolites_{file_stem}.png", group_a=g1, group_b=g2)
        plot_detection_tier(results,
                            f"{kind_label.capitalize()} detected at different rates ({g1} vs {g2})",
                            f"detection_metabolites_{file_stem}.png", group_a=g1, group_b=g2)

    plot_dir = FIG_DIR / folder_name
    if per_feature_plots:
        plot_dir.mkdir(exist_ok=True)

    for _, row in results.iterrows():
        f = row['Feature']
        va, vb = d.loc[d[group_col] == g1, f], d.loc[d[group_col] == g2, f]
        fc = fold_change_stats(va, vb, n_boot=n_boot, as_ratio=as_ratio,
                               detect_min=detect_min)
        for k in ('FC', 'log2FC', 'FC_CI_low', 'FC_CI_high', 'Median_A', 'Median_B',
                  'Median_Diff', 'Detect_A', 'Detect_B', 'N_A', 'N_B'):
            results.loc[results['Feature'] == f, k] = fc[k]
        for k, v in effect_size_stats(va, vb).items():
            results.loc[results['Feature'] == f, k] = v
        for k, v in detection_test(va, vb, detect_min=detect_min).items():
            results.loc[results['Feature'] == f, k] = v

        if not per_feature_plots:
            continue
        vals = pd.to_numeric(d[f], errors='coerce')
        pdf = pd.DataFrame({'Value': np.log1p(vals) if log_transform else vals,
                            'Group': d[group_col]}).dropna()
        if pdf.empty:
            continue
        plt.figure(figsize=(7, 6))
        ax = sns.boxplot(data=pdf, x='Group', y='Value', order=group_order,
                         palette=list(palette), showfliers=False)
        sns.stripplot(data=pdf, x='Group', y='Value', order=group_order,
                      color='black', alpha=0.5, size=4)
        ylab = 'Value (log1p)' if log_transform else 'Value'
        if log_axis:
            real = pdf.loc[pdf['Value'] > detect_min, 'Value']
            if len(real):
                ax.set_yscale('log')
                ax.set_ylim(bottom=real.min() * 0.5)
                ylab = 'Value (log scale)'
        eff = 'FC' if as_ratio else 'Δ'
        plt.title(f"{met_label(f)}\n"
                  f"rank-sum p = {row['p_value']:.3g}   FDR = {row['FDR']:.3g}\n"
                  f"{fc['FC_label']}   ({eff} = {g1} / {g2})", fontsize=11)
        plt.ylabel(ylab)
        plt.xlabel('')
        plt.tight_layout()
        clean = str(f).replace('/', '_').replace(' ', '_').replace('\n', '')
        plt.savefig(plot_dir / f"{clean}.png", dpi=200)
        plt.close()

    if per_feature_plots:
        print(f"  Saved {len(results)} boxplots to {plot_dir.name}/")

    sel = results[results['p_value'] <= panel_p_max].sort_values('p_value')
    feat_p = []
    for _, rw in sel.iterrows():
        fc, lo, hi = rw.get('FC', np.nan), rw.get('FC_CI_low', np.nan), rw.get('FC_CI_high', np.nan)
        if pd.isna(fc):
            md = rw.get('Median_Diff', np.nan)
            if not as_ratio and pd.notna(md):
                lbl = f"\u0394 median = {md:+.4g}"
            else:
                lbl = (f"FC n/a (detected {int(rw['Detect_A'])}/{int(rw['N_A'])} "
                       f"vs {int(rw['Detect_B'])}/{int(rw['N_B'])})")
        else:
            lbl = f"FC {fc:.2f}x"
            if pd.notna(lo):
                lbl += f" [{lo:.2f}-{hi:.2f}]"
        cd = rw.get('Cliffs_delta', np.nan)
        if pd.notna(cd):
            lbl += f"   \u03b4 = {cd:+.2f}"
        feat_p.append((rw['Feature'], rw['p_value'], lbl))

    note = (f"Raw rank-sum p-values, not FDR-corrected (min FDR = {results['FDR'].min():.3f}). "
            f"{len(results)} {kind_label} screened; ~{panel_p_max*len(results):.1f} would be "
            f"expected below p = {panel_p_max} by chance alone.")
    if total_signal_note and pd.notna(ratio):
        note += (f" Total ion signal differs {ratio:.2f}x between the groups "
                 f"(p={p_tot:.3g}), so fold changes carry that offset.")
    eff = 'FC' if as_ratio else '\u0394'
    plot_combined_boxplot_panel(
        d, feat_p, group_col, group_order, list(palette),
        f'{kind_label.capitalize()}: {g1} (n={n1}) vs {g2} (n={n2})\n'
        f'{len(feat_p)} {kind_label} at raw p \u2264 {panel_p_max}   ({eff} = {g1} / {g2})',
        f'boxplot_combined_{file_stem}_p{str(panel_p_max).replace("0.", "0")}.png',
        n_cols=4, footnote=note, log_transform=log_transform, log_axis=log_axis)

    top = results.sort_values('p_value').head(15)
    print(f"\n  {'feature':<40s} {'p':>10s} {'FDR':>8s} {'log2FC':>8s} {'delta':>7s}")
    for _, r in top.iterrows():
        print(f"  {met_label(r['Feature'])[:40]:<40s} {r['p_value']:>10.4g} {r['FDR']:>8.3g} "
              f"{r.get('log2FC', np.nan):>8.3g} {r.get('Cliffs_delta', np.nan):>7.2f}")
    return results


def analyze_race_boxplots(df, met_cols, **kw):
    """Black vs non-Black metabolites, cancer patients.

    Uses individual-level race, not the area-level PCT-* percentages -- those measure
    neighbourhood composition and using them for a patient-race question would be an
    ecological fallacy. Asian n=7 is too small to stand alone, so the contrast is
    Black vs everyone else.
    """
    d = df[df['Cohort'] == 'CANCER'].copy()
    if 'Race' not in d.columns:
        print("  No Race column -- skipping race analysis.")
        return None
    d['RaceGroup'] = np.where(d['Race'].astype(str).str.strip().str.upper().eq('BLACK'),
                              'Black', 'Non-Black')
    return analyze_group_metabolites(
        d, met_cols, 'RaceGroup', ['Black', 'Non-Black'],
        'Black vs Non-Black metabolite analysis (cancer patients)',
        'black_vs_nonblack', 'race_metabolites_boxplots', **kw)


def analyze_race_dnam(df, meth_cols, **kw):
    """Black vs non-Black DNA methylation markers, cancer patients.

    Same contrast as analyze_race_boxplots but with the methylation-appropriate
    settings: no fold change (values 0-6 against an eps of 2.96), a 1e-6 detection
    floor, a log-scaled raw axis, and no volcano (a log2FC axis needs a ratio).
    """
    d = df[df['Cohort'] == 'CANCER'].copy()
    if 'Race' not in d.columns:
        print("  No Race column -- skipping race DNAm analysis.")
        return None
    d['RaceGroup'] = np.where(d['Race'].astype(str).str.strip().str.upper().eq('BLACK'),
                              'Black', 'Non-Black')
    return analyze_group_metabolites(
        d, meth_cols, 'RaceGroup', ['Black', 'Non-Black'],
        'Black vs Non-Black DNA methylation analysis (cancer patients)',
        'dnam_black_vs_nonblack', 'race_dnam_boxplots',
        as_ratio=False, detect_min=METH_DETECT_MIN,
        log_transform=False, log_axis=True, make_volcano=False,
        kind_label='DNAm markers', total_signal_note=False,
        palette=('#7570b3', '#d95f02'), **kw)


def analyze_adi_boxplots(df, met_cols, tail_frac=None, **kw):
    """High vs low Area Deprivation Index, cancer patients.

    tail_frac=None       -> median split, everyone assigned to a group.
    tail_frac=0.3        -> top 30% vs bottom 30%, middle 40% dropped. Contrasting the
                            extremes raises the exposure difference between the arms,
                            which can offset the loss of n; but it discards ~40% of the
                            cohort and the resulting groups are no longer representative
                            of the full ADI range, so the two splits answer slightly
                            different questions and should be reported as such.
    """
    d = df[df['Cohort'] == 'CANCER'].copy()
    if 'ADI' not in d.columns:
        print("  No ADI column -- skipping ADI analysis.")
        return None
    adi = pd.to_numeric(d['ADI'], errors='coerce')
    n_missing = int(adi.isna().sum())

    if tail_frac is None:
        med = adi.median()
        d['ADI_Group'] = np.where(adi > med, 'High ADI', 'Low ADI')
        d = d[adi.notna()].copy()
        print(f"\n  ADI median split at {med:.1f} ({n_missing} missing, excluded)")
        title = f'High vs Low ADI metabolite analysis (median split at {med:.1f}, cancer patients)'
        stem, folder = 'high_vs_low_adi', 'adi_metabolites_boxplots'
    else:
        lo_cut, hi_cut = adi.quantile(tail_frac), adi.quantile(1 - tail_frac)
        grp = pd.Series(np.nan, index=d.index, dtype=object)
        grp[adi <= lo_cut] = 'Low ADI'
        grp[adi >= hi_cut] = 'High ADI'
        d['ADI_Group'] = grp
        n_mid = int(grp.isna().sum() - n_missing)
        d = d[grp.notna()].copy()
        pct = int(round(tail_frac * 100))
        print(f"\n  ADI tails: bottom {pct}% (<= {lo_cut:.1f}) vs top {pct}% (>= {hi_cut:.1f}); "
              f"{n_mid} middle patients dropped, {n_missing} missing")
        title = (f'High vs Low ADI metabolite analysis '
                 f'(top {pct}% >= {hi_cut:.1f} vs bottom {pct}% <= {lo_cut:.1f}, cancer patients)')
        stem = f'high_vs_low_adi_top{pct}_bottom{pct}'
        folder = f'adi_metabolites_boxplots_top{pct}_bottom{pct}'

    return analyze_group_metabolites(
        d, met_cols, 'ADI_Group', ['High ADI', 'Low ADI'], title, stem, folder,
        palette=('#b2182b', '#2166ac'), **kw)


def analyze_social_boxplots(df, social_cols, met_cols, mode_max=SOCIAL_MODE_MAX):
    print("\nRunning Social median splits for all metabolites...")
    plot_dir = FIG_DIR / 'social_boxplots_significant'
    plot_dir.mkdir(exist_ok=True)
    
    df_cases = df[df['Cohort'] == 'CANCER'].copy()
    df_controls = df[df['Cohort'] == 'No'].copy()
    
    results = []
    
    for soc in social_cols:
        if soc not in df.columns: continue
        
        # Check variance in both cohorts
        soc_cases = pd.to_numeric(df_cases[soc], errors='coerce').dropna()
        soc_controls = pd.to_numeric(df_controls[soc], errors='coerce').dropna()
        
        if len(soc_cases) < 5 or len(soc_controls) < 5: continue
        
        mode_case_pct = soc_cases.value_counts(normalize=True).max()
        mode_control_pct = soc_controls.value_counts(normalize=True).max()
        
        if mode_case_pct > mode_max or mode_control_pct > mode_max:
            continue  # Ignore deeply right-tailed/degenerate parameter distributions
            
        med_case = soc_cases.median()
        med_control = soc_controls.median()
        
        cases_high_mask = pd.to_numeric(df_cases[soc], errors='coerce') > med_case
        cases_low_mask = pd.to_numeric(df_cases[soc], errors='coerce') <= med_case
        controls_high_mask = pd.to_numeric(df_controls[soc], errors='coerce') > med_control
        controls_low_mask = pd.to_numeric(df_controls[soc], errors='coerce') <= med_control
        
        for met in met_cols:
            if met not in df.columns: continue
            
            # Cases
            val_case_high = np.log1p(pd.to_numeric(df_cases.loc[cases_high_mask, met], errors='coerce')).dropna()
            val_case_low = np.log1p(pd.to_numeric(df_cases.loc[cases_low_mask, met], errors='coerce')).dropna()
            
            p_case = np.nan
            if len(val_case_high) > 3 and len(val_case_low) > 3:
                try:
                    _, p_case = ranksums(val_case_low, val_case_high)
                except Exception:
                    pass
            
            # Controls
            val_control_high = np.log1p(pd.to_numeric(df_controls.loc[controls_high_mask, met], errors='coerce')).dropna()
            val_control_low = np.log1p(pd.to_numeric(df_controls.loc[controls_low_mask, met], errors='coerce')).dropna()
            
            p_control = np.nan
            if len(val_control_high) > 3 and len(val_control_low) > 3:
                try:
                    _, p_control = ranksums(val_control_low, val_control_high)
                except Exception:
                    pass
                
            if pd.notna(p_case) or pd.notna(p_control):
                # FC direction matches both the test above and the plot order: low / high.
                # No bootstrap CI here -- this loop covers every social x metabolite pair.
                fc_case = fold_change_stats(df_cases.loc[cases_low_mask, met],
                                            df_cases.loc[cases_high_mask, met])
                fc_ctrl = fold_change_stats(df_controls.loc[controls_low_mask, met],
                                            df_controls.loc[controls_high_mask, met])
                results.append({
                    'Social_Indicator': soc,
                    'Metabolite': met,
                    'Median_Case': med_case,
                    'Median_Control': med_control,
                    'P_Cases': p_case,
                    'P_Controls': p_control,
                    'FC_Cases_Low_over_High': fc_case['FC'],
                    'log2FC_Cases_Low_over_High': fc_case['log2FC'],
                    'FC_Controls_Low_over_High': fc_ctrl['FC'],
                    'log2FC_Controls_Low_over_High': fc_ctrl['log2FC'],
                    'FC_Label_Cases': fc_case['FC_label'],
                    'FC_Label_Controls': fc_ctrl['FC_label'],
                })
                
    res_df = pd.DataFrame(results)
    if res_df.empty:
        print("No mathematically valid variance models found for categorical social splits.")
        return
        
    # Correct P_Cases and P_Controls globally
    valid_cases = res_df.dropna(subset=['P_Cases']).copy()
    valid_controls = res_df.dropna(subset=['P_Controls']).copy()
    
    if not valid_cases.empty:
        valid_cases['FDR_Cases'] = fdrcorrection(valid_cases['P_Cases'])[1]
        res_df = res_df.merge(valid_cases[['Social_Indicator', 'Metabolite', 'FDR_Cases']], on=['Social_Indicator', 'Metabolite'], how='left')
    else:
        res_df['FDR_Cases'] = np.nan
        
    if not valid_controls.empty:
        valid_controls['FDR_Controls'] = fdrcorrection(valid_controls['P_Controls'])[1]
        res_df = res_df.merge(valid_controls[['Social_Indicator', 'Metabolite', 'FDR_Controls']], on=['Social_Indicator', 'Metabolite'], how='left')
    else:
        res_df['FDR_Controls'] = np.nan
        
    res_df['FDR_Cases'] = res_df['FDR_Cases'].fillna(1.0)
    res_df['FDR_Controls'] = res_df['FDR_Controls'].fillna(1.0)
    res_df['P_Cases'] = res_df['P_Cases'].fillna(1.0)
    res_df['P_Controls'] = res_df['P_Controls'].fillna(1.0)
    
    # Filter for hits (FDR < 0.05 on either)
    sig_mask = (res_df['FDR_Cases'] < 0.05) | (res_df['FDR_Controls'] < 0.05)
    sig_df = res_df[sig_mask].copy()
    
    print(f"Found {len(sig_df)} statistically significant social-metabolite pairs (FDR < 0.05).")
    if sig_df.empty:
        print("Falling back to Top 20 nominal P-values for plotting baseline as safety net...")
        sig_df = res_df.sort_values('P_Cases').head(20)
        
    sig_df = sig_df.sort_values(by=['FDR_Cases', 'FDR_Controls'])
        
    # Standardize loop for matplotlib rendering
    for _, row in sig_df.iterrows():
        soc = row['Social_Indicator']
        met = row['Metabolite']
        med_case = row['Median_Case']
        med_control = row['Median_Control']
        p_case = row['P_Cases']
        p_control = row['P_Controls']
        
        df_plot = []
        
        # Cases High
        m1 = pd.to_numeric(df_cases[soc], errors='coerce') > med_case
        v1 = pd.to_numeric(df_cases.loc[m1, met], errors='coerce').dropna()
        for v in v1: df_plot.append({'Group': f'Cases (> {med_case:.1f})', 'Value': np.log1p(v)})
            
        # Cases Low
        m2 = pd.to_numeric(df_cases[soc], errors='coerce') <= med_case
        v2 = pd.to_numeric(df_cases.loc[m2, met], errors='coerce').dropna()
        for v in v2: df_plot.append({'Group': f'Cases (<= {med_case:.1f})', 'Value': np.log1p(v)})
            
        # Controls High
        m3 = pd.to_numeric(df_controls[soc], errors='coerce') > med_control
        v3 = pd.to_numeric(df_controls.loc[m3, met], errors='coerce').dropna()
        for v in v3: df_plot.append({'Group': f'Controls (> {med_control:.1f})', 'Value': np.log1p(v)})
            
        # Controls Low
        m4 = pd.to_numeric(df_controls[soc], errors='coerce') <= med_control
        v4 = pd.to_numeric(df_controls.loc[m4, met], errors='coerce').dropna()
        for v in v4: df_plot.append({'Group': f'Controls (<= {med_control:.1f})', 'Value': np.log1p(v)})
            
        df_p = pd.DataFrame(df_plot)
        if df_p.empty: continue
        
        plt.figure(figsize=(10, 6))
        
        group_order = [
            f'Cases (<= {med_case:.1f})', f'Cases (> {med_case:.1f})',
            f'Controls (<= {med_control:.1f})', f'Controls (> {med_control:.1f})'
        ]
        
        sns.boxplot(data=df_p, x='Group', y='Value', 
                    palette=['#d95f02', '#e7298a', '#1b9e77', '#7570b3'], showfliers=False, order=group_order)
                    
        sns.stripplot(data=df_p, x='Group', y='Value', 
                      dodge=True, color='black', alpha=0.5, size=4, order=group_order)
        
        title_str = (f"{met} vs {soc}\n"
                     f"Cases P: {p_case:.3e} | Controls P: {p_control:.3e}\n"
                     f"Cases {row.get('FC_Label_Cases', 'FC n/a')} | "
                     f"Controls {row.get('FC_Label_Controls', 'FC n/a')}   (FC = low / high)")
        plt.title(title_str, fontsize=11)
        plt.ylabel("Value (log1p)")
        plt.xlabel("")
        plt.tight_layout()
        
        clean_name = f"{soc}_{met}".replace('/', '_').replace(' ', '_').replace('\n', '')
        
        soc_dir = plot_dir / f"{soc.replace('/', '_')}"
        soc_dir.mkdir(exist_ok=True)
        
        plt.savefig(soc_dir / f"{clean_name}.png", dpi=300)
        plt.close()
        
    print(f"Top significant plots structured inside {plot_dir}")
    return res_df

def cluster_order(corr):
    """Order rows/cols by hierarchical clustering on 1 - |rho| so co-varying
    blocks sit adjacent instead of in input order."""
    if corr.shape[0] < 3:
        return list(corr.columns)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        d = 1.0 - corr.abs().fillna(0.0).values
        np.fill_diagonal(d, 0.0)
        d = (d + d.T) / 2.0
        return [corr.columns[i] for i in leaves_list(linkage(squareform(d, checks=False),
                                                             method='average'))]
    except Exception as exc:
        print(f"  (clustering failed, keeping input order: {exc})")
        return list(corr.columns)


def spearman_fdr_matrix(data):
    """BH-FDR q-value for every unique pair, as a square DataFrame.

    q is computed once over the whole upper triangle, so subsetting a panel later
    does not change the q-values displayed in it.
    """
    cols = list(data.columns)
    k = len(cols)
    pval = np.ones((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            pair = data.iloc[:, [i, j]].dropna()
            if len(pair) < 5 or pair.iloc[:, 0].nunique() < 3 or pair.iloc[:, 1].nunique() < 3:
                p = np.nan
            else:
                p = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])[1]
            pval[i, j] = pval[j, i] = p

    iu = np.triu_indices(k, 1)
    flat = pval[iu]
    ok = ~np.isnan(flat)
    qf = np.full_like(flat, np.nan)
    if ok.sum():
        qf[ok] = fdrcorrection(flat[ok])[1]
    q = np.zeros((k, k))
    q[iu] = qf
    q[(iu[1], iu[0])] = qf
    return pd.DataFrame(q, index=cols, columns=cols)


def analyze_metabolite_correlations(df, met_cols):
    print("\nCalculating global metabolite-metabolite Spearman correlations...")
    df_met = pd.DataFrame()
    for m in met_cols:
        df_met[m] = np.log1p(pd.to_numeric(df[m], errors='coerce'))

    df_met = df_met.dropna(axis=1, how='all')

    # Calculate Spearman correlation matrix
    corr_matrix = df_met.corr(method='spearman')
    q_matrix = spearman_fdr_matrix(df_met)

    # Flatten upper triangle to get unique pairs
    pairs = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)).stack().reset_index()
    pairs.columns = ['Met_1', 'Met_2', 'Rho']
    pairs['Abs_Rho'] = pairs['Rho'].abs()
    
    top_pos_pairs = pairs.sort_values(by='Rho', ascending=False).head(20)
    top_neg_pairs = pairs.sort_values(by='Rho', ascending=True).head(20)
    top_abs_pairs = pairs.sort_values(by='Abs_Rho', ascending=False).head(20)
    
    # 1. Heatmaps for specific subsets uniquely captured from the pairs
    def plot_subset_heatmap(pairs_df, title, filename):
        unique_mets = list(set(pairs_df['Met_1'].tolist() + pairs_df['Met_2'].tolist()))
        subset_corr = corr_matrix.loc[unique_mets, unique_mets]
        order = cluster_order(subset_corr)
        subset_corr = subset_corr.loc[order, order]
        qsub_full = q_matrix.loc[order, order]

        # Symmetric matrix: draw the lower triangle only (upper half and the
        # diagonal of 1.00s are redundant). The first row and last column of that
        # triangle are structurally empty, so drop them. Label only FDR-sig cells.
        subset_corr = subset_corr.iloc[1:, :-1]
        qsub = qsub_full.iloc[1:, :-1]
        clean_labels_y = [met_label(c, maxlen=30) for c in subset_corr.index]
        clean_labels = [met_label(c, maxlen=30) for c in subset_corr.columns]
        mask = np.triu(np.ones(subset_corr.shape, dtype=bool), k=1)
        sig = (qsub.values < 0.05) & subset_corr.notna().values
        annot = np.where(sig, np.vectorize(lambda v: f"{v:.2f}")(subset_corr.fillna(0).values), "")
        n_sig, n_cells = int((sig & ~mask).sum()), int((~mask).sum())

        fig, ax = plt.subplots(figsize=(15, 13))
        sns.heatmap(subset_corr, mask=mask, annot=annot, fmt="", cmap='RdBu_r',
                    center=0, vmin=-1, vmax=1, square=True,
                    linewidths=0.5, linecolor='white', annot_kws={"size": 10},
                    xticklabels=clean_labels, yticklabels=clean_labels_y,
                    cbar_kws={'shrink': 0.45, 'label': r'Spearman $\rho$',
                              'ticks': [-1, -0.5, 0, 0.5, 1]}, ax=ax)

        ax.set_title(title, fontsize=14, pad=16, loc='left')
        ax.tick_params(left=False, bottom=False)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        plt.setp(ax.get_yticklabels(), rotation=0)
        for s in ax.spines.values():
            s.set_visible(False)
        fig.text(0.01, 0.01,
                 f"Lower triangle only (matrix is symmetric). Coefficients labelled for "
                 f"{n_sig}/{n_cells} pairs at FDR < 0.05; unlabelled cells are still coloured. "
                 f"Rows/columns ordered by hierarchical clustering on 1 - |rho|. "
                 f"Metabolites selected as the top-20 pairs, so coefficients here are "
                 f"selection-biased and should not be read as typical.",
                 fontsize=8.5, color='#444444', va='bottom', wrap=True)

        plt.tight_layout(rect=[0, 0.035, 1, 1])
        plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  {filename}: {n_sig}/{n_cells} shown pairs at FDR<0.05")
        
    plot_subset_heatmap(top_pos_pairs, "Correlation Matrix (Top 20 Positive Pairs Hub)", "heatmap_metabolite_corr_top20_positive.png")
    plot_subset_heatmap(top_neg_pairs, "Correlation Matrix (Top 20 Negative Pairs Hub)", "heatmap_metabolite_corr_top20_negative.png")
    plot_subset_heatmap(top_abs_pairs, "Correlation Matrix (Top 20 Absolute Pairs Hub)", "heatmap_metabolite_corr_top20_absolute.png")
    
    # 2. Line plots tracking literal paired distributions across patients
    def plot_line_grid(pairs_df, title, filename):
        fig, axes = plt.subplots(4, 5, figsize=(25, 18))
        axes = axes.flatten()
        
        for idx, row in pairs_df.reset_index(drop=True).iterrows():
            m1 = row['Met_1']
            m2 = row['Met_2']
            rho = row['Rho']
            
            c1 = met_label(m1, maxlen=20)
            c2 = met_label(m2, maxlen=20)
            
            # Plot chronological patient order unsorted
            df_plot = df_met[[m1, m2]].dropna()
            x_vals = range(len(df_plot))
            
            ax = axes[idx]
            ax.plot(x_vals, df_plot[m1], label='Met 1', alpha=0.8, color='tab:red', linewidth=1.5)
            ax.set_ylabel(c1, color='tab:red', fontsize=8)
            ax.tick_params(axis='y', labelcolor='tab:red', labelsize=8)
            
            # Using twin Y axes natively aligns completely different orders of magnitude without standardizing
            ax2 = ax.twinx()
            ax2.plot(x_vals, df_plot[m2], label='Met 2', alpha=0.8, color='tab:blue', linewidth=1.5)
            ax2.set_ylabel(c2, color='tab:blue', fontsize=8)
            ax2.tick_params(axis='y', labelcolor='tab:blue', labelsize=8)
            
            ax.set_title(f"Rho: {rho:.2f}\n{c1} vs {c2}", fontsize=11)
            ax.set_xticks([])
            
        plt.suptitle(title, fontsize=18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.96])
        plt.savefig(FIG_DIR / filename, dpi=300)
        plt.close()
        
    plot_line_grid(top_pos_pairs, "Top 20 Positively Coupled Metabolites (Raw Log1p Line Plots)", "line_metabolite_corr_positive.png")
    plot_line_grid(top_neg_pairs, "Top 20 Negatively Coupled Metabolites (Raw Log1p Line Plots)", "line_metabolite_corr_negative.png")
    
    print(f"Saved subset heatmaps and literal line plots to {FIG_DIR}/")
    return corr_matrix

def analyze_social_correlations(df, social_cols):
    print("\nCalculating social indicator Spearman correlations...")
    df_soc = pd.DataFrame()
    for s in social_cols:
        if s in df.columns:
            df_soc[s] = pd.to_numeric(df[s], errors='coerce')
            
    df_soc = df_soc.dropna(axis=1, how='all')
    
    corr_matrix = df_soc.corr(method='spearman')
    q_matrix = spearman_fdr_matrix(df_soc)

    order = cluster_order(corr_matrix)
    corr_matrix = corr_matrix.loc[order, order]
    q_matrix = q_matrix.loc[order, order]

    # Symmetric matrix: lower triangle only, with the structurally empty first row
    # and last column dropped. Label only FDR-significant cells so the +/-0.1 haze
    # stays visible as colour without inviting over-reading.
    corr_plot = corr_matrix.iloc[1:, :-1]
    q_plot = q_matrix.iloc[1:, :-1]
    mask = np.triu(np.ones(corr_plot.shape, dtype=bool), k=1)
    sig = (q_plot.values < 0.05) & corr_plot.notna().values
    annot = np.where(sig, np.vectorize(lambda v: f"{v:.2f}")(corr_plot.fillna(0).values), "")
    n_sig, n_cells = int((sig & ~mask).sum()), int((~mask).sum())

    fig, ax = plt.subplots(figsize=(19, 17))
    sns.heatmap(corr_plot, mask=mask, annot=annot, fmt="", cmap='RdBu_r',
                annot_kws={"size": 9}, center=0, vmin=-1, vmax=1, square=True,
                linewidths=0.5, linecolor='white',
                xticklabels=corr_plot.columns, yticklabels=corr_plot.index,
                cbar_kws={'shrink': 0.45, 'label': r'Spearman $\rho$',
                          'ticks': [-1, -0.5, 0, 0.5, 1]}, ax=ax)

    ax.set_title("Social indicators: Spearman correlation", fontsize=15, pad=16, loc='left')
    ax.tick_params(left=False, bottom=False, labelsize=11)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    plt.setp(ax.get_yticklabels(), rotation=0)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.text(0.01, 0.01,
             f"Lower triangle only (matrix is symmetric). Coefficients labelled for "
             f"{n_sig}/{n_cells} pairs at FDR < 0.05; unlabelled cells are still coloured. "
             f"Indicators ordered by hierarchical clustering on 1 - |rho|. "
             f"PCT-F and PCT-M are complementary by construction (rho = -1.00).",
             fontsize=9, color='#444444', va='bottom')

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(FIG_DIR / 'heatmap_social_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  heatmap_social_correlation.png: {n_sig}/{n_cells} pairs FDR<0.05")
    return corr_matrix

def analyze_treatments(df):
    print("Generating treatment visualizations...")
    df_cancer = df[df['Cohort'] == 'CANCER'].copy()

    # Column names
    sbrt_col = 'SBRT'
    chemo_col = 'Chemo receive  yes or none'
    immuno_col = 'Immunotherapy Receive  yes or none'
    target_col = 'Targeted therapy yes or none'
    
    # 1. Clean the columns for count plots
    def clean_treatment(series):
        # Map values to 'Yes' or 'None'
        s = series.astype(str).str.strip().str.upper()
        # If contains YES -> Yes, else None
        return np.where(s.str.contains('YES'), 'Yes', 'None')
        
    df_cancer['SBRT_clean'] = clean_treatment(df_cancer.get(sbrt_col, pd.Series(['None']*len(df_cancer))))
    df_cancer['Chemo_clean'] = clean_treatment(df_cancer.get(chemo_col, pd.Series(['None']*len(df_cancer))))
    df_cancer['Immuno_clean'] = clean_treatment(df_cancer.get(immuno_col, pd.Series(['None']*len(df_cancer))))
    df_cancer['Target_clean'] = clean_treatment(df_cancer.get(target_col, pd.Series(['None']*len(df_cancer))))

    # 2. Plot treatment counts (2x2 grid)
    fig, axes = plt.subplots(2, 2, figsize=(12, 6))
    axes = axes.flatten()
    
    sns.countplot(y='SBRT_clean', data=df_cancer, ax=axes[0], color='gray', order=['Yes', 'None'])
    axes[0].set_title('Count of SBRT')
    axes[0].set_ylabel('')
    
    sns.countplot(y='Chemo_clean', data=df_cancer, ax=axes[1], color='gray', order=['Yes', 'None'])
    axes[1].set_title('Count of Chemo receive  yes or none')
    axes[1].set_ylabel('')
    
    sns.countplot(y='Immuno_clean', data=df_cancer, ax=axes[2], color='gray', order=['Yes', 'None'])
    axes[2].set_title('Count of Immunotherapy Receive  yes or none')
    axes[2].set_ylabel('')
    
    sns.countplot(y='Target_clean', data=df_cancer, ax=axes[3], color='gray', order=['Yes', 'None'])
    axes[3].set_title('Count of Targeted therapy yes or none')
    axes[3].set_ylabel('')
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'bar_treatment_counts.png', dpi=300)
    plt.close()
    
    # 3. Treatment Combinations
    def get_combo(row):
        combo = []
        if row['SBRT_clean'] == 'Yes': combo.append('SBRT')
        if row['Chemo_clean'] == 'Yes': combo.append('Chemo')
        if row['Immuno_clean'] == 'Yes': combo.append('Immuno')
        if row['Target_clean'] == 'Yes': combo.append('Targeted')
        if not combo: return 'None'
        return '+'.join(combo)
        
    df_cancer['Treatment_Combo'] = df_cancer.apply(get_combo, axis=1)
    
    combo_counts = df_cancer['Treatment_Combo'].value_counts().reset_index()
    combo_counts.columns = ['Combination', 'Count']
    combo_counts = combo_counts.sort_values('Count', ascending=True)
    
    plt.figure(figsize=(10, 6))
    plt.barh(combo_counts['Combination'], combo_counts['Count'], color='gray')
    plt.xlabel('Number of Patients')
    plt.ylabel('Combination')
    plt.title('Distribution of All Treatment Combinations')
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'bar_treatment_combinations.png', dpi=300)
    plt.close()
    
    # 4. Treatment Stage Table
    if 'Stage' in df_cancer.columns:
        # Cross tabulation
        stage_combo = pd.crosstab(df_cancer['Treatment_Combo'], df_cancer['Stage'])
        
        # Plot as table
        fig, ax = plt.subplots(figsize=(10, 0.5 * len(stage_combo) + 1))
        ax.axis('tight')
        ax.axis('off')
        
        table = ax.table(cellText=stage_combo.values, 
                         colLabels=stage_combo.columns, 
                         rowLabels=stage_combo.index, 
                         loc='center',
                         cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        
        # Style table headers
        for key, cell in table.get_celld().items():
            if key[0] == 0 or key[1] == -1:
                cell.set_text_props(weight='bold')
                cell.set_facecolor('#EBEBEB')
            else:
                if key[0] % 2 == 0:
                    cell.set_facecolor('#F5F5F5')
                    
        plt.tight_layout()
        plt.savefig(FIG_DIR / 'table_treatment_stage.png', dpi=300, bbox_inches='tight')
        plt.close()

    # 5. Treatment Timing Stacked Bar Chart
    timing_cols = {
        'SBRT': 'Baseline_C00669_Follow Up Blood Draw before  (No) or After (yes)  SBRT',
        'Chemotherapy': 'Baseline_C00669_Follow Up Blood Draw before  (No) or After (yes)  Chemotherapy',
        'Immunotherapy': 'Baseline_C00669_Follow Up Blood Draw before  (No) or After (yes)  Immunotherapy'
    }
    
    timing_data = []
    for therapy, col in timing_cols.items():
        if col in df_cancer.columns:
            s = df_cancer[col].astype(str).str.strip().str.upper()
            before_count = s.str.contains('NO').sum()
            after_count = s.str.contains('YES').sum()
            if before_count + after_count > 0:
                timing_data.append({
                    'Therapy': therapy,
                    'Before Treatment': before_count,
                    'After Treatment': after_count,
                    'Total': before_count + after_count
                })
                
    if timing_data:
        timing_df = pd.DataFrame(timing_data)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        y_pos = np.arange(len(timing_df))
        
        # Plot before (green)
        ax.barh(y_pos, timing_df['Before Treatment'], color='#4CAF50', label='Before Treatment', height=0.6)
        
        # Plot after (blue) on top of before
        ax.barh(y_pos, timing_df['After Treatment'], left=timing_df['Before Treatment'], color='#2196F3', label='After Treatment', height=0.6)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(timing_df['Therapy'])
        ax.invert_yaxis()  # Labels read top-to-bottom
        ax.set_xlabel('Number of Patients')
        ax.set_title('Follow-Up Metabolite Measurements: Timing Relative to Treatment')
        ax.legend(loc='lower right')
        
        # Add 'n=' annotations
        for i, row in timing_df.iterrows():
            ax.text(row['Total'] + 0.5, i, f"n={row['Total']}", va='center', fontweight='bold')
            
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(FIG_DIR / 'bar_treatment_timing.png', dpi=300)
        plt.close()

def plot_treatment_timing_stats(df):
    """Two figures showing blood draw timing relative to treatment start for cancer patients."""
    print("Plotting treatment timing statistics...")
    df_c = df[df['Cohort'] == 'CANCER'].copy()

    date_cols = {
        'chemo_start': 'Start Date of Chemotherapy',
        'io_start':    'date of IO  Start',
        'sbrt_start':  'SBRT start date',
        'blood_draw':  'Date of baseline blood draw (visit)',
    }
    for key, col in date_cols.items():
        if col in df_c.columns:
            df_c[col] = pd.to_datetime(df_c[col], errors='coerce')

    draw_col  = 'Date of baseline blood draw (visit)'
    chemo_col = 'Start Date of Chemotherapy'
    io_col    = 'date of IO  Start'
    sbrt_col  = 'SBRT start date'
    chemo_rec = 'Chemo receive  yes or none'
    io_rec    = 'Immunotherapy Receive  yes or none'
    sbrt_rec  = 'SBRT'

    df_c['days_chemo'] = (df_c[draw_col] - df_c[chemo_col]).dt.days
    df_c['days_io']    = (df_c[draw_col] - df_c[io_col]).dt.days
    df_c['days_sbrt']  = (df_c[draw_col] - df_c[sbrt_col]).dt.days

    chemo_yes = df_c[chemo_rec].astype(str).str.strip().str.upper().str.contains('YES', na=False)
    io_yes    = df_c[io_rec].astype(str).str.strip().str.upper().str.contains('YES', na=False)
    sbrt_yes  = df_c[sbrt_rec].astype(str).str.strip().str.upper().str.contains('YES', na=False)

    treatments = {
        'Chemotherapy':  df_c.loc[chemo_yes, 'days_chemo'].dropna(),
        'Immunotherapy': df_c.loc[io_yes,    'days_io'].dropna(),
        'SBRT':          df_c.loc[sbrt_yes,  'days_sbrt'].dropna(),
    }

    color_post = '#d62728'
    color_pre  = '#1f77b4'
    n = len(df_c)

    # ------------------------------------------------------------------
    # Figure 1: stacked horizontal bar — pre vs post counts + medians
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 4.5))
    trt_labels = list(treatments.keys())
    n_post_list, n_pre_list, med_post_list, med_pre_list = [], [], [], []

    for d in treatments.values():
        post = d[d >= 0]; pre = d[d < 0]
        n_post_list.append(len(post))
        n_pre_list.append(len(pre))
        med_post_list.append(post.median() if len(post) else np.nan)
        med_pre_list.append(pre.median()  if len(pre)  else np.nan)

    y = np.arange(len(trt_labels))
    bar_h = 0.55

    ax.barh(y, n_post_list, height=bar_h, color=color_post,
            label='Drawn AFTER treatment start')
    ax.barh(y, n_pre_list,  height=bar_h, left=n_post_list, color=color_pre,
            label='Drawn BEFORE treatment start')

    for i, (np_, npr, mp, mpr) in enumerate(
            zip(n_post_list, n_pre_list, med_post_list, med_pre_list)):
        if np_ > 0:
            ax.text(np_/2, i, f'n={np_}\nmedian={mp:.0f}d',
                    ha='center', va='center', fontsize=9, color='white', fontweight='bold')
        if npr > 0:
            ax.text(np_ + npr/2, i, f'n={npr}\nmedian={mpr:.0f}d',
                    ha='center', va='center', fontsize=9, color='white', fontweight='bold')
        ax.text(np_ + npr + 0.3, i, f'n={np_+npr}',
                va='center', fontsize=9.5, color='black')

    ax.set_yticks(y)
    ax.set_yticklabels(trt_labels)
    ax.set_xlabel('Number of patients')
    ax.set_title(f'Blood Draw Timing Relative to Treatment Start\n({n} cancer patients)', pad=10)
    # ax.legend(loc='lower right', framealpha=0.9)
    ax.legend(fontsize=9)
    ax.set_xlim(0, max(a+b for a, b in zip(n_post_list, n_pre_list)) + 5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'bar_treatment_timing_counts.png', dpi=300, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------
    # Figure 2: strip plot — distribution of days relative to start
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.5))

    for ax, (trt, d) in zip(axes, treatments.items()):
        post = d[d >= 0]; pre = d[d < 0]
        rng = np.random.default_rng(42)
        j_post = rng.uniform(-0.15, 0.15, len(post))
        j_pre  = rng.uniform(-0.15, 0.15, len(pre))

        ax.scatter(j_post, post.values, color=color_post, alpha=0.75, s=55, zorder=3,
                   label=f'After  (n={len(post)})')
        ax.scatter(j_pre,  pre.values,  color=color_pre,  alpha=0.75, s=55, zorder=3,
                   label=f'Before (n={len(pre)})')

        if len(post):
            ax.hlines(post.median(), -0.35, 0.35, colors=color_post,
                      linewidth=2.5, linestyle='--', zorder=4,
                      label=f'Median after:  {post.median():.0f}d')
        if len(pre):
            ax.hlines(pre.median(), -0.35, 0.35, colors=color_pre,
                      linewidth=2.5, linestyle='--', zorder=4,
                      label=f'Median before: {pre.median():.0f}d')

        ax.axhline(0, color='black', linewidth=1.2, linestyle='-', zorder=2, alpha=0.6)
        ax.text(0.37, 0, ' treatment\n start', va='center', fontsize=8, color='black', alpha=0.7)
        ax.set_title(trt, fontweight='bold')
        ax.set_ylabel('Days (blood draw − treatment start)')
        ax.set_xticks([])
        ax.set_xlim(-0.6, 0.65)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.legend(fontsize=8.5, loc='upper left', framealpha=0.9)

    fig.suptitle(
        f'Distribution of Days Between Blood Draw and Treatment Start\n'
        f'({n} cancer patients; positive = drawn after treatment started)',
        fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'scatter_treatment_timing_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Treatment timing figures saved.")


def _timing_strip(ax, d, ref_label, ylabel):
    """One strip-plot panel of day offsets relative to a reference date."""
    color_post, color_pre = '#d62728', '#1f77b4'
    post, pre = d[d >= 0], d[d < 0]
    rng = np.random.default_rng(42)
    ax.scatter(rng.uniform(-0.15, 0.15, len(post)), post.values, color=color_post,
               alpha=0.75, s=55, zorder=3, label=f'After  (n={len(post)})')
    ax.scatter(rng.uniform(-0.15, 0.15, len(pre)), pre.values, color=color_pre,
               alpha=0.75, s=55, zorder=3, label=f'Before (n={len(pre)})')
    if len(post):
        ax.hlines(post.median(), -0.35, 0.35, colors=color_post, linewidth=2.5,
                  linestyle='--', zorder=4, label=f'Median after:  {post.median():.0f}d')
    if len(pre):
        ax.hlines(pre.median(), -0.35, 0.35, colors=color_pre, linewidth=2.5,
                  linestyle='--', zorder=4, label=f'Median before: {pre.median():.0f}d')
    ax.axhline(0, color='black', linewidth=1.2, zorder=2, alpha=0.6)
    ax.text(0.37, 0, f' {ref_label}', va='center', fontsize=8, color='black', alpha=0.7)
    ax.set_ylabel(ylabel)
    ax.set_xticks([])
    ax.set_xlim(-0.6, 0.65)
    for s in ('top', 'right', 'bottom'):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8.5, loc='upper left', framealpha=0.9)


def plot_followup_end_timing_stats(df):
    """Follow-up blood draw relative to END of treatment, per modality.

    Only immunotherapy has an end date, so only its panel is a true end reference and
    can be drawn as a day offset. Chemotherapy and SBRT fall back to the sheet's
    'Follow Up Blood Draw before/After' flags.

    What those flags mean is UNRESOLVED and the figure says so rather than picking a
    side. Calibrated on the IO patients who have both a flag and both dates, the flag
    matches 'draw after IO start' far more often than 'draw after IO end' -- but
    `date of IO End` has defects of its own (more end dates than start dates, and end
    dates in the future), so the alternative reading is that the flag means 'after
    treatment ended' and the end dates are unreliable. Both counts are printed on the
    figure. Resolving it needs the data owner, not more analysis.
    """
    print("\nPlotting follow-up draw vs END of treatment...")
    d = df[df['Cohort'] == 'CANCER'].copy()
    for c in [FU_DRAW_DATE_COL, 'date of IO  Start', 'date of IO End']:
        if c in d.columns:
            d[c] = pd.to_datetime(d[c], errors='coerce')
    fu = d[FU_DRAW_DATE_COL] if FU_DRAW_DATE_COL in d.columns else pd.Series(pd.NaT, index=d.index)
    n_fu = int(fu.notna().sum())
    if n_fu < 5:
        print("  Too few follow-up draws — skipping.")
        return

    # Calibration, recomputed so the figure never states a stale number.
    ios, ioe = d.get('date of IO  Start'), d.get('date of IO End')
    fl_io = d[FU_AFTER_IO_COL].astype(str).str.strip().str.upper().str.contains('YES', na=False)
    cal = fu.notna() & ios.notna() & ioe.notna()
    agree_start = int((fl_io[cal] == (fu > ios)[cal]).sum())
    agree_end = int((fl_io[cal] == (fu > ioe)[cal]).sum())
    n_cal = int(cal.sum())
    # Defects in the end-date column itself, which is why the calibration above is not
    # decisive: if some end dates are planned rather than actual, the flag could be right.
    n_end_no_start = int((ioe.notna() & ios.isna()).sum())
    n_future_end = int((ioe > pd.Timestamp.today()).sum())
    fl_ch = d[FU_AFTER_CHEMO_COL].astype(str).str.strip().str.upper().str.contains('YES', na=False) \
        if FU_AFTER_CHEMO_COL in d.columns else pd.Series(False, index=d.index)
    n_same = int((fl_io[fu.notna()] == fl_ch[fu.notna()]).sum())

    def flag_of(col):
        return d[col].astype(str).str.strip().str.upper().str.contains('YES', na=False) \
            if col in d.columns else pd.Series(False, index=d.index)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))

    # --- Chemotherapy: no end date, flag fallback ---
    for ax, (name, flag_col) in zip([axes[0], axes[2]],
                                    [('Chemotherapy', FU_AFTER_CHEMO_COL),
                                     ('SBRT', FU_AFTER_SBRT_COL)]):
        fl = flag_of(flag_col)[fu.notna()]
        counts = [int((~fl).sum()), int(fl.sum())]
        ax.bar(['Before', 'After'], counts, color=['#1f77b4', '#d62728'], alpha=0.8, width=0.55)
        for i, c in enumerate(counts):
            ax.text(i, c + 0.4, str(c), ha='center', fontsize=11, fontweight='bold')
        extra = '  (flag identical to IO in %d/%d rows)' % (n_same, int(fu.notna().sum())) \
            if name == 'Chemotherapy' else ''
        ax.set_title(f'{name}\nNO END DATE — sheet flag only (n={int(fu.notna().sum())})',
                     fontweight='bold', color='#a03030')
        ax.set_ylabel('Patients with a follow-up draw')
        ax.set_ylim(0, max(counts) * 1.25 if max(counts) else 1)
        ax.text(0.5, 0.93, '"After" = start or end? unresolved' + extra,
                transform=ax.transAxes, ha='center', fontsize=8.5, color='#a03030',
                style='italic')
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)

    # --- Immunotherapy: real end date ---
    delta = (fu - ioe).dt.days.dropna()
    if len(delta):
        _timing_strip(axes[1], delta, 'IO end', 'Days (follow-up draw − IO end)')
        axes[1].set_title(f'Immunotherapy\nTRUE END DATE (n={len(delta)})',
                          fontweight='bold', color='#20693a')
    else:
        axes[1].axis('off')

    fig.suptitle('Follow-up Blood Draw Relative to END of Treatment\n'
                 f'({n_fu} cancer patients with a follow-up draw)', fontsize=12, y=1.04)
    fig.text(0.005, -0.05,
             f'Only immunotherapy has an end date; for chemotherapy and SBRT the sheet\'s before/after flag is '
             f'the only timing information. What "After" means in that flag is UNRESOLVED. On the {n_cal} IO '
             f'patients with both a flag and both dates it matches "draw after IO START" {agree_start}/{n_cal} '
             f'times but "draw after IO END" only {agree_end}/{n_cal} — which favours a start reading. Against '
             f'that, `date of IO End` has {n_end_no_start} patients with an end date but no start date and '
             f'{n_future_end} end dates in the future, so it may record planned rather than actual completion, '
             f'in which case the flag could mean end and the dates are wrong. The two readings imply ~{agree_end} '
             f'vs ~{int(fl_io[fu.notna()].sum())} genuinely post-treatment IO patients. Resolving this needs the '
             f'data owner, not more analysis.',
             fontsize=9, color='#a03030', va='top', wrap=True)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'scatter_followup_vs_treatment_end.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  flag calibration on {n_cal} IO patients: matches after-START {agree_start}, "
          f"after-END {agree_end}")
    if len(delta):
        print(f"  IO (true end date): n={len(delta)}  after end={int((delta>=0).sum())}  "
              f"before end={int((delta<0).sum())}  median={delta.median():+.0f}d")
    for name, col in [('Chemotherapy', FU_AFTER_CHEMO_COL), ('SBRT', FU_AFTER_SBRT_COL)]:
        fl = flag_of(col)[fu.notna()]
        print(f"  {name:14s} (flag): after={int(fl.sum())}  before={int((~fl).sum())}")
    print("  Saved scatter_followup_vs_treatment_end.png")


def plot_followup_timing_stats(df):
    """Companion to plot_treatment_timing_stats, for the FOLLOW-UP blood draw.

    Panels 1-3 place the follow-up draw relative to each treatment's start; panel 4
    places it relative to IO end. Only immunotherapy has an end date recorded -- the
    workbook has no chemotherapy or SBRT end date at all -- so 'draw is after treatment
    finished' is verifiable for IO alone. That limitation is stated on the figure
    rather than left for the reader to infer from a missing panel.
    """
    print("\nPlotting follow-up blood draw timing...")
    d = df[df['Cohort'] == 'CANCER'].copy()
    for c in ['Date of baseline blood draw (visit)', FU_DRAW_DATE_COL,
              'Start Date of Chemotherapy', 'date of IO  Start', 'date of IO End',
              'SBRT start date']:
        if c in d.columns:
            d[c] = pd.to_datetime(d[c], errors='coerce')

    fu = d[FU_DRAW_DATE_COL] if FU_DRAW_DATE_COL in d.columns else pd.Series(pd.NaT, index=d.index)
    n_fu = int(fu.notna().sum())
    if n_fu < 5:
        print("  Too few follow-up draws — skipping.")
        return

    panels = [
        ('Chemotherapy start', 'Start Date of Chemotherapy', 'chemo start'),
        ('Immunotherapy start', 'date of IO  Start', 'IO start'),
        ('SBRT start', 'SBRT start date', 'SBRT start'),
        ('Immunotherapy END', 'date of IO End', 'IO end'),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(18, 5.5))
    for ax, (title, col, ref) in zip(axes, panels):
        if col not in d.columns:
            ax.axis('off')
            continue
        delta = (fu - d[col]).dt.days.dropna()
        if delta.empty:
            ax.axis('off')
            ax.set_title(f'{title}\n(no dates)', fontweight='bold')
            continue
        _timing_strip(ax, delta, ref, 'Days (follow-up draw − reference)')
        ax.set_title(f'{title}\n(n={len(delta)})', fontweight='bold')

    fig.suptitle(
        f'Distribution of Days Between FOLLOW-UP Blood Draw and Treatment Dates\n'
        f'({n_fu} cancer patients with a follow-up draw; positive = drawn after the reference date)',
        fontsize=12, y=1.03)
    fig.text(0.005, -0.02,
             'No chemotherapy or SBRT end date exists in the workbook, so only immunotherapy '
             'can be checked for "draw after treatment finished"; the other three panels are '
             'relative to treatment START.',
             fontsize=9, color='#a03030', va='top')
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'scatter_followup_timing_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Numeric summary, since the counts are the part that matters for interpretation.
    print(f"  follow-up draws: {n_fu}")
    for title, col, _ in panels:
        if col not in d.columns:
            continue
        delta = (fu - d[col]).dt.days.dropna()
        if delta.empty:
            continue
        print(f"    {title:22s} n={len(delta):2d}  after={int((delta>=0).sum()):2d}  "
              f"before={int((delta<0).sum()):2d}  median={delta.median():+.0f}d")
    print("  Saved scatter_followup_timing_distribution.png")


def plot_therapy_timeline(df, deceased_only=False):
    """Swimmer plot: one row per cancer patient (or deceased-only subset), showing
    days-from-diagnosis for blood draws, treatment starts, and survival endpoint."""
    subset_label = 'Deceased Patients' if deceased_only else 'Cancer Patients'
    print(f"\nPlotting therapy timeline ({subset_label})...")

    df_c = df[df['Cohort'] == 'CANCER'].copy()
    if deceased_only:
        df_c = df_c[df_c['Event'] == 1].copy()

    diag_col    = 'Date of Initial Diagnosis'
    draw_col    = 'Date of baseline blood draw (visit)'
    fu_draw_col = 'Baseline_C00669_Follow up Date blood was drawn'
    chemo_col   = 'Start Date of Chemotherapy'
    io_col      = 'date of IO  Start'
    sbrt_col    = 'SBRT start date'

    for col in [diag_col, draw_col, fu_draw_col, chemo_col, io_col, sbrt_col]:
        if col in df_c.columns:
            df_c[col] = pd.to_datetime(df_c[col], errors='coerce')

    def days_from_diag(col):
        return (df_c[col] - df_c[diag_col]).dt.days

    df_c['d_draw']    = days_from_diag(draw_col)    if draw_col    in df_c.columns else np.nan
    df_c['d_fu_draw'] = days_from_diag(fu_draw_col) if fu_draw_col in df_c.columns else np.nan
    df_c['d_chemo']   = days_from_diag(chemo_col)   if chemo_col   in df_c.columns else np.nan
    df_c['d_io']      = days_from_diag(io_col)       if io_col      in df_c.columns else np.nan
    df_c['d_sbrt']    = days_from_diag(sbrt_col)    if sbrt_col    in df_c.columns else np.nan

    df_c = df_c.sort_values('Time to Event', ascending=True).reset_index(drop=True)
    n = len(df_c)

    # Fewer rows → larger per-row height and font for the deceased-only plot
    if deceased_only:
        row_h    = 0.52
        label_fs = 14
    else:
        row_h    = 0.38
        label_fs = 30

    fig_h   = max(10, n * row_h)
    fig, ax = plt.subplots(figsize=(22, fig_h))

    dead_color   = 'r'
    censor_color = 'b'

    stage_colors = {
        'Stage I':   '#2166ac',
        'Stage II':  '#74add1',
        'Stage III': '#f46d43',
        'Stage IV':  '#d73027',
    }
    x_max   = df_c['Time to Event'].max()
    label_x = x_max * 1.03

    patient_labels = []
    for i, (_, row) in enumerate(df_c.iterrows()):
        tte   = row['Time to Event']
        event = row['Event']

        bar_color  = dead_color if event == 1 else censor_color
        end_marker = 'x' if event == 1 else '|'
        ax.hlines(i, 0, tte, color=bar_color, linewidth=3, alpha=0.55, zorder=1)
        ax.scatter(tte, i, marker=end_marker, color=bar_color, s=80, zorder=4, linewidths=2.0)

        if pd.notna(row['d_draw']):
            ax.scatter(row['d_draw'],    i, marker='*', color='#2ca02c', s=260, zorder=5)
        if pd.notna(row.get('d_fu_draw')) and pd.notna(row['d_fu_draw']):
            ax.scatter(row['d_fu_draw'], i, marker='*', color='#98df8a', s=200, zorder=5)
        if pd.notna(row['d_chemo']):
            ax.scatter(row['d_chemo'],   i, marker='o', color='#1f77b4', s=240, zorder=5)
        if pd.notna(row['d_io']):
            ax.scatter(row['d_io'],      i, marker='D', color='#ff7f0e', s=230, zorder=5)
        if pd.notna(row['d_sbrt']):
            ax.scatter(row['d_sbrt'],    i, marker='s', color='#9467bd', s=230, zorder=5)

        stage_str = str(row.get('Stage_orig', '')) if pd.notna(row.get('Stage_orig')) else 'Unknown'
        if stage_str in ('nan', 'None', ''): stage_str = 'Unknown'
        s_color = stage_colors.get(stage_str, '#888888')
        ax.text(label_x, i, stage_str, va='center', ha='left',
                fontsize=label_fs, color=s_color, fontweight='bold')

        pid = row.get('Patient ID', '')
        patient_labels.append(str(pid) if pd.notna(pid) and str(pid) not in ('nan', '') else f'P{i+1}')

    ax.set_yticks(range(n))
    ax.set_yticklabels(patient_labels, fontsize=label_fs)
    ax.set_xlim(0, x_max * 1.30)

    leg_ms = 12 if deceased_only else 20
    leg_fs = 11 if deceased_only else 30
    legend_els = [
        Line2D([0], [0], color=dead_color,   linewidth=2.5, label='Dead'),
        Line2D([0], [0], color=censor_color, linewidth=2.5, label='Censored (alive)'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='#2ca02c', markersize=leg_ms,     label='Baseline blood draw'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='#98df8a', markersize=leg_ms - 4, label='Follow-up blood draw'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=leg_ms,     label='Chemo start'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#ff7f0e', markersize=leg_ms,     label='IO start'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#9467bd', markersize=leg_ms,     label='SBRT start'),
        Line2D([0], [0], marker='x', color=dead_color,   markersize=leg_ms, linewidth=0, label='Death'),
        Line2D([0], [0], marker='|', color=censor_color, markersize=leg_ms, linewidth=0, label='Last follow-up'),
    ]
    ax.legend(handles=legend_els, fontsize=leg_fs, framealpha=0.9, ncol=2)

    ax.set_xlabel('Days from initial diagnosis', fontsize=label_fs)
    ax.set_ylabel('Patient ID', fontsize=label_fs)
    ax.set_title(
        f'Therapy Timeline — {n} {subset_label}\n(sorted by survival time, bottom = shortest)',
        fontsize=label_fs, fontweight='bold'
    )
    ax.tick_params(axis='x', labelsize=label_fs)
    ax.tick_params(axis='y', labelsize=label_fs)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.tight_layout()
    out = FIG_DIR / ('timeline_therapy_deceased.png' if deceased_only else 'timeline_therapy_swimmer.png')
    plt.savefig(out, dpi=120, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out.name}')


def analyze_smoking_pack_year_split(df, met_cols, n_top=10):
    """Differential metabolite analysis in CANCER patients split by Pack Year History
    (above vs at/below the median). Produces a volcano plot and top boxplots."""
    print("\n" + "="*70)
    print("SMOKING (PACK-YEAR MEDIAN SPLIT) METABOLITE ANALYSIS — CANCER PATIENTS")
    print("="*70)

    df_cancer = df[df['Cohort'] == 'CANCER'].copy()
    py = pd.to_numeric(df_cancer.get('Pack Year History'), errors='coerce')
    df_cancer = df_cancer[py.notna()].copy()
    py = py[py.notna()]
    if py.empty:
        print("  No Pack Year History values — skipping.")
        return None

    median_py = py.median()
    df_cancer['PackYear_Group'] = np.where(py > median_py, 'High Pack-Year', 'Low Pack-Year')
    n_high = int((df_cancer['PackYear_Group'] == 'High Pack-Year').sum())
    n_low  = int((df_cancer['PackYear_Group'] == 'Low Pack-Year').sum())
    print(f"  Median pack-years: {median_py:.1f}")
    print(f"  High Pack-Year (> {median_py:.1f}):  {n_high}")
    print(f"  Low Pack-Year  (<= {median_py:.1f}): {n_low}")

    if n_high < 5 or n_low < 5:
        print("  Insufficient samples per group — skipping.")
        return None

    results = rank_sum_test(df_cancer, met_cols, is_log=True,
                            group_col='PackYear_Group',
                            group1='High Pack-Year', group2='Low Pack-Year')
    if results is None or results.empty:
        print("  No rank-sum results.")
        return None

    sig = results[results['FDR'] < 0.05]
    print(f"  FDR < 0.05: {len(sig)} metabolites")

    # ---- Volcano ----
    plot_volcano(results,
                 f"Metabolites (High n={n_high} vs Low n={n_low} pack-year, cancer)",
                 "volcano_metabolites_packyear_split.png",
                 group_a='High pack-year', group_b='Low pack-year')

    # ---- Top boxplots (FDR-significant, else nominal top-N) ----
    inc = results[(results['FDR'] < 0.05) & (results['Direction'] == 'Increased')].sort_values('FDR')
    inc = inc.head(n_top) if not inc.empty else results[results['Direction'] == 'Increased'].sort_values('p_value').head(n_top)
    dec = results[(results['FDR'] < 0.05) & (results['Direction'] == 'Decreased')].sort_values('FDR')
    dec = dec.head(n_top) if not dec.empty else results[results['Direction'] == 'Decreased'].sort_values('p_value').head(n_top)

    plot_boxplots(df_cancer, inc,
                  "Top metabolites higher in high pack-year cancer patients",
                  "boxplot_top_metabolites_packyear_increased.png",
                  is_log=True, hue_col='PackYear_Group',
                  hue_order=['High Pack-Year', 'Low Pack-Year'])
    plot_boxplots(df_cancer, dec,
                  "Top metabolites lower in high pack-year cancer patients",
                  "boxplot_top_metabolites_packyear_decreased.png",
                  is_log=True, hue_col='PackYear_Group',
                  hue_order=['High Pack-Year', 'Low Pack-Year'])

    return results


def analyze_methylation_by_stage(df, meth_cols):
    """DNA methylation markers: Metastatic (Stage IV) vs Non-Metastatic (Stage I-III).
    Boxplots labeled with RAW rank-sum p-values only (no FDR)."""
    print("\n" + "="*70)
    print("DNA METHYLATION: METASTATIC (IV) vs NON-METASTATIC (I-III)")
    print("="*70)
    dstg = df[(df['Cohort'] == 'CANCER') & (df['Stage'].isin(['Metastatic', 'Non-Metastatic']))].copy()
    n_meta = int((dstg['Stage'] == 'Metastatic').sum())
    n_non  = int((dstg['Stage'] == 'Non-Metastatic').sum())
    print(f"  Metastatic (Stage IV): {n_meta} | Non-Metastatic (Stage I-III): {n_non}")

    rows = []
    for m in meth_cols:
        if m not in dstg.columns:
            continue
        a = pd.to_numeric(dstg.loc[dstg['Stage'] == 'Non-Metastatic', m], errors='coerce').dropna()
        b = pd.to_numeric(dstg.loc[dstg['Stage'] == 'Metastatic', m], errors='coerce').dropna()
        if len(a) > 3 and len(b) > 3:
            _, p = ranksums(a, b)
            direction = 'Higher in Metastatic' if b.median() > a.median() else 'Higher in Non-Metastatic'
            rows.append({'Feature': m, 'p_value': p, 'Direction': direction})
            print(f"  {str(m):<10} p={p:.4f}  {direction}")

    # Only Feature + p_value -> plot_boxplots shows p-values only (no q / log2FC)
    res = pd.DataFrame(rows)
    if res.empty:
        print("  No methylation markers testable.")
        return res
    plot_boxplots(dstg, res,
                  "DNA methylation: Metastatic (IV) vs Non-Metastatic (I-III)",
                  "boxplot_methylation_metastatic_vs_nonmetastatic.png",
                  is_log=False, hue_col='Stage', hue_order=['Metastatic', 'Non-Metastatic'])
    return res


def _stage_early_late(df):
    """Return cancer df with a StageGroup column: Early = Stage I-II, Late = Stage III-IV
    (unknown stage excluded)."""
    dfc = df[df['Cohort'] == 'CANCER'].copy()
    dfc['StageGroup'] = 'Other'
    dfc.loc[dfc['Stage_orig'].isin(['Stage I', 'Stage II']), 'StageGroup'] = 'Early'
    dfc.loc[dfc['Stage_orig'].isin(['Stage III', 'Stage IV']), 'StageGroup'] = 'Late'
    return dfc[dfc['StageGroup'].isin(['Early', 'Late'])].copy()


def analyze_early_vs_late(df, met_cols, n_top=10):
    """Differential metabolites: Early (Stage I-II) vs Late (Stage IV). Stage III excluded."""
    print("\n" + "="*70)
    print("EARLY (Stage I-II) vs LATE (Stage III-IV) METABOLITE ANALYSIS")
    print("="*70)
    grp = _stage_early_late(df)
    n_early = int((grp['StageGroup'] == 'Early').sum())
    n_late  = int((grp['StageGroup'] == 'Late').sum())
    print(f"  Early (Stage I-II): {n_early} | Late (Stage III-IV): {n_late}")
    if n_early < 5 or n_late < 5:
        print("  Insufficient samples — skipping.")
        return None

    results = rank_sum_test(grp, met_cols, is_log=True, group_col='StageGroup',
                            group1='Late', group2='Early')
    if results is None or results.empty:
        print("  No rank-sum results.")
        return None
    print(f"  FDR<0.05: {int((results['FDR'] < 0.05).sum())} | raw p<0.05: {int((results['p_value'] < 0.05).sum())}")

    plot_volcano(results,
                 f"Metabolites (Late/Stage III-IV n={n_late} vs Early/Stage I-II n={n_early})",
                 "volcano_metabolites_early_vs_late.png",
                 group_a='Late (III-IV)', group_b='Early (I-II)')

    inc = results[(results['FDR'] < 0.05) & (results['Direction'] == 'Increased')].sort_values('FDR')
    inc = inc.head(n_top) if not inc.empty else results[results['Direction'] == 'Increased'].sort_values('p_value').head(n_top)
    dec = results[(results['FDR'] < 0.05) & (results['Direction'] == 'Decreased')].sort_values('FDR')
    dec = dec.head(n_top) if not dec.empty else results[results['Direction'] == 'Decreased'].sort_values('p_value').head(n_top)

    plot_boxplots(grp, inc,
                  f"Top {n_top} metabolites increased in Late (III-IV) vs Early (I-II)",
                  "boxplot_top_metabolites_late_increased.png",
                  is_log=True, hue_col='StageGroup', hue_order=['Early', 'Late'])
    plot_boxplots(grp, dec,
                  f"Top {n_top} metabolites decreased in Late (III-IV) vs Early (I-II)",
                  "boxplot_top_metabolites_late_decreased.png",
                  is_log=True, hue_col='StageGroup', hue_order=['Early', 'Late'])
    return results


def analyze_methylation_early_vs_late(df, meth_cols):
    """DNA methylation markers: Early (Stage I-II) vs Late (Stage IV). Stage III excluded.
    Boxplots labeled with RAW rank-sum p-values only (no FDR)."""
    print("\n" + "="*70)
    print("DNA METHYLATION: EARLY (Stage I-II) vs LATE (Stage III-IV)")
    print("="*70)
    grp = _stage_early_late(df)
    n_early = int((grp['StageGroup'] == 'Early').sum())
    n_late  = int((grp['StageGroup'] == 'Late').sum())
    print(f"  Early (Stage I-II): {n_early} | Late (Stage III-IV): {n_late}")

    rows = []
    for m in meth_cols:
        if m not in grp.columns:
            continue
        a = pd.to_numeric(grp.loc[grp['StageGroup'] == 'Early', m], errors='coerce').dropna()
        b = pd.to_numeric(grp.loc[grp['StageGroup'] == 'Late', m], errors='coerce').dropna()
        if len(a) > 3 and len(b) > 3:
            _, p = ranksums(a, b)
            direction = 'Higher in Late' if b.median() > a.median() else 'Higher in Early'
            rows.append({'Feature': m, 'p_value': p, 'Direction': direction})
            print(f"  {str(m):<10} p={p:.4f}  {direction}")

    res = pd.DataFrame(rows)  # Feature + p_value only -> boxplots show p-values only
    if res.empty:
        print("  No methylation markers testable.")
        return res
    plot_boxplots(grp, res,
                  "DNA methylation: Early (I-II) vs Late (III-IV)",
                  "boxplot_methylation_early_vs_late.png",
                  is_log=False, hue_col='StageGroup', hue_order=['Early', 'Late'])
    return res


def analyze_metastatic_vs_nonmetastatic(df, met_cols, n_top=10):
    """Differential metabolite analysis: Metastatic (Stage IV) vs Non-Metastatic (Stage I–III),
    cancer patients only. Produces volcano, boxplots, and KDE distribution figures."""
    print("\n" + "="*70)
    print("METASTATIC vs NON-METASTATIC METABOLITE ANALYSIS")
    print("="*70)

    df_cancer = df[df['Cohort'] == 'CANCER'].copy()
    df_staged = df_cancer[df_cancer['Stage'].isin(['Metastatic', 'Non-Metastatic'])].copy()

    n_meta    = (df_staged['Stage'] == 'Metastatic').sum()
    n_nonmeta = (df_staged['Stage'] == 'Non-Metastatic').sum()
    n_missing = len(df_cancer) - len(df_staged)
    print(f"  Metastatic    (Stage IV):      {n_meta}")
    print(f"  Non-Metastatic (Stage I-III):  {n_nonmeta}")
    print(f"  Missing / unknown stage:       {n_missing}")

    if n_meta < 5 or n_nonmeta < 5:
        print("  Insufficient samples — skipping.")
        return None

    # ---- Rank-sum test ----
    results = rank_sum_test(df_staged, met_cols, is_log=True,
                            group_col='Stage',
                            group1='Metastatic', group2='Non-Metastatic')
    if results is None or results.empty:
        print("  No rank-sum results.")
        return None

    sig = results[results['FDR'] < 0.05]
    inc = sig[sig['Direction'] == 'Increased'].sort_values('FDR')
    dec = sig[sig['Direction'] == 'Decreased'].sort_values('FDR')
    print(f"\n  FDR < 0.05 : {len(sig)} metabolites "
          f"({len(inc)} increased in Metastatic, {len(dec)} decreased)\n")

    clean_name = met_label

    if not inc.empty:
        print("  === INCREASED IN METASTATIC (FDR < 0.05) ===")
        for _, r in inc.iterrows():
            print(f"    {clean_name(r['Feature']):<50}  FDR={r['FDR']:.4f}  p={r['p_value']:.2e}")

    if not dec.empty:
        print(f"\n  === DECREASED IN METASTATIC (FDR < 0.05) — top 20 ===")
        for _, r in dec.head(20).iterrows():
            print(f"    {clean_name(r['Feature']):<50}  FDR={r['FDR']:.4f}  p={r['p_value']:.2e}")

    # fall back to nominal top-N if fewer than n_top pass FDR
    top_inc = (inc.head(n_top) if len(inc) >= 1
               else results[results['Direction'] == 'Increased'].sort_values('p_value').head(n_top))
    top_dec = (dec.head(n_top) if len(dec) >= 1
               else results[results['Direction'] == 'Decreased'].sort_values('p_value').head(n_top))

    # ---- Volcano ----
    plot_volcano(results,
                 f"Metabolites (Metastatic n={n_meta} vs Non-Metastatic n={n_nonmeta})",
                 "volcano_metabolites_metastatic_vs_nonmetastatic.png",
                 group_a='Metastatic', group_b='Non-Metastatic')

    # ---- Boxplots ----
    plot_boxplots(df_staged, top_inc,
                  f"Top {n_top} metabolites increased in Metastatic vs Non-Metastatic",
                  "boxplot_top_metabolites_metastatic_increased.png",
                  is_log=True, hue_col='Stage',
                  hue_order=['Metastatic', 'Non-Metastatic'])

    plot_boxplots(df_staged, top_dec,
                  f"Top {n_top} metabolites decreased in Metastatic vs Non-Metastatic",
                  "boxplot_top_metabolites_metastatic_decreased.png",
                  is_log=True, hue_col='Stage',
                  hue_order=['Metastatic', 'Non-Metastatic'])

    # ---- Single-metabolite boxplot: cysteine S-sulfate (linear y, groups via legend) ----
    css = results[results['Feature'].str.contains('cysteine-S-sulfate', case=False, na=False)]
    if not css.empty:
        css_col = css.iloc[0]['Feature']
        q_val   = css.iloc[0].get('FDR', np.nan)
        l2fc    = css.iloc[0].get('log2FC', np.nan)
        cssdf = df_staged[[css_col, 'Stage']].copy()
        cssdf[css_col] = pd.to_numeric(cssdf[css_col], errors='coerce')
        cssdf = cssdf.dropna(subset=[css_col])
        cssdf['_x'] = ''
        plt.figure(figsize=(6, 6))
        ax = sns.boxplot(data=cssdf, x='_x', y=css_col, hue='Stage',
                         hue_order=['Metastatic', 'Non-Metastatic'],
                         palette={'Metastatic': '#e6550d', 'Non-Metastatic': '#3182bd'},
                         showfliers=False, width=0.6, linewidth=1.0)
        ax.set_title("Cysteine S-sulfate: Metastatic vs Non-Metastatic\n"
                     f"q={q_val:.2e},  $log_2$(FC)={l2fc:.2f}", fontsize=12)
        ax.set_ylabel("Cysteine S-sulfate")
        ax.set_xlabel("")
        ax.set_xticks([])
        ax.legend(title='Stage', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        plt.tight_layout()
        plt.savefig(FIG_DIR / "boxplot_cysteine_s_sulfate_metastatic.png", dpi=300, bbox_inches='tight')
        plt.close()

    # ---- KDE distribution plots ----
    df_meta_only    = df_staged[df_staged['Stage'] == 'Metastatic']
    df_nonmeta_only = df_staged[df_staged['Stage'] == 'Non-Metastatic']
    color_meta    = '#e6550d'
    color_nonmeta = '#3182bd'

    for top_df, direction, label in [
        (top_inc, 'Increased', 'increased_in_metastatic'),
        (top_dec, 'Decreased', 'decreased_in_metastatic'),
    ]:
        if top_df.empty:
            continue

        n = len(top_df)
        fig, axes = plt.subplots(n, 1, figsize=(9, 2.2 * n), sharex=False)
        if n == 1:
            axes = [axes]

        for i, (ax, (_, row)) in enumerate(zip(axes, top_df.iterrows())):
            feat = row['Feature']
            fdr  = row['FDR']
            name = clean_name(feat)

            v_meta    = np.log1p(pd.to_numeric(df_meta_only[feat],    errors='coerce').dropna())
            v_nonmeta = np.log1p(pd.to_numeric(df_nonmeta_only[feat], errors='coerce').dropna())

            if len(v_meta) > 2:
                sns.kdeplot(v_meta,    ax=ax, fill=True, alpha=0.45,
                            color=color_meta,    linewidth=1.5,
                            label=f'Metastatic     (n={len(v_meta)})')
            if len(v_nonmeta) > 2:
                sns.kdeplot(v_nonmeta, ax=ax, fill=True, alpha=0.45,
                            color=color_nonmeta, linewidth=1.5,
                            label=f'Non-Metastatic (n={len(v_nonmeta)})')

            ax.set_title(f'{name}   [FDR = {fdr:.2e}]', fontsize=9, loc='left', pad=3)
            ax.set_ylabel('Density', fontsize=8)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            if i == 0:
                ax.legend(fontsize=8.5, loc='upper right', framealpha=0.85)
            else:
                legend = ax.get_legend()
                if legend:
                    legend.remove()

            if i < n - 1:
                ax.set_xlabel('')
                ax.tick_params(axis='x', labelbottom=False)
            else:
                ax.set_xlabel('Value (log1p transformed)', fontsize=10)

        fig.suptitle(
            f'Top {n} metabolites {direction.lower()} in Metastatic vs Non-Metastatic  '
            f'(Metastatic n={n_meta}, Non-Metastatic n={n_nonmeta})',
            fontsize=11, y=1.005
        )
        plt.tight_layout()
        fname = f'kde_top{n_top}_metabolites_{label}.png'
        plt.savefig(FIG_DIR / fname, dpi=300, bbox_inches='tight')
        plt.close()
        print(f'  Saved {fname}')

    return results


def plot_top_km_metabolites_panel(df, met_cols, n_top=20, n_cols=4, p_max=None, filename=None):
    """Rank metabolites by log-rank p-value (median split) and plot their KM curves
    as a single panel figure with n_cols columns.

    Selection is either fixed-size or threshold-based:
      * p_max is None -> the n_top most significant metabolites.
      * p_max set     -> every metabolite with log-rank p <= p_max, however many
                         that is. These are RAW log-rank p-values, not FDR: with
                         145 metabolites screened, a p <= 0.1 cut is a hypothesis-
                         generating filter and is expected to include false
                         positives (~14 by chance alone at this threshold).
    """
    if p_max is None:
        print(f"\nPlotting top {n_top} metabolite KM curves (median split)...")
    else:
        print(f"\nPlotting metabolite KM curves with log-rank p <= {p_max} (median split)...")

    df_cancer = df[df['Cohort'] == 'CANCER'].copy()
    df_cancer = df_cancer.dropna(subset=['Time to Event', 'Event'])
    df_cancer = df_cancer[df_cancer['Time to Event'] > 0]
    n_pts = len(df_cancer)

    short_name = met_label

    # ── Score every metabolite ────────────────────────────────────────────────
    scores = []
    for feat in met_cols:
        if feat not in df_cancer.columns:
            continue
        vals = np.log1p(pd.to_numeric(df_cancer[feat], errors='coerce').clip(lower=0))
        valid = vals.notna() & df_cancer['Time to Event'].notna() & df_cancer['Event'].notna()
        v = vals[valid]
        if v.nunique() < 3 or valid.sum() < 10:
            continue
        med = v.median()
        hi = valid & (vals > med)
        lo = valid & (vals <= med)
        if hi.sum() < 4 or lo.sum() < 4:
            continue
        try:
            res = logrank_test(
                df_cancer.loc[hi, 'Time to Event'], df_cancer.loc[lo, 'Time to Event'],
                event_observed_A=df_cancer.loc[hi, 'Event'],
                event_observed_B=df_cancer.loc[lo, 'Event'],
            )
            scores.append({'Feature': feat, 'p_value': res.p_value, 'median': med})
        except Exception:
            continue

    if not scores:
        print("  No scoreable metabolites found.")
        return

    all_scored = pd.DataFrame(scores).sort_values('p_value')
    if p_max is None:
        scored = all_scored.head(n_top).reset_index(drop=True)
        print(f"  Scored {len(all_scored)} metabolites; using top {len(scored)}")
    else:
        scored = all_scored[all_scored['p_value'] <= p_max].reset_index(drop=True)
        print(f"  Scored {len(all_scored)} metabolites; "
              f"{len(scored)} with log-rank p <= {p_max}")
        if scored.empty:
            print(f"  No metabolite reaches p <= {p_max} — nothing to plot.")
            return

    # ── Build panel ───────────────────────────────────────────────────────────
    n_feats = len(scored)
    n_rows  = int(np.ceil(n_feats / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 5, n_rows * 4.2),
                             squeeze=False)
    axes_flat = axes.flatten()

    hi_color = '#d62728'
    lo_color = '#1f77b4'

    for idx, row in scored.iterrows():
        ax   = axes_flat[idx]
        feat = row['Feature']
        med  = row['p_value']   # used below as p_val

        vals = np.log1p(pd.to_numeric(df_cancer[feat], errors='coerce').clip(lower=0))
        valid = vals.notna() & df_cancer['Time to Event'].notna() & df_cancer['Event'].notna()
        median_val = vals[valid].median()

        hi = valid & (vals > median_val)
        lo = valid & (vals <= median_val)

        T_hi = df_cancer.loc[hi, 'Time to Event']
        E_hi = df_cancer.loc[hi, 'Event']
        T_lo = df_cancer.loc[lo, 'Time to Event']
        E_lo = df_cancer.loc[lo, 'Event']

        kmf = KaplanMeierFitter()
        kmf.fit(T_hi, E_hi, label=f'High (n={hi.sum()})')
        kmf.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.12, color=hi_color)

        kmf2 = KaplanMeierFitter()
        kmf2.fit(T_lo, E_lo, label=f'Low (n={lo.sum()})')
        kmf2.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.12, color=lo_color)

        # Re-run log-rank to get p-value (scored DataFrame stored it under 'p_value')
        p_val = scored.loc[idx, 'p_value']

        name = short_name(feat)
        ax.set_title(f'{name}\np = {p_val:.3e}', fontsize=9, fontweight='bold')
        ax.set_xlim(0)
        ax.set_ylim(0, 1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Axis labels only on edge panels
        col_idx = idx % n_cols
        row_idx = idx // n_cols
        ax.set_xlabel('Days' if row_idx == n_rows - 1 else '', fontsize=8)
        ax.set_ylabel('Survival probability' if col_idx == 0 else '', fontsize=8)
        ax.tick_params(labelsize=7)

        legend = ax.get_legend()
        if legend:
            legend.set_title('')
            for t in legend.get_texts():
                t.set_fontsize(7.5)

    # Hide unused axes
    for i in range(n_feats, len(axes_flat)):
        axes_flat[i].axis('off')

    if p_max is None:
        suptitle = (f'Top {n_feats} metabolites — KM curves (median split)\n'
                    f'NSCLC patients n={n_pts}')
        default_name = f'km_top{n_feats}_metabolites_panel.png'
    else:
        suptitle = (f'{n_feats} metabolites with log-rank p ≤ {p_max} — KM curves '
                    f'(median split)\n'
                    f'NSCLC patients n={n_pts}   |   raw p-values, not FDR-corrected '
                    f'({len(all_scored)} metabolites screened)')
        default_name = f'km_metabolites_logrank_p{str(p_max).replace(".", "")}_panel.png'

    fig.suptitle(suptitle, fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout()
    out = FIG_DIR / (filename or default_name)
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out.name}')

    # Print ranking table
    print(f"\n{'='*65}")
    if p_max is None:
        print(f"TOP {n_feats} METABOLITES BY LOG-RANK P-VALUE (MEDIAN SPLIT)")
    else:
        print(f"{n_feats} METABOLITES WITH LOG-RANK P <= {p_max} (MEDIAN SPLIT)")
    print(f"{'='*65}")
    for _, r in scored.iterrows():
        print(f"  {r.name+1:>2}.  p={r['p_value']:.4e}   {short_name(r['Feature'])}")
    print(f"{'='*65}\n")


def plot_top_sdoh_km_panel(df, social_cols, n_top=10, n_cols=4, filename='km_top10_sdoh_panel.png'):
    """Plot KM survival curves for the top SDOH features by log-rank p-value on median split."""
    print(f"\nPlotting top {n_top} SDOH KM curves (median split)...")
    df_cancer = df[df['Cohort'] == 'CANCER'].dropna(subset=['Time to Event', 'Event']).copy()
    df_cancer = df_cancer[df_cancer['Time to Event'] > 0]
    n_pts = len(df_cancer)

    scores = []
    for feat in social_cols:
        vals = pd.to_numeric(df_cancer[feat], errors='coerce')
        valid = vals.notna()
        if valid.sum() < 10 or vals[valid].nunique() < 3:
            continue
        med = vals[valid].median()
        hi = valid & (vals > med)
        lo = valid & (vals <= med)
        if hi.sum() < 4 or lo.sum() < 4:
            continue
        try:
            res = logrank_test(df_cancer.loc[hi, 'Time to Event'], df_cancer.loc[lo, 'Time to Event'],
                               event_observed_A=df_cancer.loc[hi, 'Event'], event_observed_B=df_cancer.loc[lo, 'Event'])
            scores.append({'Feature': feat, 'p_value': res.p_value, 'median': med})
        except Exception:
            continue

    if not scores:
        return
    scored = pd.DataFrame(scores).sort_values('p_value').head(n_top).reset_index(drop=True)
    n_feats = len(scored)
    n_rows = int(np.ceil(n_feats / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4.5, n_rows * 4), squeeze=False)
    axes_flat = axes.flatten()

    for idx, row in scored.iterrows():
        ax = axes_flat[idx]
        feat = row['Feature']
        p_val = row['p_value']
        vals = pd.to_numeric(df_cancer[feat], errors='coerce')
        valid = vals.notna()
        med = vals[valid].median()
        hi = valid & (vals > med)
        lo = valid & (vals <= med)

        kmf_hi = KaplanMeierFitter().fit(df_cancer.loc[hi, 'Time to Event'], df_cancer.loc[hi, 'Event'], label=f'High (n={hi.sum()})')
        kmf_hi.plot_survival_function(ax=ax, color='#d62728', ci_show=True, ci_alpha=0.12)
        kmf_lo = KaplanMeierFitter().fit(df_cancer.loc[lo, 'Time to Event'], df_cancer.loc[lo, 'Event'], label=f'Low (n={lo.sum()})')
        kmf_lo.plot_survival_function(ax=ax, color='#1f77b4', ci_show=True, ci_alpha=0.12)

        ax.set_title(f'{feat}\np = {p_val:.3e}', fontsize=10, fontweight='bold')
        ax.set_xlim(0)
        ax.set_ylim(0, 1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        col_i = idx % n_cols
        row_i = idx // n_cols
        ax.set_xlabel('Days' if row_i == n_rows - 1 else '', fontsize=8.5)
        ax.set_ylabel('Survival probability' if col_i == 0 else '', fontsize=8.5)
        ax.tick_params(labelsize=8)

    for i in range(n_feats, len(axes_flat)):
        axes_flat[i].axis('off')

    plt.suptitle(f'Top {n_feats} Social Determinants of Health (SDOH) — KM Curves (Median Split)\nNSCLC patients n={n_pts}',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    out = FIG_DIR / filename
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out.name}')


def plot_dnam_km_panel(df, meth_cols, filename='km_dnam_panel.png'):
    """Plot KM survival curves for all DNA methylation markers on median split."""
    print("\nPlotting DNA Methylation KM curves (median split)...")
    df_cancer = df[df['Cohort'] == 'CANCER'].dropna(subset=['Time to Event', 'Event']).copy()
    df_cancer = df_cancer[df_cancer['Time to Event'] > 0]
    n_pts = len(df_cancer)

    n_feats = len(meth_cols)
    n_cols = 3
    n_rows = int(np.ceil(n_feats / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, n_rows * 4.2), squeeze=False)
    axes_flat = axes.flatten()

    for idx, feat in enumerate(meth_cols):
        ax = axes_flat[idx]
        vals = pd.to_numeric(df_cancer[feat], errors='coerce')
        valid = vals.notna()
        med = vals[valid].median()
        hi = valid & (vals > med)
        lo = valid & (vals <= med)

        res = logrank_test(df_cancer.loc[hi, 'Time to Event'], df_cancer.loc[lo, 'Time to Event'],
                           event_observed_A=df_cancer.loc[hi, 'Event'], event_observed_B=df_cancer.loc[lo, 'Event'])

        kmf_hi = KaplanMeierFitter().fit(df_cancer.loc[hi, 'Time to Event'], df_cancer.loc[hi, 'Event'], label=f'High (n={hi.sum()})')
        kmf_hi.plot_survival_function(ax=ax, color='#d62728', ci_show=True, ci_alpha=0.12)
        kmf_lo = KaplanMeierFitter().fit(df_cancer.loc[lo, 'Time to Event'], df_cancer.loc[lo, 'Event'], label=f'Low (n={lo.sum()})')
        kmf_lo.plot_survival_function(ax=ax, color='#1f77b4', ci_show=True, ci_alpha=0.12)

        clean_meth = feat.replace('. P', '').replace('.P', '').replace('.', '').strip()
        ax.set_title(f'DNAm: {clean_meth}\np = {res.p_value:.3e}', fontsize=10.5, fontweight='bold')
        ax.set_xlim(0)
        ax.set_ylim(0, 1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        col_i = idx % n_cols
        row_i = idx // n_cols
        ax.set_xlabel('Days' if row_i == n_rows - 1 else '', fontsize=9)
        ax.set_ylabel('Survival probability' if col_i == 0 else '', fontsize=9)
        ax.tick_params(labelsize=8)

    for i in range(n_feats, len(axes_flat)):
        axes_flat[i].axis('off')

    plt.suptitle(f'DNA Methylation Markers (n={n_feats}) — KM Curves (Median Split)\nNSCLC patients n={n_pts}',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    out = FIG_DIR / filename
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out.name}')


def plot_clinical_vars_km_panel(df, filename='km_clinical_vars_panel.png'):
    """Plot KM survival curves for key clinical and demographic variables."""
    print("\nPlotting Clinical & Demographic KM curves...")
    df_cancer = df[df['Cohort'] == 'CANCER'].dropna(subset=['Time to Event', 'Event']).copy()
    df_cancer = df_cancer[df_cancer['Time to Event'] > 0]
    n_pts = len(df_cancer)

    clin_vars = [
        ('Asthma', 'Asthma (Yes vs No)', 'cat'),
        ('HTN', 'Hypertension (HTN)', 'cat'),
        ('Sex', 'Sex (Female vs Male)', 'cat'),
        ('CHF', 'Congestive Heart Failure (CHF)', 'cat'),
        ('Race', 'Race (White vs Black vs Other)', 'cat'),
        ('Stage', 'Stage (Non-Metastatic vs Metastatic)', 'cat'),
        ('Size 1 of tumor (cms)', 'Tumor Size (cms)', 'cont'),
        ('Age at the time of diagnosis', 'Age at Diagnosis', 'cont'),
        ('Pack Year History', 'Pack Year History', 'cont'),
        ('CAD', 'Coronary Artery Disease (CAD)', 'cat')
    ]

    n_cols = 4
    n_rows = int(np.ceil(len(clin_vars) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4.5, n_rows * 4), squeeze=False)
    axes_flat = axes.flatten()

    for idx, (feat, display_name, vtype) in enumerate(clin_vars):
        ax = axes_flat[idx]
        if feat not in df_cancer.columns:
            continue
        temp = df_cancer[[feat, 'Time to Event', 'Event']].dropna().copy()
        if temp.empty:
            continue

        if vtype == 'cont':
            vals = pd.to_numeric(temp[feat], errors='coerce')
            valid = vals.notna()
            med = vals[valid].median()
            hi = valid & (vals > med)
            lo = valid & (vals <= med)
            res = logrank_test(temp.loc[hi, 'Time to Event'], temp.loc[lo, 'Time to Event'],
                               event_observed_A=temp.loc[hi, 'Event'], event_observed_B=temp.loc[lo, 'Event'])
            p_val = res.p_value
            kmf_hi = KaplanMeierFitter().fit(temp.loc[hi, 'Time to Event'], temp.loc[hi, 'Event'], label=f'High (>{med:.1f}, n={hi.sum()})')
            kmf_hi.plot_survival_function(ax=ax, color='#d62728', ci_show=True, ci_alpha=0.12)
            kmf_lo = KaplanMeierFitter().fit(temp.loc[lo, 'Time to Event'], temp.loc[lo, 'Event'], label=f'Low (≤{med:.1f}, n={lo.sum()})')
            kmf_lo.plot_survival_function(ax=ax, color='#1f77b4', ci_show=True, ci_alpha=0.12)
        else:
            cats = [c for c in temp[feat].astype(str).str.strip().unique() if c not in ['nan', 'None', 'NA', 'NAN', '']]
            if len(cats) < 2:
                continue
            res = multivariate_logrank_test(temp['Time to Event'], temp[feat].astype(str), temp['Event'])
            p_val = res.p_value
            palette = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e']
            for c_idx, cat in enumerate(sorted(cats)):
                c_mask = temp[feat].astype(str).str.strip() == cat
                kmf = KaplanMeierFitter().fit(temp.loc[c_mask, 'Time to Event'], temp.loc[c_mask, 'Event'], label=f'{cat} (n={c_mask.sum()})')
                kmf.plot_survival_function(ax=ax, color=palette[c_idx % len(palette)], ci_show=True, ci_alpha=0.10)

        ax.set_title(f'{display_name}\np = {p_val:.3e}', fontsize=10, fontweight='bold')
        ax.set_xlim(0)
        ax.set_ylim(0, 1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        col_i = idx % n_cols
        row_i = idx // n_cols
        ax.set_xlabel('Days' if row_i == n_rows - 1 else '', fontsize=8.5)
        ax.set_ylabel('Survival probability' if col_i == 0 else '', fontsize=8.5)
        ax.tick_params(labelsize=8)

    for i in range(len(clin_vars), len(axes_flat)):
        axes_flat[i].axis('off')

    plt.suptitle(f'Key Clinical & Demographic Predictors — KM Curves\nNSCLC patients n={n_pts}', fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    out = FIG_DIR / filename
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  Saved {out.name}')


def analyze_pretreatment_vs_control(df, met_cols, n_top=15):
    """Cancer vs control differential analysis restricted to pre-treatment blood draws.
    A cancer patient is pre-treatment if their baseline blood draw preceded the start
    of every treatment they received (chemo, IO, SBRT).
    Produces: volcano plot, boxplots (top increased / decreased), KDE distributions.
    """
    print("\n" + "="*70)
    print("PRE-TREATMENT CANCER vs CONTROL DIFFERENTIAL METABOLITE ANALYSIS")
    print("="*70)

    draw_col  = 'Date of baseline blood draw (visit)'
    chemo_col = 'Start Date of Chemotherapy'
    io_col    = 'date of IO  Start'
    sbrt_col  = 'SBRT start date'
    chemo_rec = 'Chemo receive  yes or none'
    io_rec    = 'Immunotherapy Receive  yes or none'
    sbrt_rec  = 'SBRT'

    df_cancer = df[df['Cohort'] == 'CANCER'].copy()
    for col in [draw_col, chemo_col, io_col, sbrt_col]:
        if col in df_cancer.columns:
            df_cancer[col] = pd.to_datetime(df_cancer[col], errors='coerce')

    def received(col):
        return df_cancer[col].astype(str).str.strip().str.upper().str.contains('YES', na=False)

    chemo_yes = received(chemo_rec) if chemo_rec in df_cancer.columns else pd.Series(False, index=df_cancer.index)
    io_yes    = received(io_rec)    if io_rec    in df_cancer.columns else pd.Series(False, index=df_cancer.index)
    sbrt_yes  = received(sbrt_rec)  if sbrt_rec  in df_cancer.columns else pd.Series(False, index=df_cancer.index)

    draw = df_cancer[draw_col]
    days_chemo = (draw - df_cancer[chemo_col]).dt.days if chemo_col in df_cancer.columns else pd.Series(np.nan, index=df_cancer.index)
    days_io    = (draw - df_cancer[io_col]).dt.days    if io_col    in df_cancer.columns else pd.Series(np.nan, index=df_cancer.index)
    days_sbrt  = (draw - df_cancer[sbrt_col]).dt.days  if sbrt_col  in df_cancer.columns else pd.Series(np.nan, index=df_cancer.index)

    # Pre-treatment: for every received treatment, draw was BEFORE start (days < 0)
    post_chemo = chemo_yes & (days_chemo >= 0)
    post_io    = io_yes    & (days_io    >= 0)
    post_sbrt  = sbrt_yes  & (days_sbrt  >= 0)
    is_pre_treatment = ~(post_chemo | post_io | post_sbrt)

    df_pre = df_cancer[is_pre_treatment].copy()
    df_pre['Cohort'] = 'CANCER'
    df_ctrl = df[df['Cohort'] == 'No'].copy()

    n_pre  = len(df_pre)
    n_ctrl = len(df_ctrl)
    n_post = len(df_cancer) - n_pre
    print(f"  Pre-treatment cancer patients : {n_pre}")
    print(f"  Post-treatment cancer patients: {n_post}  (excluded)")
    print(f"  Controls                      : {n_ctrl}")

    df_combined = pd.concat([df_pre, df_ctrl], ignore_index=True)

    # Differential analysis
    results = rank_sum_test(df_combined, met_cols, is_log=True)
    if results is None or results.empty:
        print("  No results from rank-sum test.")
        return pd.DataFrame()

    n_sig = (results['FDR'] < 0.05).sum()
    n_inc = ((results['FDR'] < 0.05) & (results['Direction'] == 'Increased')).sum()
    n_dec = ((results['FDR'] < 0.05) & (results['Direction'] == 'Decreased')).sum()
    print(f"  FDR<0.05: {n_sig} metabolites  ({n_inc} higher in cancer, {n_dec} lower in cancer)")

    # Volcano plot
    plot_volcano(
        results,
        f'Pre-treatment Cancer (n={n_pre}) vs Control (n={n_ctrl})',
        'volcano_pretreatment_vs_control.png',
        group_a='Pre-treatment Cancer', group_b='Control'
    )
    plot_detection_tier(
        results,
        f'Metabolites detected at different rates: pre-treatment cancer (n={n_pre}) vs control (n={n_ctrl})',
        'detection_pretreatment_vs_control.png',
        group_a='Pre-treatment Cancer', group_b='Control'
    )

    # Boxplots — top increased
    inc = results[(results['FDR'] < 0.05) & (results['Direction'] == 'Increased')].sort_values('FDR').head(n_top)
    if inc.empty:
        inc = results[results['Direction'] == 'Increased'].sort_values('p_value').head(n_top)
        print(f"  No FDR-significant increases; using top {n_top} nominal.")
    plot_boxplots(
        df_combined, inc,
        f'Top metabolites higher in pre-treatment cancer vs control  (n_cancer={n_pre}, n_ctrl={n_ctrl})',
        'boxplot_pretreatment_increased.png',
        is_log=True, hue_col='Cohort', hue_order=['CANCER', 'No']
    )

    # Boxplots — top decreased
    dec = results[(results['FDR'] < 0.05) & (results['Direction'] == 'Decreased')].sort_values('FDR').head(n_top)
    if dec.empty:
        dec = results[results['Direction'] == 'Decreased'].sort_values('p_value').head(n_top)
        print(f"  No FDR-significant decreases; using top {n_top} nominal.")
    plot_boxplots(
        df_combined, dec,
        f'Top metabolites lower in pre-treatment cancer vs control  (n_cancer={n_pre}, n_ctrl={n_ctrl})',
        'boxplot_pretreatment_decreased.png',
        is_log=True, hue_col='Cohort', hue_order=['CANCER', 'No']
    )

    # KDE distribution plots — prefix keeps filenames distinct from full-cohort plots
    plot_metabolite_distributions(df_combined, results, n_top=10, prefix='pretreatment_')

    return results


def pair_baseline_followup(met_cols, all_cols):
    """Return [(baseline_col, followup_col), ...] by prefix-swapping the shared suffix.
    Baseline_<suffix> <-> Follow-Up_<suffix>. Only pairs present in all_cols are returned.
    """
    all_set = set(all_cols)
    pairs = []
    for b in met_cols:
        if not b.startswith('Baseline_'):
            continue
        fu = 'Follow-Up_' + b[len('Baseline_'):]
        if fu in all_set:
            pairs.append((b, fu))
    return pairs





def _paired_change_table(df_grp, pairs, min_pairs=5, eps=None):
    """Per-metabolite paired pre->post change with Wilcoxon signed-rank + BH-FDR.
    Returns a DataFrame with Feature, n_pairs, Median_Diff (log1p), log2FC, Direction, p_value, FDR.
    """
    if eps is None:
        eps = LOG2FC_EPS
    rows = []
    for b, f in pairs:
        rb = pd.to_numeric(df_grp[b], errors='coerce')
        rf = pd.to_numeric(df_grp[f], errors='coerce')
        m = rb.notna() & rf.notna()
        rb, rf = rb[m].clip(lower=0), rf[m].clip(lower=0)
        n = len(rb)
        if n == 0:
            continue
        lb, lf = np.log1p(rb), np.log1p(rf)
        diff = (lf - lb)
        med_diff = float(diff.median())
        log2fc = float(np.median(np.log2((rf + eps) / (rb + eps))))
        p = np.nan
        if n >= min_pairs and int((diff != 0).sum()) >= 1:
            try:
                p = wilcoxon(lb, lf).pvalue
            except Exception:
                p = np.nan
        rows.append({'Feature': b, 'n_pairs': n, 'Median_Diff': med_diff, 'log2FC': log2fc,
                     'Direction': 'Increased' if med_diff > 0 else 'Decreased', 'p_value': p})
    r = pd.DataFrame(rows)
    r['FDR'] = np.nan
    if not r.empty and r['p_value'].notna().any():
        valid = r['p_value'].notna()
        r.loc[valid, 'FDR'] = fdrcorrection(r.loc[valid, 'p_value'].values)[1]
    return r


def _paired_boxplot(df_grp, pairs_sel, title, filename, stats_map=None):
    """Grouped pre-vs-post boxplot (log1p) for the selected metabolites."""
    long = []
    for b, f in pairs_sel:
        lbl = met_label(b)
        # Paired design: the fold change is the median of the per-patient post/pre
        # ratios, not a ratio of group medians, so it matches the Wilcoxon above.
        if stats_map is not None and b in stats_map:
            l2 = stats_map[b]
            if pd.notna(l2):
                lbl += f"\n{2 ** l2:.2f}x"
        rb = np.log1p(pd.to_numeric(df_grp[b], errors='coerce').clip(lower=0))
        rf = np.log1p(pd.to_numeric(df_grp[f], errors='coerce').clip(lower=0))
        mask = rb.notna() & rf.notna()
        for v in rb[mask]:
            long.append({'Metabolite': lbl, 'Timepoint': 'Baseline (pre)', 'Value': v})
        for v in rf[mask]:
            long.append({'Metabolite': lbl, 'Timepoint': 'Follow-Up (post)', 'Value': v})
    if not long:
        return
    ldf = pd.DataFrame(long)
    plt.figure(figsize=(min(16, 2 + 1.2 * len(pairs_sel)), 7))
    sns.boxplot(data=ldf, x='Metabolite', y='Value', hue='Timepoint',
                hue_order=['Baseline (pre)', 'Follow-Up (post)'],
                palette={'Baseline (pre)': '#7570b3', 'Follow-Up (post)': '#d95f02'},
                showfliers=False, width=0.7, linewidth=1.0)
    plt.title(title, fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Value (log1p)')
    plt.xlabel('')
    plt.legend(title='', bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()
    plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()


def analyze_pre_post_treatment(df, met_cols, io_col='IO Responder Yes or No', min_pairs=5, n_box=8):
    """Paired BEFORE (baseline) vs AFTER (follow-up) treatment metabolite analysis, cancer patients.
    (a) paired pre->post change (any-treatment cohort); (b) same within every treatment regimen group;
    (c) association of per-patient change (delta) with survival and IO response.
    Genuinely post-treatment = follow-up draw drawn AFTER >=1 received treatment (timing flags).
    Returns (change_df, by_treatment_df, delta_io_df).
    """
    print("\n" + "="*70)
    print("BEFORE vs AFTER TREATMENT (paired) METABOLITE ANALYSIS")
    print("="*70)

    dfc = df[df['Cohort'] == 'CANCER'].copy()
    pairs = pair_baseline_followup(met_cols, df.columns)
    print(f"  Paired baseline<->follow-up metabolites available: {len(pairs)}")

    # --- post-treatment flags (follow-up draw AFTER each therapy) ---
    def flag_after(col):
        if col not in dfc.columns:
            return pd.Series(False, index=dfc.index)
        return dfc[col].astype(str).str.strip().str.upper().str.contains('YES', na=False)
    has_fu = dfc[FU_DRAW_DATE_COL].notna() if FU_DRAW_DATE_COL in dfc.columns else pd.Series(False, index=dfc.index)
    af_io, af_ch, af_sb = flag_after(FU_AFTER_IO_COL), flag_after(FU_AFTER_CHEMO_COL), flag_after(FU_AFTER_SBRT_COL)

    # Data-quality note: IO and chemo timing flags are often identical in this sheet
    if int((af_io != af_ch).sum()) == 0 and int(af_io.sum()) > 0:
        print("  [warn] IO and Chemo 'after' flags are identical for all patients — possible copy-paste in source sheet.")

    groups = {
        'Any treatment':        (has_fu & (af_io | af_ch | af_sb),            'volcano_pre_post_treatment.png'),
        'After IO':             (has_fu & af_io,                              'volcano_pre_post_after_io.png'),
        'After Chemo':          (has_fu & af_ch,                              'volcano_pre_post_after_chemo.png'),
        'After SBRT':           (has_fu & af_sb,                              'volcano_pre_post_after_sbrt.png'),
        'IO+Chemo+SBRT':        (has_fu & af_io & af_ch & af_sb,             'volcano_pre_post_regimen_io_chemo_sbrt.png'),
        'IO+Chemo (no SBRT)':   (has_fu & af_io & af_ch & ~af_sb,            'volcano_pre_post_regimen_io_chemo.png'),
        'SBRT only':            (has_fu & ~af_io & ~af_ch & af_sb,           'volcano_pre_post_regimen_sbrt_only.png'),
    }

    by_treatment = []
    sig_sets = {}
    change_df = None
    for name, (mask, vfile) in groups.items():
        grp = dfc[mask]
        n = len(grp)
        res = _paired_change_table(grp, pairs, min_pairs=min_pairs)
        n_valid = int(res['p_value'].notna().sum()) if not res.empty else 0
        n_sig = int((res['FDR'] < 0.05).sum()) if 'FDR' in res.columns else 0
        print(f"\n  [{name}] n_patients={n} | metabolites tested (n_pairs>={min_pairs}): {n_valid} | FDR<0.05: {n_sig}")
        if res.empty:
            continue
        res = res.copy(); res['Regimen'] = name
        by_treatment.append(res)
        sig_sets[name] = set(res.loc[res['FDR'] < 0.05, 'Feature']) if n_valid else set()

        if name == 'Any treatment':
            change_df = res.copy()

        # Volcano only when we have valid FDRs; else descriptive-only note
        if n_valid > 0:
            plot_volcano(res, f"Pre vs Post — {name} (n={n})", vfile,
                         group_a='Higher post', group_b='Higher pre')
            # Paired boxplots for top changers
            ranked = res[res['p_value'].notna()].sort_values('p_value')
            sel = ranked[ranked['FDR'] < 0.05].head(n_box)
            if sel.empty:
                sel = ranked.head(n_box)
            sel_pairs = [(b, 'Follow-Up_' + b[len('Baseline_'):]) for b in sel['Feature']]
            bfile = vfile.replace('volcano_', 'boxplot_top_').replace('.png', '.png')
            _paired_boxplot(grp, sel_pairs, f"Top pre->post changes — {name} (n={n})", bfile,
                            stats_map=dict(zip(sel['Feature'], sel['log2FC'])))
        else:
            print(f"    (descriptive only — no group reached n_pairs>={min_pairs}; volcano skipped)")

    by_treatment_df = pd.concat(by_treatment, ignore_index=True) if by_treatment else pd.DataFrame()

    # --- regimen overlap: what each regimen adds ---
    def _report_diff(a, b):
        if a in sig_sets and b in sig_sets:
            only_a = sig_sets[a] - sig_sets[b]
            print(f"  Changed in [{a}] but not [{b}]: {len(only_a)}"
                  + (f" e.g. {[met_label(x) for x in list(only_a)[:5]]}" if only_a else ""))
    print("\n  --- Regimen significant-set comparison (FDR<0.05) ---")
    for k, v in sig_sets.items():
        print(f"    {k}: {len(v)} significant")
    _report_diff('IO+Chemo+SBRT', 'IO+Chemo (no SBRT)')
    _report_diff('IO+Chemo (no SBRT)', 'SBRT only')

    # ---------------- (c) change vs survival & IO response ----------------
    delta_io_df = pd.DataFrame()
    post_mask = has_fu & (af_io | af_ch | af_sb)
    post = dfc[post_mask].copy()
    if change_df is not None and not change_df.empty and len(post) >= min_pairs:
        # metabolites to track: FDR<0.05 in any-treatment, else top 15 by p
        sig = change_df[change_df['FDR'] < 0.05]['Feature'].tolist()
        if not sig:
            sig = change_df[change_df['p_value'].notna()].sort_values('p_value')['Feature'].head(15).tolist()
        # per-patient delta columns
        delta_cols = []
        for b in sig:
            f = 'Follow-Up_' + b[len('Baseline_'):]
            if f not in post.columns:
                continue
            d = np.log1p(pd.to_numeric(post[f], errors='coerce').clip(lower=0)) - \
                np.log1p(pd.to_numeric(post[b], errors='coerce').clip(lower=0))
            dc = 'Delta_' + met_label(b)
            post[dc] = d
            delta_cols.append(dc)
        if delta_cols:
            post['Mean_Signed_Change'] = post[delta_cols].mean(axis=1)

            # Survival: mean signed change vs Time to Event (Spearman)
            sv = post[['Mean_Signed_Change', 'Time to Event', 'Event']].dropna(subset=['Mean_Signed_Change', 'Time to Event'])
            if len(sv) >= 5:
                rho, pv = spearmanr(sv['Mean_Signed_Change'], sv['Time to Event'])
                plt.figure(figsize=(7, 6))
                colors = sv['Event'].map({1: '#d95f02', 0: '#1b9e77'}) if 'Event' in sv.columns else '#333333'
                plt.scatter(sv['Mean_Signed_Change'], sv['Time to Event'], c=colors, alpha=0.8, s=40)
                plt.title(f"Mean metabolite change vs survival (post-tx, n={len(sv)})\nSpearman rho={rho:.2f}, p={pv:.3g}", fontsize=12)
                plt.xlabel("Mean signed change (log1p, post - pre)")
                plt.ylabel("Time to Event (days)")
                plt.tight_layout()
                plt.savefig(FIG_DIR / 'scatter_change_vs_survival.png', dpi=300, bbox_inches='tight')
                plt.close()
                print(f"\n  Change vs survival: Spearman rho={rho:.2f}, p={pv:.3g} (n={len(sv)})")

            # IO response: do metabolite CHANGES differ by responder?
            if io_col in post.columns:
                dio = post.dropna(subset=[io_col]).copy()
                dio[io_col] = dio[io_col].astype(str).str.strip().str.title()
                dio = dio[dio[io_col].isin(['Yes', 'No'])]
                n_yes = int((dio[io_col] == 'Yes').sum())
                n_no  = int((dio[io_col] == 'No').sum())
                # rank_sum_test requires >5 per group
                if n_yes > 5 and n_no > 5:
                    delta_io_df = rank_sum_test(dio, delta_cols, is_log=False,
                                                group_col=io_col, group1='Yes', group2='No')
                    if delta_io_df is not None and not delta_io_df.empty:
                        top = delta_io_df.sort_values('p_value').head(min(n_box, len(delta_io_df)))
                        plot_boxplots(dio, top,
                                      "Metabolite CHANGE (post-pre) by IO response",
                                      "boxplot_delta_metabolites_io_responder.png",
                                      is_log=False, hue_col=io_col, hue_order=['Yes', 'No'])
                        print(f"  IO-response delta test: {int((delta_io_df['FDR']<0.05).sum())} FDR<0.05 of {len(delta_io_df)}")
                else:
                    print(f"  IO-response delta test SKIPPED: too few responders in post-tx cohort "
                          f"(Yes={n_yes}, No={n_no}; need >5 each).")

    return change_df, by_treatment_df, delta_io_df


def analyze_metabolite_normalization_score(df, met_cols, pretreatment_results=None, n_top_mets=50):
    """Analysis F: Metabolite Normalization Score.
    NS(i,m) = (post_log − pre_log) * sign(ctrl_mean − pre_log) / (|ctrl_mean − pre_log| + eps)
    NS = +1 → fully normalized toward healthy; NS = 0 → no change; NS < 0 → moved away.
    Produces heatmap, scatter vs survival, KM, and bar chart by treatment.
    """
    import re

    print("\n" + "="*70)
    print("ANALYSIS F: METABOLITE NORMALIZATION SCORE")
    print("="*70)

    all_cols     = df.columns.tolist()
    kegg_pat     = re.compile(r'C\d{5}')

    # ── Identify follow-up metabolite columns (prefix Follow-Up_) ─────────────
    fu_by_kegg = {}
    for col in all_cols:
        if not col.startswith('Follow-Up_'):
            continue
        m = kegg_pat.search(col)
        if m and m.group() not in fu_by_kegg:
            fu_by_kegg[m.group()] = col

    print(f"  Follow-up metabolite columns found: {len(fu_by_kegg)}")

    paired_mets = [(b, fu_by_kegg[kegg_pat.search(b).group()])
                   for b in met_cols
                   if kegg_pat.search(b) and kegg_pat.search(b).group() in fu_by_kegg]
    print(f"  Metabolites matched baseline↔follow-up: {len(paired_mets)}")

    if not paired_mets:
        print("  No matched metabolite pairs — skipping Analysis F.")
        return pd.DataFrame()

    # ── Identify paired patients ──────────────────────────────────────────────
    df_cancer = df[df['Cohort'] == 'CANCER'].copy()
    df_ctrl   = df[df['Cohort'] == 'No'].copy()

    fu_draw_col = 'Baseline_C00669_Follow up Date blood was drawn'
    bl_draw_col = 'Date of baseline blood draw (visit)'
    for col in [fu_draw_col, bl_draw_col]:
        if col in df_cancer.columns:
            df_cancer[col] = pd.to_datetime(df_cancer[col], errors='coerce')

    # Paired = has a recorded follow-up blood draw date (proxy for having follow-up metabolites)
    if fu_draw_col in df_cancer.columns:
        has_fu = df_cancer[fu_draw_col].notna()
    else:
        # Fallback: any follow-up metabolite is non-null
        fu_sample = list(fu_by_kegg.values())[:10]
        has_fu = df_cancer[[c for c in fu_sample if c in df_cancer.columns]].apply(
            pd.to_numeric, errors='coerce').notna().any(axis=1)
    df_paired = df_cancer[has_fu].copy()
    print(f"  Paired patients (both draws): {len(df_paired)}")

    if len(df_paired) < 5:
        print("  Too few paired patients — skipping.")
        return pd.DataFrame()

    # ── Select dysregulated metabolites ──────────────────────────────────────
    if pretreatment_results is not None and not pretreatment_results.empty:
        sig_set = set(pretreatment_results[pretreatment_results['FDR'] < 0.05]['Feature'])
        sel = [(b, f) for b, f in paired_mets if b in sig_set]
        if len(sel) < 5:
            top_set = set(pretreatment_results.sort_values('p_value').head(n_top_mets)['Feature'])
            sel = [(b, f) for b, f in paired_mets if b in top_set]
        if sel:
            paired_mets = sel
    print(f"  Metabolites used in normalization score: {len(paired_mets)}")

    # ── Control means (log1p scale) ───────────────────────────────────────────
    ctrl_means = {}
    for base_col, _ in paired_mets:
        v = np.log1p(pd.to_numeric(df_ctrl[base_col], errors='coerce').clip(lower=0))
        if v.notna().sum() > 3:
            ctrl_means[base_col] = float(v.mean())

    # ── Compute per-patient normalization scores ──────────────────────────────
    eps       = 0.05
    chemo_rec = 'Chemo receive  yes or none'
    io_rec    = 'Immunotherapy Receive  yes or none'
    sbrt_rec  = 'SBRT'

    def trt_label(row):
        parts = []
        for col, name in [(chemo_rec, 'Chemo'), (io_rec, 'IO'), (sbrt_rec, 'SBRT')]:
            if col in row.index and 'YES' in str(row[col]).upper():
                parts.append(name)
        return '+'.join(parts) if parts else 'None'

    ns_rows    = []
    per_met_ns = {}  # pid -> {base_col: ns_val}

    for _, prow in df_paired.iterrows():
        pid = str(prow.get('Patient ID', '')) if pd.notna(prow.get('Patient ID', np.nan)) else 'unk'
        ns_vals = {}
        for base_col, fu_col in paired_mets:
            if base_col not in ctrl_means:
                continue
            pre = pd.to_numeric(prow.get(base_col), errors='coerce')
            fu  = pd.to_numeric(prow.get(fu_col),   errors='coerce')
            if pd.isna(pre) or pd.isna(fu):
                continue
            pre_log  = np.log1p(max(0.0, pre))
            fu_log   = np.log1p(max(0.0, fu))
            ctrl_mu  = ctrl_means[base_col]
            denom    = abs(ctrl_mu - pre_log) + eps
            ns       = (fu_log - pre_log) * np.sign(ctrl_mu - pre_log) / denom
            ns_vals[base_col] = float(np.clip(ns, -3, 3))

        if len(ns_vals) < 3:
            continue

        tte   = pd.to_numeric(prow.get('Time to Event'), errors='coerce')
        event = pd.to_numeric(prow.get('Event'),         errors='coerce')

        ns_rows.append({
            'Patient_ID':    pid,
            'Mean_NS':       float(np.mean(list(ns_vals.values()))),
            'Time_to_Event': tte,
            'Event':         event,
            'Treatment':     trt_label(prow),
            'n_mets':        len(ns_vals),
        })
        per_met_ns[pid] = ns_vals

    if not ns_rows:
        print("  No normalization scores computed.")
        return pd.DataFrame()

    ns_df = pd.DataFrame(ns_rows).dropna(subset=['Time_to_Event'])
    print(f"  Normalization scores computed for {len(ns_df)} patients")

    if len(ns_df) < 5:
        print("  Too few patients — skipping plots.")
        return ns_df

    # ── Figure 1: Heatmap ─────────────────────────────────────────────────────
    ns_sorted  = ns_df.sort_values('Time_to_Event').reset_index(drop=True)
    pid_order  = ns_sorted['Patient_ID'].tolist()
    met_order  = [b for b, _ in paired_mets if b in ctrl_means]

    heat_data = np.full((len(pid_order), len(met_order)), np.nan)
    for ri, pid in enumerate(pid_order):
        for ci, base_col in enumerate(met_order):
            if pid in per_met_ns and base_col in per_met_ns[pid]:
                heat_data[ri, ci] = per_met_ns[pid][base_col]

    short_labels = [met_label(c, maxlen=22) for c in met_order]

    fig_w = min(28, max(12, len(met_order) * 0.35))
    fig_h = min(22, max(6,  len(pid_order) * 0.30))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heat_data, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ylab = [f"{r['Patient_ID']}  {r['Time_to_Event']:.0f}d" for _, r in ns_sorted.iterrows()]
    ax.set_yticks(range(len(pid_order)))
    ax.set_yticklabels(ylab, fontsize=6)
    ax.set_xticks(range(len(met_order)))
    ax.set_xticklabels(short_labels, rotation=90, fontsize=5.5)
    ax.set_xlabel('Metabolite', fontsize=10)
    ax.set_ylabel('Patient (sorted by survival, top = shortest)', fontsize=10)
    ax.set_title(
        'Metabolite Normalization Score Heatmap\n'
        '(+1 = fully normalized toward healthy, −1 = moved away)',
        fontsize=11, fontweight='bold'
    )
    plt.colorbar(im, ax=ax, label='Normalization Score', shrink=0.6)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'heatmap_normalization_score.png', dpi=200, bbox_inches='tight')
    plt.close()
    print('  Saved heatmap_normalization_score.png')

    # ── Figure 2: Scatter — Mean NS vs Survival ───────────────────────────────
    rho, pval = spearmanr(ns_df['Mean_NS'], ns_df['Time_to_Event'])

    fig, ax = plt.subplots(figsize=(8, 6))
    event_colors = ns_df['Event'].map({1: '#d62728', 0: '#1f77b4'}).fillna('#aaaaaa')
    ax.scatter(ns_df['Mean_NS'], ns_df['Time_to_Event'],
               c=event_colors, s=70, alpha=0.85, edgecolors='k', linewidth=0.4)
    if len(ns_df) >= 4:
        m_fit, b_fit = np.polyfit(ns_df['Mean_NS'], ns_df['Time_to_Event'], 1)
        x_fit = np.linspace(ns_df['Mean_NS'].min(), ns_df['Mean_NS'].max(), 100)
        ax.plot(x_fit, m_fit * x_fit + b_fit, 'k--', linewidth=1.2, alpha=0.6)
    for _, r in ns_df.iterrows():
        ax.annotate(str(r['Patient_ID']), (r['Mean_NS'], r['Time_to_Event']),
                    fontsize=6, alpha=0.55, xytext=(2, 2), textcoords='offset points')
    legend_els = [
        mpatches.Patch(facecolor='#d62728', label='Dead'),
        mpatches.Patch(facecolor='#1f77b4', label='Censored / Alive'),
    ]
    ax.legend(handles=legend_els, fontsize=10)
    ax.set_xlabel('Mean Normalization Score', fontsize=12)
    ax.set_ylabel('Survival Time (days)', fontsize=12)
    ax.set_title(
        f'Metabolite Normalization Score vs Survival Time\n'
        f'Spearman ρ = {rho:.3f}  (p = {pval:.3e})  n = {len(ns_df)}',
        fontsize=11, fontweight='bold'
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'scatter_normalization_score_vs_survival.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Saved scatter_normalization_score_vs_survival.png  (ρ={rho:.3f}, p={pval:.3e})')

    # ── Figure 3: KM — High vs Low Normalizers ────────────────────────────────
    km_valid = ns_df.dropna(subset=['Mean_NS', 'Time_to_Event', 'Event']).copy()
    if len(km_valid) >= 8:
        med_ns = km_valid['Mean_NS'].median()
        km_valid['NS_Group'] = np.where(km_valid['Mean_NS'] > med_ns, 'High NS', 'Low NS')
        T_hi = km_valid.loc[km_valid['NS_Group'] == 'High NS', 'Time_to_Event']
        E_hi = km_valid.loc[km_valid['NS_Group'] == 'High NS', 'Event']
        T_lo = km_valid.loc[km_valid['NS_Group'] == 'Low NS',  'Time_to_Event']
        E_lo = km_valid.loc[km_valid['NS_Group'] == 'Low NS',  'Event']
        try:
            lr_p = logrank_test(T_hi, T_lo, event_observed_A=E_hi, event_observed_B=E_lo).p_value
        except Exception:
            lr_p = np.nan
        p_str = f'p = {lr_p:.3f}' if pd.notna(lr_p) else 'p = N/A'

        fig, ax = plt.subplots(figsize=(8, 6))
        KaplanMeierFitter().fit(T_hi, E_hi, label=f'High NS (n={len(T_hi)})').plot_survival_function(
            ax=ax, ci_show=True, ci_alpha=0.15, color='#2ca02c')
        KaplanMeierFitter().fit(T_lo, E_lo, label=f'Low NS (n={len(T_lo)})').plot_survival_function(
            ax=ax, ci_show=True, ci_alpha=0.15, color='#d62728')
        ax.set_title(
            f'KM: High vs Low Metabolite Normalization Score\n{p_str}  (median split, n={len(km_valid)})',
            fontsize=11, fontweight='bold'
        )
        ax.set_xlabel('Days', fontsize=11)
        ax.set_ylabel('Survival Probability', fontsize=11)
        ax.set_ylim(0, 1.05)
        ax.set_xlim(0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        plt.savefig(FIG_DIR / 'km_normalization_score.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f'  Saved km_normalization_score.png  (log-rank {p_str})')

    # ── Figure 4: Bar chart — Mean NS by Treatment ────────────────────────────
    trt_ns = (ns_df.groupby('Treatment')['Mean_NS']
              .agg(['mean', 'sem', 'count'])
              .reset_index()
              .rename(columns={'mean': 'Mean', 'sem': 'SEM', 'count': 'N'}))
    trt_ns = trt_ns[trt_ns['N'] >= 2].sort_values('Mean', ascending=False)

    if len(trt_ns) >= 2:
        fig, ax = plt.subplots(figsize=(max(6, len(trt_ns) * 1.4), 5))
        colors_bar = plt.cm.Set2(np.linspace(0, 1, len(trt_ns)))
        ax.bar(range(len(trt_ns)), trt_ns['Mean'],
                      yerr=trt_ns['SEM'], capsize=5,
                      color=colors_bar, edgecolor='k', linewidth=0.6, alpha=0.85)
        ax.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
        max_err = trt_ns['SEM'].max()
        for i, (_, row) in enumerate(trt_ns.iterrows()):
            y_top = row['Mean'] + row['SEM'] + max_err * 0.5
            ax.text(i, y_top, f"n={int(row['N'])}", ha='center', fontsize=9)
        ax.set_xticks(range(len(trt_ns)))
        ax.set_xticklabels(trt_ns['Treatment'], rotation=25, ha='right')
        ax.set_xlabel('Treatment Regimen', fontsize=11)
        ax.set_ylabel('Mean Normalization Score (± SEM)', fontsize=11)
        ax.set_title(
            'Metabolome Normalization by Treatment Regimen\n'
            '(higher = treatment returned metabolome closer to healthy)',
            fontsize=11, fontweight='bold'
        )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        plt.savefig(FIG_DIR / 'bar_normalization_score_by_treatment.png', dpi=300, bbox_inches='tight')
        plt.close()
        print('  Saved bar_normalization_score_by_treatment.png')

    print(f"\n  Mean NS: {ns_df['Mean_NS'].mean():.3f}  ±  {ns_df['Mean_NS'].std():.3f}")
    print(f"  Range:   [{ns_df['Mean_NS'].min():.3f},  {ns_df['Mean_NS'].max():.3f}]")

    return ns_df[['Patient_ID', 'Mean_NS', 'Time_to_Event', 'Event', 'Treatment', 'n_mets']]


def analyze_univariate_cox(df, met_cols, meth_cols=None, social_cols=None, n_top=20):
    """Univariate Cox PH for each feature; volcano plot + top-N table figure."""
    print("\nRunning Univariate Cox PH analysis for survival...")

    df_cox = df[df['Cohort'] == 'CANCER'].copy()
    df_cox = df_cox.dropna(subset=['Time to Event', 'Event'])
    df_cox = df_cox[(df_cox['Time to Event'] > 0)]
    n_pts = len(df_cox)
    print(f"  Using {n_pts} cancer patients with valid survival data")

    clinical_continuous = ['Age at the time of diagnosis', 'BMI', 'Pack Year History', 'Size 1 of tumor (cms)']

    feature_sets = []
    for col in (met_cols or []):
        feature_sets.append((col, 'Metabolite', True))
    for col in (meth_cols or []):
        feature_sets.append((col, 'Methylation', False))
    for col in (social_cols or []):
        feature_sets.append((col, 'Social', False))
    for col in clinical_continuous:
        if col in df_cox.columns:
            feature_sets.append((col, 'Clinical', False))

    records = []
    cph = CoxPHFitter()
    for feat, category, do_log in feature_sets:
        if feat not in df_cox.columns:
            continue
        vals = pd.to_numeric(df_cox[feat], errors='coerce')
        if do_log:
            vals = np.log1p(vals.clip(lower=0))
        tmp = pd.DataFrame({
            'duration': df_cox['Time to Event'].values,
            'event':    df_cox['Event'].values,
            feat:       vals.values,
        }).dropna()
        if len(tmp) < 10 or tmp[feat].nunique() < 3:
            continue
        # Skip near-zero variance
        if tmp[feat].std() < 1e-8:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cph.fit(tmp, duration_col='duration', event_col='event', show_progress=False)
            summary = cph.summary
            hr     = float(summary.loc[feat, 'exp(coef)'])
            ci_low = float(summary.loc[feat, 'exp(coef) lower 95%'])
            ci_hi  = float(summary.loc[feat, 'exp(coef) upper 95%'])
            pval   = float(summary.loc[feat, 'p'])
            log_hr = float(summary.loc[feat, 'coef'])
            records.append({
                'Feature':    feat,
                'Category':   category,
                'HR':         hr,
                'CI_lower':   ci_low,
                'CI_upper':   ci_hi,
                'log_HR':     log_hr,
                'p_value':    pval,
            })
        except Exception:
            continue

    if not records:
        print("  No Cox results obtained.")
        return pd.DataFrame()

    results = pd.DataFrame(records)
    results['FDR'] = fdrcorrection(results['p_value'])[1]
    results = results.sort_values('p_value').reset_index(drop=True)

    results['Display'] = results['Feature'].apply(met_label)

    n_sig = (results['p_value'] < 0.05).sum()
    n_fdr = (results['FDR'] < 0.05).sum()
    print(f"  Features tested: {len(results)}  |  p<0.05: {n_sig}  |  FDR<0.05: {n_fdr}")

    # ── Volcano plot ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    log10p = -np.log10(results['p_value'].clip(lower=1e-300))
    log_hr = results['log_HR']

    hr_thr  = np.log(VOLCANO_FC)   # 1.25 on the HR scale, matching the other volcanoes
    p_thr   = -np.log10(0.05)

    def point_color(row):
        sig = row['p_value'] < 0.05
        if sig and row['log_HR'] > hr_thr:
            return '#d62728'   # significantly high HR
        elif sig and row['log_HR'] < -hr_thr:
            return '#1f77b4'   # significantly low HR
        elif sig:
            return '#ff7f0e'   # significant but small effect
        return '#aaaaaa'       # not significant

    colors = results.apply(point_color, axis=1)
    ax.scatter(log_hr, log10p, c=colors, s=28, alpha=0.75, linewidths=0)

    ax.axhline(p_thr,   color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(hr_thr,  color='grey',  linestyle=':',  linewidth=0.8, alpha=0.6)
    ax.axvline(-hr_thr, color='grey',  linestyle=':',  linewidth=0.8, alpha=0.6)

    # Label top 15 by p-value
    top_label = results.head(15)
    for _, row in top_label.iterrows():
        ax.annotate(
            row['Display'],
            xy=(row['log_HR'], -np.log10(row['p_value'])),
            fontsize=6.5, ha='center', va='bottom',
            xytext=(0, 4), textcoords='offset points',
        )

    # Legend patches
    legend_els = [
        mpatches.Patch(facecolor='#d62728', label='High HR, p<0.05'),
        mpatches.Patch(facecolor='#1f77b4', label='Low HR, p<0.05'),
        mpatches.Patch(facecolor='#ff7f0e', label='Sig, small effect'),
        mpatches.Patch(facecolor='#aaaaaa', label='Not significant'),
    ]
    ax.legend(handles=legend_els, fontsize=9, loc='upper left')
    ax.set_xlabel('log(Hazard Ratio)', fontsize=12)
    ax.set_ylabel('-log₁₀(p-value)', fontsize=12)
    ax.set_title('Univariate Cox PH: Volcano Plot (Survival)', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'volcano_univariate_cox_survival.png', dpi=300, bbox_inches='tight')
    plt.close()
    print('  Saved volcano_univariate_cox_survival.png')

    # ── Table figure ─────────────────────────────────────────────────────────
    top_n = results.head(n_top).copy()
    top_n['Rank'] = range(1, len(top_n) + 1)
    top_n['HR (95% CI)'] = top_n.apply(
        lambda r: f"{r['HR']:.2f} ({r['CI_lower']:.2f}–{r['CI_upper']:.2f})", axis=1
    )
    top_n['p-value']  = top_n['p_value'].apply(lambda v: f"{v:.3e}")
    top_n['FDR']      = top_n['FDR'].apply(lambda v: f"{v:.3e}")
    top_n['Feature Name'] = top_n['Display']

    table_cols = ['Rank', 'Feature Name', 'Category', 'HR (95% CI)', 'p-value', 'FDR']
    table_data = top_n[table_cols].values.tolist()
    col_labels = table_cols

    row_h   = 0.38
    fig_h   = row_h * (len(table_data) + 1.5)
    fig, ax = plt.subplots(figsize=(14, max(fig_h, 4)))
    ax.axis('off')
    tbl = ax.table(
        cellText=table_data,
        colLabels=col_labels,
        cellLoc='left',
        loc='center',
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.auto_set_column_width(col=list(range(len(col_labels))))

    # Style header
    for j in range(len(col_labels)):
        cell = tbl[0, j]
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold')

    # Alternating row shading
    for i in range(1, len(table_data) + 1):
        for j in range(len(col_labels)):
            tbl[i, j].set_facecolor('#f2f2f2' if i % 2 == 0 else 'white')

    ax.set_title(
        f'Top {n_top} Features — Univariate Cox PH (Survival)\n'
        f'Cancer patients n={n_pts}',
        fontsize=11, fontweight='bold', pad=8
    )
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'table_univariate_cox_top20.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  Saved table_univariate_cox_top20.png')

    # Console summary
    print(f"\n{'='*80}")
    print(f"TOP {n_top} FEATURES — UNIVARIATE COX PH")
    print(f"{'='*80}")
    print(top_n[['Rank', 'Feature Name', 'Category', 'HR (95% CI)', 'p-value', 'FDR']].to_string(index=False))
    print(f"{'='*80}\n")

    return results


def plot_metabolite_distributions(df, met_results, n_top=10, prefix=''):
    """One figure per direction (increased / decreased): n_top rows, each row is a
    metabolite showing overlaid KDE distributions for cancer cases vs controls.
    prefix: optional filename prefix to avoid overwriting other analyses' plots."""
    print(f"Plotting top {n_top} metabolite distributions (cancer vs control)...")

    df_cancer  = df[df['Cohort'] == 'CANCER']
    df_control = df[df['Cohort'] == 'No']
    n_cancer   = len(df_cancer)
    n_control  = len(df_control)

    clean_name = met_label

    for direction, label, color_cancer, color_ctrl in [
        ('Increased', 'increased_in_cancer', '#d62728', '#1f77b4'),
        ('Decreased', 'decreased_in_cancer', '#d62728', '#1f77b4'),
    ]:
        top = (met_results[
                   (met_results['FDR'] < 0.05) &
                   (met_results['Direction'] == direction)
               ]
               .sort_values('FDR')
               .head(n_top))

        if top.empty:
            print(f"  No FDR-significant {direction.lower()} metabolites to plot.")
            continue

        n = len(top)
        fig, axes = plt.subplots(n, 1, figsize=(9, 2.2 * n), sharex=False)
        if n == 1:
            axes = [axes]

        for i, (ax, (_, row)) in enumerate(zip(axes, top.iterrows())):
            feat = row['Feature']
            fdr  = row['FDR']
            name = clean_name(feat)

            v_cancer  = np.log1p(pd.to_numeric(df_cancer[feat],  errors='coerce').dropna())
            v_control = np.log1p(pd.to_numeric(df_control[feat], errors='coerce').dropna())

            if len(v_cancer) > 2:
                sns.kdeplot(v_cancer,  ax=ax, fill=True, alpha=0.45,
                            color=color_cancer, linewidth=1.5,
                            label=f'Cancer  (n={len(v_cancer)})')
            if len(v_control) > 2:
                sns.kdeplot(v_control, ax=ax, fill=True, alpha=0.45,
                            color=color_ctrl,  linewidth=1.5,
                            label=f'Control (n={len(v_control)})')

            ax.set_title(f'{name}   [FDR = {fdr:.2e}]', fontsize=16, fontweight='bold', loc='left', pad=3)
            ax.set_ylabel('Density', fontsize=14)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            # Legend only on first row; remove from all others
            if i == 0:
                ax.legend(fontsize=12, loc='upper right', framealpha=0.85)
            else:
                legend = ax.get_legend()
                if legend:
                    legend.remove()

            # x-axis label only on last row
            if i < n - 1:
                ax.set_xlabel('')
                ax.tick_params(axis='x', labelbottom=False)
            else:
                ax.set_xlabel('Value (log1p transformed)', fontsize=12)

        direction_title = direction.lower()
        fig.suptitle(
            f'Top {n} metabolites {direction_title} in cancer vs control  '
            f'(Cancer n={n_cancer}, Control n={n_control})',
            fontsize=14, fontweight='bold', y=1.015
        )
        plt.tight_layout()
        fname = f'kde_{prefix}top{n_top}_metabolites_{label}.png'
        plt.savefig(FIG_DIR / fname, dpi=300, bbox_inches='tight')
        plt.close()
        print(f'  Saved {fname}')

# ─────────────────────────────────────────────────────────────────────────────
# Table 1 — baseline characteristics of the cancer cohort
# ─────────────────────────────────────────────────────────────────────────────

# Tumor sequencing / IHC marker columns, in workbook order. CDK12 and RET appear
# twice in the sheet; dict.fromkeys keeps the first occurrence of each.
MUTATION_COLS = list(dict.fromkeys([
    'EGFR', 'ROS1', 'KRAS', 'HRAS', 'BRAF', 'TP53', 'PIK3CA', 'ALK', 'CDKN2A', 'CDK12',
    'MET', 'GNAS R20', 'EML4-ALK', 'RET', 'STK11', 'FGFR', 'MDM2', 'GATA3', 'KEAPA',
    'SMARKA4', 'RBM10', 'CHD2', 'MLH3', 'ASXL1', 'Notch1', 'NFKB1A', 'NFKB2', 'CDKN1B',
    'KEAP1', 'ATRX', 'SMARCA4+', 'KMT2D', 'PTPRT', 'ERBB3', 'TTF1-/Napsin', 'APC',
    'CCNE1', 'MTOR', 'TOP2A', 'AKT2', 'NY-ESO-1', 'NRG1',
]))
PDL1_COLS = {'PDL-1 <1%': '< 1%', 'PDL-1 to 25%': '1–25%',
             'PDL-25 to 50%': '25–50%', 'PDL-1 >50%': '> 50%'}
COMORBIDITY_COLS = ['HTN', 'DM', 'COPD', 'Emphysema', 'Asthma', 'CAD', 'CHF', 'Thyroid']
CHEMO_AGENTS = ['Carboplatin', 'Cisplatin', 'Pemetrexed', 'Paclitaxel', 'Etoposide', 'Docetaxel',
                'Abraxane', 'Gemcitabine', 'Lurbinectedin', 'Topotecan']
IO_AGENTS = ['Pembrolizumab', 'Durvalumab', 'Nivolumab', 'Atezolizumab', 'Bevacizumab',
             'Ramucirumab', 'Ipilimumab']
TARGETED_AGENTS = ['Osimertinib', 'Erlotinib', 'Crizotinib', 'Alectinib', 'Afatinib',
                   'Lorlatinib', 'Sotorasib', 'Selpercatinib', 'Capmatinib']
SOCIAL_COLS_ALL = list(SOCIAL_LABELS.keys())

TUMOR_SIZE_IMPLAUSIBLE_CM = 15.0  # values above this are treated as data-entry errors
MUT_MIN_N = 3                     # mutations seen in fewer patients are collapsed into one row


def _get_col(df, name):
    """Return a single Series for `name`.

    The workbook has duplicated column names (CDK12, RET), so df[name] can come back
    as a DataFrame; coalesce those into one Series taking the first non-blank value.
    """
    if name not in df.columns:
        return None
    col = df[name]
    if isinstance(col, pd.DataFrame):
        out = col.iloc[:, 0].astype(object)
        blank = out.isna() | out.astype(str).str.strip().isin(['', 'nan', 'None'])
        for j in range(1, col.shape[1]):
            out = out.where(~blank, col.iloc[:, j])
            blank = out.isna() | out.astype(str).str.strip().isin(['', 'nan', 'None'])
        return out
    return col


def _clean_cat(series):
    """Blank-like values become NaN so they are counted as missing, not as a level."""
    s = series.astype(str).str.strip()
    return s.replace({'': np.nan, 'nan': np.nan, 'NaN': np.nan, 'None': np.nan,
                      'NA': np.nan, 'N/A': np.nan})


def _clean_flag(series):
    """Yes/blank flags → Yes/No. In this workbook a blank means 'not reported positive',
    which is how the reference table's EGFR No 71 / Yes 9 split arises."""
    s = series.astype(str).str.strip().str.upper()
    return pd.Series(np.where(s.str.startswith('YES'), 'Yes', 'No'), index=series.index)


def _fmt_n_pct(count, denom):
    if not denom:
        return '—'
    return f"{int(count)} ({100.0 * count / denom:.1f}%)"


def _fmt_mean_sd(series, dec=1):
    v = pd.to_numeric(series, errors='coerce').dropna()
    if v.empty:
        return '—'
    return f"{v.mean():.{dec}f} ± {v.std():.{dec}f}"


def _fmt_median_iqr(series, dec=2):
    v = pd.to_numeric(series, errors='coerce').dropna()
    if v.empty:
        return '—'
    return f"{v.median():.{dec}f} [{v.quantile(0.25):.{dec}f}–{v.quantile(0.75):.{dec}f}]"


def _rows_categorical(df, col, label, order=None, transform=None, missing_label='Not reported',
                      show_denom=True):
    """Header row '<label>, n (%)' followed by one indented row per level.

    Percentages use the non-missing denominator; the denominator is shown in the
    header whenever it is smaller than the cohort.
    """
    s = _get_col(df, col)
    if s is None:
        print(f"  [Table 1] column not found, row skipped: {col!r}")
        return []

    s = _clean_cat(s)
    if transform is not None:
        s = s.map(lambda v: v if pd.isna(v) else transform(v))

    denom = int(s.notna().sum())
    n_missing = int(s.isna().sum())
    counts = s.value_counts()

    header = (f"{label}, n (%)" if denom == len(df) or not show_denom
              else f"{label}, n (%)  [n={denom}]")
    rows = [(header, '')]

    levels = list(order) if order is not None else counts.index.tolist()
    if order is not None:
        levels += [lev for lev in counts.index if lev not in order]
    for lev in levels:
        if lev in counts.index:
            rows.append((f"    {lev}", _fmt_n_pct(counts[lev], denom)))
    if n_missing:
        rows.append((f"    {missing_label}", _fmt_n_pct(n_missing, len(df))))
    return rows


def _row_flag(df, col, label=None):
    """Single indented 'Yes' count row for a yes/blank flag column."""
    s = _get_col(df, col)
    if s is None:
        return []
    n_yes = int((_clean_flag(s) == 'Yes').sum())
    return [(f"    {label or col}", _fmt_n_pct(n_yes, len(df)))]


def _row_continuous(df, col, label, style='mean_sd', dec=1):
    s = _get_col(df, col)
    if s is None:
        print(f"  [Table 1] column not found, row skipped: {col!r}")
        return []
    v = pd.to_numeric(s, errors='coerce')
    n_missing = int(v.isna().sum())
    suffix = f"  [n={int(v.notna().sum())}]" if n_missing else ''
    val = _fmt_mean_sd(v, dec) if style == 'mean_sd' else _fmt_median_iqr(v, dec)
    return [(f"{label}{suffix}", val)]


def build_table_one(df, met_cols=None, filename_stem='table_baseline_cancer'):
    """Publication-style Table 1 for the cancer analysis cohort.

    Summarises the cleaned cancer cohort (NSCLC; SCLC already excluded by
    DataLoader.clean_data) across demographics, survival follow-up, tumor
    characteristics, smoking, comorbidities, treatment/response, molecular and
    epigenetic markers, social indicators, and assay availability.

    Writes <filename_stem>.csv / .md / .png into FIG_DIR and returns the tidy frame.
    """
    d = df[df['Cohort'] == 'CANCER'].copy()
    n = len(d)
    print(f"\nBuilding Table 1 for the cancer analysis cohort (n={n})...")

    blocks = []  # list of (block_name, [(characteristic, value), ...])

    # ── Demographics ─────────────────────────────────────────────────────────
    rows = []
    rows += _row_continuous(d, 'Age at the time of diagnosis', 'Age at diagnosis, years — mean ± SD')
    rows += _row_continuous(d, 'Age at the time of diagnosis', '    median [IQR]',
                            style='median_iqr', dec=1)
    rows += _rows_categorical(d, 'Sex', 'Sex', order=['Female', 'Male'])
    rows += _rows_categorical(d, 'Race', 'Race', order=['Black', 'White', 'Asian'])
    rows += _row_continuous(d, 'BMI', 'BMI, kg/m² — mean ± SD')
    blocks.append(('Demographics', rows))

    # ── Vital status and follow-up ───────────────────────────────────────────
    event = pd.to_numeric(_get_col(d, 'Event'), errors='coerce')
    tte = pd.to_numeric(_get_col(d, 'Time to Event'), errors='coerce')
    n_dead, n_alive = int((event == 1).sum()), int((event == 0).sum())
    rows = [
        ('Vital status, n (%)', ''),
        ('    Alive (censored)', _fmt_n_pct(n_alive, n)),
        ('    Dead (event observed)', _fmt_n_pct(n_dead, n)),
        ('Time to event, days — median [IQR]', _fmt_median_iqr(tte, dec=0)),
        ('    range', f"{tte.min():.0f}–{tte.max():.0f}" if tte.notna().any() else '—'),
        ('Follow-up in survivors, days — median [IQR]', _fmt_median_iqr(tte[event == 0], dec=0)),
        ('Time to death in decedents, days — median [IQR]', _fmt_median_iqr(tte[event == 1], dec=0)),
    ]
    diag = pd.to_datetime(_get_col(d, 'Date of Initial Diagnosis'), errors='coerce')
    draw = pd.to_datetime(_get_col(d, 'Date of baseline blood draw (visit)'), errors='coerce')
    rows.append(('Baseline draw relative to diagnosis, days — median [IQR]',
                 _fmt_median_iqr((draw - diag).dt.days, dec=0)))
    if diag.notna().any():
        rows.append(('Year of diagnosis, range', f"{diag.min():%Y}–{diag.max():%Y}"))
    blocks.append(('Vital status and follow-up', rows))

    # ── Tumor characteristics ────────────────────────────────────────────────
    rows = _rows_categorical(d, 'Histological Classification', 'Histological classification')
    rows += _rows_categorical(d, 'Stage_orig', 'Stage (AJCC clinical)',
                              order=['Stage I', 'Stage II', 'Stage III', 'Stage IV'])
    up = lambda v: str(v).upper().replace('NO', 'N0')  # noqa: E731 — 'NO' is a typo for N0 in 1 row
    rows += _rows_categorical(d, 'T', 'Tumor (T)', order=['T1', 'T2', 'T3', 'T4', 'TX'],
                              transform=lambda v: str(v).upper())
    rows += _rows_categorical(d, 'N', 'Nodes (N)', order=['N0', 'N1', 'N2', 'N3', 'NX'],
                              transform=up)
    rows += _rows_categorical(d, 'M', 'Metastasis (M)', order=['M0', 'M1', 'MX'],
                              transform=lambda v: str(v).upper())

    size = pd.to_numeric(_get_col(d, 'Size 1 of tumor (cms)'), errors='coerce')
    size_ok = size[size <= TUMOR_SIZE_IMPLAUSIBLE_CM]
    n_implausible = int((size > TUMOR_SIZE_IMPLAUSIBLE_CM).sum())
    rows.append(('Size of tumor, cm — mean ± SD', _fmt_mean_sd(size, dec=2)))
    if n_implausible:
        rows.append((f'    excluding {n_implausible} value(s) > {TUMOR_SIZE_IMPLAUSIBLE_CM:.0f} cm',
                     _fmt_mean_sd(size_ok, dec=2)))
    rows.append(('    median [IQR]', _fmt_median_iqr(size, dec=2)))
    rows += _rows_categorical(d, 'Primary Tumor Location', 'Primary tumor location',
                              transform=lambda v: str(v).upper())
    blocks.append(('Tumor characteristics', rows))

    # ── Smoking exposure ─────────────────────────────────────────────────────
    rows = _rows_categorical(d, 'Smoking Status', 'Smoking status',
                             order=['Never', 'Former', 'Current'])
    packyr = pd.to_numeric(_get_col(d, 'Pack Year History'), errors='coerce')
    smoke = _clean_cat(_get_col(d, 'Smoking Status'))
    rows.append(('Pack-year history — mean ± SD', _fmt_mean_sd(packyr)))
    ever = packyr[smoke.isin(['Former', 'Current'])]
    rows.append((f'    ever-smokers only [n={int(ever.notna().sum())}]', _fmt_mean_sd(ever)))
    rows.append(('    median [IQR]', _fmt_median_iqr(packyr, dec=1)))
    blocks.append(('Smoking exposure', rows))

    # ── Comorbidities ────────────────────────────────────────────────────────
    rows = [('Comorbidities, n (%)', '')]
    for c in COMORBIDITY_COLS:
        rows += _row_flag(d, c)
    blocks.append(('Comorbidities', rows))

    # ── Treatment and response ───────────────────────────────────────────────
    tx_cols = [('SBRT', 'SBRT'),
               ('Chemo receive  yes or none', 'Chemotherapy'),
               ('Immunotherapy Receive  yes or none', 'Immunotherapy'),
               ('Targeted therapy yes or none', 'Targeted therapy')]
    tx_given = {lab: (_clean_flag(_get_col(d, col)) == 'Yes') for col, lab in tx_cols}
    n_modalities = pd.concat(list(tx_given.values()), axis=1).sum(axis=1)

    rows = [('Treatment received, n (%)', '')]
    for _, lab in tx_cols:
        rows.append((f"    {lab}", _fmt_n_pct(int(tx_given[lab].sum()), n)))
    rows.append(('    No treatment recorded', _fmt_n_pct(int((n_modalities == 0).sum()), n)))

    # Every observed modality combination, so the mutually-exclusive groups sum to n
    short = {'SBRT': 'SBRT', 'Chemotherapy': 'Chemo', 'Immunotherapy': 'IO',
             'Targeted therapy': 'Targeted'}
    combo = pd.Series(
        [' + '.join(short[lab] for lab in short if tx_given[lab].iloc[i]) or 'None'
         for i in range(n)], index=d.index)
    rows.append(('Treatment combination, n (%)', ''))
    for lab, cnt in combo.value_counts().items():
        rows.append((f"    {lab}", _fmt_n_pct(cnt, n)))
    rows.append(('Number of modalities, n (%)', ''))
    for k in range(int(n_modalities.max()) + 1):
        cnt = int((n_modalities == k).sum())
        if cnt:
            rows.append((f"    {k}", _fmt_n_pct(cnt, n)))

    # Individual agents, denominators restricted to the patients given that modality
    for lab, agents in [('Chemotherapy', CHEMO_AGENTS),
                        ('Immunotherapy', IO_AGENTS),
                        ('Targeted therapy', TARGETED_AGENTS)]:
        given = tx_given[lab]
        counts = {a: int((_clean_flag(_get_col(d, a)) == 'Yes').sum())
                  for a in agents if _get_col(d, a) is not None}
        counts = {a: v for a, v in counts.items() if v}
        if not counts:
            continue
        rows.append((f'{lab} agents, n (%) of {int(given.sum())} treated', ''))
        for a, v in sorted(counts.items(), key=lambda kv: -kv[1]):
            rows.append((f"    {a}", _fmt_n_pct(v, int(given.sum()))))

    rows += _rows_categorical(d, 'ECOG status at the end', 'ECOG performance status (last recorded)',
                              transform=lambda v: str(v).split('.')[0],
                              missing_label='Not recorded')

    # Response is only meaningful among the IO-treated, so use that subgroup as the denominator
    io_given = tx_given['Immunotherapy']
    n_io = int(io_given.sum())
    rows += _rows_categorical(d[io_given], 'IO Responder Yes or No',
                              f'Immunotherapy response (IO-treated only, n={n_io})',
                              order=['Yes', 'No'],
                              transform=lambda v: str(v).strip().title(),
                              missing_label='Not recorded', show_denom=False)
    rows.append(('    Not applicable (no immunotherapy)', _fmt_n_pct(n - n_io, n)))
    blocks.append(('Treatment and response', rows))

    # ── Molecular alterations ────────────────────────────────────────────────
    mut_counts = {}
    for c in MUTATION_COLS:
        s = _get_col(d, c)
        if s is not None:
            mut_counts[c] = int((_clean_flag(s) == 'Yes').sum())
    frequent = {k: v for k, v in mut_counts.items() if v >= MUT_MIN_N}
    rare = {k: v for k, v in mut_counts.items() if 0 < v < MUT_MIN_N}

    # A patient with no positive call anywhere was most likely never sequenced, so report
    # that group explicitly instead of implying all n patients were tested.
    any_call = pd.concat([_clean_flag(_get_col(d, c)) == 'Yes' for c in mut_counts],
                         axis=1).any(axis=1)
    n_called = int(any_call.sum())
    rows = [('Tumor sequencing / IHC, n (%)', ''),
            ('    ≥1 positive marker call', _fmt_n_pct(n_called, n)),
            ('    No positive call (testing status unknown)', _fmt_n_pct(n - n_called, n)),
            (f'Markers detected in ≥{MUT_MIN_N} patients, n (%) of {n}', '')]
    for c, cnt in sorted(frequent.items(), key=lambda kv: -kv[1]):
        rows.append((f"    {c}", _fmt_n_pct(cnt, n)))
    if rare:
        rows.append((f"    Other markers, each in <{MUT_MIN_N} patients",
                     f"{len(rare)} markers (listed below)"))

    # The four PD-L1 bands are one ordinal variable, not four independent flags
    pdl1_flags = {lab: (_clean_flag(_get_col(d, c)) == 'Yes')
                  for c, lab in PDL1_COLS.items() if _get_col(d, c) is not None}
    n_bands = pd.concat(list(pdl1_flags.values()), axis=1).sum(axis=1)
    n_tested = int((n_bands >= 1).sum())
    n_ambiguous = int((n_bands > 1).sum())
    rows.append((f'PD-L1 tumor proportion score, n (%) of {n_tested} tested', ''))
    for lab, mask in pdl1_flags.items():
        rows.append((f"    {lab}", _fmt_n_pct(int(mask.sum()), n_tested)))
    if n_ambiguous:
        rows.append(('    (more than one band flagged)', _fmt_n_pct(n_ambiguous, n_tested)))
    rows.append(('    No PD-L1 band recorded', _fmt_n_pct(n - n_tested, n)))
    blocks.append(('Molecular alterations', rows))

    # ── Epigenetic markers ───────────────────────────────────────────────────
    rows = [('DNA methylation (β), median [IQR]', '')]
    for c in METH_COLS_ALL:
        s = _get_col(d, c)
        if s is not None:
            rows.append((f"    {METH_LABELS.get(c, c)}", _fmt_median_iqr(s, dec=3)))
    blocks.append(('Epigenetic markers', rows))

    # ── Social indicators ────────────────────────────────────────────────────
    rows = [('Neighborhood-level indices, median [IQR]', '')]
    for c in SOCIAL_COLS_ALL:
        s = _get_col(d, c)
        if s is None:
            continue
        v = pd.to_numeric(s, errors='coerce')
        dec = 0 if c == 'INC' else 2
        lab = f"    {c} — {SOCIAL_LABELS[c]}"
        if v.isna().any():
            lab += f"  [n={int(v.notna().sum())}]"
        rows.append((lab, _fmt_median_iqr(v, dec=dec)))
    blocks.append(('Social indicators', rows))

    # ── Assay availability ───────────────────────────────────────────────────
    exclude_keywords = ['Date', 'Follow Up', 'Draw', 'Immunotherapy', 'Chemotherapy', 'SBRT']
    all_baseline = [c for c in df.columns
                    if str(c).startswith('Baseline_')
                    and not any(k in str(c) for k in exclude_keywords)]
    # The pipeline only keeps KEGG-annotated features (the 'Baseline_C…' prefix)
    kegg_baseline = [c for c in all_baseline if str(c).startswith('Baseline_C')]
    if met_cols is None:
        sub = d[kegg_baseline].apply(pd.to_numeric, errors='coerce')
        met_cols = [c for c in kegg_baseline if (sub[c] == 0).mean() <= pct_zeros]
    fu_cols = [c for c in df.columns if str(c).startswith('Follow-Up_')]
    n_fu = int(d[fu_cols].apply(pd.to_numeric, errors='coerce').notna().any(axis=1).sum()) if fu_cols else 0
    immune_present = [c for c in IMMUNE_COLS if _get_col(d, c) is not None]
    n_immune = int(pd.concat([pd.to_numeric(_get_col(d, c), errors='coerce')
                              for c in immune_present], axis=1).notna().any(axis=1).sum()) if immune_present else 0

    rows = [
        ('Baseline metabolite features quantified', f"{len(all_baseline)}"),
        ('    KEGG-annotated (used by the pipeline)', f"{len(kegg_baseline)}"),
        (f'    retained after ≤{pct_zeros:.0%} zero filter', f"{len(met_cols)}"),
        ('Patients with a paired follow-up draw', _fmt_n_pct(n_fu, n)),
        ('Patients with DNAm marker panel', _fmt_n_pct(int(pd.concat(
            [pd.to_numeric(_get_col(d, c), errors='coerce') for c in METH_COLS_ALL],
            axis=1).notna().any(axis=1).sum()), n)),
        ('Patients with tumor-infiltrate panel', _fmt_n_pct(n_immune, n)),
        ('Social indices available', f"{len([c for c in SOCIAL_COLS_ALL if _get_col(d, c) is not None])} of 30"),
    ]
    blocks.append(('Assay availability', rows))

    # ── Tidy frame ───────────────────────────────────────────────────────────
    tidy = pd.DataFrame(
        [(blk, char, val) for blk, rws in blocks for char, val in rws],
        columns=['Block', 'Characteristic', 'Value'],
    )

    notes = _table_one_notes(d, n_implausible, packyr, smoke, n_immune)
    if rare:
        notes.append(f"Markers positive in <{MUT_MIN_N} patients: "
                     + ', '.join(f"{k} ({v})" for k, v in sorted(rare.items())) + ".")

    # ── Outputs ──────────────────────────────────────────────────────────────
    title = f"Table 1. Baseline characteristics of the cancer cohort (N={n})"
    csv_path = FIG_DIR / f'{filename_stem}.csv'
    tidy.to_csv(csv_path, index=False)

    md_path = FIG_DIR / f'{filename_stem}.md'
    with open(md_path, 'w') as fh:
        fh.write(f"# {title}\n\n")
        fh.write("NSCLC cancer patients after the exclusions applied in `DataLoader.clean_data`.\n\n")
        fh.write("| Characteristic | Value |\n|---|---|\n")
        for blk, rws in blocks:
            fh.write(f"| **{blk}** | |\n")
            for char, val in rws:
                fh.write(f"| {char.replace('    ', '&nbsp;&nbsp;&nbsp;&nbsp;')} | {val} |\n")
        fh.write("\n## Notes\n\n")
        for note in notes:
            fh.write(f"- {note}\n")

    _render_table_one_png(blocks, title, notes, FIG_DIR / f'{filename_stem}.png')

    print(f"  Saved {csv_path.name}, {md_path.name}, {filename_stem}.png  ({len(tidy)} rows)")
    print("\n" + "=" * 80)
    print(title.upper())
    print("=" * 80)
    print(tidy[['Characteristic', 'Value']].to_markdown(index=False))
    print("-" * 80)
    for note in notes:
        print(f"  * {note}")
    print("=" * 80 + "\n")

    return tidy


def _table_one_notes(d, n_implausible, packyr, smoke, n_immune):
    """Footnotes, including the data-quality problems found in the workbook."""
    n = len(d)
    notes = [
        "Cohort: cancer patients after the exclusions in `DataLoader.clean_data` "
        "(SCLC excluded; baseline metabolite values only).",
        "Percentages use the non-missing denominator; a reduced denominator is shown as [n=…] "
        "in the row label.",
        "Molecular and treatment flags are recorded as 'yes' or blank in the workbook, so a blank "
        "is counted as 'not reported positive' rather than as missing.",
    ]

    # Data-quality issues found while validating against the reference table
    dq = []
    if n_implausible:
        size = pd.to_numeric(_get_col(d, 'Size 1 of tumor (cms)'), errors='coerce')
        biggest = size.max()
        nxt = size[size <= TUMOR_SIZE_IMPLAUSIBLE_CM].max()
        dq.append(f"tumor size of {biggest:.1f} cm in {n_implausible} patient(s) "
                  f"(next largest {nxt:.1f} cm) — likely a decimal-point data-entry error")

    m = _clean_cat(_get_col(d, 'M')).map(lambda v: v if pd.isna(v) else str(v).upper())
    stg = _clean_cat(_get_col(d, 'Stage_orig'))
    bad_m1 = int(((m == 'M1') & (stg != 'Stage IV')).sum())
    bad_m0 = int(((m == 'M0') & (stg == 'Stage IV')).sum())
    if bad_m1 or bad_m0:
        dq.append(f"{bad_m1} patient(s) coded M1 but not staged IV, and {bad_m0} coded M0 but "
                  f"staged IV — AJCC M1 implies Stage IV, so these rows disagree with themselves")

    zero_py = int(((packyr == 0) & smoke.isin(['Former', 'Current'])).sum())
    if zero_py:
        dq.append(f"{zero_py} ever-smoker(s) with a pack-year history of 0")

    n_typo = int((_get_col(d, 'N').astype(str).str.strip().str.upper() == 'NO').sum())
    if n_typo:
        dq.append(f"{n_typo} nodal stage recorded as 'NO' (letter O) — counted as N0")

    io_given = _clean_flag(_get_col(d, 'Immunotherapy Receive  yes or none')) == 'Yes'
    io_resp = _clean_cat(_get_col(d, 'IO Responder Yes or No')).notna()
    resp_no_io, io_no_resp = int((io_resp & ~io_given).sum()), int((io_given & ~io_resp).sum())
    if resp_no_io or io_no_resp:
        dq.append(f"{resp_no_io} patient(s) have an immunotherapy response recorded without being "
                  f"flagged as IO-treated (excluded from the response row above, which is "
                  f"restricted to the IO-treated group), and {io_no_resp} IO-treated patient(s) "
                  f"have no response recorded")

    if dq:
        notes.append("Data-quality issues present in the workbook and reported here as-is, "
                     "without correction: " + "; ".join(dq) + ".")

    notes.append(
        "The workbook records no explicit 'tested' indicator for sequencing or PD-L1, so per-marker "
        "percentages are over the whole cohort and cannot separate 'tested negative' from "
        "'never tested'. Patients with no positive call at all are reported as their own row."
    )
    notes.append(
        f"Tumor-infiltrate percentages (macrophages, B/NK/T cells, CD4, CD8) are measured in only "
        f"{n_immune} of {n} patients and are therefore omitted from the table body."
    )
    notes.append("TTF1/Napsin and NY-ESO-1 are immunohistochemistry markers rather than sequenced "
                 "mutations, but share the same yes/blank encoding.")
    return notes


def _render_table_one_png(blocks, title, notes, out_path, n_panels=3):
    """Render Table 1 as a multi-panel figure, matching the styling used by
    analyze_univariate_cox's table figure."""
    # Keep blocks intact but balance panel heights: contiguous partition of the blocks
    # that minimises the tallest panel (binary search on the height cap).
    sizes = [len(rws) + 1 for _, rws in blocks]  # +1 for the block header row

    def _panels_needed(cap):
        cnt, cur = 1, 0
        for s in sizes:
            if cur and cur + s > cap:
                cnt, cur = cnt + 1, s
            else:
                cur += s
        return cnt

    lo, hi = max(sizes), sum(sizes)
    while lo < hi:
        mid = (lo + hi) // 2
        if _panels_needed(mid) <= n_panels:
            hi = mid
        else:
            lo = mid + 1

    panels, cur = [], []
    for (blk, rws), s in zip(blocks, sizes):
        if cur and len(cur) + s > lo:
            panels.append(cur)
            cur = []
        cur.append(('__BLOCK__', blk, ''))
        cur.extend(('row', char, val) for char, val in rws)
    if cur:
        panels.append(cur)

    max_rows = max(len(p) for p in panels)
    fig_h = 0.235 * max_rows + 1.6
    fig, axes = plt.subplots(1, len(panels), figsize=(9.5 * len(panels), fig_h))

    for ax, panel in zip(np.atleast_1d(axes), panels):
        ax.axis('off')
        cell_text = [[c.replace('    ', '      '), v] for _, c, v in panel]
        cell_text += [['', '']] * (max_rows - len(panel))  # pad so panels align

        tbl = ax.table(cellText=cell_text, colLabels=['Characteristic', 'Value'],
                       cellLoc='left', loc='upper center', colWidths=[0.72, 0.28])
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)
        tbl.scale(1, 1.18)

        for j in range(2):
            hdr = tbl[0, j]
            hdr.set_facecolor('#2c3e50')
            hdr.set_text_props(color='white', fontweight='bold')

        for i, (kind, _, _) in enumerate(panel, start=1):
            for j in range(2):
                cell = tbl[i, j]
                cell.set_edgecolor('#d9d9d9')
                if kind == '__BLOCK__':
                    cell.set_facecolor('#dfe6ec')
                    cell.set_text_props(fontweight='bold')
                else:
                    cell.set_facecolor('#f7f7f7' if i % 2 == 0 else 'white')
        for i in range(len(panel) + 1, max_rows + 1):
            for j in range(2):
                tbl[i, j].set_edgecolor('white')
                tbl[i, j].set_facecolor('white')

    fig.suptitle(title, fontsize=13, fontweight='bold', y=0.995)
    footnote = '\n'.join(f"• {note}" for note in notes)
    fig.text(0.01, 0.005, footnote, fontsize=7.5, va='bottom', ha='left', wrap=True)
    plt.tight_layout(rect=[0, 0.055, 1, 0.985])
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Probabilistic Quotient Normalization (PQN)
# ─────────────────────────────────────────────────────────────────────────────

def pqn_normalize(X: pd.DataFrame, reference_mask: np.ndarray) -> tuple[pd.DataFrame, pd.Series]:
    """Probabilistic quotient normalization.

    Each sample is divided by the median ratio of its detected metabolites to a
    reference profile. Zeros are non-detections, not measurements, so they are
    excluded from the quotient rather than counted as a ratio of zero.
    """
    Xp = X.mask(X <= 0)
    reference = Xp.loc[reference_mask].median(axis=0)
    quotients = Xp.div(reference, axis=1)
    dilution = quotients.median(axis=1).replace(0, np.nan).fillna(1.0)
    return X.div(dilution, axis=0), dilution


def apply_pqn_normalization(df):
    """Apply Probabilistic Quotient Normalization (PQN) to all baseline and follow-up metabolites.

    Reference profile r_j: median intensity across healthy controls (Cohort == 'No') for detected values (> 0).
    Sample dilution factor d_i: median quotient across detected metabolites for that sample.
    Each sample's metabolite profile is divided by its dilution factor.
    Returns: (df_pqn, dilution_baseline, dilution_followup)
    """
    df = df.copy()
    is_ctrl = (df['Cohort'] == 'No').to_numpy()

    # 1. Baseline metabolites
    exclude_kw = ['Date', 'Follow Up', 'Draw', 'Immunotherapy', 'Chemotherapy', 'SBRT']
    b_cols = [c for c in df.columns if str(c).startswith('Baseline_C') and not any(k in str(c) for k in exclude_kw)]

    if b_cols:
        X_base = df[b_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
        Xp_base = X_base.mask(X_base <= 0)
        ref_base = Xp_base.loc[is_ctrl].median(axis=0)

        quotients_base = Xp_base.div(ref_base, axis=1)
        dilution_base = quotients_base.median(axis=1).replace(0, np.nan).fillna(1.0)

        df[b_cols] = X_base.div(dilution_base, axis=0)
        print(f"\n[PQN Normalization] Applied to {len(b_cols)} baseline metabolites across {len(df)} samples.")
        print(f"  Fitted dilution factor d_i: Controls median = {dilution_base[is_ctrl].median():.3f}, "
              f"Cancer median = {dilution_base[~is_ctrl].median():.3f} (ratio = {dilution_base[is_ctrl].median() / dilution_base[~is_ctrl].median():.2f}x)")
    else:
        dilution_base = pd.Series(1.0, index=df.index)

    # 2. Follow-up metabolites (paired to baseline control reference)
    pairs = pair_baseline_followup(b_cols, df.columns)
    if pairs and b_cols:
        b_paired = [p[0] for p in pairs]
        fu_paired = [p[1] for p in pairs]

        X_fu = df[fu_paired].apply(pd.to_numeric, errors='coerce').fillna(0.0)
        Xp_fu = X_fu.mask(X_fu <= 0)
        ref_fu = ref_base[b_paired].values

        quotients_fu = Xp_fu.div(ref_fu, axis=1)
        dilution_fu = quotients_fu.median(axis=1).replace(0, np.nan).fillna(1.0)

        df[fu_paired] = X_fu.div(dilution_fu, axis=0)
        has_fu = (X_fu > 0).sum(axis=1) > 5
        print(f"  [PQN Follow-up] Normalized {len(fu_paired)} paired follow-up metabolites.")
        print(f"  {has_fu.sum()} samples with follow-up draws, median dilution = {dilution_fu[has_fu].median():.3f}\n")
    else:
        dilution_fu = pd.Series(1.0, index=df.index)

    return df, dilution_base, dilution_fu


def plot_pqn_dilution_distribution(dilution_base, is_cancer, filename='dilution_factor_distribution.png'):
    """Plot distribution of fitted PQN dilution factors for cancer cases vs controls."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(np.log2(dilution_base[is_cancer]), bins=20, alpha=0.6, color='#d62728', label='Cancer')
    ax.hist(np.log2(dilution_base[~is_cancer]), bins=20, alpha=0.6, color='#1f77b4', label='Control')
    ax.axvline(0, color='black', linestyle='--', lw=1, alpha=0.7, label='Reference (log2=0)')
    ax.set_xlabel(r'$\log_2$(fitted dilution factor)', fontsize=12)
    ax.set_ylabel('Sample count', fontsize=12)
    ax.set_title('Distribution of Fitted PQN Dilution Factors\n(Cases vs Controls)', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


def prepare_metabolites(df, use_pqn=USE_PQN):
    """Normalize or retain raw metabolite values based on the use_pqn flag."""
    global LOG2FC_EPS
    if use_pqn:
        LOG2FC_EPS = LOG2FC_EPS_PQN
        print("\n" + "="*70)
        print("METABOLITE NORMALIZATION: PQN (Probabilistic Quotient Normalization)")
        print("="*70)
        df, dilution_base, dilution_fu = apply_pqn_normalization(df)
        plot_pqn_dilution_distribution(dilution_base, (df['Cohort'] == 'CANCER').to_numpy())
        return df, dilution_base, dilution_fu
    else:
        LOG2FC_EPS = LOG2FC_EPS_RAW
        print("\n" + "="*70)
        print("METABOLITE NORMALIZATION: RAW PEAK INTENSITIES (PQN DISABLED)")
        print("="*70)
        return df, pd.Series(1.0, index=df.index), pd.Series(1.0, index=df.index)


def main(use_pqn=USE_PQN):
    loader = DataLoader(FILE_PATH)
    _dc.provenance_banner(FILE_PATH)
    df_raw = loader.load_data()
    df_raw, _audit = _dc.apply_corrections(
        df_raw, audit_path=OUTPUT_DIR / 'data_corrections_oct15.csv')
    df_clean = loader.clean_data(df_raw)
    
    cases_count = len(df_clean[df_clean['Cohort'] == 'CANCER'])
    controls_count = len(df_clean[df_clean['Cohort'] == 'No'])
    print("\n" + "="*50)
    print("FINAL COHORT SIZES (After all exclusions)")
    print("="*50)
    print(f"Total valid CANCER cases: {cases_count}")
    print(f"Total valid CONTROLS:     {controls_count}")
    print(f"Total samples processed:  {len(df_clean)}")
    print("="*50 + "\n")

    # Metabolite Normalization: PQN or Raw
    df_clean, dilution_base, dilution_fu = prepare_metabolites(df_clean, use_pqn=use_pqn)

    # 2. Get Biomarker columns
    met_cols, meth_cols, social_cols, met_cols_cc = loader.get_feature_groups(df_clean)

    # 0. Table 1: baseline characteristics of the cancer cohort
    build_table_one(df_clean, met_cols=met_cols)

    # 1. Clinical Data Exploration Plot
    plot_clinical_distribution(df_clean)

    # 1.5. Treatment Data Visualizations
    analyze_treatments(df_clean)
    plot_treatment_timing_stats(df_clean)
    plot_followup_timing_stats(df_clean)
    plot_therapy_timeline(df_clean)
    plot_therapy_timeline(df_clean, deceased_only=True)

    # 3. Wilcoxon Tests (Metabolites - log transformed)
    print("Testing Metabolites (cancer vs control: no sparsity filter)...")
    met_results = rank_sum_test(df_clean, met_cols_cc, is_log=True)
    
    # Increased in cases
    met_inc = met_results[(met_results['FDR'] < 0.05) & (met_results['Direction'] == 'Increased')]
    met_inc = met_inc.sort_values('FDR').head(15)
    
    # Decreased in cases
    met_dec = met_results[(met_results['FDR'] < 0.05) & (met_results['Direction'] == 'Decreased')]
    met_dec = met_dec.sort_values('FDR').head(15)
    # set_trace()
    
    plot_boxplots(df_clean, met_inc, 
                  "Top 15 metabolites with significantly increased levels in cases (FDR < 0.05)",
                  "boxplot_top_metabolites_increased.png", is_log=True, hue_col='Cohort', hue_order=['CANCER', 'No'])
                  
    plot_boxplots(df_clean, met_dec, 
                  "Top 15 metabolites with significantly decreased levels in cases (FDR < 0.05)",
                  "boxplot_top_metabolites_decreased.png", is_log=True, hue_col='Cohort', hue_order=['CANCER', 'No'])
                  
    plot_volcano(met_results, "Metabolites (Cancer vs Control)", "volcano_metabolites_cancer_vs_control.png",
                 group_a='Cancer', group_b='Control')
    plot_detection_tier(met_results, "Metabolites detected at different rates (NSCLC vs Control)",
                        "detection_metabolites_cancer_vs_control.png",
                        group_a='Cancer', group_b='Control')
    plot_metabolite_distributions(df_clean, met_results, n_top=10)
    # also a cross-cohort comparison, so it uses the unfiltered list
    pretreatment_results = analyze_pretreatment_vs_control(df_clean, met_cols_cc, n_top=15)
    analyze_metastatic_vs_nonmetastatic(df_clean, met_cols, n_top=10)
    analyze_methylation_by_stage(df_clean, meth_cols)
    analyze_early_vs_late(df_clean, met_cols, n_top=10)
    analyze_methylation_early_vs_late(df_clean, meth_cols)
    analyze_smoking_pack_year_split(df_clean, met_cols, n_top=10)
    plot_top_km_metabolites_panel(df_clean, met_cols, n_top=10, n_cols=4, filename='km_top10_metabolites_panel.png')
    plot_top_km_metabolites_panel(df_clean, met_cols, n_top=20, n_cols=4, filename='km_top20_metabolites_panel.png')
    # Threshold-based companion: every metabolite reaching log-rank p <= 0.1
    plot_top_km_metabolites_panel(df_clean, met_cols, n_cols=4, p_max=0.1)
    plot_top_sdoh_km_panel(df_clean, social_cols, n_top=10, n_cols=4, filename='km_top10_sdoh_panel.png')
    plot_dnam_km_panel(df_clean, meth_cols, filename='km_dnam_panel.png')
    plot_clinical_vars_km_panel(df_clean, filename='km_clinical_vars_panel.png')

    # 3.5 Univariate Cox PH survival analysis
    cox_results = analyze_univariate_cox(df_clean, met_cols, meth_cols, social_cols, n_top=20)

    # 3.6 Analysis F: Metabolite Normalization Score
    ns_results = analyze_metabolite_normalization_score(
        df_clean, met_cols, pretreatment_results=pretreatment_results
    )

    # 3.7 Before vs after treatment (paired) metabolite analysis
    prepost_change, prepost_by_treatment, prepost_delta_io = analyze_pre_post_treatment(df_clean, met_cols)

    # 4. Wilcoxon Tests (Methylation)
    print("Testing Methylation...")
    meth_results = rank_sum_test(df_clean, meth_cols, is_log=False)
    
    meth_sig = meth_results[(meth_results['FDR'] < 0.05)]
    meth_top = meth_sig.sort_values('FDR').head(15)
    
    plot_boxplots(df_clean, meth_top,
                  "DNAm differences (Cancer vs Control)",
                  "boxplot_top_methylation_diffs.png", is_log=False, hue_col='Cohort', hue_order=['CANCER', 'No'])
                  
    # Methylation splitting isolated across the Cancer Survival median
    print("Running Methylation split over survival outcomes...")
    df_dead = df_clean[(df_clean['Cohort'] == 'CANCER') & (df_clean['Event'] == 1)].copy()
    if not df_dead.empty:
        med_surv = df_dead['Time to Event'].median()
        df_dead['Survival_Group'] = np.where(df_dead['Time to Event'] > med_surv, 'High Survival', 'Low Survival')
        
        # Metabolites split over survival outcomes for Volcano plot
        print("Running Metabolites split over survival outcomes for Volcano plot...")
        met_surv_res = rank_sum_test(df_dead, met_cols, is_log=True, group_col='Survival_Group', group1='High Survival', group2='Low Survival')
        plot_volcano(met_surv_res, "Metabolites (High vs Low Survival)", "volcano_metabolites_survival_split.png")
        
        meth_surv_res = rank_sum_test(df_dead, meth_cols, is_log=False, group_col='Survival_Group', group1='Low Survival', group2='High Survival')
        
        # Pull top markers predicting strict mortality timeline
        meth_surv_sig = meth_surv_res[meth_surv_res['p_value'] < 0.05]
        if meth_surv_sig.empty:
            meth_surv_sig = meth_surv_res.sort_values('p_value').head(15)
        else:
            meth_surv_sig = meth_surv_sig.sort_values('p_value').head(15)
            
        plot_boxplots(df_dead, meth_surv_sig,
                      "DNAm differences in Dead Patients (Low vs High Survival)",
                      "boxplot_top_methylation_survival_diffs.png", is_log=False, hue_col='Survival_Group', hue_order=['Low Survival', 'High Survival'])
                      
        # Social Indicators splitting isolated across the Cancer Survival median
        print("Running Social Indicators split over survival outcomes...")
        social_surv_res = rank_sum_test(df_dead, social_cols, is_log=False, group_col='Survival_Group', group1='Low Survival', group2='High Survival')
        
        # Pull top markers
        soc_surv_sig = social_surv_res.sort_values('p_value').head(15)
            
        plot_boxplots(df_dead, soc_surv_sig,
                      "Social Indicators in Dead Patients (Low vs High Survival)",
                      "boxplot_top_social_indicators_survival_diffs.png", is_log=False, standardize=True, hue_col='Survival_Group', hue_order=['Low Survival', 'High Survival'])
                  
    # 5. KM Plot generation
    print("Generating KM Plots for all features...")
    res_met = generate_km_plots(df_clean, met_cols, is_categorical=False)
    res_meth = generate_km_plots(df_clean, meth_cols, is_categorical=False)
    res_soc = generate_km_plots(df_clean, social_cols, is_categorical=False)
    
    clinical_continuous = ['Age at the time of diagnosis', 'BMI', 'Pack Year History', 'Size 1 of tumor (cms)']
    clinical_categorical = ['Sex', 'Race', 'Smoking Status', 'Stage', 'Histological Classification', 'Primary Tumor Location']
    comorb_cols = ['HTN', 'CHF', 'CAD', 'DM', 'Thyroid', 'Emphysema', 'Asthma', 'COPD']
    
    res_clin_cont = generate_km_plots(df_clean, clinical_continuous, is_categorical=False)
    res_clin_cat = generate_km_plots(df_clean, clinical_categorical + comorb_cols, is_categorical=True)
                  
    # Combine and print top 20
    all_res = pd.concat([res_met, res_meth, res_soc, res_clin_cont, res_clin_cat], ignore_index=True)
    if not all_res.empty:
        all_res = all_res.dropna(subset=['P_Value'])
        if len(all_res) > 0:
            all_res['FDR'] = fdrcorrection(all_res['P_Value'])[1]
            top_20 = all_res.sort_values('P_Value').head(20)
            print("\n" + "="*80)
            print("TOP 20 FEATURES CORRELATED WITH SURVIVAL (Log-Rank Test)")
            print("="*80)
            print(top_20.to_markdown(index=False))
            print("="*80 + "\n")
            
    # 7. Survival split box plots
    print("\nRunning Boxplot splits for survival...")
    res_surv_box = analyze_survival_boxplots(df_clean, met_cols + meth_cols + clinical_continuous,
                                             met_cols=met_cols, meth_cols=meth_cols, panel_p_max=0.1)

    # 7.5 Race (Black vs non-Black) and ADI (high vs low) metabolite differences
    res_race = analyze_race_boxplots(df_clean, met_cols)
    res_adi  = analyze_adi_boxplots(df_clean, met_cols)

    # 8. Boxplot splits for social parameters significance filters
    res_soc_box = analyze_social_boxplots(df_clean, social_cols, met_cols, SOCIAL_MODE_MAX)
                  
    # 9. Metabolite Correlations
    met_corr_matrix = analyze_metabolite_correlations(df_clean, met_cols)
                  
    # 10. Social Correlations
    soc_corr_matrix = analyze_social_correlations(df_clean, social_cols)

    # 10.2. Social x Metabolite cross-correlations (cases and controls)
    print("\nSocial x metabolite cross-correlations...")
    analyze_social_metabolite_correlations(df_clean, social_cols, met_cols, cohort='CANCER')
    analyze_social_metabolite_correlations(df_clean, social_cols, met_cols, cohort='No')


    # 10.5. IO Responder Differential Analysis
    print("\nRunning IO Responder Split (Yes vs No) on Metabolites...")
    io_col = 'IO Responder Yes or No'
    if io_col in df_clean.columns:
        df_io = df_clean.dropna(subset=[io_col]).copy()
        df_io[io_col] = df_io[io_col].astype(str).str.strip().str.title()
        
        io_met_results = rank_sum_test(df_io, met_cols, is_log=True, group_col=io_col, group1='Yes', group2='No')
        
        if io_met_results is not None and not io_met_results.empty:
            io_met_sig = io_met_results.sort_values('p_value').head(15)
            plot_boxplots(df_io, io_met_sig,
                          "Top 15 Metabolites: IO Responder (Yes vs No)",
                          "boxplot_top_metabolites_io_responder.png", is_log=True, 
                          hue_col=io_col, hue_order=['Yes', 'No'])
                          
    # 11. Excel Export
    excel_path = FIG_DIR / 'all_statistical_tests.xlsx'
    print(f"\nWriting all statistical test results to {excel_path} ...")
    with pd.ExcelWriter(excel_path) as writer:
        if met_results is not None and not met_results.empty:
            met_results.to_excel(writer, sheet_name='Diff_Metabolites', index=False)
        if meth_results is not None and not meth_results.empty:
            meth_results.to_excel(writer, sheet_name='Diff_Methylation', index=False)
        if all_res is not None and not all_res.empty:
            all_res.to_excel(writer, sheet_name='KM_LogRank_Overall', index=False)
        if res_surv_box is not None and not res_surv_box.empty:
            res_surv_box.to_excel(writer, sheet_name='Survival_Boxplot_Stats', index=False)
        if res_soc_box is not None and not res_soc_box.empty:
            res_soc_box.to_excel(writer, sheet_name='Social_Boxplot_Stats', index=False)
        if res_race is not None and not res_race.empty:
            res_race.to_excel(writer, sheet_name='Race_Metabolites', index=False)
        if res_adi is not None and not res_adi.empty:
            res_adi.to_excel(writer, sheet_name='ADI_Metabolites', index=False)
        if met_corr_matrix is not None and not met_corr_matrix.empty:
            met_corr_matrix.to_excel(writer, sheet_name='Metabolite_Correlations', index=True)
        if soc_corr_matrix is not None and not soc_corr_matrix.empty:
            soc_corr_matrix.to_excel(writer, sheet_name='Social_Correlations', index=True)
        if 'met_surv_res' in locals() and not met_surv_res.empty:
            met_surv_res.to_excel(writer, sheet_name='Surv_Diff_Metabolites', index=False)
        if 'meth_surv_res' in locals() and not meth_surv_res.empty:
            meth_surv_res.to_excel(writer, sheet_name='Surv_Diff_Methylation', index=False)
        if 'social_surv_res' in locals() and not social_surv_res.empty:
            social_surv_res.to_excel(writer, sheet_name='Surv_Diff_Social', index=False)
        if 'io_met_results' in locals() and io_met_results is not None and not io_met_results.empty:
            io_met_results.to_excel(writer, sheet_name='IO_Responder_Diffs', index=False)
        if cox_results is not None and not cox_results.empty:
            cox_results.to_excel(writer, sheet_name='Univariate_Cox_Survival', index=False)
        if 'ns_results' in locals() and ns_results is not None and not ns_results.empty:
            ns_results.to_excel(writer, sheet_name='Normalization_Score', index=False)
        if 'prepost_change' in locals() and prepost_change is not None and not prepost_change.empty:
            prepost_change.to_excel(writer, sheet_name='PrePost_Change', index=False)
        if 'prepost_by_treatment' in locals() and prepost_by_treatment is not None and not prepost_by_treatment.empty:
            prepost_by_treatment.to_excel(writer, sheet_name='PrePost_By_Treatment', index=False)
        if 'prepost_delta_io' in locals() and prepost_delta_io is not None and not prepost_delta_io.empty:
            prepost_delta_io.to_excel(writer, sheet_name='PrePost_Delta_IO', index=False)
            
    print(f"\nDone! All plots and statistics saved securely to {FIG_DIR}.")

def table_one_only(use_pqn=USE_PQN):
    """Regenerate just Table 1, without running the full (slow) analysis pipeline."""
    loader = DataLoader(FILE_PATH)
    _dc.provenance_banner(FILE_PATH)
    _raw, _audit = _dc.apply_corrections(
        loader.load_data(), audit_path=OUTPUT_DIR / 'data_corrections_oct15.csv')
    df_clean = loader.clean_data(_raw)
    df_clean, _, _ = prepare_metabolites(df_clean, use_pqn=use_pqn)
    met_cols, _, _, _ = loader.get_feature_groups(df_clean)
    build_table_one(df_clean, met_cols=met_cols)


if __name__ == "__main__":
    import sys
    use_pqn = USE_PQN
    if '--no-pqn' in sys.argv or '--raw' in sys.argv:
        use_pqn = False
    elif '--pqn' in sys.argv:
        use_pqn = True

    if '--table-one' in sys.argv:
        table_one_only(use_pqn=use_pqn)
    else:
        main(use_pqn=use_pqn)
