!pip install pydeseq2

%%writefile config.py
RAW_COUNTS_PATH   = "data/raw_counts.csv"   # genes x samples raw count matrix (input)
SAMPLE_INFO_PATH  = "data/sample_info.csv"  # metadata: which sample is which condition/timepoint/replicate
OUTPUT_DIR        = "outputs"               # every step writes its output .pkl here

# --- Column names expected in sample_info.csv ---
# YOU NEED TO MAKE THIS FILE ON YOUR OWN - YOU CAN SEE SAMPL_INFO.CSV FOR REFERENCE HOW IT SHOUDL LOOK LIKE
# THIS FILE DEFINES WHAT EACH COLUMN MEANS ? LIKE WHAT DOES IT MEAN FOR A COLUMN NAMED amp_14_t_0= ampiclin 1/4 x mic is used, t=0 sec 

SAMPLE_ID_COL   = "sample_id"
TIMEPOINT_COL   = "timepoint"
CONDITION_COL   = "condition"
REPLICATE_COL   = "replicate"

#NOW WE WILL INTRODUCE THRESHOLDS FOR FILTERING UNWANTED NOISY ERRORFUL SAMPLES 

MIN_COUNT_PER_GENE     = 10
MIN_SAMPLES_EXPRESSED  = 2
MIN_TOTAL_LIBRARY_SIZE = 1e5
REPLICATE_CORR_THRESHOLD = 0.85

NORMALIZATION_METHOD = "deseq2"   
LOG_PSEUDOCOUNT = 1               #global variable set to 1 so that if counts is 0, their log doesn't go to minus infinity
EXPORT_FOR_DYNGENIE3 = True       #Acts like a toggle switch, and the whole downstream of this part can be controlled from here

# -------------------------- LOADING DATA -------------------------------------------------
# PURPOSE: load raw counts + sample metadata from disk,
# verify the two files actually refer to the same samples,
# print a quick summary, and cache both as pickles so the
# next step doesn't need to re-parse CSVs.

%%writefile step1_load_data.py
import pandas as pd
import os
import config

# DEFINED TO READ THE RAW COUNTS FILE & SHOW ERROR IF THE PATH OF THE UPLOADED RAW COUNTS- FILE IS NOT CORRECT
def load_raw_counts(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Could not find raw counts file at '{path}'.")
    return pd.read_csv(path, index_col=0)    #By default, pd.read_csv() creates a numerical row index (0, 1, 2, ...) 
                                             # and treats your gene names (like Gene_A, Gene_B) as a standard data column.
                                             # INDEX_COL=0 TELLS PANDA THAT FIRST ROW IS GENE ID NOT RAW COUNTS;

# DEFINED TO READ SAMPLE INFO FILE & SHOW ERROR IF THE PATH OF THE UPLOADED SAMPLE INFO- FILE IS NOT CORRECT
def load_sample_info(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Could not find sample info file at '{path}'.")
    return pd.read_csv(path)

#CHECKING IF THE LIST WE HAVE UPLOADDED CORRECTLY FITS WITH EACH OTHER;
#EXAMPLE- IF SAMPLE_A IS WRITTEN IN RAWCOUNTS FILE SAME MUST BE IN SAMPLE INFO FILE WITHOUT ANY EROR OR TYPO
#CREATING SETS TO KEEP EXPRESSION VALUES AND OTHER COLUMNS SEPERATE IN SETS TO DO MATHEMATICAL OPERATIONS 
def validate_alignment(counts: pd.DataFrame, sample_info: pd.DataFrame):
    count_samples = set(counts.columns)                         # A SET WITH NAME OF ALL SAMPLES WHICH IS COLUMN HEADER FROM RAWCOUNTS FILE
    meta_samples = set(sample_info[config.SAMPLE_ID_COL])       # A SET WITH NAME OF ALL SAMPLES WHICH IS FIRST ROW FROM SAMPLEINFO FILE 
    missing_in_meta = count_samples - meta_samples              # COMPARISON THROUGH SET SUBSTRACTION
    missing_in_counts = meta_samples - count_samples
    problems = []                                               # IF ANY MISMATCH IS FOUND STOP THE PIPELINE 
    if missing_in_meta: problems.append(f"  - Missing from sample_info.csv: {sorted(missing_in_meta)}")
    if missing_in_counts: problems.append(f"  - Missing from count matrix: {sorted(missing_in_counts)}")
    if problems: raise ValueError("Sample mismatch:\n" + "\n".join(problems))   #TO PRINT THE PROBLEM IN OUTPUT 

    # RE-ORDERING STEPS--> sample columns in our count matrix (counts) and the rows in our metadata (sample_info) must be in the EXACT SAME ORDER
    # Explicitly rename 'index' back to config's SAMPLE_ID_COL after reindexing
    sample_info = sample_info.set_index(config.SAMPLE_ID_COL).loc[counts.columns]
    sample_info = sample_info.reset_index().rename(columns={"index": config.SAMPLE_ID_COL})

    print(f"✓ All {len(count_samples)} samples match.")
    return sample_info

# DEFINING FUNVCN To PRINT THE OUTPUT AND SHOW (HOW MANY RETAINED / FROM HOW MANY WERE AVAILABLE) & OTHER THINGS TOO
def summarize_input(counts: pd.DataFrame, sample_info: pd.DataFrame):
    print("\n── Input summary ──")
    print(f"Genes:               {counts.shape[0]}")       # NO OF ROWS=GENES
    print(f"Samples:             {counts.shape[1]}")       # NO OF COLUMNS= SAMPLES
    print(f"Timepoints:           {sorted(sample_info[config.TIMEPOINT_COL].unique())}")
    print(f"Conditions:           {sample_info[config.CONDITION_COL].unique().tolist()}\n")

# CALLING FUNCN TO PRINT THE SUMMARY
if __name__ == "__main__":
    counts = load_raw_counts(config.RAW_COUNTS_PATH)              # FOR LOADING both files from the paths defined in config.py
    sample_info = load_sample_info(config.SAMPLE_INFO_PATH)       
    sample_info = validate_alignment(counts, sample_info)
    summarize_input(counts, sample_info)                          #CALLING THE FUNCN DEFINED ABOVE FOR SUMMARY
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)                 # MAKING sure the outputs/ folder exists before writing into it.
    counts.to_pickle(f"{config.OUTPUT_DIR}/01_raw_counts.pkl")          # TURNING THE FILE TO PICKLE (PRESERVES DATATYPES USED IN FILE)
    sample_info.to_pickle(f"{config.OUTPUT_DIR}/01_sample_info.pkl")


# -------------------------- QUALITY CONTROL -------------------------------------------------

# PURPOSE: three independent QC checks —
#   1. flag samples that were sequenced too shallowly
#   2. drop genes that are essentially never expressed (noise, not signal)
#   3. flag replicate pairs that don't correlate well (bad replicate / mislabeled sample) 

%%writefile step2_quality_control.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import config
sns.set_style("whitegrid")

def check_library_sizes(counts: pd.DataFrame) -> pd.Series:
    lib_sizes = counts.sum(axis=0)
    flagged = lib_sizes[lib_sizes < config.MIN_TOTAL_LIBRARY_SIZE]   # MIN_TOTAL_LIBRARY_SIZE as defined at the top as Global VAribale;
    if len(flagged) > 0:                                            } 
        for s, v in flagged.items():                                } # TO PRINT WHICH SAMPLES HAVE LOW DEPTH READ
            print(f"  ⚠ Low depth -> {s}: {v:.0f} reads")          }
    else:                                            
        print(f"✓ All samples have enough reads.")   
    return lib_sizes

# CHECKING FOR LOW COUNT GENES;
# (counts >= MIN_COUNT_PER_GENE) --> True wherever a gene has at least 10 reads(as defined at the top) in that particular sample
# >= MIN_SAMPLES_EXPRESSED -> True/False mask per gene: "is this gene expressed in enough samples to be worth keeping?"
def filter_low_count_genes(counts: pd.DataFrame) -> pd.DataFrame:
    expressed_mask = (counts >= config.MIN_COUNT_PER_GENE).sum(axis=1) >= config.MIN_SAMPLES_EXPRESSED # sum(axis=1)--> sums each row
    filtered = counts.loc[expressed_mask]
    print(f"✓ Kept {filtered.shape[0]} / {counts.shape[0]} genes.")
    return filtered

#CHCEKING FOR REPLICATE DUPLICATES IF PRESENT TO BE REMOVED
# COMPARING WITH PERASON CORRELATION - BUT IT IS HIGHLY SKEWED RIGHT, SO FIRST LOG TRANSFROM VALUES
def check_replicate_correlation(counts: pd.DataFrame, sample_info: pd.DataFrame):
    log_counts = np.log2(counts + 1)                      # LOG TRANSFFROMING compresses the skewness for better comparison
    corr_matrix = log_counts.corr(method="pearson")      # PAIRWISE CORRELATION bw EACH COLUMN(SAMPLES)--> sample x sample correlation matrix
    grouped = sample_info.groupby([config.TIMEPOINT_COL, config.CONDITION_COL])[config.SAMPLE_ID_COL].apply(list)
    for group_key, samples in grouped.items():
        if len(samples) < 2: continue
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                r = corr_matrix.loc[samples[i], samples[j]]
                if r < config.REPLICATE_CORR_THRESHOLD:
                    print(f"  ⚠ Low corr -> {group_key}: {samples[i]} vs {samples[j]} (r = {r:.3f})")
    return corr_matrix

if __name__ == "__main__":
    counts = pd.read_pickle(f"{config.OUTPUT_DIR}/01_raw_counts.pkl")
    sample_info = pd.read_pickle(f"{config.OUTPUT_DIR}/01_sample_info.pkl")
    print("\n── Step 2: Quality Control ──")
    lib_sizes = check_library_sizes(counts)
    counts_filtered = filter_low_count_genes(counts)
    corr_matrix = check_replicate_correlation(counts_filtered, sample_info)
    counts_filtered.to_pickle(f"{config.OUTPUT_DIR}/02_filtered_counts.pkl")

%%writefile step3_normalization.py
import pandas as pd
import numpy as np
import config

def normalize_deseq2(counts: pd.DataFrame, sample_info: pd.DataFrame) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    counts_t = counts.T.astype(int)
    meta = sample_info.set_index(config.SAMPLE_ID_COL).loc[counts_t.index]
    design_col = config.CONDITION_COL
    if meta[design_col].nunique() < 2: design_col = config.TIMEPOINT_COL; meta[design_col] = meta[design_col].astype(str)
    dds = DeseqDataSet(counts=counts_t, metadata=meta, design=f"~{design_col}", inference=DefaultInference(), quiet=True)
    dds.fit_size_factors()
    normalized = counts_t.div(dds.obs["size_factors"], axis=0).T
    print("✓ DESeq2 normalisation complete.")
    return normalized

# for printing the result if this step
if __name__ == "__main__":
    counts = pd.read_pickle(f"{config.OUTPUT_DIR}/02_filtered_counts.pkl")
    sample_info = pd.read_pickle(f"{config.OUTPUT_DIR}/01_sample_info.pkl")
    print("\n── Step 3: Normalisation ──")
    if config.NORMALIZATION_METHOD == "deseq2":
        normalized = normalize_deseq2(counts, sample_info)
    normalized.to_pickle(f"{config.OUTPUT_DIR}/03_normalized_counts.pkl")

%%writefile step4_log_transform_and_explore.py
import pandas as pd
import numpy as np
import config

if __name__ == "__main__":
    normalized = pd.read_pickle(f"{config.OUTPUT_DIR}/03_normalized_counts.pkl")
    print("\n── Step 4: Log Transform ──")
    log_data = np.log2(normalized + config.LOG_PSEUDOCOUNT)
    print(f"✓ Log2 transform applied.")
    log_data.to_pickle(f"{config.OUTPUT_DIR}/04_log_normalized_counts.pkl")

%%writefile step5_export_for_dyngenie3.py
import pandas as pd
import numpy as np
import pickle
import config

def export_per_condition(log_data: pd.DataFrame, sample_info: pd.DataFrame):
    meta = sample_info.set_index(config.SAMPLE_ID_COL)
    conditions = meta[config.CONDITION_COL].unique()
    for cond in conditions:
        cond_samples = meta[meta[config.CONDITION_COL] == cond].index
        cond_data = log_data[cond_samples]
        cond_meta = meta.loc[cond_samples]
        avg_by_timepoint = {}
        for tp in sorted(cond_meta[config.TIMEPOINT_COL].unique()):
            tp_samples = cond_meta[cond_meta[config.TIMEPOINT_COL] == tp].index
            avg_by_timepoint[tp] = cond_data[tp_samples].mean(axis=1)
        avg_matrix = pd.DataFrame(avg_by_timepoint).T

        # Strip out to plain numpy/python types -- no pandas objects, no dtype version issues
        time_points = np.array(sorted(cond_meta[config.TIMEPOINT_COL].unique()), dtype=float)
        gene_names = [str(g) for g in avg_matrix.columns]
        TS_data = [np.array(avg_matrix.values, dtype=float)]
        decay_rates = None

        with open(f"{config.OUTPUT_DIR}/05_dynGENIE3_input_{cond}.pkl", "wb") as f:
            pickle.dump((TS_data, time_points, decay_rates, gene_names), f, protocol=4)

    print("✓ Reshaped matrices exported for dynGENIE3.")

if __name__ == "__main__":
    log_data = pd.read_pickle(f"{config.OUTPUT_DIR}/04_log_normalized_counts.pkl")
    sample_info = pd.read_pickle(f"{config.OUTPUT_DIR}/01_sample_info.pkl")
    print("\n── Step 5: Export for dynGENIE3 ──")
    export_per_condition(log_data, sample_info)

## Upload Your Own Data

This pipeline expects two files inside a `data/` folder:

1. **`raw_counts.csv`** — a raw count matrix, **genes as rows, samples as columns**. The first column must be the gene ID (this becomes the row index).

2. **`sample_info.csv`** — one row per sample, with these columns:
   - `sample_id` — must exactly match the column names used in `raw_counts.csv`
   - `timepoint` — numeric (e.g. 0, 1, 4, 8, 24)
   - `condition` — e.g. `control`, `antibiotic`
   - `replicate` — replicate number (e.g. 1, 2)

Run the cell below, then use the file picker to upload both CSVs (any names are fine — you'll be asked which file is which). They'll be copied into `data/raw_counts.csv` and `data/sample_info.csv` automatically.

import os
import shutil
from google.colab import files

os.makedirs("data", exist_ok=True)

print("Select your RAW COUNTS file (genes x samples):")
uploaded_counts = files.upload()
counts_filename = list(uploaded_counts.keys())[0]
shutil.move(counts_filename, "data/raw_counts.csv")

print("\nSelect your SAMPLE INFO file (sample_id, timepoint, condition, replicate):")
uploaded_info = files.upload()
info_filename = list(uploaded_info.keys())[0]
shutil.move(info_filename, "data/sample_info.csv")

print("\n✓ Files saved to data/raw_counts.csv and data/sample_info.csv")

%%writefile run_pipeline.py
import subprocess
import sys
import os
import config

STEPS = [
    "step1_load_data.py",
    "step2_quality_control.py",
    "step3_normalization.py",
    "step4_log_transform_and_explore.py",
    "step5_export_for_dyngenie3.py",
]

if __name__ == "__main__":
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)
    print("=" * 45)
    print(" RUNNING FULL DATA PREPARATION PIPELINE")
    print("=" * 45)
    for step_script in STEPS:
        result = subprocess.run([sys.executable, step_script])
        if result.returncode != 0:
            print(f"\n✗ Pipeline stopped — '{step_script}' failed.")
            sys.exit(1)
    print("=" * 45)
    print(" PIPELINE COMPLETE. Find outputs in 'outputs/'")
    print("=" * 45)

!python run_pipeline.py
