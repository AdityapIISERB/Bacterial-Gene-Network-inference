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

MIN_COUNT_PER_GENE     = 10          # edgeR default=10; doi:10.12688/f1000research.8987.2
MIN_SAMPLES_EXPRESSED  = 2           #Although we are using pydeseq here, but edgeR automatically sets it to 2
MIN_TOTAL_LIBRARY_SIZE = 5e6         #Pearson Correlation, 
REPLICATE_CORR_THRESHOLD = 0.90     

NORMALIZATION_METHOD = "deseq2"   
LOG_PSEUDOCOUNT = 1               #global variable set to 1 so that if counts is 0, their log doesn't go to minus infinity
EXPORT_FOR_DYNGENIE3 = True       #Acts like a toggle switch, and the whole downstream of this part can be controlled from here

# -------------------------- READING DATA ----------------------------------------------------------
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


# -------------------------- QUALITY CONTROL ---------------------------------------------------------

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
    if len(flagged) > 0:                                             
        for s, v in flagged.items():                                 # TO PRINT WHICH SAMPLES HAVE LOW DEPTH READ
            print(f"  ⚠ Low depth -> {s}: {v:.0f} reads")          
    else:                                            
        print(f"✓ All samples have enough reads.")   
    return lib_sizes

# CHECKING FOR LOW COUNT GENES;
# (counts >= MIN_COUNT_PER_GENE) --> True wherever a gene has at least 10 reads(as defined at the top) in that particular sample
# >= MIN_SAMPLES_EXPRESSED -> True/False mask per gene: "is this gene expressed in enough samples to be worth keeping?"
# Creates a True/False DataFrame for each cell of the exact same dimensions as counts. Every entry with > 10 reads becomes True, while entries with < 10 reads become False. 
# fro each cell true=1 & false=0; then sums each row and tells how many genes is > MIN_SAMPLES_EXPRESSED (global variable)
def filter_low_count_genes(counts: pd.DataFrame) -> pd.DataFrame:
    expressed_mask = (counts >= config.MIN_COUNT_PER_GENE).sum(axis=1) >= config.MIN_SAMPLES_EXPRESSED    #sum(axis=1)--> sums each row
    filtered = counts.loc[expressed_mask]          # filtered only retains where expressed masks is TRUE
    print(f"✓ Kept {filtered.shape[0]} / {counts.shape[0]} genes.")
    return filtered

# CHECKING FOR BAD OR DUPLICATE REPLICATES IF PRESENT TO BE REMOVED
# PEARSON CORRELATION tells - checks if two biological replicates behave identically;  
# compare bw 2 replicates; LOG_rep_1 vs LOG_rep_2--> r=1(if one goes up, another goes up too), r=-1(opposite), r=0(no pattern)
# COMPARING WITH PERASON CORRELATION - BUT IT IS HIGHLY SKEWED RIGHT, SO FIRST LOG TRANSFORM VALUES
# is this gene expressed in enough samples to be worth keeping ?
# ACTUALLY REMOVES one sample from each low-correlating replicate pair. Between two disagreeing replicates we can't know for certain
# WHICH one is "wrong" from correlation alone -- but the one with the lower total library size is the more likely technical failure, so
# that's the tie-breaker used here. This runs AFTER filter_low_depth_samples, so by this point most genuinely broken (very-low-depth) 
# samples are already gone; this catches replicate pairs that still disagree for reasons other than raw depth.

def check_and_filter_replicates(counts: pd.DataFrame, sample_info: pd.DataFrame,
                                 lib_sizes: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    log_counts = np.log2(counts + 1)                      # LOG TRANSFFROMING compresses the skewness for better comparison
    corr_matrix = log_counts.corr(method="pearson")      # PAIRWISE CORRELATION bw EACH COLUMN(SAMPLES)--> sample x sample correlation matrix
    #grouped on the basis of condition and timepoint
    grouped = sample_info.groupby([config.TIMEPOINT_COL, config.CONDITION_COL])[config.SAMPLE_ID_COL].apply(list)
    to_drop = set()
    for group_key, samples in grouped.items():
        if len(samples) < 2: continue               #to check if there is atleast two samples present to compare otherwise skips(continue)
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                s1, s2 = samples[i], samples[j]
                if s1 in to_drop or s2 in to_drop:
                    continue                        # already dropping one of this pair 
                r = corr_matrix.loc[s1, s2]
                if r < config.REPLICATE_CORR_THRESHOLD:         # THIS CHECKS IF r<0.90(REPLICATE_CORR_THRESHOLD)(defined in config at top) those 2 samples are replicates;
                    worse = s1 if lib_sizes.get(s1, 0) < lib_sizes.get(s2, 0) else s2     #Second check for lib sizes to decide which replicate to drop
                    print(f"  ⚠ Low corr -> {group_key}: {s1} vs {s2} (r = {r:.3f}) — dropping {worse}")
                    to_drop.add(worse)
    if to_drop:
        counts = counts.drop(columns=list(to_drop))
        sample_info = sample_info[~sample_info[config.SAMPLE_ID_COL].isin(to_drop)].reset_index(drop=True)
        # recompute on the now-cleaned data so the returned matrix reflects what's ACTUALLY in the final dataset, not the pre-drop version
        corr_matrix = np.log2(counts + 1).corr(method="pearson")
    return counts, sample_info, corr_matrix

#TO PRINT THE OUTPUT OF THIS STEP:
if __name__ == "__main__":
    counts = pd.read_pickle(f"{config.OUTPUT_DIR}/01_raw_counts.pkl")
    sample_info = pd.read_pickle(f"{config.OUTPUT_DIR}/01_sample_info.pkl")
    print("\n── Step 2: Quality Control ──")
    lib_sizes = check_library_sizes(counts)                                     #check1- seq depth
    counts_filtered = filter_low_count_genes(counts)                            #check2- low count genes
    counts_filtered, sample_info, corr_matrix = check_and_filter_replicates(counts_filtered, sample_info, lib_sizes)     #check3- bad replicates (now actually applies the drop)
    counts_filtered.to_pickle(f"{config.OUTPUT_DIR}/02_filtered_counts.pkl")      # converting back to pickle
    sample_info.to_pickle(f"{config.OUTPUT_DIR}/02_filtered_sample_info.pkl")     # cleaned sample_info, so downstream steps use it too


# -------------------------- NORMALIZATION --------------------------------------------------------------
# PYDESEQ2 NEEDS GENES AS COLUMNS & SAMPLES AS ROWS [OUR EXACT OPPPOSITE], so we will transpose our file
# Say, SampleA (10M reads) and SampleB (5M reads) both passed quality control, Sample A has TWICE as many total reads as Sample B.
# Twice not because expressed highly, but sequenced twice deeply
# lets say , sampleA = 200 count [10M total depth]    {depth will be calculated by adding counts of all genes for sampleA}
#          , sampleB = 100 count [5M total depth]
# Actually both are expressed equally 200/10Million =100/5Million= 0.02%  --> CPM MEHOD(will not be used here)
# MEDIAN OF RATIOS METHOD for Normalization is used:
  #1- pseudo-reference = Calculate geometric mean of all expression values of a gene across all samples 
  #2- Gene Ratios = divide each gene expression with its pseudo-ref calculated
  #3- Size FActor = median of all those ratios is taken for a sample
  # NORMALIZED COUNT = ( RAW COUNTS/SIZE FACTOR )
##  IT WORKS BECAUSE GEOMETRIC MEANS IS BETTER FOR CALCULATION WHAT FOLD EXPRESSION IS CHANGED

%%writefile step3_normalization.py
import pandas as pd
import numpy as np
import config

def normalize_deseq2(counts: pd.DataFrame, sample_info: pd.DataFrame) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    counts_t = counts.T.astype(int)     # T.astype = transpose our row&column of row count matrix
    meta = sample_info.set_index(config.SAMPLE_ID_COL).loc[counts_t.index]   
    design_col = config.CONDITION_COL   # CONDITION as the factor is used DESeq2 uses to fit dispersion/size models
    if meta[design_col].nunique() < 2: design_col = config.TIMEPOINT_COL; meta[design_col] = meta[design_col].astype(str)
    dds = DeseqDataSet(counts=counts_t, metadata=meta, design=f"~{design_col}", inference=DefaultInference(), quiet=True)
    dds.fit_size_factors()
    normalized = counts_t.div(dds.obs["size_factors"], axis=0).T
    print("✓ DESeq2 normalisation complete.")
    return normalized

# for printing the result of this step
if __name__ == "__main__":
    counts = pd.read_pickle(f"{config.OUTPUT_DIR}/02_filtered_counts.pkl")
    sample_info = pd.read_pickle(f"{config.OUTPUT_DIR}/02_filtered_sample_info.pkl")
    print("\n── Step 3: Normalisation ──")
    if config.NORMALIZATION_METHOD == "deseq2":
        normalized = normalize_deseq2(counts, sample_info)
    normalized.to_pickle(f"{config.OUTPUT_DIR}/03_normalized_counts.pkl")
#---------------------------BOX PLOT COMPARISON --------------------------------------------------------
# PURPOSE: produce the two QC plots-  a before/after boxplot of log-scaled counts per sample
# compare RAW counts (post-filtering, pre-normalization)
# against DESeq2-NORMALIZED counts, to visually confirm that normalization actually aligned the sample distributions.
# ---------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import config

sns.set_style("whitegrid")

def plot_boxplot_comparison(raw_counts: pd.DataFrame, norm_counts: pd.DataFrame,
                             out_path: str = "outputs/qc_boxplot_comparison.png"):
    """
    Draws two side-by-side boxplots: log2(counts+1) per sample,
    BEFORE normalization (left) and AFTER normalization (right).

    What to look for: on the left, box medians will sit at different
    heights across samples (different sequencing depth / composition).
    On the right, medians should line up much more closely — that
    alignment is the visual proof that normalization worked.
    """
    # log2(x + 1): pseudocount of 1 avoids log2(0) = -inf. This is a
    # DISPLAY-ONLY transform — it does not touch step4's actual
    # log-transformed pipeline output, it's just so the boxplot isn't
    # crushed by a handful of extremely highly-expressed genes.
    log_raw = np.log2(raw_counts + 1)
    log_norm = np.log2(norm_counts + 1)

    fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

    # sns.boxplot expects "long-form" data by default when given a
    # DataFrame directly it treats each column as one box — exactly
    # what we want (one box per sample).
    sns.boxplot(data=log_raw, ax=axes[0], color="lightcoral")
    axes[0].set_title("Before Normalization (raw counts)")
    axes[0].set_ylabel("log2(count + 1)")
    axes[0].tick_params(axis='x', rotation=90)

    sns.boxplot(data=log_norm, ax=axes[1], color="mediumseagreen")
    axes[1].set_title("After Normalization (DESeq2)")
    axes[1].tick_params(axis='x', rotation=90)

    fig.suptitle("QC: Count Distribution Before vs. After Normalization", fontsize=14)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved boxplot comparison -> {out_path}")

# -------------------------- LOG TRANSFORM --------------------------------------------------------------
## log2(x + pseudocount). The pseudocount(1)(from config.py) is added becaus#e log2(0) is (-infinity)
# genes with a normalize count of 0 in some sample would otherwise break the transform.

%%writefile step4_log_transform_and_explore.py
import pandas as pd
import numpy as np
import config

if __name__ == "__main__":
    normalized = pd.read_pickle(f"{config.OUTPUT_DIR}/03_normalized_counts.pkl")
    sample_info = pd.read_pickle(f"{config.OUTPUT_DIR}/02_filtered_sample_info.pkl")
    print("\n── Step 4: Log Transform ──")
    log_data = np.log2(normalized + config.LOG_PSEUDOCOUNT)
    print(f"✓ Log2 transform applied.")
    log_data.to_pickle(f"{config.OUTPUT_DIR}/04_log_normalized_counts.pkl")

    print("\n── Step 4b: Skewness Diagnostic Plot ──")


# -------------------------- Plotting Skewness --------------------------------------------------------------
# Here We are plotting skewness of data [how many genes (density) have what count]
# Y-axis= density; X-axis= count
# For 2 samples - amplicllin 1/4 x MIC at t=0 rep1; & amplicllin 1/16 x MIC at t=0 rep1
 
    sample_1 = "amp_1_4_t0_rep1"
    sample_2 = "amp_1_16_t0_rep1"
 
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 9)) #For side by side, we use axes
 
    # Plot 1: Top-Left (raw/normalized, sample 1)
    axes[0, 0].hist(normalized[sample_1], bins=50, color='darkgreen')
    axes[0, 0].set_xlabel('Read Count')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title(f'Genome-wide Count Profile ({sample_1})')
 
    # Plot 2: Top-Right (raw/normalized, sample 2)
    axes[0, 1].hist(normalized[sample_2], bins=50, color='darkgreen')
    axes[0, 1].set_xlabel('Read Count')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title(f'Genome-wide Count Profile ({sample_2})')
 
    # Plot 3: Bottom-Left (log-transformed, sample 1)
    axes[1, 0].hist(log_data[sample_1], bins=50, color='darkgreen')
    axes[1, 0].set_xlabel('Log2(Read Count + 1)')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title(f'Log2 Transformed ({sample_1})')
 
    # Plot 4: Bottom-Right (log-transformed, sample 2)
    axes[1, 1].hist(log_data[sample_2], bins=50, color='darkgreen')
    axes[1, 1].set_xlabel('Log2(Read Count + 1)')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title(f'Log2 Transformed ({sample_2})')
 
    plt.tight_layout()
    plt.savefig(f"{config.OUTPUT_DIR}/04b_raw_vs_log_histogram.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved histogram comparison -> {config.OUTPUT_DIR}/04b_raw_vs_log_histogram.png")

# ---------------------------------------PCA PLOT----------------------------------------------
    print("\n── Step 4c: PCA Plot (samples) ──")

    from sklearn.decomposition import PCA
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D

    expr_matrix = log_data.T                     # rows = samples, columns = genes (what PCA expects)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(expr_matrix.values)   # fit + project samples onto PC1 & PC2
    pct_var = pca.explained_variance_ratio_ * 100 # how much of the total variance each PC captures

    pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=expr_matrix.index)
    pca_df = pca_df.merge(
        sample_info.set_index(config.SAMPLE_ID_COL)[[config.CONDITION_COL, config.TIMEPOINT_COL]],
        left_index=True, right_index=True
    )

    # one color family (colormap) per condition; gradient within a condition = timepoint
    # (light -> dark = early -> late). Update these keys if your condition names differ.
    cmaps = {
        "amp_1_4":  LinearSegmentedColormap.from_list("amp14", ["#c6dbef", "#08306b"]),  # blues
        "amp_1_16": LinearSegmentedColormap.from_list("amp16", ["#c7e9c0", "#00441b"]),  # greens
        "cipro":    LinearSegmentedColormap.from_list("cipro", ["#fdd0a2", "#7f2704"]),  # oranges
    }
    tp_min = pca_df[config.TIMEPOINT_COL].min()
    tp_max = pca_df[config.TIMEPOINT_COL].max()

    fig, ax = plt.subplots(figsize=(9, 7))
    for cond, cmap in cmaps.items():
        subset = pca_df[pca_df[config.CONDITION_COL] == cond]
        if subset.empty:
            continue
        norm_tp = (subset[config.TIMEPOINT_COL] - tp_min) / (tp_max - tp_min)
        colors = cmap(norm_tp.values)
        ax.scatter(subset["PC1"], subset["PC2"], c=colors, s=110,
                   edgecolor="black", linewidth=0.6, label=cond, zorder=3)

    ax.set_xlabel(f"PC1 ({pct_var[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pct_var[1]:.1f}% variance)")
    ax.set_title("PCA of Samples")
    ax.grid(alpha=0.3)

    cond_legend = [Line2D([0], [0], marker="o", color="w", markerfacecolor=cmaps[c](0.7),
                          markeredgecolor="black", markersize=10, label=c) for c in cmaps]
    leg1 = ax.legend(handles=cond_legend, title=config.CONDITION_COL, loc="upper left", bbox_to_anchor=(1.02, 1))
    ax.add_artist(leg1)

    gray_cmap = LinearSegmentedColormap.from_list("gray", ["#dddddd", "#222222"])
    tp_vals = sorted(pca_df[config.TIMEPOINT_COL].unique())
    tp_legend = [Line2D([0], [0], marker="o", color="w",
                        markerfacecolor=gray_cmap((t - tp_min) / (tp_max - tp_min)),
                        markeredgecolor="black", markersize=9, label=f"t={t}") for t in tp_vals]
    ax.legend(handles=tp_legend, title="Timepoint\n(light→dark = early→late)", loc="lower left", bbox_to_anchor=(1.02, 0))

    plt.tight_layout()
    plt.savefig(f"{config.OUTPUT_DIR}/04c_pca_plot.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved PCA plot -> {config.OUTPUT_DIR}/04c_pca_plot.png")


# -------------------------- EXPORT READY FOR dynGenie3 --------------------------------------------------------------
#dynGENIE3 expects input structured in a specific way:
  #1 - Diff files for diff conditions --> one file per condn;
  #2- dynGENIE3 models time-series dynamics & expects one single expression value per gene at each timepoint;
    # In this step- FOr 2 replicates in same condn, same timepoint--> we average it and make it a single representative of that point;
  #3- expects time series data as Numpy array rather than pandas data frames that we have used 
  #4- expects Rows= timepoints % Columns= genes  [will be done via Transpose]
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
        #cond_samples is now= ['S1_ctrl', 'S2_ctrl...'], or column present with the header sample_id  
        cond_data = log_data[cond_samples]           #picks only those columns which have names extracted based on conditons
        cond_meta = meta.loc[cond_samples]
        avg_by_timepoint = {}
        for tp in sorted(cond_meta[config.TIMEPOINT_COL].unique()):  #Iterates through every unique timepoint in sorted order t=0, 1, 2, 4...
            tp_samples = cond_meta[cond_meta[config.TIMEPOINT_COL] == tp].index  
            avg_by_timepoint[tp] = cond_data[tp_samples].mean(axis=1)  #averaging out replicates of same condn + same timpepoint 
        avg_matrix = pd.DataFrame(avg_by_timepoint).T                  #Transpose for making it --> row=timepoints & col=genes

        # To convert into Numpy 
        time_points = np.array(sorted(cond_meta[config.TIMEPOINT_COL].unique()), dtype=float)
        gene_names = [str(g) for g in avg_matrix.columns]
        TS_data = [np.array(avg_matrix.values, dtype=float)]
        decay_rates = None

        #For converting into PICKLE:
        with open(f"{config.OUTPUT_DIR}/05_dynGENIE3_input_{cond}.pkl", "wb") as f:
            pickle.dump((TS_data, time_points, decay_rates, gene_names), f, protocol=4)

    print("✓ Reshaped matrices exported for dynGENIE3.")

if __name__ == "__main__":
    log_data = pd.read_pickle(f"{config.OUTPUT_DIR}/04_log_normalized_counts.pkl")
    sample_info = pd.read_pickle(f"{config.OUTPUT_DIR}/02_filtered_sample_info.pkl")
    print("\n── Step 5: Export for dynGENIE3 ──")
    export_per_condition(log_data, sample_info)

# -------------------------- UPLOADING DATA --------------------------------------------------------------
## Upload Your Own Data
#This pipeline expects two files inside a `data/` folder:

# 1. "raw_counts.csv" — a raw count matrix, **genes as rows, samples as columns**. The first column must be the gene ID (this becomes the row index).

# 2. "sample_info.csv" — one row per sample, with these columns: YOU HAVE TO MAKE THIS FILE MANUALLY
  #- sample_id — must exactly match the column names used in `raw_counts.csv`
  #- timepoint — numeric (e.g. 0, 1, 4, 8, 24)
  #- condition — e.g. control, antibiotic
  #- replicate — replicate number (e.g. 1, 2)

# RUN FOLLOWING CODE, then use the file picker to upload both CSVs (any names are fine — you'll be asked which file is which). 
# They'll be copied into "data/raw_counts.csv" and "data/sample_info.csv" automatically.

import os
import shutil
from google.colab import files

# Create the data/ folder if it doesn't already exist (exist_ok=True prevents an error if it's already there).
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


# -------------------------- CORE EXECUTION SCRIPT --------------------------------------------------------------
# PURPOSE: seperated script — runs all 5 steps in order as separate subprocesses, stopping immediately if any step fails

%%writefile run_pipeline.py
import subprocess
import sys
import os
import config

# The exact order matters: each step reads the previous step'spickled output (e.g. step2 needs 01_raw_counts.pkl from step1).
STEPS = [
    "step1_load_data.py",
    "step2_quality_control.py",
    "step3_normalization.py",
    "step4_log_transform_and_explore.py",
    "step5_export_for_dyngenie3.py",
]

if __name__ == "__main__":
    # Make sure outputs/ exists before any step tries to write into it.
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)

    print("=" * 45)
    print(" RUNNING FULL DATA PREPARATION PIPELINE")
    print("=" * 45)

    for step_script in STEPS:
        # sys.executable -> guarantees we use the SAME python interpreter that's running the previous code;
        # Each step is run as its own subprocess, no chance of failing other part of code
        result = subprocess.run([sys.executable, step_script])

        # subprocess.run's returncode is 0 on success, non-zero if the script raised an uncaught exception (like the FileNotFoundError)
        # Checking this immediately prevents other step from running on top of a step that never actually finished;
        if result.returncode != 0:
            print(f"\n✗ Pipeline stopped — '{step_script}' failed.")
            sys.exit(1)

    print("=" * 45)
    print(" PIPELINE COMPLETE. Find outputs in 'outputs/'")
    print("=" * 45)
# -------------------------- RUN EXECUTION SCRIPT --------------------------------------------------------------
!python run_pipeline.py
