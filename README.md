# Dynamic Transcriptomic States and GRN Rewiring in Antibiotic Resistance Evolution

Summer research project
Supervisor: Dr. Garima Rani and Dr. Neetika

## Research Question

Can dynamic transcriptomic states and gene regulatory network (GRN) rewiring — rather than mutations alone — predict and redirect the evolution of antibiotic resistance in bacteria?

Most resistance-evolution studies focus on genetic mutations as the primary driver. This project instead asks whether **regulatory rewiring** (changes in which genes control which other genes, over time, under antibiotic stress) can itself act as a predictive signal or even redirect resistance outcomes.

## Approach

The project treats gene expression as a **time series** under two conditions — `control` and `antibiotic` — and infers a directed gene regulatory network for each condition. Comparing the two networks (edges gained/lost, hub gene shifts, bifurcation points) is the core analysis for detecting rewiring.

Computational workflow:

```
RNA-seq raw counts
        │
        ▼
Preprocessing (QC → DESeq2 normalization → log2 transform)
        │
        ▼
Per-condition time-series export
        │
        ├──► dynGENIE3 (control)     ──► ranked regulatory links (control)
        └──► dynGENIE3 (antibiotic)  ──► ranked regulatory links (antibiotic)
                        │
                        ▼
        RegulonDB validation of inferred edges
                        │
                        ▼
        Hub gene centrality analysis
                        │
                        ▼
        DREM bifurcation analysis (divergence points between conditions)
                        │
                        ▼
        Network visualization & control-vs-antibiotic comparison
```

## Data Format Requirements

**`data/raw_counts.csv`** — raw count matrix, genes as rows, samples as columns. First column is the gene ID.

**`data/sample_info.csv`** — one row per sample, with columns:
| column | description |
|---|---|
| `sample_id` | must exactly match a column name in `raw_counts.csv` |
| `timepoint` | numeric (e.g. 0, 1, 4, 8, 24) |
| `condition` | e.g. `control`, `antibiotic` |
| `replicate` | replicate number |

## Pipeline Steps

### Step 1 — Synthetic Data Generation

`generate_test_data.py` produces a small synthetic RNA-seq count matrix (`data/raw_counts.csv`) and matching sample metadata (`data/sample_info.csv`) for two conditions (`control`, `antibiotic`) across several timepoints, with replicates. This exists purely to **sanity-check the pipeline end to end** before running it on real data — it seeds a few low-depth samples and low-correlation replicates on purpose so that Step 2's QC checks have something to flag. Run it with:

```bash
python generate_test_data.py
```

Skip this step once you're working with your own real sequencing data — real files just need to be placed at `data/raw_counts.csv` and `data/sample_info.csv` in the same format (see Data Format Requirements below).

### Step 2 — Data Preprocessing

`run_pipeline.py` chains together the full preprocessing sequence:

1. **Load & validate** (`step1_load_data.py`) — reads the raw counts and sample metadata, checks that every sample in the count matrix has matching metadata (and vice versa), and reorders metadata to align with the count matrix column order.
2. **Quality control** (`step2_quality_control.py`) — flags samples with low total library size, filters out genes with negligible counts across samples, and flags replicate pairs with low expression correlation (possible outliers or mislabeled samples).
3. **Normalization** (`step3_normalization.py`) — computes DESeq2 size factors (via `pydeseq2`) and normalizes counts to correct for differences in sequencing depth between samples.
4. **Log transform** (`step4_log_transform_and_explore.py`) — applies `log2(x + pseudocount)` to stabilize variance and bring expression values onto a more model-friendly scale.
5. **Export for dynGENIE3** (`step5_export_for_dyngenie3.py`) — reshapes the log-normalized data into one matrix per condition, averaging replicates at each timepoint, producing `outputs/05_dynGENIE3_input_<condition>_averaged.csv` files ready for network inference.

Run the whole sequence with:

```bash
python run_pipeline.py
```
**Output will be** - `TS_data_antibiotic.pkl` and `TS_data_control.pkl`



### Step 3 — Run the dynGENIE3 Pipeline

Each condition's exported CSV is converted into the `(TS_data, time_points, decay_rates, gene_names)` format dynGENIE3 expects, then run through the inference function **separately per condition** — dynGENIE3 infers one network from one time-series matrix, so `control` and `antibiotic` are processed independently:

```python
VIM, alphas, prediction_score, stability_score, treeEstimators = dynGENIE3(
    TS_data,
    time_points,
    alpha='from_data',   # estimate decay rates from the data itself
    gene_names=gene_names,
    tree_method='RF',
    ntrees=1000,
    nthreads=1
)
```

Each run produces a ranked list of regulatory links (`regulatory_network_ranking_<condition>.txt`) — the weighted edges of the inferred network. Comparing the two ranked lists (edges gained/lost, hub gene shifts) is how regulatory rewiring under antibiotic stress is detected.

## Key References

- Marku & Pancaldi (2023) — normalization and preprocessing rationale for dynamic GRN inference
- Huynh-Thu et al. — GENIE3
- Huynh-Thu & Geurts (2018) — dynGENIE3
- Ernst et al. — DREM
- Esquerre et al. (2014) — mRNA decay rate estimates
- Love et al. (2014) — DESeq2
- Robinson & Oshlack (2010) — TMM normalization
