# Bacterial-Gene-Network-inference
# Dynamic Transcriptomic States and GRN Rewiring in Antibiotic Resistance Evolution

Summer research project — IISER Bhopal
Supervisor: Dr. Garima Rani · Collaborator: Neetika

## Research Question

Can dynamic transcriptomic states and gene regulatory network (GRN) rewiring — rather than mutations alone — predict and redirect the evolution of antibiotic resistance in bacteria?

Most resistance-evolution studies focus on genetic mutations as the primary driver. This project instead asks whether **regulatory rewiring** (changes in which genes control which other genes, over time, under antibiotic stress) can itself act as a predictive signal or even a redirectable lever for resistance outcomes.

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

## Repository Structure

```
├── Data_Preprocessing.ipynb          # End-to-end preprocessing pipeline (Colab notebook)
├── config.py                         # Central config: paths, thresholds, normalization method
├── step1_load_data.py                # Load + validate raw counts and sample metadata
├── step2_quality_control.py          # Library size checks, low-count gene filtering, replicate correlation
├── step3_normalization.py            # DESeq2-based size-factor normalization (pydeseq2)
├── step4_log_transform_and_explore.py# log2(x + pseudocount) transform
├── step5_export_for_dyngenie3.py     # Reshape + average by condition/timepoint for dynGENIE3
├── run_pipeline.py                   # Runs steps 1–5 end to end
├── data/                             # User-supplied input (not tracked in repo)
│   ├── raw_counts.csv                # genes (rows) x samples (columns)
│   └── sample_info.csv               # sample_id, timepoint, condition, replicate
└── outputs/                          # Pipeline outputs (pickled intermediates + dynGENIE3-ready CSVs)
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

## Usage

```bash
# 1. Place your own raw_counts.csv and sample_info.csv inside data/
# 2. Run the full preprocessing pipeline
python run_pipeline.py

# 3. Feed the per-condition CSVs into dynGENIE3 separately
#    (outputs/05_dynGENIE3_input_control_averaged.csv)
#    (outputs/05_dynGENIE3_input_antibiotic_averaged.csv)
```

Each condition is run through dynGENIE3 **separately** — the tool infers one network per expression matrix, so control and antibiotic networks are generated independently and compared afterward.

## Status

- [x] Modular preprocessing pipeline built and tested (QC → normalization → log transform → export)
- [x] Bug fix: sample metadata realignment after `reset_index()` in `step1_load_data.py`
- [x] Verified pipeline runs end-to-end on real uploaded data (not just synthetic test data)
- [x] dynGENIE3 inference run per condition, producing ranked regulatory link lists
- [x] Preliminary network visualization (`networkx`) with hub-gene-weighted node sizing
- [ ] RegulonDB validation of inferred edges
- [ ] Hub gene centrality analysis (formal)
- [ ] DREM bifurcation analysis
- [ ] Control vs. antibiotic network comparison (shared-layout visualization)

## Key References

- Marku & Pancaldi (2023) — normalization and preprocessing rationale for dynamic GRN inference
- Huynh-Thu et al. — GENIE3
- Huynh-Thu & Geurts (2018) — dynGENIE3
- Ernst et al. — DREM
- Esquerre et al. (2014) — mRNA decay rate estimates
- Love et al. (2014) — DESeq2
- Robinson & Oshlack (2010) — TMM normalization

## Notes

This is an active work-in-progress research pipeline; interfaces (e.g. `config.py` parameters) may change as the reference-paper review for methodology justification is finalized.
