# HIVphyloSeq

All analyses used publicly available sequences downloaded on Feb. 24, 2025, from the [LANL HIV Sequence Database](https://www.hiv.lanl.gov/). See [`sequence_selection.md`](sequence_selection.md) for the full filtering workflow.

---

## Pipeline Overview

```
LANL sequences
    │
    ▼
Sequence filtering (→ 584 full-length subtype B genomes)
    │
    ▼
Extract 7 subregions (env, gag-pol, protease, RT, integrase, vif, nef)
    │
    ▼
Multiple sequence alignment (DECIPHER)
    │
    ▼
Maximum-likelihood phylogeny (RAxML-NG)
    │
    ├─► ClusterPicker ──┐
    │                   ├─► V-measure concordance (compare_cluster_full.R)
    └─► HIV-TRACE ──────┘               │
                                        ▼
                                vmeasure_results.tsv
```

---

## Software Requirements

| Tool | Version | Notes |
|------|---------|-------|
| R | ≥ 4.5 | |
| DECIPHER | Bioconductor | alignment |
| Biostrings | Bioconductor | FASTA I/O |
| clevr | CRAN | V-measure |
| dplyr, readr, stringr, tibble, jsonlite | CRAN | data handling |
| RAxML-NG | 1.2.2 | phylogeny |
| ClusterPicker | 1.2.5 | clustering |
| HIV-TRACE | — | clustering |

All tools were run inside a Singularity/Apptainer container:
```bash
singularity shell /path/to/multiple-alignment-container.sif
```

---

## Step 1 — Multiple Sequence Alignment (DECIPHER)

Alignments were generated using `DECIPHER::AlignSeqs()` with iterative refinement (5 iterations, 3 refinements, 20 processors). Run once per region FASTA file:

```r
library(DECIPHER)
library(Biostrings)

sequences <- readDNAStringSet("path/to/region_sequences.fasta")
aligned   <- AlignSeqs(sequences, iterations = 5, refinements = 3, processors = 20)
writeXStringSet(aligned, "DECIPHER_aligned_region_sequences.fasta")
```

---

## Step 2 — Maximum-Likelihood Phylogeny (RAxML-NG)

Run once per region. The example below is for `Env_CDS`; repeat for each of the 8 regions (7 subregions + full-length), adjusting `--msa` and `--prefix` accordingly.

```bash
raxml-ng --all \
    --msa    path/to/phy_file/DECIPHER_aligned_Env_CDS_sequences.phy \
    --model  GTR+G \
    --prefix path/to/tree/Env_CDS/bootstrap_1000/Env_CDS_mytree_bs \
    --seed   123 \
    --bs-trees 1000 \
    --threads 16
```

**Model**: GTR + Gamma-distributed rate heterogeneity  
**Bootstrap**: 1,000 nonparametric replicates; support values mapped onto the best-scoring ML tree (`.raxml.support` file used downstream)

---

## Step 3 — Transmission Cluster Identification

### ClusterPicker

Run once per region × bootstrap threshold × genetic distance threshold combination. Bootstrap thresholds evaluated: 70, 90, 99. Genetic distance thresholds evaluated: 0.5%–4.5% in 0.5% steps.

| Argument | Value | Description |
|----------|-------|-------------|
| 1 | `DECIPHER_aligned_Env_CDS_sequences.fasta` | DECIPHER-aligned FASTA file |
| 2 | `Env_CDS_mytree_bs.raxml.support` | RAxML-NG bootstrap support tree |
| 3 | `90` | Initial bootstrap support threshold (%) for cluster detection |
| 4 | `90` | Final bootstrap support threshold (%) — set equal to argument 3 |
| 5 | `0.03` | Maximum pairwise genetic distance within a cluster (3%) |
| 6 | `2` | Minimum cluster size (at least 2 sequences) |
| 7 | `gap` | Distance method: gap positions are ignored in distance calculation |

```bash
java -jar ClusterPicker_1.2.5.jar \
    DECIPHER_aligned_Env_CDS_sequences.fasta \
    Env_CDS_mytree_bs.raxml.support \
    90 90 \
    0.03 \
    2 \
    gap
```

Repeat for each region, adjusting the FASTA file, tree file, bootstrap threshold (90 → 70 or 99), and genetic distance (0.03 → 0.005–0.045) accordingly.

### HIV-TRACE

Run once per region × genetic distance threshold combination. Genetic distance thresholds evaluated: 0.5%–4.5% in 0.5% steps. Use the appropriate HXB2 reference per region (`HXB2_pol` for gag-pol, protease, RT, integrase, and full-length; `HXB2_env` for env; `HXB2_vif` for vif; `HXB2_nef` for nef).

| Flag | Value | Description |
|------|-------|-------------|
| `-i` | `DECIPHER_aligned_Env_CDS_sequences.fasta` | DECIPHER-aligned FASTA file |
| `-a` | `0.015` | Ambiguity threshold: sequences with >1.5% ambiguous bases are excluded |
| `-r` | `HXB2_env` | HXB2 reference name for the genomic region |
| `-t` | `0.03` | Genetic distance threshold for cluster membership (3%) |
| `-m` | `200` | Minimum nucleotide overlap required between two sequences |
| `-g` | `0.9` | Minimum overlap fraction (90% of the shorter sequence must overlap) |
| `-o` | `Env_CDS_hivtrace_GD030.json` | Output JSON file |

```bash
hivtrace \
    -i DECIPHER_aligned_Env_CDS_sequences.fasta \
    -a 0.015 \
    -r HXB2_env \
    -t 0.03 \
    -m 200 \
    -g 0.9 \
    -o Env_CDS_hivtrace_GD030.json
```

Repeat for each region, adjusting `-i`, `-r`, `-t`, and `-o` accordingly.

---

## Step 4 — V-measure Concordance Analysis

Compares both ClusterPicker and HIV-TRACE subregion clusters against the full-length genome reference (ClusterPicker, GD = 3.0%, BS = 90%) using the `clevr` R package. Produces a single output table covering all regions, methods, and thresholds.

```bash
Rscript compare_cluster_full.R regions.txt
```

Output: `vmeasure_results.tsv` — one row per region × method × GD threshold × bootstrap combination.

Edit the path variables at the top of `compare_cluster_full.R` before running.

---

## File Structure

```
HIVphyloSeq/
├── README.md                  # this file
├── sequence_selection.md      # LANL filtering criteria and cohort description
└── compare_cluster_full.R     # V-measure concordance analysis (clevr R package)
```
