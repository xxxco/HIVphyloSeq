# Sequence Selection

Sequences were obtained from the [Los Alamos National Laboratory (LANL) HIV Sequence Database](https://www.hiv.lanl.gov/). The filtering pipeline below reduced the full LANL dataset to the final analysis cohort of 584 full-length HIV-1 subtype B genomes.

## Filtering Steps

| Step | Criterion | Remaining (n) | Excluded (n) |
|------|-----------|---------------|--------------|
| Full LANL dataset | — | 2,727,128 | — |
| Subtype filter | HIV-1 subtype B only | 1,361,074 | 1,366,054 |
| Geographic filter | Germany or Belgium only | 60,943 | 1,300,131 |
| Temporal filter | Sampling year 2015–2019 | 8,441 | 52,502 |
| Quality control | Passed LANL QC annotations | 6,481 | 1,960 |
| Length filter | Sequence length > 7,000 bp | 1,495 | 4,986 |
| One per patient | Earliest sequence per patient retained; ties broken by lower accession number | 731 | 764 |
| Complete gag coverage | Sequences fully spanning the gag region | 598 | 133 |
| Additional QC | Excluded sequences with issues relative to HXB2 reference | 584 | 14 |

**Final analysis set: n = 584**

## Rationale for Criteria

- **Subtype B**: Restricts analysis to a single subtype to reduce confounding in phylogenetic inference; subtype B predominates in Germany and Belgium.
- **Germany / Belgium**: Small geographic area increases the likelihood of clustering among sequences.
- **2015–2019**: Avoids molecular-clock bias over longer periods; excludes the period when SARS-CoV-2 may have affected HIV diagnosis and sampling patterns.
- **LANL QC**: Removes sequences flagged for excessive ambiguous nucleotides, frameshifts, hypermutation, or poor HXB2 alignment — artifacts that distort phylogenetic inference and genetic distance estimation.
- **Length > 7,000 bp**: Retains near-complete genomes suitable for full-genome analysis.
- **One per patient**: Prevents over-representation of individuals with multiple longitudinal sequences.
- **Complete gag coverage**: Required for consistent subregion extraction across all downstream analyses.
- **CRF exclusion**: Recombinant strains (CRFs) were excluded to reduce confounding in clustering analyses.

## Subregions Analyzed

Seven commonly deposited HIV-1 genomic subregions were extracted from the final 584 sequences using their mapped HXB2 coordinates:

| Subregion | HXB2 Start | HXB2 End | Length (bp) | Final N |
|-----------|-----------|---------|-------------|---------|
| Protease | 2,253 | 2,549 | 297 | 584 |
| Reverse Transcriptase (RT) | 2,550 | 4,229 | 1,680 | 584 |
| Integrase | 4,230 | 5,096 | 867 | 584 |
| Vif | 5,041 | 5,619 | 579 | 584 |
| Nef | 8,797 | 9,417 | 621 | 579 |
| Env | 6,225 | 8,795 | 2,571 | 582 |
| Gag-Pol | 790 | 5,096 | 4,307 | 584 |

These regions correspond to genomic portions frequently sequenced for drug resistance testing, vaccine studies, and immune-escape analyses, and represent the most commonly deposited subregions in the LANL HIV Sequence Database.
