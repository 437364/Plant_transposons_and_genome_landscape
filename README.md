# Plant Transposons and Genome Landscape

Analysis code and data for the doctoral dissertation *"Plant Transposons and Genome Landscape"* (Marie Krátká, Masaryk University, Faculty of Science).

This repository contains the bioinformatics pipelines, statistical analyses, and visualization code used to investigate how centromere architecture shapes LTR retrotransposon dynamics across 68 plant species (37 holocentric, 31 monocentric).

## Overview

The central question is whether holocentric chromosomes — which lack localized centromeres and associated pericentromeric heterochromatin — expose transposable elements to more uniform recombination-based removal compared to monocentric chromosomes. We test this through:

1. **Comparative TE annotation** across 68 species using DANTE-LTR and RepeatMasker
2. **Spatial distribution analysis** of TE families along chromosomes using DTW clustering
3. **Statistical modeling** of solo LTR ratios, LTR identity, TE density, and Gypsy/Copia proportions using GLMMs
4. **Epigenetic landscape comparison** between *Luzula sylvatica* (holocentric) and *Juncus effusus* (monocentric)
5. **RT lineage phylogenetic analysis** comparing TE family structure between species (analysis by P. Jedlička)

## Repository Structure

```
├── data/                           # Metadata and small reference files
│   ├── metadata.tsv                # Species, accessions, centromere architecture, QC metrics
│   ├── accessions_genome.tsv       # NCBI FTP links for genome downloads
│   └── pruned_tree.tre             # Pruned angiosperm phylogeny
│
├── transposon_annotation/          # TE detection and annotation pipeline
│   ├── rules/                      # Snakemake rules (numbered by execution order)
│   ├── scripts/                    # Helper scripts called by the pipeline
│   └── notebooks/                  # Data aggregation, windowing, visualization
│
├── dtw_clustering/                 # DTW-based spatial distribution analysis
│   └── notebooks/                  # Clustering of TE density, LTR identity, and solo ratio profiles
│
├── glmm_analysis/                  # Generalized linear mixed models
│   └── glmm_analysis_v2.md         # R analysis (rendered from .Rmd) with full results
│
├── epigenetic_data/                # ChIP-seq and methylation analysis (Lsyl vs Jeff)
│   ├── notebooks/                  # Correlation matrices, density comparisons
│   └── scripts/                    # TE counting in genomic windows
│
├── rt_lineage_clustering/          # RT sequence clustering and phylogenetics
│   └── scripts/                    # CD-HIT, phylogenetic analysis, visualization
│
└── envs/                           # Conda environment specs and container references
    ├── *.yml                       # Conda environment exports
    ├── r_package_versions.csv      # R environment for GLMM analysis
    └── docker_urls                 # Docker container sources for Singularity
```

## Pipeline

### 1. Transposon Annotation (`transposon_annotation/`)

Snakemake-based pipeline run on Metacentrum HPC. Rules are numbered by execution order:

| Rule | Step |
|------|------|
| `01_download_genome.smk` | Download genomes from NCBI, extract chromosomes |
| `02_busco.smk` | Assembly quality assessment (BUSCO) |
| `03_gaps.smk` | Identify assembly gaps per chromosome |
| `04_dante_ltr.smk` | Detect LTR retrotransposons (DANTE + DANTE-LTR) |
| `05_prepare_library.smk` | Build per-species TE libraries, extract LTR coordinates |
| `06_find_elements.smk` | RepeatMasker annotation, filtering, solo LTR detection |
| `06.1_chunk_large_chromosomes.smk` | Split large chromosomes for RepeatMasker (>100 Mb) - alternative for failed snakemake jobs |
| `06.2_chunked_repeatmasker.smk` | Run RepeatMasker on chunks - alternative for failed snakemake jobs|

Key scripts:
- `find_LTR.py` — Extracts LTR boundary coordinates from DANTE-LTR output (replaces LTR_retriever's BLAST-based approach)
- `solo_finder.pl` — Identifies solo LTRs from RepeatMasker output (from LTR_retriever, Shujun Ou - https://github.com/oushujun/LTR_retriever)
- `filter_repeatmasker_records.py` — Filters RM records with inconsistent library sequence sizes
- `join_chunks.py` — Merges chunked RepeatMasker results with overlap resolution

Notebooks:
- `create_element_database.ipynb` — Aggregates all elements into master tables (`intact_elements_dante.tsv`, `partial_elements_dante.tsv`, `solo_elements.tsv`)
- `make_windows.ipynb` — Creates per-species sliding genomic windows (based on percentage or fixed size) with element counts, solo ratios, and LTR identity per window
- `prep_table_for_glmm.ipynb` — Joins element, chromosome, window, and species metadata into a single table for GLMM analysis
- `te_abundance_analysis.ipynb` — TE family frequency tables, LTR identity distributions
- `solo_plots.ipynb` — Solo ratio and LTR identity visualization across species
- `species_characteristics.ipynb` — Species overview plots (genome size, chromosome number, QC)
- `prune_tree.ipynb` — Prunes published angiosperm phylogeny to study species
- `pygenometracks.ipynb` — Prepares BigWig tracks for karyogram-style visualization
- `assembly_qc.ipynb` — Compiles BUSCO and gap statistics

### 2. DTW Clustering (`dtw_clustering/`)

Dynamic Time Warping analysis of TE spatial distributions along chromosomes, comparing holocentric and monocentric species.

- `dtw_clusters_full_length.ipynb` — Clustering of proportional TE family distributions (element counts per positional bin)
- `dtw_clusters_LTR_identity.ipynb` — Clustering of mean LTR identity profiles
- `dtw_clusters_solo_ratio.ipynb` — Clustering of solo ratio profiles

All notebooks use Ward linkage on DTW distance matrices, with chi-square tests for association between cluster membership and centromere architecture.

### 3. GLMM Analysis (`glmm_analysis/`)

Beta-family GLMMs testing the effects of centromere architecture and chromosome size on:
- Solo LTR ratio (proxy for NAHR-mediated TE removal)
- LTR identity (proxy for element age)
- TE density
- Gypsy proportion
- Ty3-gypsy CRM proportion

Models fitted with `glmmTMB` using random effects `(1 | Family/Genus/Species)` and heterogeneous dispersion `~Centromere_architecture`. Two nested models per response: additive, and interaction (architecture × log chromosome size)

**Input:** `data/df_for_linear_models.tsv` (generated by `prep_table_for_glmm.ipynb`; not included due to size — regenerate from the master element tables).

### 4. Epigenetic Data (`epigenetic_data/`)

Comparison of epigenetic landscapes between *L. sylvatica* (holocentric) and *J. effusus* (monocentric) using ChIP-seq (CenH3, H3K9me2, H3K4me3) and DNA methylation (CpG, CHG, CHH) data.

- `Lsyl_local_patterns.ipynb` / `Jeff_local_patterns.ipynb` — 100 kb window analysis, Spearman correlation matrices
- `element_density_Lsyl_Jeff.ipynb` — TE density comparisons (10 Mb windows), centromeric element proportions
- `epigenetic_workflow.md` — Complete documentation of data processing (ChIP-seq, Bismark methylation, Helixer gene annotation, TideCluster satellite annotation, deepTools metaplots)

See `epigenetic_workflow.md` for the full processing pipeline.

### 5. RT Lineage Clustering (`rt_lineage_clustering/`)

Phylogenetic analysis of reverse transcriptase sequences to compare TE family structure between *L. sylvatica* and *J. effusus*. Scripts developed by Pavel Jedlička.

Pipeline: RT FASTA → `cd_hit_rts_run_i99_80.sh` (CD-HIT clustering at multiple identity thresholds) → `TE_pipeline.py` (centroid extraction, MAFFT alignment, pairwise distances, heatmaps) → `TE_phylo_analysis.R` (FastTree phylogenies, MPD/MNTD with picante, permutation tests) → `summarize_TE_families.py` + `summary_fig.py` (cross-family summary and visualization).

## Data Availability

### Included in this repository
- `metadata.tsv` — Species metadata (68 species, centromere architecture, accessions, QC metrics)
- `accessions_genome.tsv` — NCBI FTP links for downloading all genome assemblies
- `pruned_tree.tre` — Pruned phylogenetic tree

### External data sources
- Genome assemblies: NCBI GenBank/RefSeq (accessions in `accessions_genome.tsv`)
- *L. sylvatica* ChIP-seq and EM-seq: Mata-Sucre, Krátká et al. (2024) *Nature Communications* ([doi:10.1038/s41467-024-53944-5](https://doi.org/10.1038/s41467-024-53944-5))
- *J. effusus* ChIP-seq: Dias et al. (2024) *The Plant Journal* ([doi:10.1111/tpj.16712](https://doi.org/10.1111/tpj.16712))
- Angiosperm phylogeny: source tree pruned to study species (see `prune_tree.ipynb`)

## Software and Environments

Conda environment specifications are in `envs/`. Singularity container sources are in `envs/docker_urls`.

| Tool | Environment / Container | Purpose |
|------|------------------------|---------|
| DANTE, DANTE-LTR | `dante_ltr.yml` | LTR retrotransposon detection |
| RepeatMasker | Metacentrum module | Similarity-based TE annotation |
| BUSCO | `busco.yml` | Assembly quality assessment |
| samtools, seqtk | Metacentrum modules | Sequence manipulation |
| Bismark 0.19.1 | `docker://argrosso/bismark:0.19.1` | Bisulfite alignment and methylation calling |
| MACS3 | `docker://parvsachdeva/macs3:latest` | ChIP-seq peak calling |
| Helixer v0.3.3 | `docker://gglyptodon/helixer-docker:helixer_v0.3.3_cuda_11.8.0-cudnn8` | Gene annotation |
| deepTools | `chipseq.yml` | ChIP-seq signal processing, metaplots |
| TideCluster | `tidecluster.yml` | Satellite repeat annotation |
| CD-HIT 4.8.1 | Metacentrum module | RT sequence clustering |
| MAFFT, FastTree, trimAl | Metacentrum modules | Alignment and phylogenetics |
| R (glmmTMB, DHARMa, picante, ape) | local R installation | Statistical modeling, phylogenetic analysis |
| Python (dtaidistance, seaborn, sklearn) | `python_tools.yml` / `merge_chunks.yml` | DTW clustering, visualization |
| pyGenomeTracks | `pygenometracks.yml` | Karyogram-style visualization |
| ete3 | `ete3.yml` | Phylogenetic tree manipulation |

## Computational Environment

All heavy computation (TE annotation, RepeatMasker, ChIP-seq processing) was performed on [MetaCentrum](https://metavo.metacentrum.cz/) (Czech national grid infrastructure). Downstream analyses (statistics, visualization) were run locally.