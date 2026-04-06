# Epigenetic Data Workflow

Documentation of ChIP-seq, DNA methylation, gene/satellite annotation, and downstream analyses for *Luzula sylvatica* (Lsyl, holocentric) and *Juncus effusus* (Jeff, monocentric).

---

## 1. Data Sources

### Luzula sylvatica
- **ChIP-seq & methyl-seq:** Processed by collaborators for the Mata-Sucre, Krátká et al. (2024) *Nature Communications* paper.
  - `CpG_OT_Lsyl_methyl-seq_bismark_bt2_pe.deduplicated.bigwig`
  - `CHG_OT_Lsyl_methyl-seq_bismark_bt2_pe.deduplicated.bigwig`
  - `CHH_OT_Lsyl_methyl-seq_bismark_bt2_pe.deduplicated.bigwig`
- **Published pipeline:** Raw reads trimmed with cutadapt (Q > 20, adapters, 150 bp, single-end), mapped with bowtie2 (default), duplicates removed and filtered (`samtools view -q 10 -F 4`), ChIP vs input log2 ratio normalized by RPKM via `bamCompare` (deepTools).
- **Methylation:** Bismark pipeline, coverage files converted to BigWig.
- **Reference genome:** Chromosome-scale assembly from Mata-Sucre et al. (2024).
- **CenH3 peaks:** `CenH3_150k_merged_units.bed` — MACS3 broad peak calling (same pipeline as Jeff below), merged at 150 kb distance. Reused from the published analysis.

### Juncus effusus
- **ChIP-seq & methyl-seq:** Was previously published in Dias et al 2024, raw sequencing data received from André Marques. Processed on Metacentrum.
- **Reference genome:** GCA_027726005.1, unplaced scaffolds excluded
- **Marks used:** CenH3, H3K9me2, H3K4me3. H3K27me3 was also processed but excluded from the analysis because this mark is not available for Lsyl.

---

## 2. Juncus ChIP-seq Processing

**Location:** `data/chip-seq_juncus/`

### 2.1 Environment setup
```bash
mamba create -n chipseq -c bioconda -c conda-forge cutadapt bowtie2 samtools deeptools
```

### 2.2 Read trimming (cutadapt)
```bash
cutadapt -q 20 -a AGATCGGAAGAGC --length 150 -m 20 -j 8 \
  -o trimmed/${basename}_trimmed.fastq.gz $file
```
Parameters: quality threshold Q > 20, Illumina universal adapter, max length 150 bp, min length 20 bp.

### 2.3 Merge technical replicates
```bash
cat trimmed/5560_${sample}_run761_trimmed.fastq.gz \
    trimmed/5560_${sample}_run762_trimmed.fastq.gz \
    > trimmed_merged/5560_${sample}_merged.fastq.gz
```
Samples A–N merged across runs 761 and 762.

### 2.4 Alignment and filtering
```bash
bowtie2 -p 16 -q --very-sensitive-local \
  -x ${REF} -U $file \
  | samtools view -@ 4 -h -b -S - \
  | samtools sort -@ 4 -o aligned_corrected/${basename}.sorted.bam -
samtools index aligned_corrected/${basename}.sorted.bam
```
Reference index: `data/local_patterns/Jeff_ref/Juncus_effusus_bt2`

### 2.5 Individual coverage BigWigs
```bash
bamCoverage -b $f \
  -o bigwig/${basename}.log2.RPKM.cov.bw \
  -p 8 --ignoreDuplicates --normalizeUsing RPKM \
  --binSize 20 --smoothLength 60 --centerReads --extendReads 220
```

### 2.6 Merge biological replicates and compute ChIP/input ratio
```bash
samtools merge -@ 8 ${track}_merged.sorted.bam ${track}_rep1.sorted.bam ${track}_rep2.sorted.bam

bamCompare \
  -b1 ${track}_merged.sorted.bam \
  -b2 control_merged.sorted.bam \
  -o bigwig/${track}_merged-input.log2.RPKM.comp.bw \
  -p 8 --ignoreDuplicates --normalizeUsing RPKM \
  --binSize 20 --smoothLength 60 --centerReads --extendReads 220 \
  --scaleFactorsMethod None
```
Output: log2(ChIP/input) BigWig tracks for CenH3, H3K9me2, H3K4me3 (and H3K27me3, not used downstream).

### 2.7 Peak calling (CenH3)
```bash
# MACS3 broad peaks
macs3 callpeak -t cenh3_merged.sorted.bam -c control_merged.sorted.bam \
  -f BAM -g $GENOME_SIZE --broad --min-length 1000 --broad-cutoff 0.1 \
  -n cenh3_merged --outdir peaks/macs3

# EPIC2
epic2 -t cenh3_merged.sorted.bam -c control_merged.sorted.bam --chromsizes data/local_patterns/Jeff_ref/chrom.sizes -m 0 --output peaks/epic2/cenh3_merged_peaks.bed

# Consensus: intersect MACS3 and EPIC2, merge at 50 kb
bedtools merge -d 50000 -i peaks/macs3/cenh3_merged_peaks.broadPeak | bedtools sort > peaks/macs_merged.bed
bedtools merge -d 50000 -i peaks/epic2/cenh3_merged_peaks.bed | bedtools sort > peaks/epic_merged.bed
bedtools intersect -a peaks/macs_merged.bed -b peaks/epic_merged.bed -u > peaks/macs_in_epic.bed
bedtools intersect -a peaks/epic_merged.bed -b peaks/macs_merged.bed -u > peaks/epic_in_macs.bed
cat peaks/macs_in_epic.bed peaks/epic_in_macs.bed | bedtools sort | bedtools merge > peaks/consensus_peaks.bed
```

---

## 3. Juncus Methylation Processing

**Location:** `data/methylation_juncus/`

### 3.1 Bismark genome preparation
```bash
singularity exec bismark_0.19.1.sif bismark_genome_preparation data/local_patterns/Jeff_ref/
```

### 3.2 Alignment (paired-end)
```bash
singularity exec ${SINGULARITY_IMAGE} bismark \
  --genome ${GENOME_DIR} --bowtie2 --multicore 4 \
  --output_dir aligned --temp_dir tmp \
  -1 raw/run713_R1.fastq.gz -2 raw/run713_R2.fastq.gz
```
Two runs (713, 714) aligned separately.

### 3.3 Merge and deduplicate
```bash
samtools merge -@ 8 aligned/Juncus_merged.bam aligned/run713_*.bam aligned/run714_*.bam
samtools sort -n -@ 8 -o aligned/Juncus_merged.sorted.bam aligned/Juncus_merged.bam
singularity exec ${SINGULARITY_IMAGE} deduplicate_bismark \
  --paired --output_dir deduplicated --bam aligned/Juncus_merged.sorted.bam
```

### 3.4 Methylation extraction
```bash
singularity exec ${SINGULARITY_IMAGE} bismark_methylation_extractor \
  --paired-end --comprehensive --bedGraph --counts --multicore 4 \
  --buffer_size 10G --output methylation_calls \
  deduplicated/Juncus_merged.sorted.deduplicated.bam
```

### 3.5 Context-specific bedGraph and BigWig conversion
```bash
# For each context (CpG, CHG, CHH):
singularity exec ${SINGULARITY_IMAGE} bismark2bedGraph --CX --buffer_size 10G \
  -o ${context}_context.bedGraph ${context}_context_Juncus_merged.sorted.deduplicated.txt

# Sort and convert to BigWig
gunzip -c ${context}_context.bedGraph.gz | \
  awk 'NR==1 && /^track/ {print; next} {print | "sort -k1,1 -k2,2n"}' \
  > ${context}_context.sorted.bedGraph
bedGraphToBigWig ${context}_context.sorted.bedGraph chrom.sizes bigwigs/${context}_context.bw
```

---

## 4. Gene and Satellite Annotations

### 4.1 Gene annotation (Helixer)
```bash
singularity exec --nv $IMG Helixer.py \
  --fasta-path $INPUT --lineage land_plant \
  --gff-output-path data/local_patterns/{species}/{species}_helixer.gff3 \
  --model-filepath land_plant_v0.3_m_0200.h5 \
  --subsequence-length 213840
```
Run on GPU node. Model: `land_plant_v0.3_m_0200.h5`. Applied to both Lsyl and Jeff genomes.

### 4.2 Satellite annotation (TideCluster/TideHunter)
```bash
TideCluster.py run_all -f ${genome_fasta} -pr data/local_patterns/{species}/{species} -c $PBS_NUM_PPN
```
Applied to both Lsyl and Jeff genomes.

## 5. Centromeric vs Non-centromeric Element Classification

Elements classified by overlap with CenH3 peaks:
```bash
# Example for Lsyl solo LTRs:
bedtools intersect -a solo_LTRs.gff3 -b CenH3_150k_merged_units.bed -wa > solo_centromeric.gff3

# Full elements split similarly for both species:
bedtools intersect -a full_dante.gff3 -b ${peaks_bed} -wa > full_centromeric.gff3
bedtools intersect -a full_dante.gff3 -b ${peaks_bed} -v  > full_noncentromeric.gff3
```
Peak files used: `CenH3_150k_merged_units.bed` (Lsyl, from published analysis, MACS3+EPIC2 merged at 150 kb) and `consensus_peaks.bed` (Jeff, MACS3+EPIC2 consensus merged at 50 kb).

---

## 6. Window-based Correlation Analysis

**Notebooks:** `Lsyl_local_patterns.ipynb`, `Jeff_local_patterns.ipynb`

### 6.1 Create 100 kb windows
```bash
bedtools makewindows -g chrom.sizes -w 100000 \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3}' > bins.bed
```

### 6.2 TE elements → BED → BigWig → per-window summary
Within the notebooks: elements exported as center-point BED files, converted to BigWig via `bedtools genomecov -bga` + `bedGraphToBigWig`, then summarized per window.

### 6.3 Epigenetic marks → per-window summary
```bash
# ChIP-seq (Jeff example):
bigWigAverageOverBed data/chip-seq_juncus/bigwig/cenh3_merged-input.log2.RPKM.comp.bw bins.bed CenH3_average.tab
bigWigAverageOverBed data/chip-seq_juncus/bigwig/h3k9me2_merged-input.log2.RPKM.comp.bw bins.bed H3K9_average.tab
bigWigAverageOverBed data/chip-seq_juncus/bigwig/h3k4me3_merged-input.log2.RPKM.comp.bw bins.bed H3K4_average.tab

# Methylation (Jeff example):
bigWigAverageOverBed data/methylation_juncus/bigwigs/CpG_context.bw bins.bed CpG.tab
bigWigAverageOverBed data/methylation_juncus/bigwigs/CHG_context.bw bins.bed CHG.tab
bigWigAverageOverBed data/methylation_juncus/bigwigs/CHH_context.bw bins.bed CHH.tab
```
Same approach for Lsyl using published BigWig tracks.

### 6.4 Spearman correlation matrix
All `.tab` files merged into a single table per species (`merged_bin_info.tab`). Pairwise Spearman correlation with Benjamini-Hochberg FDR correction (threshold: p ≤ 0.01). Non-significant cells masked white.

**Features:** CenH3, H3K9me2, H3K4me3, CpG, CHG, CHH, full (intact element count), partial, solo, LTR identity (mean per window), gene density, satellite density.

**Output figures:** `{species}_spearman_heatmap.svg`, `{species}_correlation_heatmap_subset.svg`

---

## 7. Metaplots and Heatmaps (deepTools)

**Tool:** `computeMatrix scale-regions` + `plotProfile` / `plotHeatmap`

All metaplots use: `-m 4000 -b 2000 -a 2000` (4 kb body, 2 kb flanks).

### 7.1 Full/partial/solo element profiles
```bash
# Epigenetic marks over TE categories:
computeMatrix scale-regions \
  -S CenH3.bw H3K9.bw H3K4.bw \
  -R full_dante.bed partial_dante.bed solo_LTRs.bed \
  -o ${species}_epigenetic_marks.matrix.gz -m 4000 -b 2000 -a 2000

# Methylation over TE categories:
computeMatrix scale-regions \
  -S CpG.bw CHG.bw CHH.bw \
  -R full_dante.bed partial_dante.bed solo_LTRs.bed \
  -o ${species}_methylation_marks.matrix.gz -m 4000 -b 2000 -a 2000
```

### 7.2 Centromeric vs non-centromeric elements (Angela and Athila)
```bash
computeMatrix scale-regions \
  -S ${epigenetic_bws} \
  -R angela_centromeric.bed angela_noncentromeric.bed \
     athila_centromeric.bed athila_noncentromeric.bed \
  -o ${species}_centromeric_epigenetic_marks.matrix.gz -m 4000 -b 2000 -a 2000
```
Same for methylation tracks. Both species.

### 7.3 TE family-specific profiles
```bash
computeMatrix scale-regions \
  -S CpG.bw CHG.bw CHH.bw CenH3.bw H3K9.bw H3K4.bw \
  -R full_Angela.bed full_Athila.bed full_Tork.bed full_Ivana.bed full_random.bed \
  -o ${species}_dominant_families_full.matrix.gz -m 4000 -b 2000 -a 2000
```
Random control: `bedtools shuffle -i full_dante.bed -g genome.fai > full_random.bed`


### 7.4 Visualization parameters

**Profile plots:**
- Histone marks: `--yMax 1 1 1`
- Methylation: `--yMax 100 85 25`, `--yMin 60 20 0`
- Colors: `"#034bca" "#ff2dc2" "#ffa443"` (marks) or `"#005083" "#483cfc" "#eb24f4" "#ff4986"` (quartiles)

**Heatmaps:**
- Custom colormap: 16-step gradient from white through purple/olive to black
- Histone marks: `--zMin -2.5`, `--zMax 1`
- Methylation: context-specific z-ranges
- `--sortRegions descend --sortUsing mean`

---