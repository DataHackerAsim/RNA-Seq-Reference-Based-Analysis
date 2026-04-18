# Methodology — Complete Reproducible RNA-Seq Analysis Pipeline

This document provides a comprehensive, step-by-step methodology for performing reference-based RNA-Seq differential gene expression analysis using the Galaxy platform (usegalaxy.eu). The pipeline covers the entire analytical workflow from raw sequencing data through functional enrichment of differentially expressed genes. All tool versions, parameters, and their biological or statistical justifications are documented to ensure full reproducibility.

---

## Table of Contents

1. [Overview of the Analytical Framework](#1-overview-of-the-analytical-framework)
2. [Platform Setup and Data Acquisition](#2-platform-setup-and-data-acquisition)
3. [Quality Assessment of Raw Sequencing Reads](#3-quality-assessment-of-raw-sequencing-reads)
4. [Adapter and Quality Trimming](#4-adapter-and-quality-trimming)
5. [Splice-Aware Read Alignment](#5-splice-aware-read-alignment)
6. [Gene-Level Read Quantification](#6-gene-level-read-quantification)
7. [Differential Gene Expression Analysis](#7-differential-gene-expression-analysis)
8. [Visualization of Expression Patterns](#8-visualization-of-expression-patterns)
9. [Functional Enrichment Analysis](#9-functional-enrichment-analysis)
10. [Summary of Parameters and Thresholds](#10-summary-of-parameters-and-thresholds)
11. [Troubleshooting Reference](#11-troubleshooting-reference)

---

## 1. Overview of the Analytical Framework

### 1.1 Experimental Context

The dataset used in this analysis originates from a *Drosophila melanogaster* RNA-Seq experiment investigating the transcriptomic effects of a treatment condition. It comprises 7 biological samples in total: 4 untreated controls (GSM461176, GSM461177, GSM461178, GSM461182) and 3 treated samples (GSM461179, GSM461180, GSM461181). A subset of these samples were sequenced as paired-end (PE) reads while others were sequenced as single-end (SE) reads, which introduces a technical confounding variable that must be accounted for during statistical modeling.

### 1.2 Pipeline Architecture

The analysis follows a well-established reference-based RNA-Seq workflow that can be decomposed into four major phases:

**Phase I — Preprocessing:** Raw FASTQ reads undergo quality assessment (Falco/MultiQC), followed by adapter removal and quality trimming (Cutadapt). This phase ensures that only high-confidence bases enter the alignment stage.

**Phase II — Alignment and Quantification:** Trimmed reads are aligned to the *Drosophila melanogaster* reference genome (dm6 assembly) using the splice-aware aligner RNA STAR, which is specifically designed to handle the exon-exon junctions inherent in mRNA-derived reads. Aligned reads are then assigned to genomic features (genes) using featureCounts to produce a raw count matrix.

**Phase III — Statistical Inference:** The raw count matrix is subjected to differential expression testing using DESeq2, which applies variance-stabilizing normalization, dispersion estimation via an empirical Bayes shrinkage approach, and hypothesis testing via the Wald test. A multi-factor design matrix is employed to separate the biological effect of treatment from the technical effect of sequencing type.

**Phase IV — Biological Interpretation:** Statistically significant differentially expressed genes (DEGs) are visualized using hierarchical clustering heatmaps and subsequently subjected to Gene Ontology (GO) and KEGG pathway enrichment analysis using goseq, which corrects for gene length bias — a known confound in RNA-Seq enrichment studies.

### 1.3 Platform and Reproducibility

All analyses were conducted on the European Galaxy server (https://usegalaxy.eu), a freely accessible, browser-based bioinformatics platform that requires no local software installation. Galaxy provides full provenance tracking through its History system, ensuring that every tool invocation, parameter setting, and intermediate output is recorded automatically. Two separate Galaxy Histories were used to organize the workflow:

- **History 1** (`RNA-Seq Reference-Based Analysis`): Encompasses all steps from data import through read counting (Phases I and II).
- **History 2** (`RNA-Seq DESeq2 Analysis`): Encompasses differential expression analysis through functional enrichment (Phases III and IV), using pre-computed count files from all 7 samples to ensure adequate statistical power.

---

## 2. Platform Setup and Data Acquisition

### 2.1 History Initialization

A new Galaxy History was created and named `RNA-Seq Reference-Based Analysis` to serve as the organizational container for all upstream processing steps. In Galaxy, a History functions analogously to a project directory, storing all input datasets, tool configurations, and generated outputs in chronological order.

### 2.2 Sequencing Data Import

Subsampled paired-end FASTQ files (~5 MB each, reduced from ~1.5 GB originals) were imported directly from the Zenodo data repository. Subsampling accelerates the tutorial workflow without altering the analytical logic. The following four files were fetched via the Galaxy Upload utility using the Paste/Fetch Data option:

| File | Sample | Condition | Read Direction |
|------|--------|-----------|----------------|
| `GSM461177_1_subsampled.fastqsanger` | GSM461177 | Untreated | Forward (R1) |
| `GSM461177_2_subsampled.fastqsanger` | GSM461177 | Untreated | Reverse (R2) |
| `GSM461180_1_subsampled.fastqsanger` | GSM461180 | Treated | Forward (R1) |
| `GSM461180_2_subsampled.fastqsanger` | GSM461180 | Treated | Reverse (R2) |

**Source URLs (Zenodo Record 6457007):**

```
https://zenodo.org/record/6457007/files/GSM461177_1_subsampled.fastqsanger
https://zenodo.org/record/6457007/files/GSM461177_2_subsampled.fastqsanger
https://zenodo.org/record/6457007/files/GSM461180_1_subsampled.fastqsanger
https://zenodo.org/record/6457007/files/GSM461180_2_subsampled.fastqsanger
```

### 2.3 Datatype Verification

Each uploaded file was manually inspected to confirm that Galaxy assigned the correct datatype of `fastqsanger`. This distinction is critical because Galaxy uses datatype metadata to determine tool compatibility. The `fastqsanger` designation specifies that quality scores are encoded using the Sanger/Illumina 1.8+ Phred+33 encoding scheme, which is the standard for all modern Illumina sequencing platforms. Files incorrectly typed as generic `fastq` were manually corrected via the Edit Attributes panel (pencil icon → Datatype tab → set to `fastqsanger` → Save).

### 2.4 Paired-End Collection Construction

Forward and reverse reads for each sample were organized into a paired dataset collection using Galaxy's collection builder. This data structure preserves the mate-pair relationship between R1 and R2 files, ensuring that all downstream tools (Cutadapt, RNA STAR) process both mates together and maintain read pairing integrity throughout the pipeline.

**Procedure:**
1. All 4 FASTQ files were selected using the checkbox interface in the History panel.
2. The "Build List of Dataset Pairs" option was selected from the batch operations dropdown.
3. Galaxy auto-paired files based on the `_1` / `_2` suffix convention, producing two pairs: GSM461177 (untreated) and GSM461180 (treated).
4. The resulting collection was named `2 PE fastqs`.

---

## 3. Quality Assessment of Raw Sequencing Reads

### 3.1 Rationale

Quality control (QC) of raw sequencing data is an essential first step in any next-generation sequencing analysis pipeline. Illumina sequencing chemistry can introduce several types of artifacts, including progressive quality degradation toward the 3' end of reads (due to incomplete fluorescent signal termination), residual adapter sequences from library preparation, PCR duplicates causing overrepresented sequences, and GC content biases. Identifying these issues before alignment allows for informed parameterization of the trimming and filtering steps that follow.

### 3.2 Collection Flattening

Galaxy's MultiQC tool cannot directly ingest paired dataset collections, as it expects a flat (non-nested) list of inputs. The `Flatten collection` tool was therefore applied to the `2 PE fastqs` collection to convert the hierarchical paired structure into a simple linear list of 4 individual FASTQ files, while preserving all dataset metadata.

**Tool:** Flatten collection  
**Input:** `2 PE fastqs`

### 3.3 Per-File Quality Reporting with Falco

Falco (version 1.2.4+galaxy0), a high-performance drop-in replacement for the widely used FastQC tool, was run on the flattened collection to generate individual quality reports for each of the 4 FASTQ files. Falco produces identical output formats to FastQC but executes significantly faster due to its C++ implementation, making it preferable for large-scale or time-sensitive analyses.

**Tool:** Falco (version 1.2.4+galaxy0)  
**Input:** Output of Flatten collection (applied as a collection)  
**Outputs:** RawData (text-format statistics) and Webpage (interactive HTML reports)

### 3.4 Quality Metrics Evaluated

The following metrics were examined in each Falco report:

| Metric | Acceptable Range | Implication of Deviation |
|--------|-----------------|--------------------------|
| Per-base sequence quality | Median Phred ≥ 20 across all positions | Low-quality bases at read termini necessitate 3' trimming |
| Per-sequence quality scores | Unimodal peak at Phred ≥ 28 | Bimodal distribution suggests a subset of poor-quality reads |
| Per-base sequence content | Uniform proportions after position ~12 | Early-position bias is normal (random hexamer priming); late-position bias suggests contamination |
| GC content distribution | Matches theoretical species distribution | Deviation suggests contamination or extreme library bias |
| Sequence duplication levels | Low duplication expected for RNA-Seq | High duplication may indicate low library complexity or PCR over-amplification |
| Adapter content | Minimal or absent | Significant adapter presence requires explicit adapter trimming |
| Overrepresented sequences | None flagged | Flagged sequences may indicate adapter read-through or rRNA contamination |

### 3.5 Aggregated Quality Summary with MultiQC

All individual Falco reports were aggregated into a single interactive HTML dashboard using MultiQC (version 1.27+galaxy4). This enables simultaneous comparison of quality metrics across all samples, facilitating the identification of sample-specific anomalies or batch effects.

**Tool:** MultiQC (version 1.27+galaxy4)  
**Parameters:**
- Software generating logs: FastQC (Falco outputs are FastQC-compatible)
- FastQC output type: Raw data
- Input: Falco RawData collection

---

## 4. Adapter and Quality Trimming

### 4.1 Rationale

Based on QC assessment, reads require trimming to remove low-quality bases that could lead to misalignment or false variant calls, and to discard reads that become too short to map unambiguously after trimming. Cutadapt is a widely validated tool for this purpose and offers native support for paired-end data, which is critical because trimming must be coordinated across both mates — if one mate is discarded, the other must also be removed to maintain pairing.

### 4.2 Trimming Parameters and Justification

**Tool:** Cutadapt (version 5.2+galaxy0)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Mode | Paired-end Collection | Preserves mate-pair relationships; if one read is discarded, its mate is also removed |
| Input Collection | `2 PE fastqs` | Original paired collection (not the flattened version) |
| Quality cutoff (R1) | 20 | Removes trailing bases with Phred quality score < 20, corresponding to a base call error probability > 1%. This threshold is a standard community convention balancing data retention and accuracy |
| Minimum length (R1) | 20 | Discards reads shorter than 20 bp after trimming, as very short reads map ambiguously to multiple genomic loci and inflate multi-mapping rates |
| Additional outputs | Report enabled | Generates per-adapter statistics for downstream QC aggregation |

**Outputs:**
- Trimmed paired FASTQ reads (used as input for alignment)
- Trimming report (used for MultiQC aggregation)

### 4.3 Post-Trimming Quality Verification

A second MultiQC aggregation was performed on the Cutadapt report outputs to quantify the impact of trimming.

**Tool:** MultiQC (version 1.27+galaxy4)  
**Parameters:**
- Software generating logs: Cutadapt/Trim Galore!
- Input: Cutadapt Report collection

**Metrics assessed:**
- Total number and percentage of read pairs removed (expected to be low for high-quality libraries)
- Total base pairs trimmed from R1 vs. R2 (R2 typically shows more trimming due to the chemistry of Illumina sequencing, where the second read generally has lower quality)

---

## 5. Splice-Aware Read Alignment

### 5.1 Rationale for Splice-Aware Alignment

RNA-Seq reads originate from mature mRNA molecules, which have undergone post-transcriptional splicing to remove intronic sequences. When these reads are aligned back to the genomic reference (which retains introns), a substantial fraction of reads will span exon-exon junctions. Standard genomic aligners (e.g., BWA, Bowtie2) are not designed to handle such split alignments and would either fail to map junction-spanning reads or produce incorrect alignments. RNA STAR (Spliced Transcripts Alignment to a Reference) uses a two-pass alignment strategy that first identifies splice junctions from the data and then re-maps reads using these discovered junctions, achieving both high sensitivity and high speed.

```
Genome:     ═══EXON 1═══════[------INTRON------]═══════EXON 2═══
mRNA:       ═══EXON 1═══════════════════════════════════EXON 2═══
Read:                   ◄━━━━━━ read spans junction ━━━━━━►
```

### 5.2 Genome Annotation Import

The gene annotation file in GTF format was imported from Zenodo. This file defines the genomic coordinates of all known genes, transcripts, and exons in the *Drosophila melanogaster* BDGP6.32 genome assembly (UCSC dm6), and is required by both RNA STAR (to guide splice junction identification) and featureCounts (to assign reads to genes).

**Source:**
```
https://zenodo.org/record/6457007/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
```

Datatype was verified as `gtf` (not `gff` or `gff3`), as GTF is the required input format for RNA STAR.

### 5.3 RNA STAR Alignment Parameters

**Tool:** RNA STAR (version 2.7.11b+galaxy0)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Read type | Paired-end (as collection) | Maintains mate-pair information during alignment |
| Input reads | Cutadapt trimmed reads collection | Quality-filtered reads from the previous step |
| Genome source | Built-in index | Uses pre-indexed dm6 genome on the Galaxy server, avoiding the need for manual index construction |
| Genome annotation | Use genome reference without builtin gene-model but provide a GTF | Allows use of a specific annotation version (BDGP6.32.109) rather than the server's default annotation, ensuring consistency |
| Reference genome | Fly (*Drosophila melanogaster*): dm6 Full | UCSC dm6 assembly corresponding to Ensembl BDGP6 |
| Gene model GTF | `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz` | Matching annotation version |
| sjdbOverhang (junction overhang) | 36 | Calculated as read length minus 1 (37 − 1 = 36). This parameter specifies the length of the genomic sequence flanking each splice junction used for constructing the splice junction database. Setting it to the maximum read length minus one maximizes alignment sensitivity at junction boundaries |
| Per-gene output | GeneCounts | Produces a raw gene-level count matrix using the STAR internal counting algorithm, providing a secondary quantification to cross-validate against featureCounts |
| Coverage output | Bedgraph format | Generates strand-specific genome coverage tracks for visualization in genome browsers |

**Outputs:**
- `mapped.bam` — Primary alignment file containing all read-to-genome mappings with CIGAR strings encoding splice junctions
- `log` — Alignment summary statistics (mapping rates, multi-mapping rates, chimeric reads)
- `reads per gene` — STAR-internal gene counts
- `splice junctions.bed` — Coordinates of all detected splice junctions with supporting read counts
- Coverage bedgraph files — Per-strand genome coverage for visualization

### 5.4 Alignment Quality Assessment

MultiQC was used to aggregate RNA STAR log files and evaluate alignment performance across samples.

**Tool:** MultiQC (version 1.27+galaxy4)  
**Parameters:**
- Software generating logs: STAR
- Output type: Log
- Input: RNA STAR log collection

**Benchmark thresholds:**

| Metric | Expected Range | Interpretation |
|--------|---------------|----------------|
| Uniquely mapped reads | > 75% | Indicates high-quality library and appropriate reference genome |
| Multi-mapped reads | < 10% | Excessive multi-mapping suggests repetitive element contamination or inappropriate reference |
| Unmapped reads (too short) | < 5% | High values suggest over-aggressive trimming |
| Unmapped reads (other) | < 5% | High values suggest reference genome mismatch or contamination with non-target organisms |

---

## 6. Gene-Level Read Quantification

### 6.1 Rationale

Following alignment, each read has been assigned a genomic position. The quantification step counts how many reads overlap with each annotated gene, producing a raw count matrix where rows represent genes and columns represent samples. These raw counts serve as the input for statistical differential expression analysis. featureCounts was selected for this purpose due to its computational efficiency and well-documented behavior with paired-end RNA-Seq data.

### 6.2 featureCounts Parameters

**Tool:** featureCounts (version 2.1.1+galaxy0)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Alignment input | RNA STAR `mapped.bam` collection | BAM files from the alignment step |
| Strand specificity | Unstranded | The library preparation protocol did not preserve strand information; specifying this correctly prevents systematic undercounting that would occur if a stranded protocol were incorrectly assumed |
| Annotation source | GTF file in history | `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz` — same annotation used during alignment to ensure coordinate consistency |
| Feature type | exon | Counts reads overlapping exonic regions only, which is the standard approach for gene-level quantification in RNA-Seq. Intronic reads typically represent pre-mRNA or genomic DNA contamination |
| Gene identifier attribute | gene_id | The GTF attribute used to aggregate exon-level counts to the gene level |
| Output format | Gene-ID "\t" read-count | Produces a two-column tabular format directly compatible with DESeq2 and MultiQC |
| Gene-length file | Yes | Generates a file containing the effective length of each gene (sum of non-overlapping exonic bases), which is required by goseq for length-bias correction during enrichment analysis |
| Paired-end handling | Count fragments (pairs as 1) | Each mate-pair is counted once per fragment rather than twice per read, preventing artificial inflation of counts |
| Minimum mapping quality | 10 | Excludes reads with a mapping quality score below 10 (corresponding to > 10% probability of incorrect alignment). This removes ambiguously mapped reads that could introduce noise into the count matrix |

**Outputs:**
- Counts — Gene-level raw count matrix
- Summary — Read assignment statistics (assigned, unassigned due to ambiguity/no feature/low quality)
- Feature lengths — Gene length file for downstream enrichment analysis

### 6.3 Quantification Quality Assessment

MultiQC was applied to the featureCounts summary output to evaluate quantification efficiency.

**Tool:** MultiQC (version 1.27+galaxy4)  
**Parameters:**
- Software generating logs: featureCounts
- Input: featureCounts Summary collection

**Key metrics:**
- Percentage of reads successfully assigned to genes (target: > 60%)
- Percentage of unassigned reads categorized by reason (no feature, ambiguity, low mapping quality, etc.)

---

## 7. Differential Gene Expression Analysis

### 7.1 Rationale for Using DESeq2

Raw read counts cannot be directly compared between samples for differential expression due to several systematic biases: differences in sequencing depth (total library size) between samples, gene length effects (longer genes accumulate more reads), and RNA composition effects (a few highly expressed genes can skew normalization). DESeq2 addresses all of these through its median-of-ratios normalization method and employs a negative binomial generalized linear model (GLM) framework that is well-suited to the count-based, overdispersed nature of RNA-Seq data. Dispersion estimates are stabilized using an empirical Bayes shrinkage procedure, which is particularly valuable when sample sizes are small (as in this experiment), by borrowing information across genes with similar mean expression levels.

### 7.2 History Setup and Full Dataset Import

A second Galaxy History (`RNA-Seq DESeq2 Analysis`) was created for the differential expression and downstream analysis phases. Pre-computed featureCounts output files for all 7 samples were imported from Zenodo to ensure adequate statistical power. Using all 7 samples (rather than only the 2 processed from raw FASTQ in History 1) provides sufficient biological replication for robust statistical inference — DESeq2 requires a minimum of 2 replicates per condition, and additional replicates improve sensitivity for detecting true biological effects.

**Imported count files:**

| File | Sample | Condition | Sequencing |
|------|--------|-----------|------------|
| `GSM461176_untreat_single_featureCounts.counts` | GSM461176 | Untreated | SE |
| `GSM461177_untreat_paired_featureCounts.counts` | GSM461177 | Untreated | PE |
| `GSM461178_untreat_paired_featureCounts.counts` | GSM461178 | Untreated | PE |
| `GSM461179_treat_single_featureCounts.counts` | GSM461179 | Treated | SE |
| `GSM461180_treat_paired_featureCounts.counts` | GSM461180 | Treated | PE |
| `GSM461181_treat_paired_featureCounts.counts` | GSM461181 | Treated | PE |
| `GSM461182_untreat_single_featureCounts.counts` | GSM461182 | Untreated | SE |

### 7.3 DESeq2 Multi-Factor Model Design

A two-factor design matrix was specified to model both the biological variable of interest (Treatment) and the known technical confounding variable (Sequencing type). Including the sequencing type as a covariate allows DESeq2 to partition the variance attributable to SE vs. PE library preparation from the variance due to biological treatment, yielding more accurate estimates of the treatment effect.

**Tool:** DESeq2 (version 2.11.40.8+galaxy0)

**Factor 1 — Treatment (variable of interest):**

| Level | Samples |
|-------|---------|
| Treated | GSM461179, GSM461180, GSM461181 |
| Untreated | GSM461176, GSM461177, GSM461178, GSM461182 |

**Factor 2 — Sequencing (technical covariate):**

| Level | Samples |
|-------|---------|
| PE | GSM461177, GSM461178, GSM461180, GSM461181 |
| SE | GSM461176, GSM461179, GSM461182 |

**Additional parameters:**

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Files have header | Yes | Count files include a header row with column names |
| Input data type | Count data (featureCounts/HTSeq) | Raw integer counts, not pre-normalized values |
| Beta prior | Yes | Applies log-fold-change shrinkage to stabilize estimates for genes with low counts or high dispersion, reducing the influence of noisy estimates on downstream ranking |
| Generate visualizations | Yes | Produces diagnostic plots including PCA, sample distance heatmap, MA plot, and dispersion estimates |
| Output normalized counts | Yes | Provides size-factor-normalized count values for downstream visualization |

**Outputs:**
- DESeq2 result file — Tabular file containing, for each gene: base mean expression, log2 fold change, standard error, Wald test statistic, raw p-value, and Benjamini-Hochberg adjusted p-value (padj)
- Diagnostic plots — PCA plot (sample clustering), MA plot (fold change vs. mean expression), dispersion plot, and sample-to-sample distance heatmap
- Normalized counts — Size-factor-normalized expression values for all genes across all samples

### 7.4 Result Annotation

The DESeq2 output contains only Ensembl gene IDs without descriptive gene names or genomic coordinates. To facilitate biological interpretation, the Annotate DESeq2/DEXSeq tool was used to append gene symbols, chromosomal locations, and strand information from the GTF annotation file.

**Tool:** Annotate DESeq2/DEXSeq output tables (version 1.1.0+galaxy1)  
**Parameters:**
- Input: DESeq2 result file
- File type: DESeq2/edgeR/limma
- Reference annotation: `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`

A header row was then prepended using the Concatenate datasets tool to produce the final annotated result table with the following columns:

```
GeneID    Base mean    log2(FC)    StdErr    Wald-Stats    P-value    P-adj    Chromosome    Start    End    Strand    Feature    Gene symbol
```

### 7.5 Identification of Differentially Expressed Genes

Two sequential filtering steps were applied to identify genes exhibiting statistically significant and biologically meaningful expression changes.

**Filter 1 — Statistical significance:**

Tool: Filter data on any column using simple expressions  
Condition: `c7 < 0.05` (adjusted p-value < 0.05)  
Header lines skipped: 1

This filter retains only genes with a Benjamini-Hochberg adjusted p-value below 0.05, corresponding to a false discovery rate (FDR) of 5%. The BH correction is applied to control for the multiple testing problem inherent in simultaneously testing thousands of genes.

**Filter 2 — Biological significance (effect size):**

Tool: Filter data on any column using simple expressions  
Condition: `abs(c3) > 1` (absolute log2 fold change > 1)  
Header lines skipped: 1

This filter retains only genes with an absolute log2 fold change exceeding 1, which corresponds to a minimum 2-fold change in expression between treated and untreated conditions. This effect size threshold eliminates statistically significant but biologically trivial changes.

**Result:** 113 genes passed both filters and were designated as the final set of differentially expressed genes (DEGs).

---

## 8. Visualization of Expression Patterns

### 8.1 Extraction of Normalized Counts for DEGs

To visualize expression patterns of the 113 DEGs, their normalized count values were extracted from the full normalized count matrix using a table join operation.

**Step 1 — Join:**
Tool: Join two Datasets side by side on a specified field  
Parameters:
- Dataset 1: Normalized counts (DESeq2 output)
- Join column: Column 1 (GeneID)
- Dataset 2: Filtered DEG list (113 genes)
- Match column: Column 1 (GeneID)
- Keep non-joining lines: No
- Keep header: Yes

**Step 2 — Column selection:**
Tool: Cut columns from a table  
Parameters: Columns c1–c8 (GeneID + 7 sample columns), tab-delimited  
Output renamed: `Normalized counts for the most differentially expressed genes`

### 8.2 Heatmap of Log2-Transformed Normalized Counts

A heatmap was generated to visualize the absolute expression levels of all 113 DEGs across all 7 samples, with hierarchical clustering applied to both rows (genes) and columns (samples) to reveal co-expression patterns and sample groupings.

**Tool:** heatmap2 (version 3.2.0+galaxy1)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Input | Normalized counts for DEGs | 113 genes × 7 samples |
| Data transformation | Log2(value + 1) | Log transformation compresses the dynamic range of count data and approximates a normal distribution, making expression differences between lowly and highly expressed genes visually comparable. The pseudocount of 1 prevents undefined values for genes with zero counts |
| Clustering | Enabled | Hierarchical clustering (default: complete linkage with Euclidean distance) groups genes with similar expression profiles and samples with similar transcriptomic signatures |
| Labeling | Columns only | Gene names are omitted from row labels to prevent visual clutter with 113 genes |
| Color scheme | 2-color gradient | Provides intuitive low-to-high expression visualization |

### 8.3 Heatmap of Z-Score-Normalized Expression

A second heatmap was generated using Z-score normalization (row-wise standardization) to emphasize relative differences in expression across samples rather than absolute expression levels. This transformation is particularly useful for comparing genes that differ greatly in their basal expression levels on a common scale.

**Tool:** heatmap2 (version 3.2.0+galaxy1)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Input | Normalized counts for DEGs | Same dataset as above |
| Data transformation | None (plot as-is) | Z-scores are computed internally by the tool |
| Z-score computation | Compute on rows | For each gene (row), values are transformed as: Z = (x − μ) / σ, where μ is the gene's mean expression and σ is its standard deviation across samples. This centers each gene at zero and scales by variance |
| Clustering | Enabled | Groups genes and samples by relative expression similarity |
| Labeling | Columns only | Consistent with the previous heatmap |
| Color scheme | 3-color gradient | Diverging color scheme (e.g., green–black–red) visually distinguishes below-average (negative Z), average (zero), and above-average (positive Z) expression |

---

## 9. Functional Enrichment Analysis

### 9.1 Rationale

Identifying 113 DEGs provides a list of individual genes affected by treatment, but biological insight requires understanding the higher-order functional themes. Functional enrichment analysis tests whether specific biological categories (gene functions, pathways, cellular locations) are overrepresented among the DEGs compared to what would be expected by chance. The goseq package was selected because it corrects for gene length bias — a systematic artifact in RNA-Seq where longer genes are more likely to be detected as differentially expressed simply because they accumulate more reads, which inflates their statistical power.

### 9.2 Preparation of Input Files

#### File 1 — Gene-Level Differential Expression Status

A binary indicator file was created mapping each gene to its differential expression status (TRUE if padj < 0.05, FALSE otherwise). This file serves as the foreground/background partition for enrichment testing.

**Step 1 — Boolean computation:**
Tool: Compute on rows  
Input: DESeq2 result file  
Expression: `bool(float(c7) < 0.05)`  
Mode: Append new column  
Autodetect types: No  
Replacement for errors: False

**Step 2 — Column extraction:**
Tool: Cut columns  
Columns: c1, c8 (GeneID and boolean indicator)

**Step 3 — Case standardization:**
Tool: Change Case  
Column: c1 → Upper case  
Output renamed: `Gene IDs and differential expression`

Case conversion to uppercase ensures consistent matching with the goseq internal annotation database.

#### File 2 — Gene Length Information

Gene lengths (summed exonic length per gene) were obtained from the featureCounts `Feature lengths` output generated during Step 6. This file was similarly case-standardized.

Tool: Change Case  
Column: c1 → Upper case  
Output renamed: `Gene IDs and length`

### 9.3 Gene Ontology (GO) Enrichment Analysis

Gene Ontology provides a structured vocabulary for describing gene function across three independent hierarchies:

- **Biological Process (BP):** The larger biological program a gene contributes to (e.g., apoptotic process, immune response)
- **Molecular Function (MF):** The specific biochemical activity of the gene product (e.g., kinase activity, DNA binding)
- **Cellular Component (CC):** The subcellular location where the gene product functions (e.g., nucleus, plasma membrane)

**Tool:** goseq (version 1.50.0+galaxy0)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| DE gene file | Gene IDs and differential expression | Binary DE status for all tested genes |
| Gene length file | Gene IDs and length | Required for the Probability Weighting Function that corrects for length bias |
| Gene categories | Get categories (automatic retrieval) | Downloads GO annotations directly from the goseq database |
| Genome | dm6 (Fruit fly) | Matches the reference genome used throughout the analysis |
| Gene ID format | Ensembl Gene ID | Matches the identifier format in the DESeq2 output |
| Selected categories | GO: Biological Process, GO: Molecular Function, GO: Cellular Component | Comprehensive analysis across all three GO domains |
| Top GO terms plot | Yes | Visual summary of most significantly enriched terms |
| Extract DE genes per category | Yes | Lists which specific DEGs belong to each enriched term |

**Outputs:**
- Ranked GO term list with over-representation p-values (Wallenius approximation) and FDR-adjusted values
- Bar plot of top enriched GO terms
- Gene-category mapping showing which DEGs belong to each significantly enriched term

### 9.4 KEGG Pathway Enrichment Analysis

Complementing the GO analysis, KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway enrichment was performed to identify curated metabolic and signaling pathways that are overrepresented among the DEGs. KEGG provides a different organizational framework from GO, focusing on molecular interaction networks and reaction pathways rather than individual gene functions.

**Tool:** goseq (version 1.50.0+galaxy0)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| DE gene file | Gene IDs and differential expression | Same as GO analysis |
| Gene length file | Gene IDs and length | Same as GO analysis |
| Gene categories | Get categories | Automatic retrieval |
| Genome | dm6 | Consistent with upstream analyses |
| Gene ID format | Ensembl Gene ID | Consistent identifier format |
| Selected categories | KEGG only | Pathway-level analysis |
| Top GO terms plot | No | KEGG results are better interpreted through the KEGG pathway mapper |
| Extract DE genes per category | Yes | Identifies which DEGs contribute to each enriched pathway |

**Outputs:**
- Ranked KEGG pathway list with over-representation statistics
- Gene-pathway mapping for significantly enriched pathways

---

## 10. Summary of Parameters and Thresholds

The following table consolidates all critical parameters and their justifications across the entire pipeline:

| Tool | Parameter | Value | Biological/Statistical Rationale |
|------|-----------|-------|----------------------------------|
| Cutadapt | Quality cutoff | 20 | Phred 20 = 1% base call error rate; standard community threshold |
| Cutadapt | Minimum read length | 20 bp | Prevents ambiguous multi-mapping of ultrashort reads |
| RNA STAR | sjdbOverhang | 36 | Optimized as read length − 1 (37 − 1) for maximum junction sensitivity |
| RNA STAR | Genome | dm6 | UCSC assembly matching Ensembl BDGP6 for *D. melanogaster* |
| featureCounts | Strandedness | Unstranded | Matches the unstranded library preparation protocol |
| featureCounts | Feature type | exon | Standard for gene-level mRNA quantification |
| featureCounts | Min MAPQ | 10 | Filters ambiguously mapped reads (>10% error probability) |
| featureCounts | Paired counting | Fragment | Counts each insert once, not each mate separately |
| DESeq2 | Design formula | ~Sequencing + Treatment | Two-factor model separating technical and biological variance |
| DESeq2 | Beta prior | Enabled | Shrinks noisy fold-change estimates toward zero |
| DESeq2 | padj threshold | 0.05 | 5% FDR; standard for exploratory RNA-Seq studies |
| DESeq2 | \|log2FC\| threshold | > 1 | Minimum 2-fold change for biological relevance |
| goseq | Length correction | Enabled | Corrects for RNA-Seq gene length detection bias |
| goseq | Method | Wallenius | Non-central hypergeometric approximation accounting for length bias |

---

## 11. Troubleshooting Reference

| Issue | Likely Cause | Resolution |
|-------|-------------|------------|
| Galaxy account not activated | Confirmation email in spam | Check spam/junk folder and follow activation link |
| Files remain in upload queue | Server load or account not verified | Verify account activation; retry during off-peak hours |
| Datatype shown as `fastq` instead of `fastqsanger` | Auto-detection failure | Manually set datatype via Edit Attributes → Datatype → `fastqsanger` → Save |
| dm6 genome not visible in RNA STAR dropdown | List not fully loaded | Scroll through the complete dropdown or type "dm6" in the search field |
| Paired collection cannot be renamed | Galaxy UI limitation | Rename after collection is built, not during construction |
| featureCounts Feature lengths file absent | Option not enabled during run | Re-run featureCounts with "Create gene-length file" set to Yes |
| DESeq2 fails with factor level error | Sample misassignment between factor levels | Verify that each sample appears in exactly one level per factor |
| goseq gene ID mismatch | Case mismatch between files | Ensure all gene IDs are converted to uppercase using the Change Case tool |

---

## References

- Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15–21.
- Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.
- Liao, Y., Smyth, G.K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923–930.
- Young, M.D., Wakefield, M.J., Smyth, G.K., & Oshlack, A. (2010). Gene ontology analysis for RNA-seq: accounting for selection bias. *Genome Biology*, 11(2), R14.
- Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10–12.
- de Sena Brandine, G., & Smith, A.D. (2021). Falco: high-speed FastQC emulation for quality control of sequencing data. *F1000Research*, 8, 1874.
- Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048.
- The Galaxy Community. (2022). The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2022 update. *Nucleic Acids Research*, 50(W1), W13–W21.

---

*This methodology document accompanies the RNA-Seq reference-based analysis assignment. All analyses were performed on the European Galaxy server (https://usegalaxy.eu) following the Galaxy Training Network reference-based RNA-Seq tutorial.*
