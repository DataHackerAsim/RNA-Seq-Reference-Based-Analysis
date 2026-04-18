# Methodology: Reference-Based RNA-Seq Analysis in Galaxy

## 1. Overview

This methodology describes a complete reference-based RNA-Seq workflow carried out on **Galaxy** using `usegalaxy.eu`. The pipeline covers data import, quality control, read trimming, genome alignment, quantification, differential expression analysis, visualization, and functional enrichment. The workflow is designed to be fully reproducible and to provide a transparent record of every major analysis decision.

The analysis is organized into two Galaxy histories:

- **History 1: RNA-Seq Reference-Based Analysis**  
  Used for data import, quality control, trimming, alignment, and read counting.
- **History 2: RNA-Seq DESeq2 Analysis**  
  Used for differential expression analysis, annotation, visualization, and enrichment.

Wherever possible, the same input datasets, annotations, and thresholds were used consistently across the workflow to ensure reproducibility and biological interpretability.

---

## 2. Scientific Rationale

RNA-Seq measures transcript abundance by sequencing cDNA derived from RNA molecules. Because the reads originate from processed transcripts rather than continuous genomic DNA, analysis requires tools that can:

1. assess read quality and remove technical artifacts;
2. align reads across splice junctions;
3. summarize read evidence at the gene level;
4. model count data statistically to identify differentially expressed genes;
5. interpret the resulting gene list in a biological context.

This methodology follows a standard reference-based RNA-Seq strategy and uses well-established tools in the Galaxy ecosystem to ensure both accessibility and rigor.

---

## 3. Input Data and Reference Resources

### 3.1 Sequencing data

Paired-end FASTQ files were used for the reference-based mapping workflow. To keep the workflow efficient while preserving the logic of the analysis, subsampled files were used for the alignment-based steps. Full-sized files are also available through Zenodo for the same dataset.

The following FASTQ files were imported for the mapping workflow:

- `GSM461177_1_subsampled.fastqsanger`
- `GSM461177_2_subsampled.fastqsanger`
- `GSM461180_1_subsampled.fastqsanger`
- `GSM461180_2_subsampled.fastqsanger`

For differential expression analysis, pre-computed gene count files corresponding to seven samples were used to represent the full biological comparison more comprehensively.

### 3.2 Reference genome and annotation

The analysis used the **Drosophila melanogaster** reference genome (**dm6**) together with the matching GTF annotation file:

- `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`

The annotation file was used consistently for genome-guided spliced alignment, gene-level counting, annotation of DESeq2 results, and downstream functional interpretation.

---

## 4. Workflow Summary

1. Import and verify FASTQ input data.
2. Group paired-end reads into collections.
3. Assess read quality with Falco and summarize with MultiQC.
4. Trim low-quality bases and filter short reads using Cutadapt.
5. Align reads to the dm6 genome with RNA STAR.
6. Count reads per gene using featureCounts.
7. Run DESeq2 on pre-computed count tables with a two-factor design.
8. Annotate differential expression results.
9. Filter significant genes using adjusted p-value and fold-change thresholds.
10. Visualize expression patterns with heatmaps.
11. Perform GO and KEGG enrichment analysis using goseq.

---

## 5. Galaxy Setup and Data Organization

### 5.1 Create a new history

A dedicated Galaxy history was created for the mapping workflow and named:

**RNA-Seq Reference-Based Analysis**

A second history was created for DESeq2-based statistical analysis:

**RNA-Seq DESeq2 Analysis**

Using separate histories helps keep the workflow organized, reduces the risk of accidental overwriting, and makes the analysis easier to review and reproduce.

### 5.2 Import data

FASTQ and annotation files were uploaded to Galaxy using the **Upload Data** function and the **Paste/Fetch Data** option. After upload, each dataset was checked to confirm that:

- the file completed successfully and turned green;
- the datatype was correctly assigned;
- paired files were correctly matched.

### 5.3 Datatype verification

Each FASTQ dataset was checked to ensure it was recognized as `fastqsanger`. If a file was incorrectly identified as generic `fastq`, the datatype was corrected manually. This step is important because downstream tools expect quality scores to be in the correct encoding.

### 5.4 Build paired collections

The paired-end FASTQ files were combined into a paired collection so that each forward read was matched with its corresponding reverse read. The resulting collection was named:

**2 PE fastqs**

This organization is essential for paired-end trimming and mapping because tools must process read pairs together to preserve read pairing information throughout the workflow.

---

## 6. Quality Control

### 6.1 Purpose of quality control

Raw sequencing reads can contain technical or biological artifacts that affect downstream analysis. Common issues include:

- low-quality base calls toward read ends;
- adapter contamination;
- uneven base composition;
- abnormal GC content;
- overrepresented sequences;
- unexpectedly short or truncated reads.

Performing quality control before trimming allows assessment of whether the data require preprocessing and whether the overall sequencing run appears technically sound.

### 6.2 Flatten paired collections

Some Galaxy reporting tools work more easily with simple dataset lists than with paired collections. Therefore, the paired collection was first flattened into a single list before running per-sample quality assessment.

### 6.3 Falco quality assessment

The flattened FASTQ datasets were analyzed using **Falco** to generate quality metrics for each file.

**Purpose:**  
To inspect the per-base and per-read quality profile of each sample.

**Typical outputs included:**

- HTML quality reports for visual inspection;
- raw report tables for summary aggregation.

### 6.4 MultiQC aggregation

The Falco outputs were summarized using **MultiQC** to produce a single consolidated report across all samples.

**Purpose:**  
To provide a sample-by-sample overview of sequencing quality in one report.

### 6.5 Quality interpretation criteria

The following metrics were inspected carefully:

| Metric | Expected interpretation | Concern if abnormal |
|---|---|---|
| Per-base quality | High quality across most positions | Need for trimming or data quality concerns |
| Adapter content | Minimal or absent | Adapter removal required |
| Per-sequence quality | Majority of reads at good quality | Presence of low-quality reads |
| GC content | Consistent with the organism and library | Possible contamination or bias |
| Sequence length distribution | Consistent with input data | Possible trimming or sequencing issues |

Quality control is not only a descriptive step; it informs the trimming strategy and helps validate that downstream results are built on reliable input reads.

---

## 7. Read Trimming and Filtering

### 7.1 Purpose of trimming

Trimming removes technical sequence content and low-quality bases that could interfere with mapping. In this workflow, trimming was used to:

- remove low-quality bases from read ends;
- discard reads that became too short after trimming;
- improve overall alignment accuracy;
- reduce the impact of sequencing artifacts.

### 7.2 Cutadapt trimming strategy

**Tool:** Cutadapt

**Approach:** Paired-end trimming using the paired collection

**Key settings used:**

- quality cutoff for R1: **20**
- minimum length: **20**
- report generation: enabled

### 7.3 Rationale for the trimming parameters

A quality cutoff of 20 corresponds to a commonly used conservative threshold for removing low-confidence bases. A minimum read length of 20 ensures that only reads with sufficient informative sequence are retained for reliable mapping.

For paired-end data, the paired read relationship must be preserved. Reads that become too short after trimming are removed together so that both members of the pair remain synchronized for downstream alignment.

### 7.4 MultiQC on trimming reports

The Cutadapt reports were summarized using MultiQC to assess:

- the number of reads trimmed;
- the amount of sequence removed from each read;
- whether trimming was extensive or only minor.

This step provides an audit trail showing how much the preprocessing changed the raw sequencing data.

---

## 8. Spliced Read Alignment

### 8.1 Why spliced alignment is required

RNA-Seq reads often span exon–exon junctions because mature mRNA has undergone splicing. Standard genomic aligners are not sufficient because they do not explicitly model introns. RNA STAR is a splice-aware aligner designed specifically for this purpose.

### 8.2 Genome and annotation setup

The dm6 genome reference was used together with the provided GTF annotation. The GTF file informs STAR where exons and gene structures are located, enabling more accurate splice-junction-aware mapping and downstream gene counting.

### 8.3 RNA STAR alignment

**Tool:** RNA STAR

**Input:** Paired-end trimmed reads from Cutadapt

**Reference configuration:** Built-in dm6 genome reference with external GTF annotation

**Important settings used:**

- paired-end mode: enabled
- gene annotation file: provided GTF annotation
- junction overhang: **36**
- gene counts output: enabled
- coverage output: enabled in bedgraph format

### 8.4 Rationale for the STAR parameter choice

The junction overhang was set to **36**, which matches the read length minus one for the dataset used. This allows STAR to detect splice junctions reliably while accommodating reads that span exon boundaries.

The use of gene-count output from STAR provides an additional read-based summary that can be used for inspection and quality assurance.

### 8.5 STAR output interpretation

The main alignment outputs included:

- BAM alignment files;
- alignment log files;
- gene-level read count summaries;
- splice junction outputs;
- coverage tracks in bedgraph format.

### 8.6 STAR quality assessment with MultiQC

STAR log files were summarized using MultiQC. The main alignment metrics reviewed were:

- percentage of uniquely mapped reads;
- percentage of multi-mapped reads;
- percentage of unmapped reads;
- overall alignment consistency across samples.

These metrics help determine whether the library quality, trimming strategy, and reference selection were appropriate.

---

## 9. Read Counting

### 9.1 Purpose of gene counting

After alignment, the next step is to convert mapped reads into a gene-by-sample count matrix. This matrix is the starting point for statistical analysis of differential expression.

### 9.2 featureCounts quantification

**Tool:** featureCounts

**Input:** STAR-aligned BAM files

**Annotation file:** dm6 GTF annotation

**Key settings used:**

- feature type: **exon**
- gene identifier: **gene_id**
- paired-end handling: count each read pair as one fragment
- minimum mapping quality: **10**
- strand specificity: **unstranded**
- gene-length output: enabled

### 9.3 Rationale for counting settings

Counting exon features is appropriate for gene-level RNA-Seq quantification because mature transcript abundance is reflected in exonic read coverage. Using `gene_id` groups exon counts into gene-level totals.

A mapping quality threshold of 10 excludes weak or ambiguous alignments. The library was treated as unstranded because the sequencing design did not require stranded counting.

### 9.4 featureCounts summary review

The summary output was reviewed to assess:

- the percentage of reads assigned to genes;
- the percentage of unassigned reads;
- whether the assignment pattern suggested mapping or annotation issues.

A high assignment rate generally indicates that the reference annotation and alignment strategy were appropriate.

---

## 10. Differential Expression Analysis

### 10.1 Why use DESeq2

Raw counts cannot be compared directly across samples because of differences in sequencing depth and library composition. DESeq2 addresses these issues by:

- estimating size factors for normalization;
- modeling count variability with negative binomial statistics;
- testing whether observed expression differences are likely to be real rather than due to sampling noise.

### 10.2 Input data

For the statistical analysis, seven pre-computed count files were imported into a separate Galaxy history. These corresponded to the available biological samples and provided a broader basis for differential expression testing than the smaller alignment demonstration set.

### 10.3 Model design

A two-factor DESeq2 design was used to account for both biological condition and sequencing type.

**Factor 1: Treatment**
- treated
- untreated

**Factor 2: Sequencing**
- paired-end
- single-end

Including sequencing type in the model helps prevent technical library design from confounding the biological treatment effect.

### 10.4 DESeq2 settings

**Tool:** DESeq2

**Input type:** Count data from featureCounts

**Important options enabled:**

- headers present in input files: yes
- beta priors: enabled
- diagnostic plots: enabled
- normalized counts output: enabled

### 10.5 Statistical output

DESeq2 generated:

- gene-level test statistics;
- raw and adjusted p-values;
- log2 fold changes;
- normalized counts;
- diagnostic plots such as PCA, MA plot, dispersion estimates, and sample clustering visualizations.

### 10.6 Interpretation criteria

Genes were considered differentially expressed when they met both of the following thresholds:

- adjusted p-value < **0.05**
- absolute log2 fold change > **1**

This combined criterion balances statistical confidence and biological magnitude of change.

### 10.7 Why adjusted p-values matter

Thousands of genes are tested simultaneously in RNA-Seq analysis. Adjusted p-values control the false discovery rate and reduce the likelihood of identifying genes as significant purely by chance.

### 10.8 Annotation of DESeq2 results

The DESeq2 output table was annotated using the GTF reference so that gene identifiers could be linked to gene symbols and genomic locations. This made the results easier to interpret, report, and visualize.

The annotated table included fields such as:

- gene identifier;
- base mean;
- log2 fold change;
- standard error;
- test statistic;
- p-value;
- adjusted p-value;
- chromosome;
- genomic start and end;
- strand;
- gene symbol.

---

## 11. Identification of Significant Differentially Expressed Genes

The annotated DESeq2 results were filtered in two stages.

### 11.1 Significance filter

The first filter retained genes with adjusted p-value below 0.05.

### 11.2 Effect-size filter

The second filter retained genes with absolute log2 fold change greater than 1.

### 11.3 Final DEG set

After applying both filters, the resulting gene set represented the most robust differential expression signal in the dataset. This set was used for downstream visualization and enrichment analysis.

Using both statistical significance and effect size is important because some genes may be statistically significant but show only small expression differences that are unlikely to be biologically meaningful.

---

## 12. Visualization of Differential Expression

### 12.1 Purpose of visualization

Visualization helps interpret the structure of the DEG set and compare expression patterns across samples. It also provides a useful check for sample grouping, clustering behavior, and consistency of the selected genes.

### 12.2 Extraction of normalized counts

Normalized counts for the significant DEGs were extracted from the DESeq2 output and joined with the filtered DEG list so that only genes meeting the significance criteria were retained for plotting.

### 12.3 Heatmap of log-transformed normalized counts

A heatmap was generated from the normalized count matrix after applying a `log2(value + 1)` transformation.

**Purpose:**  
To reduce the influence of very large count values and make patterns across genes and samples easier to visualize.

### 12.4 Heatmap with z-score scaling

A second heatmap was generated using row-wise z-score standardization.

**Purpose:**  
To compare relative expression changes within each gene, independent of the gene’s absolute expression level.

**Interpretation:**

- higher-than-average expression appears as positive z-scores;
- lower-than-average expression appears as negative z-scores;
- clustering groups genes and samples with similar expression profiles.

### 12.5 Why use both heatmap types

The log-transformed heatmap preserves approximate magnitude differences, while the z-score heatmap emphasizes relative expression patterns. Together, they provide complementary views of the same DEG set.

---

## 13. Functional Enrichment Analysis

### 13.1 Purpose of enrichment analysis

Once a biologically meaningful DEG set has been identified, enrichment analysis asks what kinds of functions, processes, and pathways are over-represented among those genes. This step moves the analysis from a gene list to a biological interpretation.

### 13.2 Input preparation for goseq

To perform enrichment analysis correctly, two types of inputs were prepared:

1. a list of genes and their differential expression status;
2. a gene length file.

Gene length correction is important because longer genes are more likely to accumulate reads and can therefore be overrepresented among DEGs for technical reasons if not accounted for.

### 13.3 GO enrichment

**Tool:** goseq

**Purpose:**  
To identify over-represented Gene Ontology categories among the DEGs.

The following GO domains were analyzed:

- **Biological Process (BP)**  
  Processes and pathways the gene contributes to.
- **Molecular Function (MF)**  
  Molecular activity performed by the gene product.
- **Cellular Component (CC)**  
  The cellular location where the gene product acts.

The output included ranked enriched categories and a plot of the most over-represented terms.

### 13.4 KEGG pathway enrichment

**Tool:** goseq

**Purpose:**  
To identify KEGG pathways that are over-represented among the DEGs.

This analysis helps place the results into broader cellular and metabolic context by showing whether the significant genes cluster in known biological pathways.

### 13.5 Interpretation of enrichment results

Enrichment results should be interpreted as hypotheses about biological themes associated with the experiment. A significant term does not prove causality, but it does highlight processes that are disproportionately represented in the DEG list and therefore worthy of further investigation.

---

## 14. Parameter Justification Summary

| Tool | Parameter | Value | Justification |
|---|---:|---|---|
| Cutadapt | Quality cutoff | 20 | Removes low-confidence bases |
| Cutadapt | Minimum length | 20 | Discards fragments too short for reliable mapping |
| RNA STAR | Junction overhang | 36 | Matches read length minus one |
| featureCounts | Feature type | exon | Appropriate for gene-level RNA-Seq quantification |
| featureCounts | Minimum MAPQ | 10 | Excludes weak alignments |
| featureCounts | Strand specificity | Unstranded | Matches the library design |
| DESeq2 | Adjusted p-value cutoff | 0.05 | Controls false discovery rate |
| DESeq2 | log2FC cutoff | 1 | Requires biologically meaningful change |

---

## 15. Quality Control and Interpretation Standards

To support a strong final report, the following standards were used throughout the analysis:

- input files were verified before analysis;
- quality reports were reviewed before trimming;
- trimming reports were checked to confirm preprocessing behaved as expected;
- alignment logs were summarized to assess mapping performance;
- count summaries were checked to confirm gene-level assignment;
- DESeq2 results were filtered using both significance and effect-size criteria;
- visualizations and enrichment analyses were based only on the final DEG set.

This layered validation reduces the risk of drawing conclusions from poor-quality input or weak statistical signal.

---

## 16. Reproducibility Notes

This workflow was designed for reproducibility and clarity. A reader should be able to reproduce the analysis by following the same steps in Galaxy using the same reference genome, annotation file, and parameter choices.

To ensure reproducibility, the following practices were followed:

- organized histories were used for different analysis stages;
- input files were clearly named;
- paired-end reads were handled as pairs throughout preprocessing and alignment;
- parameter choices were explicitly recorded;
- reference genome and annotation were kept consistent across steps;
- statistical thresholds were stated in advance and applied consistently.

---

## 17. Limitations

Although the workflow is robust and standard for reference-based RNA-Seq analysis, several limitations should be acknowledged:

- subsampled reads were used for the mapping demonstration, which may not reflect the full depth of the original experiment;
- count files used for DESeq2 were pre-computed, so the statistical analysis does not re-establish the alignment stage from the full raw FASTQ files;
- enrichment results depend on the completeness and correctness of the annotation database;
- biological interpretation is strongest when paired with experimental metadata and validation data.

These limitations do not undermine the workflow, but they should be acknowledged in a polished report.

---

## 18. Troubleshooting Guide

| Issue | Likely cause | Recommended action |
|---|---|---|
| File does not upload correctly | Browser or activation issue | Recheck account setup and upload status |
| Datatype appears wrong | Galaxy auto-detected incorrectly | Manually set datatype to `fastqsanger` |
| Paired files do not match | Collection built incorrectly | Rebuild the paired collection carefully |
| STAR alignment appears poor | Poor read quality or incorrect reference | Review QC and confirm genome/annotation |
| featureCounts assigns few reads | Annotation mismatch or wrong strand option | Check GTF and library orientation |
| DESeq2 output is confusing | Input factors not defined correctly | Verify the design matrix and sample groups |
| Enrichment output is sparse | DEG list too small or too strict | Review filtering thresholds carefully |

---

## 19. Concluding Statement

This RNA-Seq methodology provides a complete, transparent, and reproducible framework for transcriptomic analysis in Galaxy. By combining robust quality control, splice-aware alignment, gene-level quantification, statistically sound differential expression testing, and biologically informative enrichment analysis, the workflow produces results that are both technically defensible and scientifically meaningful.

---

## 20. Core Workflow Tools Used

- **Falco** — read quality assessment  
- **MultiQC** — report aggregation  
- **Cutadapt** — trimming and filtering  
- **RNA STAR** — spliced alignment  
- **featureCounts** — gene-level read counting  
- **DESeq2** — differential expression analysis  
- **Annotate DESeq2/DEXSeq output tables** — result annotation  
- **heatmap2** — visualization  
- **goseq** — GO and KEGG enrichment

---

*Prepared for Galaxy-based reference RNA-Seq analysis on Drosophila melanogaster (dm6).*
