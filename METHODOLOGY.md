# 📋 Methodology — Step-by-Step Analysis Guide

This document describes every step of the RNA-Seq analysis 
in detail, including all tool parameters. Anyone can 
reproduce this analysis by following these steps exactly 
on [Galaxy](https://usegalaxy.eu).

---

## 📖 Table of Contents
1. [Data Import & Setup](#step-1--data-import--setup)
2. [Quality Control](#step-2--quality-control)
3. [Read Trimming](#step-3--read-trimming)
4. [Read Mapping](#step-4--read-mapping-rna-star)
5. [Read Counting](#step-5--read-counting-featurecounts)
6. [Differential Expression](#step-6--differential-expression-deseq2)
7. [Visualization](#step-7--visualization)
8. [Functional Enrichment](#step-8--functional-enrichment)

---

## 🔰 Before You Start

### What is Galaxy?
Galaxy is a free, web-based platform for bioinformatics analysis.
You do not need to install any software. Everything runs in your 
browser. Go to [usegalaxy.eu](https://usegalaxy.eu) and register.

### What is a History?
In Galaxy, a **History** is like a project folder. It stores all 
your input files, tool outputs and results. We use 2 histories:
- **History 1**: `RNA-Seq Reference-Based Analysis` → QC to Counting
- **History 2**: `RNA-Seq DESeq2 Analysis` → DESeq2 to Enrichment

---

## STEP 1 — Data Import & Setup

### 1.1 Create New History
```
1. Log in to usegalaxy.eu
2. Click "+" icon at top right of History panel
3. Name it: RNA-Seq Reference-Based Analysis
4. Press Enter
```

### 1.2 Import FASTQ Files

> 💡 We use subsampled files (~5MB) for faster processing.
> Full files (~1.5GB each) are also available on Zenodo.
```
1. Click Upload button (top left panel)
2. Click "Paste/Fetch Data"
3. Paste these 4 URLs:

https://zenodo.org/record/6457007/files/GSM461177_1_subsampled.fastqsanger
https://zenodo.org/record/6457007/files/GSM461177_2_subsampled.fastqsanger
https://zenodo.org/record/6457007/files/GSM461180_1_subsampled.fastqsanger
https://zenodo.org/record/6457007/files/GSM461180_2_subsampled.fastqsanger

4. Click Start → Close
5. Wait for all 4 files to turn green
```

### 1.3 Check Datatype
```
For each file:
1. Click the file → click pencil icon (Edit attributes)
2. Click "Datatype" tab
3. Confirm it shows: fastqsanger
4. If it shows "fastq" → change to "fastqsanger" → Save
```

### 1.4 Create Paired Collection
```
Purpose: Group forward (_1) and reverse (_2) reads 
         for each sample into pairs

1. Click checkbox icon at top of History panel
2. Select all 4 FASTQ files
3. Click "4 of 4 selected" dropdown
4. Choose "Build List of Dataset Pairs"
5. Confirm auto-pairing:
   - GSM461177_1 ↔ GSM461177_2 (untreated pair)
   - GSM461180_1 ↔ GSM461180_2 (treated pair)
6. Collection name: 2 PE fastqs
7. Click Build
```

---

## STEP 2 — Quality Control

### Why do Quality Control?
Raw sequencing data can contain:
- **Low quality bases** → incorrect nucleotide calls
- **Adapter sequences** → artificial sequences from library prep
- **Overrepresented sequences** → contamination
These issues can affect downstream analysis if not addressed.

### 2.1 Flatten Collection
```
Purpose: MultiQC cannot process paired collections directly,
         so we first convert to a simple list

Tool: Flatten collection
Parameter:
  - Input Collection: 2 PE fastqs
Click: Run Tool
```

### 2.2 Run Falco (Quality Check)
```
Purpose: Generate quality report for each FASTQ file

Tool: Falco (version 1.2.4+galaxy0)
Parameter:
  - Raw read data: [collection icon] Output of Flatten collection
Click: Run Tool

Outputs:
  - RawData (txt files) → used for MultiQC input
  - Webpage (html files) → visual reports per sample
```

### 2.3 What to Look For in Falco Reports
| Metric | Good Range | Action if Bad |
|--------|-----------|---------------|
| Per base quality | Phred score > 20 | Trim low quality ends |
| Adapter content | None/minimal | Use Cutadapt to remove |
| Per sequence quality | Peak at high scores | Filter low quality reads |
| GC content | Matches expected | Investigate if abnormal |
| Read length | Consistent | Note for STAR parameter |

### 2.4 Run MultiQC (Combine Falco Reports)
```
Purpose: Combine all Falco reports into one summary report

Tool: MultiQC (version 1.27+galaxy4)
Parameters:
  - Which tool generated logs?: FastQC
    (Falco is a drop-in replacement for FastQC)
  - FastQC output type: Raw data
  - FastQC output: [collection] Falco on collection: RawData
Click: Run Tool

Output: One combined HTML report for all samples
```

---

## STEP 3 — Read Trimming

### Why Trim?
After quality control, we remove:
- **Low quality bases** at read ends (Phred < 20)
- **Short reads** that are too short to map reliably
- **Adapter sequences** if present

### 3.1 Run Cutadapt
```
Tool: Cutadapt (version 5.2+galaxy0)
Parameters:
  - Single or Paired-end?: Paired-end Collection
  - Paired Collection: 2 PE fastqs
  
  Under "Other Read Trimming Options":
  - Quality cutoff (R1): 20
  
  Under "Read Filtering Options":
  - Minimum length (R1): 20
  
  Under "Additional outputs":
  - ✅ Report (per-adapter statistics)
  
Click: Run Tool

Outputs:
  - Reads (trimmed FASTQ pairs) → used for mapping
  - Report (statistics) → used for MultiQC
```

> 💡 **Why run once for paired-end?**
> Cutadapt processes both reads in a pair together. If one 
> read becomes too short after trimming, it removes BOTH 
> reads to maintain pairing. This is important — unpaired 
> reads cause problems in downstream analysis.

### 3.2 Run MultiQC on Cutadapt Reports
```
Tool: MultiQC (version 1.27+galaxy4)
Parameters:
  - Which tool generated logs?: Cutadapt/Trim Galore!
  - Output of Cutadapt: [collection] Cutadapt on collection: Report
Click: Run Tool

Check in report:
  - How many reads were removed?
  - How many basepairs were trimmed from R1 vs R2?
```

---

## STEP 4 — Read Mapping (RNA STAR)

### Why is Mapping Special for RNA-Seq?
In RNA-Seq, reads come from **processed mRNA** (with introns removed).
When we map back to the genome (which has introns), the mapper must 
handle **splice junctions** — places where reads span two exons.
RNA STAR is a **splice-aware** mapper designed for this purpose.
```
Genome: ===EXON1=====[INTRON]=====EXON2===
mRNA:   ===EXON1========================EXON2===
Read:              ←read spans junction→
```

### 4.1 Import GTF Annotation File
```
Purpose: The GTF file tells STAR where genes and exons are 
         located in the genome

1. Click Upload → Paste/Fetch Data
2. Paste:
https://zenodo.org/record/6457007/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
3. Click Start → Close
4. Verify datatype is gtf or gtf.gz (not gff)
```

### 4.2 Run RNA STAR
```
Tool: RNA STAR (version 2.7.11b+galaxy0)
Parameters:
  - Single or paired-end?: Paired-end (as collection)
  - RNA-Seq paired reads: [collection] Cutadapt on collection: Reads
  - Custom or built-in genome?: Use a built-in index
  - Reference genome with/without annotation?: 
    use genome reference without builtin gene-model but provide a gtf
  - Select reference genome: Fly (Drosophila melanogaster): dm6 Full
  - Gene model (gtf) file: Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
  - Length of genomic sequence around junctions: 36
    (= read length - 1; our reads are 37bp)
  - Per gene output: Per gene read counts (GeneCounts)
  - Compute coverage: Yes in bedgraph format

Click: Run Tool

Outputs:
  - mapped.bam → aligned reads (main output)
  - log → alignment statistics
  - reads per gene → gene counts
  - splice junctions.bed → detected splice sites
  - Coverage files (bedgraph) → strand coverage
```

### 4.3 Run MultiQC on STAR Logs
```
Tool: MultiQC (version 1.27+galaxy4)
Parameters:
  - Which tool generated logs?: STAR
  - Type of STAR output?: Log
  - STAR log output: [collection] RNA STAR on collection: log
Click: Run Tool

What to check:
  - Uniquely mapped reads: should be >75%
  - Multi-mapped reads: should be <10%
  - Unmapped reads: should be <10%
```

---

## STEP 5 — Read Counting (featureCounts)

### Why Count Reads?
After mapping, we know WHERE each read came from in the genome.
Now we need to count HOW MANY reads came from each gene.
This count represents the **expression level** of each gene.

### 5.1 Run featureCounts
```
Tool: featureCounts (version 2.1.1+galaxy0)
Parameters:
  - Alignment file: [collection] RNA STAR on collection: mapped.bam
  - Specify strand information: Unstranded
  - Gene annotation file: A GFF/GTF file in your history
  - Gene annotation file: Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
  - GFF feature type filter: exon
  - GFF gene identifier: gene_id
  - Output format: Gene-ID "\t" read-count (MultiQC/DESeq2 compatible)
  - Create gene-length file: Yes
  - Does input have read pairs?: Yes, paired-end count as 1 fragment
  
  Under "Read filtering options":
  - Minimum mapping quality: 10

Click: Run Tool

Outputs:
  - Counts → gene expression count matrix
  - Summary → assignment statistics
  - Feature lengths → gene length file (needed for goseq)
```

### 5.2 Run MultiQC on featureCounts
```
Tool: MultiQC (version 1.27+galaxy4)
Parameters:
  - Which tool generated logs?: featureCounts
  - Output of FeatureCounts: [collection] featureCounts: Summary
Click: Run Tool

What to check:
  - % reads assigned to genes: ideally >60%
  - % unassigned reads: should be low
```

---

## STEP 6 — Differential Expression (DESeq2)

### Why Use DESeq2?
Raw read counts cannot be directly compared between samples because:
- Samples have different **sequencing depths** (total read numbers)
- Longer genes naturally get **more reads**
- There may be differences in **library composition**

DESeq2 handles all these normalizations and uses statistical 
models to find genes that are truly differentially expressed.

### 6.1 Create New History
```
1. Click "+" icon at top right of History panel
2. Name it: RNA-Seq DESeq2 Analysis
3. Press Enter
```

### 6.2 Import 7 Count Files
```
1. Click Upload → Paste/Fetch Data
2. Paste all 7 URLs:

https://zenodo.org/record/6457007/files/GSM461176_untreat_single_featureCounts.counts
https://zenodo.org/record/6457007/files/GSM461177_untreat_paired_featureCounts.counts
https://zenodo.org/record/6457007/files/GSM461178_untreat_paired_featureCounts.counts
https://zenodo.org/record/6457007/files/GSM461179_treat_single_featureCounts.counts
https://zenodo.org/record/6457007/files/GSM461180_treat_paired_featureCounts.counts
https://zenodo.org/record/6457007/files/GSM461181_treat_paired_featureCounts.counts
https://zenodo.org/record/6457007/files/GSM461182_untreat_single_featureCounts.counts

3. Click Start → Close
4. Wait for all 7 files to turn green
```

> 💡 **Why use pre-computed count files?**
> The full analysis uses 7 samples (3 treated + 4 untreated) for 
> statistical power. Processing all 7 from FASTQ would take too long, 
> so pre-computed counts are provided for the DESeq2 step.

### 6.3 Run DESeq2
```
Tool: DESeq2 (version 2.11.40.8+galaxy0)
Parameters:
  - how: Select datasets per level

  FACTOR 1:
  - Factor name: Treatment
  - Level 1 name: treated
  - Count files for treated: GSM461179, GSM461180, GSM461181
  - Level 2 name: untreated
  - Count files for untreated: GSM461176, GSM461177, GSM461178, GSM461182

  Click "+ Insert Factor" to add second factor:

  FACTOR 2:
  - Factor name: Sequencing
  - Level 1 name: PE
  - Count files for PE: GSM461177, GSM461178, GSM461180, GSM461181
  - Level 2 name: SE
  - Count files for SE: GSM461176, GSM461179, GSM461182

  Other settings:
  - Files have header?: Yes
  - Choice of input data: Count data (featureCounts/HTSeq)
  - Use beta priors: Yes

  Output options:
  - ✅ Generate plots for visualizing results
  - ✅ Output normalised counts

Click: Run Tool

Outputs:
  - Result file → DESeq2 statistics for all genes
  - Plots → PCA, heatmap, MA plot, dispersion estimates
  - Normalized counts → normalized expression values
```

> 💡 **Why 2 factors?**
> We include "Sequencing" as a second factor because some samples 
> are paired-end and others are single-end. Including this in the 
> model prevents it from confounding our treatment effect results.

### 6.4 Annotate DESeq2 Results
```
Purpose: Add gene names and locations to the results table

Step 1 - Import GTF:
  Upload: https://zenodo.org/record/6457007/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz

Step 2 - Run Annotate tool:
  Tool: Annotate DESeq2/DEXSeq output tables (version 1.1.0+galaxy1)
  Parameters:
    - Tabular output of DESeq2: DESeq2 result file
    - Input file type: DESeq2/edgeR/limma
    - Reference annotation: Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
  Click: Run Tool

Step 3 - Add column headers:
  Create a new tabular file with this header (tab-separated):
  GeneID  Base mean  log2(FC)  StdErr  Wald-Stats  P-value  P-adj  Chromosome  Start  End  Strand  Feature  Gene symbol

  Tool: Concatenate datasets
  Parameters:
    - First: header file
    - Second: Annotate output
  Rename output: Annotated DESeq2 results
```

### 6.5 Filter Differentially Expressed Genes
```
Step 1 - Filter by adjusted p-value:
  Tool: Filter data on any column using simple expressions
  Parameters:
    - Filter: Annotated DESeq2 results
    - Condition: c7<0.05
    - Header lines to skip: 1
  Rename: Genes with significant adj p-value

  → This keeps only statistically significant results
  → Adjusted p-value corrects for multiple testing

Step 2 - Filter by fold change:
  Tool: Filter data on any column using simple expressions
  Parameters:
    - Filter: Genes with significant adj p-value
    - Condition: abs(c3)>1
    - Header lines to skip: 1
  Rename: Genes with significant adj p-value & abs(log2(FC)) > 1

  → abs(log2FC) > 1 means fold change > 2 or < 0.5
  → This keeps only biologically meaningful changes

Result: 113 genes pass both filters
```

---

## STEP 7 — Visualization

### 7.1 Extract Normalized Counts for DEGs
```
Step 1 - Join datasets:
  Tool: Join two Datasets side by side on a specified field
  Parameters:
    - Join: Normalized counts file (from DESeq2)
    - Using column: Column 1
    - With: Genes with significant adj p-value & abs(log2FC) > 1
    - And column: Column 1
    - Keep non-joining lines: No
    - Keep header: Yes
  Click: Run Tool

Step 2 - Cut columns:
  Tool: Cut columns from a table
  Parameters:
    - Cut columns: c1-c8
    - Delimited by: Tab
    - From: output of Join tool
  Rename: Normalized counts for the most differentially expressed genes
```

### 7.2 Heatmap of Normalized Counts
```
Purpose: Visualize expression patterns of top DEGs across all samples

Tool: heatmap2 (version 3.2.0+galaxy1)
Parameters:
  - Input: Normalized counts for the most differentially expressed genes
  - Data transformation: Log2(value+1) transform my data
  - Enable data clustering: Yes
  - Labeling: Label columns and not rows
  - Colormap: Gradient with 2 colors
Click: Run Tool

Interpretation:
  - Rows = genes, Columns = samples
  - Colors show relative expression level
  - Clustering groups similar genes/samples together
```

### 7.3 Heatmap of Z-scores
```
Purpose: Show how far each value is from the gene's mean
         (removes differences in absolute expression levels)

Tool: heatmap2 (version 3.2.0+galaxy1)
Parameters:
  - Input: Normalized counts for the most differentially expressed genes
  - Data transformation: Plot the data as it is
  - Compute z-scores prior to clustering: Compute on rows
  - Enable data clustering: Yes
  - Labeling: Label columns and not rows
  - Colormap: Gradient with 3 colors
Click: Run Tool

Interpretation:
  - Red = higher than average expression
  - Green = lower than average expression
  - Z-score = (value - mean) / standard deviation
```

---

## STEP 8 — Functional Enrichment

### Why Enrichment Analysis?
After finding 113 DEGs, we want to understand:
- What **biological processes** do these genes control?
- Which **molecular pathways** are affected?
- Is there a pattern to the genes that change?

### 8.1 Prepare Input Files for goseq

#### File 1 — Gene IDs and Differential Expression
```
Step 1 - Compute boolean column:
  Tool: Compute on rows
  Parameters:
    - Input: DESeq2 result file
    - Expression: bool(float(c7)<0.05)
    - Mode: Append
    - Autodetect column types: No
    - Replacement value: False
  Click: Run Tool

Step 2 - Cut columns:
  Tool: Cut columns from a table
  Parameters:
    - Cut columns: c1,c8
    - From: output of Compute tool
  Click: Run Tool

Step 3 - Change case:
  Tool: Change Case
  Parameters:
    - From: output of Cut tool
    - Change case of columns: c1
    - To: Upper case
  Rename: Gene IDs and differential expression
```

#### File 2 — Gene Length File
```
Step 1 - Re-run featureCounts (in History 1) with:
         "Create gene-length file: Yes"
         to generate Feature lengths output

Step 2 - Change case:
  Tool: Change Case
  Parameters:
    - From: featureCounts Feature lengths
    - Change case of columns: c1
    - To: Upper case
  Rename: Gene IDs and length
```

### 8.2 GO Enrichment Analysis
```
Purpose: Find Gene Ontology terms enriched in our DEGs
         GO categories:
         - BP = Biological Process (what process the gene is involved in)
         - MF = Molecular Function (what the gene does at molecular level)
         - CC = Cellular Component (where the gene acts in the cell)

Tool: goseq (version 1.50.0+galaxy0)
Parameters:
  - Differentially expressed genes: Gene IDs and differential expression
  - Gene lengths file: Gene IDs and length
  - Gene categories: Get categories
  - Genome: Fruit fly (dm6)
  - Gene ID format: Ensembl Gene ID
  - Categories: 
    ✅ GO: Cellular Component
    ✅ GO: Biological Process
    ✅ GO: Molecular Function
  
  Output Options:
  - Output Top GO terms plot?: Yes
  - Extract DE genes for categories?: Yes

Click: Run Tool

Outputs:
  - Ranked category list → all GO terms with statistics
  - Top over-represented GO terms plot → visual summary
  - DE genes for categories → which DEGs are in each GO term
```

### 8.3 KEGG Pathway Analysis
```
Purpose: Find KEGG metabolic/signaling pathways enriched in our DEGs
         KEGG maps molecular interactions and reactions in cells

Tool: goseq (version 1.50.0+galaxy0)
Parameters:
  - Differentially expressed genes: Gene IDs and differential expression
  - Gene lengths file: Gene IDs and length
  - Gene categories: Get categories
  - Genome: Fruit fly (dm6)
  - Gene ID format: Ensembl Gene ID
  - Categories: ✅ KEGG only

  Output Options:
  - Output Top GO terms plot?: No
  - Extract DE genes for categories?: Yes

Click: Run Tool

Outputs:
  - Ranked category list → KEGG pathways with statistics
  - DE genes for categories → which DEGs are in each pathway
```

---

## 📊 Parameters Summary Table

| Tool | Key Parameter | Value | Reason |
|------|--------------|-------|--------|
| Cutadapt | Quality cutoff | 20 | Remove Phred<20 bases |
| Cutadapt | Min length | 20 | Remove very short reads |
| RNA STAR | Junction overhang | 36 | Read length (37) - 1 |
| featureCounts | Strand | Unstranded | Library type |
| featureCounts | Feature type | exon | Count exonic reads |
| featureCounts | Min mapq | 10 | Remove poor alignments |
| DESeq2 | Beta prior | Yes | Stabilize fold changes |
| DESeq2 | padj cutoff | 0.05 | 5% false discovery rate |
| DESeq2 | log2FC cutoff | 1 | Fold change > 2 |

---

## ⚠️ Common Issues & Solutions

| Issue | Solution |
|-------|----------|
| Account not activated | Check email spam folder |
| Files not uploading | Verify account activation |
| Wrong datatype | Change to fastqsanger manually |
| dm6 genome not found | Scroll down in dropdown list |
| Collection not editable | Rename after building |
| Feature lengths missing | Re-run featureCounts with "Create gene-length file: Yes" |

---

*For questions about this analysis, refer to the 
[Galaxy Training Network](https://training.galaxyproject.org) 
tutorial on Reference-based RNA-Seq data analysis.*
