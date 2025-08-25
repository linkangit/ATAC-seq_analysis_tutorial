# Complete ATAC-seq Analysis and Visualization Tutorial

## Table of Contents
1. [Introduction to ATAC-seq](#introduction)
2. [Prerequisites and Setup](#prerequisites)
3. [Data Quality Control](#quality-control)
4. [Read Alignment](#alignment)
5. [Post-alignment Processing](#post-alignment)
6. [Peak Calling](#peak-calling)
7. [Quality Assessment](#quality-assessment)
8. [Differential Accessibility Analysis](#differential-analysis)
9. [Visualization](#visualization)
10. [Advanced Analysis](#advanced-analysis)
11. [Troubleshooting](#troubleshooting)

## 1. Introduction to ATAC-seq {#introduction}

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a method for studying chromatin accessibility genome-wide. The technique uses a hyperactive Tn5 transposase to insert sequencing adapters into accessible chromatin regions, allowing identification of open chromatin and transcription factor binding sites.

### Key Concepts:
- **Accessible chromatin**: Regions where DNA is not tightly wrapped around histones
- **Tn5 transposase**: Enzyme that fragments DNA and inserts adapters simultaneously
- **Insert size distribution**: Characteristic pattern reflecting nucleosome positioning
- **TSS enrichment**: Measure of signal quality around transcription start sites

## 2. Prerequisites and Setup {#prerequisites}

### Required Software
```bash
# Install conda/mamba package manager first
# Then install required tools:

# Create ATAC-seq environment
conda create -n atacseq -c bioconda -c conda-forge \
    fastqc multiqc trimmomatic \
    bowtie2 samtools bedtools \
    macs3 deeptools \
    r-base r-deseq2 r-diffbind \
    python=3.9

conda activate atacseq

# Additional R packages (run in R)
install.packages(c("ggplot2", "dplyr", "GenomicRanges", 
                   "ChIPseeker", "clusterProfiler"))
```

### Directory Structure
```bash
mkdir atac_analysis
cd atac_analysis
mkdir -p {raw_data,fastqc,trimmed,aligned,peaks,figures,results}
```

### Reference Genome Setup
```bash
# Download reference genome (example: human hg38)
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Build bowtie2 index
bowtie2-build hg38.fa hg38_index

# Download blacklist regions
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
```

## 3. Data Quality Control {#quality-control}

### Initial FastQC Analysis
```bash
# Run FastQC on raw reads
fastqc raw_data/*.fastq.gz -o fastqc/

# Aggregate results with MultiQC
multiqc fastqc/ -o fastqc/
```

### Adapter Trimming
```bash
# Trim adapters using Trimmomatic
for sample in raw_data/*_R1.fastq.gz; do
    base=$(basename $sample _R1.fastq.gz)
    
    trimmomatic PE -phred33 \
        raw_data/${base}_R1.fastq.gz raw_data/${base}_R2.fastq.gz \
        trimmed/${base}_R1_paired.fastq.gz trimmed/${base}_R1_unpaired.fastq.gz \
        trimmed/${base}_R2_paired.fastq.gz trimmed/${base}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Quality control after trimming
fastqc trimmed/*_paired.fastq.gz -o fastqc/
```

## 4. Read Alignment {#alignment}

### Bowtie2 Alignment
```bash
# Align reads with bowtie2
for sample in trimmed/*_R1_paired.fastq.gz; do
    base=$(basename $sample _R1_paired.fastq.gz)
    
    bowtie2 -x hg38_index \
        -1 trimmed/${base}_R1_paired.fastq.gz \
        -2 trimmed/${base}_R2_paired.fastq.gz \
        -p 8 --very-sensitive --no-discordant --no-mixed \
        | samtools sort -@ 8 -o aligned/${base}.bam
    
    # Index BAM files
    samtools index aligned/${base}.bam
done
```

### Alignment Statistics
```bash
# Generate alignment statistics
for bam in aligned/*.bam; do
    base=$(basename $bam .bam)
    samtools flagstat $bam > aligned/${base}_flagstat.txt
done

# Aggregate with MultiQC
multiqc aligned/ -o aligned/
```

## 5. Post-alignment Processing {#post-alignment}

### Remove Low-Quality and Mitochondrial Reads
```bash
for bam in aligned/*.bam; do
    base=$(basename $bam .bam)
    
    # Filter reads: remove unmapped, low quality, mitochondrial
    samtools view -@ 8 -h -f 2 -F 1804 -q 30 $bam | \
    grep -v chrM | \
    samtools sort -@ 8 -o aligned/${base}_filtered.bam
    
    # Index filtered BAM
    samtools index aligned/${base}_filtered.bam
done
```

### Remove PCR Duplicates
```bash
# Using sambamba or Picard
for bam in aligned/*_filtered.bam; do
    base=$(basename $bam _filtered.bam)
    
    # Mark and remove duplicates
    samtools markdup -r -@ 8 $bam aligned/${base}_final.bam
    samtools index aligned/${base}_final.bam
done
```

### Shift Reads for Tn5 Insertion Sites
```bash
# Create script to shift reads
cat > shift_reads.py << 'EOF'
#!/usr/bin/env python3
import pysam
import sys

def shift_read(read):
    if read.is_reverse:
        read.reference_start = read.reference_end - 1
    else:
        read.reference_start = read.reference_start + 4
    return read

input_bam = sys.argv[1]
output_bam = sys.argv[2]

with pysam.AlignmentFile(input_bam, "rb") as infile:
    with pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
        for read in infile:
            if not read.is_unmapped:
                shifted_read = shift_read(read)
                outfile.write(shifted_read)
EOF

# Apply shifting
for bam in aligned/*_final.bam; do
    base=$(basename $bam _final.bam)
    python shift_reads.py $bam aligned/${base}_shifted.bam
    samtools index aligned/${base}_shifted.bam
done
```

## 6. Peak Calling {#peak-calling}

### MACS3 Peak Calling
```bash
# Call peaks using MACS3
for bam in aligned/*_shifted.bam; do
    base=$(basename $bam _shifted.bam)
    
    macs3 callpeak -t $bam \
        -f BAMPE -n $base \
        -g hs -p 0.01 \
        --shift -75 --extsize 150 \
        --nomodel --call-summits \
        --outdir peaks/
done
```

### Filter Peaks
```bash
# Remove peaks in blacklist regions
for peak in peaks/*.narrowPeak; do
    base=$(basename $peak .narrowPeak)
    
    bedtools intersect -v \
        -a $peak \
        -b hg38-blacklist.v2.bed \
        > peaks/${base}_filtered.narrowPeak
done
```

### Create Consensus Peak Set
```bash
# Merge peaks from all samples
cat peaks/*_filtered.narrowPeak | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > peaks/consensus_peaks.bed

# Count reads in consensus peaks for each sample
for bam in aligned/*_shifted.bam; do
    base=$(basename $bam _shifted.bam)
    bedtools coverage -a peaks/consensus_peaks.bed -b $bam \
        > peaks/${base}_coverage.txt
done
```

## 7. Quality Assessment {#quality-assessment}

### TSS Enrichment Analysis
```bash
# Download TSS regions
wget -O tss.bed "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refTSS.txt.gz"
gunzip tss.bed.gz

# Calculate TSS enrichment using deepTools
computeMatrix reference-point \
    --referencePoint TSS \
    -b 2000 -a 2000 \
    -R tss.bed \
    -S aligned/*_shifted.bam \
    --skipZeros \
    -o matrix_tss.gz

plotProfile -m matrix_tss.gz \
    -out figures/tss_profile.pdf \
    --numPlotsPerRow 2
```

### Insert Size Distribution
```bash
# Calculate insert size distribution
for bam in aligned/*_final.bam; do
    base=$(basename $bam _final.bam)
    
    samtools view -f 2 $bam | \
    awk '{print sqrt($9^2)}' | \
    sort -n | uniq -c > aligned/${base}_insert_sizes.txt
done

# Plot insert size distribution (R script)
cat > plot_insert_sizes.R << 'EOF'
library(ggplot2)
library(dplyr)

files <- list.files("aligned", pattern = "*insert_sizes.txt", full.names = TRUE)
all_data <- data.frame()

for (file in files) {
    sample_name <- gsub("_insert_sizes.txt", "", basename(file))
    data <- read.table(file, col.names = c("count", "size"))
    data$sample <- sample_name
    all_data <- rbind(all_data, data)
}

p <- ggplot(all_data, aes(x = size, y = count, color = sample)) +
    geom_line() +
    xlim(0, 1000) +
    labs(title = "Insert Size Distribution",
         x = "Insert Size (bp)",
         y = "Count") +
    theme_minimal()

ggsave("figures/insert_size_distribution.pdf", p, width = 10, height = 6)
EOF

Rscript plot_insert_sizes.R
```

### FRiP Score (Fraction of Reads in Peaks)
```bash
# Calculate FRiP scores
for bam in aligned/*_shifted.bam; do
    base=$(basename $bam _shifted.bam)
    
    total_reads=$(samtools view -c $bam)
    reads_in_peaks=$(bedtools intersect -u -a $bam -b peaks/consensus_peaks.bed | samtools view -c)
    
    frip=$(echo "scale=3; $reads_in_peaks / $total_reads" | bc)
    echo -e "${base}\t${frip}" >> results/frip_scores.txt
done
```

## 8. Differential Accessibility Analysis {#differential-analysis}

### Create Count Matrix
```bash
# Create count matrix for differential analysis
cat > create_count_matrix.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import glob

# Read consensus peaks
peaks = pd.read_csv('peaks/consensus_peaks.bed', sep='\t', 
                   names=['chr', 'start', 'end'])
peaks['peak_id'] = [f"peak_{i}" for i in range(len(peaks))]

# Read coverage files
coverage_files = glob.glob('peaks/*_coverage.txt')
count_matrix = peaks[['peak_id']].copy()

for file in coverage_files:
    sample_name = file.split('/')[-1].replace('_coverage.txt', '')
    coverage = pd.read_csv(file, sep='\t', 
                          names=['chr', 'start', 'end', 'name', 
                                'score', 'strand', 'count'])
    count_matrix[sample_name] = coverage['count']

count_matrix.to_csv('results/count_matrix.csv', index=False)
EOF

python create_count_matrix.py
```

### DESeq2 Analysis
```bash
cat > differential_analysis.R << 'EOF'
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Read count matrix
counts <- read.csv("results/count_matrix.csv", row.names = 1)

# Create sample metadata (modify according to your experimental design)
samples <- colnames(counts)
condition <- c(rep("Control", 3), rep("Treatment", 3))  # Example design
colData <- data.frame(condition = factor(condition), row.names = samples)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = colData,
                             design = ~ condition)

# Filter low count peaks
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "Treatment", "Control"))
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)

# Save results
write.csv(res_df, "results/differential_peaks.csv")

# Volcano plot
volcano_data <- res_df
volcano_data$significant <- abs(volcano_data$log2FoldChange) > 1 & 
                           volcano_data$padj < 0.05

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = "Volcano Plot - Differential Accessibility",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_minimal()

ggsave("figures/volcano_plot.pdf", p_volcano, width = 8, height = 6)

# PCA plot
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    labs(title = "PCA Plot") +
    theme_minimal()

ggsave("figures/pca_plot.pdf", p_pca, width = 8, height = 6)
EOF

Rscript differential_analysis.R
```

## 9. Visualization {#visualization}

### Heatmaps of Peak Regions
```bash
# Create heatmap around peak centers
computeMatrix scale-regions \
    -S aligned/*_shifted.bam \
    -R peaks/consensus_peaks.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 2000 \
    --regionBodyLength 1000 \
    --skipZeros \
    -o matrix_peaks.gz

plotHeatmap -m matrix_peaks.gz \
    -out figures/peak_heatmap.pdf \
    --colorMap Blues \
    --missingDataColor white
```

### Genome Browser Tracks
```bash
# Convert BAM to BigWig for visualization
for bam in aligned/*_shifted.bam; do
    base=$(basename $bam _shifted.bam)
    
    bamCoverage -b $bam \
        -o figures/${base}.bw \
        --binSize 10 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398  # hg38 effective genome size
done
```

### Peak Annotation
```bash
cat > annotate_peaks.R << 'EOF'
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read peaks
peaks <- readPeakFile("peaks/consensus_peaks.bed")

# Annotate peaks
peakAnno <- annotatePeak(peaks, tssRegion = c(-2000, 2000),
                        TxDb = txdb, annoDb = "org.Hs.eg.db")

# Plot annotation
p1 <- plotAnnoPie(peakAnno)
ggsave("figures/peak_annotation_pie.pdf", p1, width = 8, height = 6)

p2 <- plotAnnoBar(peakAnno)
ggsave("figures/peak_annotation_bar.pdf", p2, width = 8, height = 6)

# Distance to TSS
p3 <- plotDistToTSS(peakAnno, title = "Distribution of peaks relative to TSS")
ggsave("figures/distance_to_tss.pdf", p3, width = 8, height = 6)

# Save annotated peaks
write.csv(as.data.frame(peakAnno), "results/annotated_peaks.csv")
EOF

Rscript annotate_peaks.R
```

## 10. Advanced Analysis {#advanced-analysis}

### Motif Analysis
```bash
# Install HOMER if not available
# Then perform motif analysis on differential peaks

# Extract sequences for upregulated peaks
awk '$3 > 1 && $6 < 0.05 {print $7}' results/differential_peaks.csv | \
    tail -n +2 > results/upregulated_peaks.txt

# Get peak coordinates for upregulated peaks
grep -f results/upregulated_peaks.txt peaks/consensus_peaks.bed \
    > results/upregulated_peaks.bed

# Motif finding (requires HOMER)
# findMotifsGenome.pl results/upregulated_peaks.bed hg38 \
#     results/motif_analysis/ -size 200 -mask
```

### Footprinting Analysis
```bash
# TOBIAS footprinting analysis (if installed)
# TOBIAS ATACorrect --bam aligned/*_shifted.bam \
#     --genome hg38.fa \
#     --peaks peaks/consensus_peaks.bed \
#     --outdir footprints/ \
#     --cores 8
```

### Pathway Enrichment
```bash
cat > pathway_analysis.R << 'EOF'
library(clusterProfiler)
library(org.Hs.eg.db)

# Read differential peaks
diff_peaks <- read.csv("results/differential_peaks.csv")
anno_peaks <- read.csv("results/annotated_peaks.csv")

# Get upregulated peaks
up_peaks <- diff_peaks[diff_peaks$log2FoldChange > 1 & 
                      diff_peaks$padj < 0.05, ]

# Match with gene annotations
up_genes <- anno_peaks[anno_peaks$peak_id %in% up_peaks$peak_id, ]
gene_list <- unique(up_genes$geneId[!is.na(up_genes$geneId)])

# GO enrichment
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Plot results
p_go <- dotplot(ego, showCategory = 20)
ggsave("figures/go_enrichment.pdf", p_go, width = 10, height = 8)

# KEGG pathway analysis
kegg <- enrichKEGG(gene = gene_list,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)

p_kegg <- dotplot(kegg, showCategory = 20)
ggsave("figures/kegg_enrichment.pdf", p_kegg, width = 10, height = 8)
EOF

Rscript pathway_analysis.R
```

## 11. Troubleshooting {#troubleshooting}

### Common Issues and Solutions

**Low Library Complexity**
- Check for excessive PCR duplication
- Verify Tn5 concentration and incubation conditions
- Consider starting material quality

**Poor TSS Enrichment**
- Check for mitochondrial contamination
- Verify proper read shifting
- Consider nuclei isolation quality

**High Background Signal**
- Implement stricter quality filtering
- Remove blacklist regions
- Check for over-tagmentation

**Low Peak Numbers**
- Adjust MACS3 parameters (try -p 0.05 instead of 0.01)
- Check sequencing depth
- Verify proper peak filtering

### Quality Control Benchmarks

**Good ATAC-seq Library:**
- TSS enrichment > 7
- FRiP score > 0.3
- Insert size shows nucleosome pattern
- Duplication rate < 25%
- Mitochondrial reads < 10%

### Memory and Performance Tips

```bash
# For large datasets, consider using more memory-efficient tools
# Use sambamba instead of samtools for some operations
# Parallelize processing where possible
# Use appropriate thread numbers based on available cores

# Example with increased memory
export JAVA_OPTS="-Xmx32G"
```

## Summary

This tutorial covers the complete ATAC-seq analysis pipeline from raw reads to biological insights. Key steps include quality control, alignment, peak calling, differential analysis, and visualization. Always validate results with appropriate quality control metrics and consider the biological context when interpreting findings.

Remember to:
1. Document your analysis parameters
2. Keep track of software versions
3. Validate results with independent methods when possible
4. Consider batch effects in multi-sample studies
5. Integrate with other genomic datasets for comprehensive insights

For specific troubleshooting or advanced applications, consult the documentation of individual tools and consider the unique aspects of your experimental design.
