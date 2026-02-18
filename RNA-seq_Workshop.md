# BIOC 600 Workshop 1: RNA-seq

## I. Gaining Access to the High-Performance Computing (HPC) Clusters

### Accessing via a secure shell (ssh)
1. Open the terminal (MacOS) or command prompt (PC)
2. Type `ssh<username>@bioc600d.calculquebec.cloud`
3. Type `<password>` (note: the password keystrokes won't appear in the terminal window as a security measure)

### Accessing via JupyterHub
<INSERT>

## I. Downloading Published Genomic Datasets
The National Center for Biotechnology Information (NCBI) hosts the Gene Expression Omnibus (GEO), a data repository that allows people to freely access published genomics datasets.

## II. Assessing the Quality of FastQ Files

```
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --output=fastqc.out

# load modules
module load fastqc/0.11.9

# perform fastqc on fastq files
fastqc chow_rep1_rna_R1.fastq.gz -t 4 -o ./fastqc/raw
fastqc chow_rep1_rna_R2.fastq.gz -t 4 -o ./fastqc/raw
fastqc chow_rep2_rna_R1.fastq.gz -t 4 -o ./fastqc/raw
fastqc chow_rep2_rna_R2.fastq.gz -t 4 -o ./fastqc/raw
fastqc hfd_rep1_rna_R1.fastq.gz -t 4 -o ./fastqc/raw
fastqc hfd_rep1_rna_R2.fastq.gz -t 4 -o ./fastqc/raw
fastqc hfd_rep2_rna_R1.fastq.gz -t 4 -o ./fastqc/raw
fastqc hfd_rep2_rna_R2.fastq.gz -t 4 -o ./fastqc/raw
```

However, writing bash scripts like this is repetitive and thus inefficient. The question we should constantly be asking ourselves is "Is there a way for me to structure my scripts so that I have to do the least amount of work possible?"

```
#!/bin/bash
#SBATCH --job-name=fastqc_array
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --output=fastqc_%A_%a.out
#SBATCH --array=1-4

# ------------------------------
# Directories and Files
# ------------------------------

SAMPLES_FILE=samples.txt
OUTPUT_DIR=./fastqc/raw

# Make directories if they don't exist
mkdir -p ${OUTPUT_DIR}

# -----------------------------
# Setup SLURM Array
# -----------------------------

# Read the sample name for this array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES_FILE})

# ------------------------------
# Load Modules
# ------------------------------

module load fastqc/0.11.9

# ------------------------------
# Perform FastQC on the FastQ files
# ------------------------------

fastqc ${SAMPLE}_R1.fastq.gz -t 4 -o ${OUTPUT_DIR}
fastqc ${SAMPLE}_R2.fastq.gz -t 4 -o ${OUTPUT_DIR}
```

## III. Trimming your FastQ Files (Optional)

```
#!/bin/bash
#SBATCH --job-name=trimming_pe
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=4
#SBATCH --output=trimming_%A_%a.out
#SBATCH --array=1-4

# ------------------------------
# Directories and Files
# ------------------------------
SAMPLES_FILE=samples.txt
ADAPTER_FILE=TruSeq3-PE.fa
OUTPUT_DIR=./fastq/trim

# Make directories if they don't exist
mkdir -p ${OUTPUT_DIR}

# -----------------------------
# Setup SLURM Array
# -----------------------------

# Read the sample name for this array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES_FILE})

# ------------------------------
# Load Modules
# ------------------------------
module load trimmomatic/0.39

# ------------------------------
# Perform Trimming with Trimmomatic
# ------------------------------

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
     -threads 4 \
     -phred33 \
     -trimlog ${SAMPLE}_trim_log \
     ${SAMPLE}_R1.fastq.gz \
     ${SAMPLE}_R2.fastq.gz \
     ${SAMPLE}_R1_trimmed.fastq.gz \
     ${SAMPLE}_R2_trimmed.fastq.gz \
     ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:True \
     LEADING:3 \
     TRAILING:3 \
     SLIDINGWINDOW:4:15 \
     MINLEN:36
```

## IV. Aligning Reads to a Reference Genome
```
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --output=align_%A_%a.out
#SBATCH --error=align_%A_%a.err
#SBATCH --array=1-4

# ------------------------------
# Directories and Files
# ------------------------------

SAMPLES_FILE=samples.txt
REF_GENOME_DIR=genomes/mm10.fa
OUTPUT_DIR=./align

# Make directories if they don't exist
mkdir -p ${OUTPUT_DIR}

# -----------------------------
# Setup SLURM Array
# -----------------------------

# Read the sample name for this array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES_FILE})

# ------------------------------
# Load Modules
# ------------------------------
module load star/2.7.11b

# ------------------------------
# Perform Alignment with STAR
# ------------------------------

STAR --runMode alignReads \
     --runThreadN 4 \
     --genomeDir ${REF_GENOME_DIR} \
     --readFilesIn ${SAMPLE}_R1_trimmed.fastq.gz ${SAMPLE}_R2_trimmed.fastq.gz \
     --readFilesCommand zcat \
     --sjdbOverhang 99 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard \
     --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE}_
```

## V. Visualizing Your Aligned Sequencing Files
```
#!/bin/bash
#SBATCH --job-name=index
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --output=index.out

# load modules
module load samtools/1.17

# index bam files
samtools index *_Aligned.sortedByCoord.out.bam
```

```
#!/bin/bash
#SBATCH --job-name=bigwig
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --output=bigwig_%A_%a.out
#SBATCH --array=1-4

# -----------------------------
# Directories and Files
# -----------------------------

SAMPLES_FILE=samples.txt

# -----------------------------
# Setup Array
# -----------------------------

# Read the sample name for this array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES_FILE})

# -----------------------------
# Modules
# -----------------------------

module load python/3.5.6

# -----------------------------
# Activate Virtual Environment
# -----------------------------

source ~/virtual_envs/deepTools_env/bin/activate

# -----------------------------
# Convert BAM to bigWig
# -----------------------------

bamCoverage -b ${SAMPLE}_Aligned.sortedByCoord.out.bam \
            --normalizeUsing CPM \
            -bs 20 \
            --smoothLength 60 \
            --numberOfProcessors 4 \
            -of bigwig \
            -o ${SAMPLE}.bw

# Deactivate Virtual Environment
deactivate
```

## VI. Counting Your Aligned Reads
```
#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --output="counts_progress.txt"

# load modules
module load StdEnv/2020 gcc/9.3.0 subread/2.0.3

# count reads
featureCounts -T 4 \
              -s 2 \
              -a <genome_gtf> \
              -o <output_file> \
              <input_bam>
````

## VII. Analyzing Differentially Expressed Genes (DEGs)
```
#R

```
## VIII. Visualizing DEGs
```
#R

```

## IX. Gene Ontology: What Pathways Are These DEGs Involved In?
```
#R
library(clusterProfiler)
library(enrichplot)
library(org.Mus.eg.db)

data <- read.table()
colnames(data)[1] <- "gene_names"

head(data)
# we want the log2 fold change
original_gene_list <- data$log2FoldChange

# name the vector
names(original_gene_list) <- data$gene_names

# omit any NA values
gene_list <- na.omit(original_gene_list)

# sort the lsit in decreasing order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing=TRUE)

# extract significant results (padj < 0.05)
sig_genes_df <- subset(data, padj < 0.05)

# from significant results, we want to filter on log2 fold change
sig_genes_df <- subset(sig_genes_df, log2FoldChange >= 2)

# name the vector
genes <- sig_genes_df$log2FoldChange
names(genes) <- sig_genes_df$gene_names

# omit NA values
genes <- na.omit(genes)

# filter on min log2 fold change (log2FoldChange >= 2)
genes <- names(genes)[abs(genes) >= 2]

# GO plot
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = org.Mus.eg.db,
                      keyType = "ENSEMBL",
                      readable = TRUE,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05
)

dotplot(go_enrich)
```
## X. Resources
### Software Documentation
- FastQC
- Trimmomatic
- STAR
- featureCounts
- DESeq2

### Recommened Readings

### Free Programming Courses/Tutorials

### Free Online Tools

