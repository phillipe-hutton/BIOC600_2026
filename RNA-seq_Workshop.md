# BIOC 600 Workshop 1: RNA-seq

## I. Downloading Published Genomic Datasets
The [National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/) hosts the [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/), a data repository that allows people to freely access published genomics datasets.

Let's download RNA-seq data from this paper:

- Xia, H., Dufour, C.R., Medkour, Y. et al. Hepatocyte FBXW7-dependent activity of nutrient-sensing nuclear receptors controls systemic energy homeostasis and NASH progression in male mice. Nat Commun 14, 6982 (2023). [https://doi.org/10.1038/s41467-023-42785-3](https://doi.org/10.1038/s41467-023-42785-3)

We'll be using this specific dataset to compare the effects of feeding a high fat diet (HFD) to a mouse: [GSE205846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%20GSE205846)

You will see a page that looks like this:

![GEO Entry #1](https://github.com/user-attachments/assets/ff9ffa4a-3c14-41f8-b0a1-e9e71f9e22e3)

Scroll down to the bottom of the page and copy the BioProject accession number:

![GEO Entry #2](https://github.com/user-attachments/assets/a9f013bd-3a06-43d8-9175-693e41aeeb38)

We can enter this into [SRA-Explorer](https://sra-explorer.info/) to select the samples we want to analyze:

![SRA Explorer](https://github.com/user-attachments/assets/7fd87a10-ce34-4f97-8518-d96b767e8c62)

We'll be using the following samples:

![SRA Explorer Samples](https://github.com/user-attachments/assets/ede2c89b-77e1-4359-b9dd-32a588a91096)

Select the "Add 4 to collection" button and then select the cart with 4 saved datasets. From the "FastQ Downloads" tab, open the "Bash script for downloading FastQ files" section and copy the script.

![Download Bash script](https://github.com/user-attachments/assets/d0bcfcd6-22e9-4429-a276-a1e9e25dda38)


Open VSCode and open a new file. Select Text File.

![VSCode Welcome](https://github.com/user-attachments/assets/1b3798ca-60a9-44e6-9ba6-bf9d20b0c27a)

In the new Text File, click "select a language" and type in "Shell Script" and select it. Paste in the copied script from SRA-Explorer:

```
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/005/SRR19621805/SRR19621805_1.fastq.gz -o SRR19621805_GSM6234437_Liver_high-fat_diet_WT_biol_rep_2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/005/SRR19621805/SRR19621805_2.fastq.gz -o SRR19621805_GSM6234437_Liver_high-fat_diet_WT_biol_rep_2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/014/SRR19621814/SRR19621814_1.fastq.gz -o SRR19621814_GSM6234428_Liver_chow_diet_WT_biol_rep_1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/014/SRR19621814/SRR19621814_2.fastq.gz -o SRR19621814_GSM6234428_Liver_chow_diet_WT_biol_rep_1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/006/SRR19621806/SRR19621806_1.fastq.gz -o SRR19621806_GSM6234436_Liver_high-fat_diet_WT_biol_rep_1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/006/SRR19621806/SRR19621806_2.fastq.gz -o SRR19621806_GSM6234436_Liver_high-fat_diet_WT_biol_rep_1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/013/SRR19621813/SRR19621813_1.fastq.gz -o SRR19621813_GSM6234429_Liver_chow_diet_WT_biol_rep_2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/013/SRR19621813/SRR19621813_2.fastq.gz -o SRR19621813_GSM6234429_Liver_chow_diet_WT_biol_rep_2_Mus_musculus_RNA-Seq_2.fastq.gz
```

You will notice that you have 8 files instead of 4. This is because this experiment was done with paired-end sequencing rather than single-end. Therefore, you will have to process them as pairs.

Let's clean this script up so that our files are easier to read. I also find that curl tends to timeout when you try to download too many fastq files, so we will alter this script to be able to use wget which is more reliable:

```
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/014/SRR19621814/SRR19621814_1.fastq.gz -O chow_rep1_rna_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/014/SRR19621814/SRR19621814_2.fastq.gz -O chow_rep1_rna_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/013/SRR19621813/SRR19621813_1.fastq.gz -O chow_rep2_rna_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/013/SRR19621813/SRR19621813_2.fastq.gz -O chow_rep2_rna_R2.fastq.gz

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/006/SRR19621806/SRR19621806_1.fastq.gz -O hfd_rep1_rna_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/006/SRR19621806/SRR19621806_2.fastq.gz -O hfd_rep1_rna_R2.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/005/SRR19621805/SRR19621805_1.fastq.gz -O hfd_rep2_rna_R1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/005/SRR19621805/SRR19621805_2.fastq.gz -O hfd_rep2_rna_R2.fastq.gz
```

- Replace `curl` with `wget`
- Add the `-c`/`--continue` subcommand which allows it to resume getting a partially-downloaded file
- Replace the lowercase `-o` with the uppercase `-O` to write the downloaded file to our newly named file

Login to the HPC using the ssh in the terminal (click Terminal > New Terminal). Once you login, you can see what's in your home directory:

```
[phutton@login1 ~]$ ls
projects  scratch
```

Change the working directory to the scratch directory and we will create new directories that we will be storing our files. It's easier to do this at the beginning rather than trying to move several large files from one directory to another.

```
[phutton@login1 ~]$ cd scratch/
[phutton@login1 scratch]$ mkdir -p hfd-rna/{align,bigwig,counts,fastq/{raw,trim,fastqc},logs,scripts}
```

Set your working directory to the raw directory nested in the fastq directory, then copy and paste our modified bash script into your terminal and hit Enter. It will take a while for this to run.

## II. Trimming and Assessing the Quality of your FastQ Files

Once you have downloaded your FastQ files, it is good practice to trim the adapters and to assess the quality of your raw sequencing files. There's no point in analyzing a poorly sequenced experiment! One of the tools we can use to do this is called fastp. Fastp is able to autodetect the presence of adapters in your FastQ files, trim them out, and give you an html report file of the quality of the files.

```
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=1:00:00
#SBATCH --mem=10G

# Load modules
module load fastp/1.0.1

fastp -i /home/phutton/scratch/hfd-rna/fastq/raw/chow_rep1_rna_R1.fastq.gz \
     -I /home/phutton/scratch/hfd-rna/fastq/raw/chow_rep1_rna_R2.fastq.gz \
     -o /home/phutton/scratch/hfd-rna/fastq/trim/chow_rep1_rna_R1_trimmed.fastq.gz \
     -O /home/phutton/scratch/hfd-rna/fastq/trim/chow_rep1_rna_R2_trimmed.fastq.gz \
     -h /home/phutton/scratch/hfd-rna/fastq/trim/chow_rep1_rna_fastp.html \
     -j /home/phutton/scratch/hfd-rna/fastq/trim/chow_rep1_rna_fastp.json

fastp -i /home/phutton/scratch/hfd-rna/fastq/raw/chow_rep2_rna_R1.fastq.gz \
     -I /home/phutton/scratch/hfd-rna/fastq/raw/chow_rep2_rna_R2.fastq.gz \
     -o /home/phutton/scratch/hfd-rna/fastq/trim/chow_rep2_rna_R1_trimmed.fastq.gz \
     -O /home/phutton/scratch/hfd-rna/fastq/trim/chow_rep2_rna_R2_trimmed.fastq.gz \
     -h /home/phutton/scratch/hfd-rna/fastq/trim/chow_rep2_rna_fastp.html \
     -j /home/phutton/scratch/hfd-rna/fastq/trim/chow_rep2_rna_fastp.json

fastp -i /home/phutton/scratch/hfd-rna/fastq/raw/hfd_rep1_rna_R1.fastq.gz \
     -I /home/phutton/scratch/hfd-rna/fastq/raw/hfd_rep1_rna_R2.fastq.gz \
     -o /home/phutton/scratch/hfd-rna/fastq/trim/hfd_rep1_rna_R1_trimmed.fastq.gz \
     -O /home/phutton/scratch/hfd-rna/fastq/trim/hfd_rep1_rna_R2_trimmed.fastq.gz \
     -h /home/phutton/scratch/hfd-rna/fastq/trim/hfd_rep1_rna_fastp.html \
     -j /home/phutton/scratch/hfd-rna/fastq/trim/hfd_rep1_rna_fastp.json

fastp -i /home/phutton/scratch/hfd-rna/fastq/raw/hfd_rep2_rna_R1.fastq.gz \
     -I /home/phutton/scratch/hfd-rna/fastq/raw/hfd_rep2_rna_R2.fastq.gz \
     -o /home/phutton/scratch/hfd-rna/fastq/trim/hfd_rep2_rna_R1_trimmed.fastq.gz \
     -O /home/phutton/scratch/hfd-rna/fastq/trim/hfd_rep2_rna_R2_trimmed.fastq.gz \
     -h /home/phutton/scratch/hfd-rna/fastq/trim/hfd_rep2_rna_fastp.html \
     -j /home/phutton/scratch/hfd-rna/fastq/trim/hfd_rep2_rna_fastp.json
```

- `-i` is the **input** file for **read1** of the paired-end sequencing files
- `I` is the **input** file for **read2** of the paired-end sequencing files
- `-o` is the **output** file for **read1** of the paired-end sequencing files
- `-O` is the **output** file for **read2** of the paired-end sequencing files
- `-h` is the output file for the html report file
- `-j` is the output file for the json report file

As you can see, writing bash scripts like this is repetitive and becomes quite inefficient — especially when you start dealing with a large number of samples. Therefore, we should be asking ourselves "*Is there a way for me to structure my scripts so that I have to do the least amount of work possible?*" There are a few ways we can improve this script:

- We can create a job array that will allow us to send each sample as its own job to the job scheduler so that we don't have to wait for each command to be run sequentially
- We can create variables that will reduce the number of times we have to modify our scripts
- We can redirect the standard output and the standard error to a log directory in case our script isn't working to help with troubleshooting

We can rewrite our script to look something like this to implement these changes:

```
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --output=/home/phutton/scratch/hfd-rna/logs/fastp_%A_%a.out
#SBATCH --error=/home/phutton/scratch/hfd-rna/logs/fastp_%A_%a.err
#SBATCH --array=1-4

# -----------------------------
# Directories and Files
# -----------------------------
SAMPLES_FILE=~/scratch/hfd-rna/samples.txt
INPUT_DIR=~/scratch/hfd-rna/fastq/raw
OUTPUT_DIR=~/scratch/hfd-rna/fastq/trim

# Make directories if they don't exist
mkdir -p "${OUTPUT_DIR}"

# -----------------------------
# Setup Array
# -----------------------------

# Read the sample name for this array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES_FILE})

# -----------------------------
# Modules
# -----------------------------

module load fastp/1.0.1

# -----------------------------
# Define Input and Output Files
# -----------------------------

INPUT_R1="${INPUT_DIR}/${SAMPLE}_R1.fastq.gz"
INPUT_R2="${INPUT_DIR}/${SAMPLE}_R2.fastq.gz"
OUTPUT_R1="${OUTPUT_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
OUTPUT_R2="${OUTPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
HTML_REPORT="${OUTPUT_DIR}/${SAMPLE}_fastp.html"
JSON_REPORT="${OUTPUT_DIR}/${SAMPLE}_fastp.json"

# -----------------------------
# Run fastp
# -----------------------------

fastp -i "${INPUT_R1}" \
-I "${INPUT_R2}" \
-o "${OUTPUT_R1}" \
-O "${OUTPUT_R2}" \
-h "${HTML_REPORT}" \
-j "${JSON_REPORT}"
```

We will create a samples.txt file using the `nano` text editor. Create this txt file in your *hfd-rna* directory: `nano samples.txt`

```
chow_rep1_rna
chow_rep2_rna
hfd_rep1_rna
hfd_rep2_rna
```

Create your trimming bash script: `nano fastp.sh` and paste in your fastp script.

To submit the job to the Slurm job scheduler, type in `sbatch fastp.sh` and hit enter. You will be able to track your progress using the `sq` command.

If you run into any issues with a job, you can enter `sacct -j <job_id> --format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode,Elapsed,MaxRSS` and you will receive the error codes for the submitted jobs to get an idea of what when wrong. You should also look at your logs directory to see the .out and .err files for your corresponding jobs for further troubleshooting.

To view the report, login to the HPC cluster using JupyterHub, go to the trim directory, and open the html files. The report will look like this:

![Fastp Report: Summary](https://github.com/user-attachments/assets/d8d0afa0-375c-4b8b-8d4a-ab0b38eed7a8)

![Fastp Report: Adapters](https://github.com/user-attachments/assets/f9fc8031-ead5-4425-bc00-2b65f4a92916)

![Fastp Report: Insert Size](https://github.com/user-attachments/assets/ef66dadd-603c-48b8-95c9-cec756809084)

![Fastp Report: Filtering Stats](https://github.com/user-attachments/assets/697601f3-4c17-494a-91e2-8f51487dba0c)







## III. Aligning Reads to a Reference Genome

Before we can align reads to a reference genome, we must first index our reference genome so that the STAR software can work with it. **For the sake of time, we have already prepared an indexed genome for the mouse mm10 reference genome that you can use**.

If you require a reference genome for a difference organism, you can do the following:

Go to the UCSC Genome Browser

```
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --cpus-per-task=8 # Total memory will be 8 CPUs * 7G = 56G —— star is memory intensive
#SBATCH --output=/home/phutton/scratch/hfd-rna/logs/align_%A_%a.out
#SBATCH --error=/home/phutton/scratch/hfd-rna/logs/align_%A_%a.err
#SBATCH --array=1-4%2 # limit to 2 concurrent jobs to manage memory usage

# ------------------------------
# Directories and Files
# ------------------------------

SAMPLES_FILE=~/scratch/hfd-rna/samples.txt
REF_GENOME_DIR=~/projects/def-sponsor00/phutton/genomes/star_mm10 # keep as this to access my indexed genome
INPUT_DIR=~/scratch/hfd-rna/fastq/trim
OUTPUT_DIR=~/scratch/hfd-rna/align

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
     --runThreadN 8 \
     --genomeDir ${REF_GENOME_DIR} \
     --readFilesIn ${INPUT_DIR}/${SAMPLE}_R1_trimmed.fastq.gz ${INPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz \
     --readFilesCommand zcat \
     --sjdbOverhang 99 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard \
     --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE}_
```

## IV. Visualizing Your Aligned Sequencing Files


```
#!/bin/bash
#SBATCH --job-name=index
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --output=/home/phutton/scratch/hfd-rna/logs/index.out

# load modules
module load samtools/1.22.1

# index bam files
samtools index -M ~/scratch/hfd-rna/align/*_Aligned.sortedByCoord.out.bam
```

```
#!/bin/bash
#SBATCH --job-name=bigwig
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --output=/home/phutton/scratch/hfd-rna/logs/bigwig_%A_%a.out
#SBATCH --error=/home/phutton/scratch/hfd-rna/logs/bigwig_%A_%a.err
#SBATCH --array=1-4

# -----------------------------
# Directories and Files
# -----------------------------
SAMPLES_FILE=~/scratch/hfd-rna/samples.txt
INPUT_DIR=~/scratch/hfd-rna/align
OUTPUT_DIR=~/scratch/hfd-rna/bigwig
VIRT_ENV=~/projects/def-sponsor00/phutton/virtual_envs/deepTools_env # path to my virtual environment with deepTools installed

# -----------------------------
# Setup Array
# -----------------------------

# Read the sample name for this array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLES_FILE})

# -----------------------------
# Modules
# -----------------------------

module load python/3.11.5

# -----------------------------
# Activate Virtual Environment
# -----------------------------

source ${VIRT_ENV}/bin/activate

# -----------------------------
# Convert BAM to bigWig
# -----------------------------

bamCoverage -b ${INPUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
            --normalizeUsing CPM \
            -bs 20 \
            --smoothLength 60 \
            --numberOfProcessors 4 \
            -of bigwig \
            -o ${OUTPUT_DIR}/${SAMPLE}.bw

# deactivate virtual environment
deactivate
```

## V. Counting Your Aligned Reads
```
#!/bin/bash
#SBATCH --job-name=FeatureCounts
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --output=/home/phutton/scratch/hfd-rna/logs/counts.out
#SBATCH --error=/home/phutton/scratch/hfd-rna/logs/counts.err

# -----------------------------
# Directories and Files
# -----------------------------
ANNOTATION_DIR=~/projects/def-sponsor00/phutton/genomes/mm10.ensGene.gtf # keep as this to access my annotation file
INPUT_DIR=~/scratch/hfd-rna/align
OUTPUT_DIR=~/scratch/hfd-rna/counts

# -----------------------------
# Load Modules
# -----------------------------

module load StdEnv/2020 gcc/12.3 subread/2.0.6

featureCounts -T 4 \
              -p \
              --countReadPairs \
              -B \
              -C \
              -s 0 \
              -t exon \
              -g gene_id \
              -a ${ANNOTATION_DIR} \
              -o ${OUTPUT_DIR}/counts.txt \
              ${INPUT_DIR}/*_Aligned.sortedByCoord.out.bam
````

We can now clean this up using R

```
#R

# Install packages
install.packages("tidyverse")
install.packages("tibble")
install.packages("BiocManager")
BiocManager::install("AnnotationDbi")
BiocManager::install("DESeq2")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pheatmap")

# Load libraries
library(tidyverse)
library(tibble)
library(BiocManager)
library(AnnotationDbi)
library(DESeq2)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(pheatmap)
```

Next we will clean up the counts file

```
# Read featureCounts output
raw_counts <- read.table(
  "counts.txt",
  header = TRUE,
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Convert to data frame
raw_counts <- as.data.frame(raw_counts)

# Keep Geneid + count columns, rename for clarity
countMatrix <- raw_counts %>%
  dplyr::select(Geneid, contains("chow_rep"), contains("hfd_rep")) %>%
  dplyr::rename(
    chow_rep1 = "/home/phutton/scratch/hfd-rna/align/chow_rep1_rna_Aligned.sortedByCoord.out.bam",
    chow_rep2 = "/home/phutton/scratch/hfd-rna/align/chow_rep2_rna_Aligned.sortedByCoord.out.bam",
    hfd_rep1 = "/home/phutton/scratch/hfd-rna/align/hfd_rep1_rna_Aligned.sortedByCoord.out.bam",
    hfd_rep2 = "/home/phutton/scratch/hfd-rna/align/hfd_rep2_rna_Aligned.sortedByCoord.out.bam"
  )

# set Geneid as row names and remove the Geneid column
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix %>% dplyr::select(-Geneid)

# Save countMatrix as a csv
write.csv(countMatrix, "countMatrix.csv")
```


## VI. Analyzing Differentially Expressed Genes (DEGs)
```
#R

```
## VII. Visualizing DEGs
```
#R

```

## VIII. Gene Ontology: What Pathways Are These DEGs Involved In?
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
## IX. Resources
### Software Documentation
- FastQC
- Trimmomatic
- STAR
- featureCounts
- DESeq2

### Recommened Readings

### Free Programming Courses/Tutorials

### Free Online Tools

