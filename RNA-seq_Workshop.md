# BIOC 600 Workshop 1: RNA-seq

![RNA-seq Workflow](https://github.com/user-attachments/assets/c285ea45-de4a-4d06-86e2-d84f65383deb)

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
[phutton@login1 scratch]$ mkdir -p hfd-rna/{align,bigwig,counts,fastq/{raw,trim,},logs,scripts}
```

Set your working directory to the raw directory nested in the fastq directory, then copy and paste our modified bash script into your terminal and hit Enter. It will take a while for this to run.

## II. Trimming and Assessing the Quality of your FastQ Files

Once you have downloaded your FastQ files, it is good practice to trim the adapters and to assess the quality of your raw sequencing files. There's no point in analyzing a poorly sequenced experiment! One of the tools we can use to do this is called fastp. Fastp is able to autodetect the presence of adapters in your FastQ files, trim them out, and give you an html report file of the quality of the files.

We can write the scripts as follows:

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

If you run into any issues with a job, you can enter 

`sacct -j <job_id> --format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode,Elapsed,MaxRSS` 

and you will receive the error codes for the submitted jobs to get an idea of what when wrong. You should also look at your logs directory to see the .out and .err files for your corresponding jobs for further troubleshooting.

To view the report, login to the HPC cluster using JupyterHub, go to the trim directory, and open the html files. The report will look like this:

![Fastp Report: Summary](https://github.com/user-attachments/assets/d8d0afa0-375c-4b8b-8d4a-ab0b38eed7a8)

![Fastp Report: Adapters](https://github.com/user-attachments/assets/f9fc8031-ead5-4425-bc00-2b65f4a92916)

![Fastp Report: Insert Size](https://github.com/user-attachments/assets/ef66dadd-603c-48b8-95c9-cec756809084)

![Fastp Report: Filtering Stats](https://github.com/user-attachments/assets/697601f3-4c17-494a-91e2-8f51487dba0c)

Once you have validated that your sequencing files are of good quality, we can go ahead and align the reads to a reference genome.

## III. Aligning Reads to a Reference Genome

We will be using the STAR (<ins>S</ins>pliced <ins>T</ins>ranscripts <ins>A</ins>lignment to a <ins>R</ins>eference) aligner to map our RNA-seq reads to the mouse mm10 reference genome. Before we can align reads to a reference genome, we must first index our reference genome so that the STAR software can work with it. **For the sake of time, we have already prepared an indexed genome for the mouse mm10 reference genome that you can use**.

If you require a reference genome for a difference organism for this workshop, let us know and we can prepare this for you.

To run the aligner, we will use this script:

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
- `--runMode` tells STAR we are aligning reads to a reference genome
- `--runThreadN 8` tells STAR that we will be using 8 cpus per task
- `--genomeDir` specifies the pathway of our indexed genome directory
- `--readFilesIn` specifies the pathway to our input fastq files 
- `--readFilesCommand zcat` specifies that our files are gzipped
- `--sjdbOverhang 99` specifies the length of the sequence around the annotated splice junction (ideal length is 99)
- `--outSAMtype BAM SortedByCoordinate` directs the output files to be in the binary alignment map (BAM) format
- `--outSAMunmapped Within` tells STAR to keep unmapped reads but specifies that they will not be used
- `--outSAMattributes Standard` indicates standard setting for SAM attributes
- `--outFileNamePrefix` specifies the output folder and the prefix you want your aligned files to start with

After the alignment has finished, you should have 5 output files:

1. Log.Final.out: a summary file of your alignment
2. Aligned.SortedByCoordinant.bam: your aligned reads
3. Log.out: a progress output file
4. Log.progress.out: a detailed progress file 
5. SJ.out.tab: a splice junction output file

You can look inside your Log.Final.out files to get an idea as to how well your reads mapped to the reference genome. We ideally want to see around 60-75% of uniquely mapped reads to be mapped to the reference genome. Lower than this is a red flag that points to issues in our alignment. There are additional steps you can take to do QC on your RNA-seq data such as running [RNA-SeQC](https://github.com/getzlab/rnaseqc) that allow you to get a better understanding of the biases in your sequencing data.

Finally, a great way to assess our data is to visualize it on a genome viewer.

## IV. Visualizing Your Aligned Sequencing Files

We will use IGV to visualize our aligned RNA-seq data. IGV can accept bam files if they're indexed, so it's good practice to index your bam files:

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

However, bam files are huge and not fun to transfer to your local hard drive. So, we will show some mercy for our local storage space by generating bigWig files from our bam files. bigWig files are much smaller and load nicely onto IGV. To generate bigWigs, we will be using the bamCoverage command from deepTools.

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
INPUT_DIR=~/projects/def-sponsor00/phutton/hfd-rna/align
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

Once your job has finished running, we can download our bigWig files from JupyterHub and load them into IGV. This is a great way to see how your treatment impacts the expression of different genes. This also allows to to combine RNA-seq data with ChIP-seq tracks to see where transcription factors of interest are binding as well as changes in histone marks.

![IGV Tracks](https://github.com/user-attachments/assets/7569dfa3-b132-4795-84b2-a83e7a039f06)

However, going through the entire genome gene by gene is impratical. To identify genes whose expression is different in our treatment vs our control condition, we will need to first count our uniquely mapped reads and test whether the change in expression is statistically significant.

## V. Counting Your Aligned Reads

We will be using featureCounts to count our uniquely mapped reads.

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
INPUT_DIR=~/projects/def-sponsor00/phutton/hfd-rna/align
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

- `-T` number of threads for the command
- `-p` specifies that input data is paired-end
- `countReadPairs` specifies that read pairs will be counted instead of reads (should be included for paired-end data)
- `-B` directs featureCounts to only count fragments with both end successfully aligned
- `-C` directs featureCounts to **NOT** count chimeric fragments (fragments that have their two ends aligned to different chromosomes)
- `-s` specifies that the data is unstranded (0), stranded (1), or reversely stranded (2)
- `-t` specifies the feature type being counted (exons in our case)
- `-g` specifies the attribute type used to group features (e.g. exons) into meta-features (e.g. genes) when gtf annotation is provided (default = gene_id)
- `-a` specifies the pathway to the gene annotation gtf file
- `-o` name of the output file

This will generate a counts.txt file with the counts assigned to an Ensembl Gene ID. To make this file useable, we will clean it up using R

First, make sure you have these packages installed and loaded:

```
#R

# Install packages
install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("AnnotationDbi")
BiocManager::install("DESeq2")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pheatmap")

# Load libraries
library(tidyverse)
library(BiocManager)
library(AnnotationDbi)
library(DESeq2)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(pheatmap)
```

Now we can clean up the counts file

```
#R

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

Next, we will set up our data to determine which which genes are differentially expressed.

## VI. Analyzing Differentially Expressed Genes (DEGs)

To obtain our list of DEGs, we will be using DESeq2. To set up the experiment for DESeq2, we will first create a meta data file that provides information about our experimental design.

```
#R

#Create sample metadata
meta <- data.frame(
  sample = c("chow_rep1", "chow_rep2", "hfd_rep1", "hfd_rep2"),
  condition = c("chow", "chow", "hfd", "hfd"),
  row.names = "sample"
)

# Save meta as a csv
write.csv(meta, "meta.csv")
```

Next, we will validate that our count matrix and our meta data is set up properly for DESeq2 to work.

```
#R

# Validation of count data and meta data for DESeq2
class(countMatrix)
class(meta)
all(colnames(countMatrix) %in% rownames(meta))
all(colnames(countMatrix) == rownames(meta))
head(countMatrix)
head(meta)
```

If done properly, you should get this in response:

```
[1] "data.frame"
> class(meta)
[1] "data.frame"
> all(colnames(countMatrix) %in% rownames(meta))
[1] TRUE
> all(colnames(countMatrix) == rownames(meta))
[1] TRUE
> head(countMatrix)
chow_rep1 chow_rep2 hfd_rep1 hfd_rep2
ENSMUSG00000102693         0         0        0        0
ENSMUSG00000064842         0         0        0        0
ENSMUSG00000051951         0         0        1        5
ENSMUSG00000102851         0         0        0        0
ENSMUSG00000103377         0         0        0        0
ENSMUSG00000104017         0         0        0        0
> head(meta)
condition
chow_rep1      chow
chow_rep2      chow
hfd_rep1        hfd
hfd_rep2        hfd
```

Before proceeding with DESeq2, we should first ensure that our controls are listed first in the factor levels for our treatment to make sure our comparison is being done properly. 

```
#R

# Ensure our controls are first in the factor levels for condition
meta$condition <- factor(meta$condition, levels = c("chow", "hfd"))
```

Now we can run DESeq2:

```
#R

# Create DESeqDataSet object from our count matrix and run DESeq2
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = meta, design = ~ condition)
dds <- DESeq(dds)
```

To obtain the name of our DESeq2 comparison, we can obtain the results name:

```
#R

# get results names
resultsNames(dds)
```

Which will give us this:

```
[1] "Intercept"        "condition_hfd_vs_chow"
```

And now we can extract the results of this comparison:

```
#R

# Extract our DESeq2 results
res <- results(dds, name = "condition_hfd_vs_chow")

# Convert to a data frame
res <- as.data.frame(res)
```

Now we must filter our results so that we only keep statistically significant values

```
#R

# Filter for significant DEGs (padj < 0.05 and FoldChange > |1.5|) and remove NAs
sig_res <- res %>%
  filter(!is.na(padj), padj < 0.05,
       !is.na(log2FoldChange), abs(log2FoldChange) > log2(1.5))
```

To be able to understand our significant results, we must convert the Ensembl Gene IDs to something that we can actually read and understand. We will use AnnotationDbi to help us with this using the org.Mm.eg.db mouse gene database.

```
#R

# Add gene symbols to the significant results
sig_res$gene_symbol <- mapIds(org.Mm.eg.db,
                              keys = rownames(sig_res),
                              column = "SYMBOL",
                              keytype = "ENSEMBL")

# Merge and create list of significant DEGs
sig_res <- sig_res %>%
  rownames_to_column(var = "ensemblID") %>%
  dplyr::select(ensemblID, gene_symbol, everything()) %>%
  filter(!is.na(gene_symbol)) %>%
  arrange(desc(log2FoldChange))

# Save significant DEGs to a csv
write.csv(sig_res, "significant_DEGs.csv")
```

Now we have an annotated list of DEGs that we can analyze.

## VII. Visualizing DEGs

First, let's see how many DEGs are upregulated vs downregulated under a high fat diet:

```
#R

# Select upregulated DEGs
upreg_DEGs <- sig_res %>%
  filter(log2FoldChange > 0)

# Select downregulated DEGs
downreg_DEGs <- sig_res %>%
  filter(log2FoldChange < 0)

# Create summary data to generate bar plot
deg_summary <- data.frame(
  category = c("Upregulated", "Downregulated"),
  count = c(nrow(upreg_DEGs), nrow(downreg_DEGs))
)

# plot upregulated vs downregulated DEGs as a bar plot
ggplot(deg_summary, aes(x = category, y = count, fill = category)) +
  geom_col(alpha = 0.7) +
  geom_text(aes(label = count), 
            vjust = -0.5,
            family = "Helvetica",
            size = 5) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 800)) +
  labs(title = "Differentially Expressed Genes in HFD vs Chow",
       x = "",
       y = "# of DEGs") +
  scale_fill_manual(values = c("Upregulated" = "firebrick", 
                               "Downregulated" = "steelblue")) +
  theme(
    plot.title = element_text(family = "Helvetica", size = 16),
    axis.title = element_text(family = "Helvetica", size = 12),
    axis.text.x = element_text(family = "Helvetica", size = 12, angle = 45, hjust = 1),
    axis.line.y = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none",
    panel.background = element_blank()
  )
```

This will give us a bar plot that looks like this:

![HFD DEGs Bar Plot](https://github.com/user-attachments/assets/0640f78f-02ad-4d1f-b61a-2939c84f78be)

We can also generate a volcano plot:

```
#R

# Volcano plot of all DEGs
plot_data <- res

plot_data$gene_symbol <- mapIds(org.Mm.eg.db,
                              keys = rownames(plot_data),
                              column = "SYMBOL",
                              keytype = "ENSEMBL")

EnhancedVolcano(
  plot_data,
  lab = plot_data$gene_symbol,
  x = 'log2FoldChange',
  xlab = bquote(~Log[2]~ "(FC)"),
  y = 'padj',
  ylab = bquote(~-Log[10]~ "(Padj)"),
  legendLabels = c("NS", "Log2(FC)", "Padj", "Padj & Log2(FC)"),
  pCutoff = 0.05,
  FCcutoff = 1.5,
  legendPosition = 'right',
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  labSize = 3.0,
  colAlpha = 0.7,
  max.overlaps = 20,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  caption = NULL,
  title = "Differentially Expressed Genes in HFD vs Chow",
)
```

Which will give us this:

![Volcano_Plot](https://github.com/user-attachments/assets/576a5eb1-e5ee-4037-a7d3-608b01d06548)

We can also use our list of DEGs to figure out what biological pathways they're involved in.

## VIII. Gene Ontology: What Pathways Are These DEGs Involved In?

Let's take our list of DEGs and see what biological processes they're associated with. We can use the [Enrichr](https://maayanlab.cloud/Enrichr/) webtool to access several Gene Set Enrichment Analyses and Gene Ontology Pathways. We will paste the gene names of our DEGs and look at MSigDB Hallmark 2020.

![Enrichr Home Page](https://github.com/user-attachments/assets/e900b87b-e149-4889-bfc5-de592f241cfd)

![Enrichr MSigDB](https://github.com/user-attachments/assets/3a85d3ca-a0b7-4102-9bcf-c71159fe4575)

We can export the entries to a table and take a closer look:

![MSigDB File](https://github.com/user-attachments/assets/3beb7e33-fd7f-4bb0-8a80-d49fd57dbdcb)

Adding the file into our working directory, we can create a bar plot to show the top 10 hits we get back ordered by p-value:

```
#R

# Read table and select top 10 significant pathways based on adjusted p-value
msigdb <- read.table("MSigDB_Hallmark_2020_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
msigdb <- as_tibble(msigdb)
msigdb <- msigdb %>%
  arrange(`P-value`) %>%
  dplyr::select(`Term`,`P-value`,`Adjusted P-value`,`Genes`) %>%
  slice_head(n = 10)

# Mutate to add -log10(P-value) and split genes into lists
msigdb <- msigdb %>%
  mutate(`-log10(P-value)` = -log10(`P-value`)) %>%
  mutate(gene_list = strsplit(Genes, ";")) %>%
  dplyr::select(-Genes)
  
# Make bar plot of top 10 significant pathways using -log10(P-value) as the x-axis and Term as the y-axis
ggplot(msigdb, aes(x = `-log10(P-value)`, y = reorder(Term, `-log10(P-value)`), fill = `-log10(P-value)`)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  labs(title = "MSigDB HallMark 2020",
       x = "-log10(P-value)",
       y = "") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 25)) +
  theme(
    plot.title = element_text(family = "Helvetica", size = 16),
    axis.title = element_text(family = "Helvetica", size = 12),
    axis.text = element_text(family = "Helvetica", size = 10),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none",
    panel.background = element_blank()
  )
```

Which gives us this:

![MSigDB_BarPlot](https://github.com/user-attachments/assets/7ea1e1c6-aa0a-41d5-91c8-d0a4a8a968b3)


Let's take a closer look at the genes enriched for "*Epithelial Mesenchymal Transition*"

```
#R

# Extract EMT DEGs List
emt_genes <- msigdb %>%
  filter(Term == "Epithelial Mesenchymal Transition") %>%
  unnest("gene_list") %>%
  dplyr::select(gene_list) %>%
  as.data.frame()

# Extract from Significant DEGs List
emt_genes <- emt_genes %>%
  mutate(gene_list = str_to_title(gene_list)) # Convert to title case
emt_genes <- filter(sig_res, gene_symbol %in% emt_genes$gene_list)
```

One way that we can visualize how the expression of this set of genes is altered under a high fat diet is to plot a heatmap using pheatmap. To do this, we will use the Regularized Logarithmic (rlog) transformation from DESeq2. DESeq2 also has the Variance Stabilizing Transformation (vst) which is better for larger datasets (anything larger than 50 genes). Since we are working with 46 genes, we will opt for the rlog transformation since it looks nicer.

```
#R

# Prepare matrix data for heatmap
genes_to_plot <- emt_genes$ensemblID
emt_degs <- rlog(dds[genes_to_plot,], blind = TRUE)
mat <- assay(emt_degs)[genes_to_plot,]
gene_symbols_vec <- setNames(emt_genes$gene_symbol, emt_genes$ensemblID)
matched_symbols <- gene_symbols_vec[rownames(mat)]
matched_symbols[is.na(matched_symbols)] <- rownames(mat)[is.na(matched_symbols)]
rownames(mat) <- matched_symbols
mat_z <- t(scale(t(mat)))
annotation <- data.frame(Condition = c("Chow", "Chow", "HFD", "HFD"))
rownames(annotation) <- colnames(mat_z)

# plot heatmap
pheatmap(mat_z,
         annotation_col = annotation,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 6,
         border_color = NA,
         scale = "none")

```

Plotting the heatmap will give us something that looks like this:

![Pheatmap_EMT](https://github.com/user-attachments/assets/63fb2778-d9bc-4c53-b837-f743057ac064)

You now have all the tools and resources to be able to analyze RNA-seq data and plot some useful figures to better understand your data.

Good luck with the rest of the workshop!

## IX. Resources

If you are interested in learning more about the tools you worked with during this workshop, if you want to dive deeper into bioinformatics, or if you find yourself struggling, I've linked below resources that I hope you will find useful.

### Software Documentation
- [Fastp](https://github.com/OpenGene/fastp)
- [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [Samtools](https://www.htslib.org/)
- [deepTools](https://deeptools.readthedocs.io/en/latest/index.html)
- [featureCounts](https://subread.sourceforge.net/featureCounts.html)
- [Analyzing RNA-seq data with DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

### Free Programming Courses/Tutorials
- [McGill Initiative in Computational Medicine](https://www.mcgill.ca/micm/)
- [Codecademy: Learn R](https://www.codecademy.com/enrolled/courses/learn-r)
- [Codecademy: Python for Programmers](https://www.codecademy.com/enrolled/courses/python-for-programmers)
- [Codecademy: Getting Started with Python for Data Science](https://www.codecademy.com/enrolled/courses/getting-started-with-python-for-data-science)

### Free Online Tools
- [SRA-Explorer](https://sra-explorer.info/): Search and download fastq files
- [Enrichr](https://maayanlab.cloud/Enrichr/): Look up Gene Set Enrichment and Gene Ontologies for a list of genes
- [UCSC Genome Browser](https://genome.ucsc.edu/): Online genome browser and a resource for downloading reference genomes and gene annotation files
- [Paletton](https://paletton.com/#uid=1000u0kllllaFw0g0qFqFg0w0aF): Guides the selection of a colour palette for designing figures
- [Bioart](https://bioart.niaid.nih.gov/): Scientific illustrations provided for free by the NIH for any use.

