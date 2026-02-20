# Welcome to the BIOC 600 Bioinformatics Workshop!

Within this GitHub repository, you will find all the information and resources you will need to successfully complete the workshop sections of BIOC 600.

The first workshop will cover RNA-seq analysis. By the end of this workshop, you will be able to:

- Work in the unix environmenet to write bash scripts and submit jobs on a high performance computing (HPC) cluster
- Analyze and perform introductory level data wrangling using R
- Download raw sequencing data from public data repositories
- Align raw sequencing files from an RNA-seq experiment to a reference genome
- Count the number of aligned reads for a gene and visualize your normalized reads on a genome browser
- Identify genes that are differentially expressed between your experimental condition and your control condition
- Determine what biological processes/molecular pathways your differentially expressed genes are involved in

The second workshop will cover ChIP-seq analysis. By the end of this workshop, you will be able to:

- Align raw sequencing files from a ChIP-seq experiment to a reference genome
- Filter mapped reads to remove PCR duplicates and low quality/unmapped reads
- Visualize ChIP peaks on a genome browswer
- Identify peaks using the MACS2 callpeaks function
- Use Homer to identify motifs for transcription factors that are enriched in your ChIP peaks

Between the first and second workshop, you will select and download either a ChIP-seq or an RNA-seq dataset that is relevant to your area of research. This should correspond to a complete experiment (e.g. a ChIP and a control, RNA-seq in two different conditions). It can be data from your own research lab, in which case you will need to be able to obtain relevant fastq files from your research lab, or from a published manuscript provided it has been deposited on the GEO database.

Between the second and third workshop, you will conduct analysis of your project with assistance as needed from the TAs. On the last workshop day you will give a ten minute presentation explaining:

1. What biological problem the paper (or your lab) is working on.
2. The details of the experiment you analyzed
3. The results

The grade will be based 50% on the presentation itself (how well you explained scientific background, why you conducted analysis you did) and 50% on completion and quality of the bioinformatic analysis.

To get started, you will need to download:

- [Visual Studio Code](https://code.visualstudio.com/Download)
- [R Studio and R](https://posit.co/download/rstudio-desktop/)

Optionally, you can create a GitHub account which will allow you to have GitHub CoPilot be available on VSCode and R Studio. There is a paid plan, but the free plan should be more than sufficient for your coding tasks. GitHub CoPilot functions as an AI autocomplete for writing scripts and it will significantly increase your speed for writing scripts. To enable GitHub CoPilot, do the following steps:

- VS Code: Add the "GitHub Copilot Chat" Extension and sign in to your GitHub account to enable the extension.
- R Studio: Open RStudio and navigate to Tools > Global Options > Copilot to enable Copilot.

![github-copilot-privacy-settings](<img width="945" height="423" alt="github-copilot-privacy-settings" src="https://github.com/user-attachments/assets/1256e71a-3eb5-4cc8-ba9b-678e816265b4" />
)

I encourage you to read through the workshop material in advance so that you are prepared for the day of the workshop. We will go in the following order:

1. Bash, Unix, and HPC
2. R and R Studio
3. RNA-seq Analysis
4. ChIP-seq Analysis

Focus on 1 and 2 at the beginning as this is the content that people tend to struggle with at the beginning. We look forward to seeing you at the first workshop!

Phillipe Hutton & Ishtiaque Hossain

BIOC 600 TAs
