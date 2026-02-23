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

To get started, you will need to download the following Integrated Development Environments (IDEs):

- [Visual Studio Code](https://code.visualstudio.com/Download)
- [R Studio and R](https://posit.co/download/rstudio-desktop/) (Download both R and R Studio)

Optionally, you can create a GitHub account which will allow you to have GitHub CoPilot be available on VSCode and R Studio. **Do not get the paid plan**—the free plan should be more than sufficient for your coding tasks. GitHub CoPilot functions as an AI autocomplete for writing scripts and it will significantly increase your speed for writing scripts. To enable GitHub CoPilot, do the following steps:

- VS Code: Add the "GitHub Copilot Chat" Extension and sign in to your GitHub account to enable the extension.
- R Studio: Open RStudio and navigate to Tools > Global Options > Copilot and sign in to enable GitHub Copilot.

Note that the usage of Copilot, as with other AI tools, should be done with caution. Check your privacy settings and be cautious about entering sensitive data into IDEs with Copilot enabled. I recommend setting your privacy settings as shown in the image below to ensure that your code won't be used for training AI models.

---
![github-copilot-privacy-settings](https://github.com/user-attachments/assets/7f859dc9-6f87-4c44-9202-faa8d65fe753)
---

Additionally, I advise that this tool should be used as an advanced "autocomplete" like when you text your friends: It's handy when it gives you the response that you're looking for—but it won't always give you the response that you actually want to use. So always check your scrips to make sure there's no errors in it.

I encourage you to read through the workshop material in advance so that you are prepared for the day of the workshop. We will go in the following order:

1. Bash, Unix, and HPC
2. R and R Studio
3. RNA-seq Analysis
4. ChIP-seq Analysis

Focus on 1 and 2 at the beginning as this is the content that people tend to struggle with at first. We look forward to seeing you at the first workshop!

Phillipe Hutton & Ishtiaque Hossain

BIOC 600 TAs
