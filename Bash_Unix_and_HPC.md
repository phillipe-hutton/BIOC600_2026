# Bash, Unix, and HPC

## High Performance Computing
To be able to analyze genomic datasets, we need a few things:

1. A computer: to work with the sequencing data
2. Storage: to store the sequencing data
3. Memory: to be able to make calculations and process the sequencing data

Given that the human genome contains approximately 3 billion base pairs and the typical RNA-seq or ChIP-seq experiment requires around 20X coverage, this can easily be about 10 Gb for each sample you sequence. In other words: your laptop does not have the hardrive space nor the processing power to be able to analyze these experiments. You will need to access high performance computing (HPC) clusters that have the infrastructure to work with these large datasets. As graduate students of a Canadian univeristy, you can request access to the HPC clusters managed by the Digital Research Alliance of Canada, which are provided free of charge for academics doing research in Canada. If you see yourself working a lot with genomics datasets, you should ask your supervisor about sponsoring an account for you.

For this workshop, you will be working on a temporary virtual cluster provided by Calcul Québec.

There are two ways to login to a HPC:

1. Using a secure shell (ssh) (recommended)
2. Using JupyterHub via a web browser

## Bash and Unix
The **<ins>B</ins>ourne <ins>A</ins>gain <ins>S</ins>hell** (Bash) is a Unix shell that allows you to work with the operating system of a computer. HPC clusters run almost exclusively on Linux (which is similar to Unix). So, to be able to work with HPC computers, you will need to learn some Bash.

Here is a cheat sheet of common Bash commands and what they do:

| Command | Function |
| ------- | -------- |
| cd pathway_to_directory | Changes the working directory |
| cd .. | Go back one directory |
| ls | List the contents of a directory |
| ls -all | List all the contents of a directory + info |
| nano filename | Open a new or existing file and edit it |
| less filename | Open a file but you can't edit it |
| cp | Copy paste a file |
| mv | Move a file from one directory to another |
| wget | Download a file from a web link |
| bash | Run a bash script |
| sbatch | Submit a script on a job scheduler |
| module avail | See what software is available to use |
| module load | Load software |
| module purge | Unloads all software |
