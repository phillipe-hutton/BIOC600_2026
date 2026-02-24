# Bash, Unix, and HPC

## High Performance Computing
To be able to analyze genomic datasets, we need a computer that has enough:

1. Storage: to store the sequencing data
2. Memory: to be able to make calculations and process the sequencing data

Given that the human genome contains approximately 3 billion base pairs and the typical RNA-seq or ChIP-seq experiment requires around 20X coverage, this can easily be about 10 Gb for each sample you sequence. In other words: your laptop does not have the hardrive space nor the computational power to be able to process data of this size. You will therefore need to access a high performance computing (HPC) cluster that has the infrastructure to work with these large datasets. As graduate students of a Canadian univeristy, you can request access to the HPC clusters managed by the [Digital Research Alliance of Canada](https://www.alliancecan.ca/en), which are provided free of charge for academics doing research in Canada. If you see yourself working a lot with genomics datasets, you should ask your supervisor about sponsoring an account for you.

For this workshop, you will be working on a temporary virtual cluster provided by Calcul Québec.

There are two ways to login to a HPC:

### Using a secure shell (ssh) (recommended)

1. Open your terminal (Linux or macOS) or command prompt (Windows)
2. Type in `ssh username@c@int.bioc600.calculquebec.cloud`
3. Enter your password

### Using JupyterHub via a web browser

1. Click this link: [https://jupyter.bioc600.calculquebec.cloud/hub/login?next=%2Fhub%2F](https://jupyter.bioc600.calculquebec.cloud/hub/login?next=%2Fhub%2F)
2. Enter your account credentials

## Bash and Unix
The **<ins>B</ins>ourne <ins>A</ins>gain <ins>S</ins>hell** (Bash) is a Unix shell that allows you to work with the operating system of a computer. HPC clusters run almost exclusively on Linux (which is similar to Unix). So, to be able to work with HPC computers, you will need to learn some Bash.

Here is a cheat sheet of common Bash commands and what they do:

| Command | Function |
| ------- | -------- |
| mkdir | Make a new directory |
| mkdir -p | Make a new directory including parent directory |
| rm | Delete file |
| rm -r | Delete directory including contents |
| cd | Changes the working directory |
| cd .. | Go back one directory |
| ls | List the contents of a directory |
| ls -all | List all the contents of a directory + info |
| nano | Open a new or existing file and edit it |
| less | Open a file but you can't edit it |
| cp | Copy paste a file |
| mv | Move a file from one directory to another |
| wget | Download a file from a web link |
| gzip | Compress a single file |
| gunzip | Decompress a single file |
| chmod u+x | Give the user execute permission |
| chmod +x | Give everybody execute permission |
| which | Use to find what version software is loaded |
| grep | Find in files |
| sed | Replace in files |
| bash | Run a bash script |
| sbatch | Submit a script on a job scheduler |
| module avail | See what software is available to use |
| module load | Load software |
| module purge | Unloads all software |

For a deeper dive into bash, you can read the following resources:
- [Bash Tutorial](https://www.w3schools.com/bash/index.php)
- [Bash cheat sheet: Top 25 commands and creating custom commands](https://www.educative.io/blog/bash-shell-command-cheat-sheet)
