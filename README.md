# COMP383 Computational Biology MiniProject Assignment

## Needed to run python script
### Languages:
Both Python3 and R must be installed for use in terminal
### Command Line Programs:
You will need Kallisto, Sleuth, Bowtie2, SPAdes, and BLAST+. Download links below.
- [Download Kallisto](https://pachterlab.github.io/kallisto/download)
- [Download Sleuth](https://pachterlab.github.io/sleuth/download)
- [Download Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Download SPAdes](https://cab.spbu.ru/software/spades/)
- [Download BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

## How to use the code:
Upon cloning the repo, simply use python3 to run the script with the following command in the linux terminal:

`python3 miniProject_amazher.py`

You will be prompted to enter either ***test*** or ***full*** depending on whether you want to run the test dataset or the full dataset (respectively). The test dataset should finish in a couple minutes; the full dataset will take 2-3 hours. *Note*: if something other than ***test*** or ***full*** is entered, the program will terminate and an error will be written to the generated `miniProject.log file`. If you wish to rerun, delete the created `/miniProject_Abdullah_Mazher` directory before trying again.

As the script runs, it will generate all output and downloaded files inside a directory named `/miniProject_Abdullah_Mazher`. The user is not required to move anything to or from this directory - everything will be handled by the code. This includes the test data - when the repo is cloned, the test data is a directory named `/test_data` and the script will automatically access that directory and copy files to the generated `/miniProject_Abdullah_Mazher` directory as needed.
