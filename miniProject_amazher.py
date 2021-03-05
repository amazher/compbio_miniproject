''' Python wrapper for CompBio 383 mini-project. This script will automate the retrieval, assembly, and analysis of HCMV transcriptomes
    using programs such as kallisto, sleuth, Bowtie2, SPAdes, and blast in the command line/terminal. '''

import os #this module will give access to OS commands within the script
import shutil #will use to move and copy test files
import sys #will use to exit program if user enters bad arg
from Bio import Entrez #use to access genbank
from Bio import SeqIO #use to read entrez handle and parse files
from Bio.Blast import NCBIWWW #use to access BLAST

# create empty dir to store all output, move into it, create empty log file
os.system('mkdir miniProject_Abdullah_Mazher')
os.chdir('miniProject_Abdullah_Mazher')
os.system('touch miniProject.log')

# prompt user to enter whether they would like to run test data or the full data
prompt = input('Which set of data would you like to run?\nPlease enter either:\n' +
               '"test" to run the test dataset, or\n"full" to run the full dataset.\n')

# if the user wants to run the test dataset, move the test data from the test_data folder to the miniProject_Abdullah_Mazher folder, and
# rename it to what the files would have been called if they were full dataset files (from fastq-dump) so the name is consistent below
if prompt == 'test':
    # go back up one directory (so script can see test_data folder)
    os.chdir('..') 
    # use shutil to copy the testfiles to the miniproject dir and rename the testfiles as if they were the result of fastq-dump of full files
    shutil.copyfile('test_files/SRR5660030.1_1_test.fastq', 'miniProject_Abdullah_Mazher/SRR5660030.1_1.fastq')
    shutil.copyfile('test_files/SRR5660030.1_2_test.fastq', 'miniProject_Abdullah_Mazher/SRR5660030.1_2.fastq')
    shutil.copyfile('test_files/SRR5660033.1_1_test.fastq', 'miniProject_Abdullah_Mazher/SRR5660033.1_1.fastq')
    shutil.copyfile('test_files/SRR5660033.1_2_test.fastq', 'miniProject_Abdullah_Mazher/SRR5660033.1_2.fastq')
    shutil.copyfile('test_files/SRR5660044.1_1_test.fastq', 'miniProject_Abdullah_Mazher/SRR5660044.1_1.fastq')
    shutil.copyfile('test_files/SRR5660044.1_2_test.fastq', 'miniProject_Abdullah_Mazher/SRR5660044.1_2.fastq')
    shutil.copyfile('test_files/SRR5660045.1_1_test.fastq', 'miniProject_Abdullah_Mazher/SRR5660045.1_1.fastq')
    shutil.copyfile('test_files/SRR5660045.1_2_test.fastq', 'miniProject_Abdullah_Mazher/SRR5660045.1_2.fastq')
    # go back to the miniProject folder and continue with script below
    os.chdir('miniProject_Abdullah_Mazher')

# if the user wants to run the full dataset, download full data from ncbi and continue with script
elif prompt == 'full':
    # retrieve transcriptomes from patient donors from SRA using wget
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1')
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1')
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1')
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1')

    # convert retrieved transcriptomes to paired-end fastq files using fastq-dump
    os.system('fastq-dump -I --split-files SRR5660030.1')
    os.system('fastq-dump -I --split-files SRR5660033.1')
    os.system('fastq-dump -I --split-files SRR5660044.1')
    os.system('fastq-dump -I --split-files SRR5660045.1')

# if the user enters something other than 'test' or 'full', exit the program
else:
    log_file = open('miniProject.log', 'w')
    # also write an error message to the log
    log_file.write('Program terminated due to unrecognized arg. Please enter either "test" or "full" when prompted in future runs.\n' +
                   'Delete the miniProject_Abdullah_Mazher directory before retrying.')
    log_file.close()
    sys.exit('Error: Program Terminated. Please only enter either "test" or "full".') #end script

''' 2 '''
cds_file = open('HCMV_transcriptome.fasta', 'w') #will write all CDS to this file
log_file = open('miniProject.log', 'w') #log file

# use biopython (Entrez efetch) to retrieve genome of HCMV (EF999921) in genbank format from NCBI
Entrez.email = 'amazher@luc.edu'
handle = Entrez.efetch(db='nucleotide', id='EF999921', rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')

# count number of CDS and write to log, append all CDS in fasta format to cds_file
CDS_count = 0 #counter for CDS
for feature in record.features: #check every feature to see if it is a CDS
    if feature.type == 'CDS':
        CDS_count += 1 #increment CDS counter
        #write the CDS to cds_file in fasta format (header has protein id, the record id, and the count of the CDS)
        cds_file.write('>' + feature.qualifiers['protein_id'][0] + '_' + str(record.id) + '_CDS' + str(CDS_count) + '\n' +
                       str(feature.extract(record.seq)) + '\n')
log_file.write('The HCMV genome (EF999921) has ' + str(CDS_count) + ' CDS.' + '\n\n') #write cds count to log
cds_file.close()

# build kallisto index using HCMV_transcriptome.fasta
os.system('time kallisto index -i HCMV_index.idx HCMV_transcriptome.fasta')

''' 3 '''
# quantify TPM of each CDS in transcriptome using kallisto
os.system('mkdir kallisto_results') #dir to store kallisto results
kallisto_short = 'time kallisto quant -i HCMV_index.idx -b 30 -t 2 -o kallisto_results/' #to shorten command below:
os.system(kallisto_short + 'SRR5660030.1 SRR5660030.1_1.fastq SRR5660030.1_2.fastq')
os.system(kallisto_short + 'SRR5660033.1 SRR5660033.1_1.fastq SRR5660033.1_2.fastq')
os.system(kallisto_short + 'SRR5660044.1 SRR5660044.1_1.fastq SRR5660044.1_2.fastq')
os.system(kallisto_short + 'SRR5660045.1 SRR5660045.1_1.fastq SRR5660045.1_2.fastq')

# make file containing table of sample/condition/path for sleuth to use
sleuth_table = open('sleuth_table.txt', 'w')
sleuth_table.write('sample condition path\n')
sleuth_table.write('SRR5660030.1 2dpi kallisto_results/SRR5660030.1\n')
sleuth_table.write('SRR5660033.1 6dpi kallisto_results/SRR5660033.1\n')
sleuth_table.write('SRR5660044.1 2dpi kallisto_results/SRR5660044.1\n')
sleuth_table.write('SRR5660045.1 6dpi kallisto_results/SRR5660045.1\n')
sleuth_table.close()

# make short R script which will be used to find differentially expressed genes between the 2 timepoints using sleuth
sleuth_script = open('sleuth_script.R', 'w')
#noteToSelf: each '' ends with \n to create new line, + is visual marker of diff line
sleuth_script.write('library(sleuth)\n' + 'table <- read.table("sleuth_table.txt", header=TRUE, stringsAsFactors=FALSE)\n' +
                    'so <- sleuth_prep(table)\n' + 'so <- sleuth_fit(so, ~condition, "full")\n' +
                    'so <- sleuth_fit(so, ~1, "reduced")\n' + 'so <- sleuth_lrt(so, "reduced", "full")\n' +
                    'library(dplyr)\n' + 's_table <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE)\n' +
                    'sleuth_significant <- dplyr::filter(s_table, qval <= 0.05) %>% dplyr::arrange(pval)\n' +
                    'write.table(sleuth_significant, file="fdr05_results.txt", quote=FALSE, row.names=FALSE)\n')
sleuth_script.close()

# run the sleuth R script in terminal using Rscript:
os.system('Rscript sleuth_script.R')

# read in the output from sleuth, parse wanted info and write to log
sleuth_output = open('fdr05_results.txt', 'r')
lines_list = [] #will store each line as a list (split at ' ') inside this list
for line in sleuth_output.readlines():
    lines_list.append(line.split(' ')) #for every line, make a list split at empty space

#write the headers (from the 0th index in lines_list):
log_file.write(lines_list[0][0] +'\t'+ lines_list[0][3] +'\t'+ lines_list[0][1] +'\t'+ lines_list[0][2] +'\n')
#loop through rest of lines from lines_list and write targetid(0), teststat(3), pval(1), qval(2)
for i in range(1, len(lines_list)):
    log_file.write(lines_list[i][0] + '\t' + lines_list[i][3] + '\t' + lines_list[i][1] + '\t' + lines_list[i][2] + '\n')
log_file.write('\n')
sleuth_output.close()

''' 4 '''
# use entrez efetch to get the full HCMV genome sequence to make bowtie index
handle = Entrez.efetch(db='nucleotide', id='EF999921', rettype='fasta', retmode='text')
record = SeqIO.read(handle, 'fasta')

HCMV_complete = open('HCMV_complete_seq.fasta', 'w') #will write complete seq to this file
HCMV_complete.write('>' + record.description + '\n' + str(record.seq)) #write the genome seq in fasta format to file
HCMV_complete.close()

# make the bowtie index using the previously created fasta file (HCMV_complete_seq.fasta)
os.system('bowtie2-build HCMV_complete_seq.fasta HCMV_genome')

# run bowtie on each fastq pair, output only reads that align at least once (--al-conc)
bowtie = 'bowtie2 --quiet -x HCMV_genome '
os.system(bowtie+'-1 SRR5660030.1_1.fastq -2 SRR5660030.1_2.fastq -S SRR5660030.1_map.sam --al-conc SRR5660030.1_mapped.fastq')
os.system(bowtie+'-1 SRR5660033.1_1.fastq -2 SRR5660033.1_2.fastq -S SRR5660033.1_map.sam --al-conc SRR5660033.1_mapped.fastq')
os.system(bowtie+'-1 SRR5660044.1_1.fastq -2 SRR5660044.1_2.fastq -S SRR5660044.1_map.sam --al-conc SRR5660044.1_mapped.fastq')
os.system(bowtie+'-1 SRR5660045.1_1.fastq -2 SRR5660045.1_2.fastq -S SRR5660045.1_map.sam --al-conc SRR5660045.1_mapped.fastq')

#first 2 are lists that just hold the file names for use in loop below, third list holds the donor associated with file
before_files = ['SRR5660030.1_1.fastq', 'SRR5660033.1_1.fastq', 'SRR5660044.1_1.fastq', 'SRR5660045.1_1.fastq']
after_files = ['SRR5660030.1_mapped.1.fastq', 'SRR5660033.1_mapped.1.fastq', 'SRR5660044.1_mapped.1.fastq', 'SRR5660045.1_mapped.1.fastq']
donors_list = ['Donor 1 (2dpi)', 'Donor 1 (6dpi)', 'Donor 3 (2dpi)', 'Donor 3 (6dpi)']

# write the number of reads in each transcriptome before and after bowtie2 mapping to log using for loop and lists above
for i in range(len(before_files)):
    before_file = open(before_files[i], 'r')
    after_file = open(after_files[i], 'r')
    before_count = 0
    after_count = 0
    for record in SeqIO.parse(before_file, 'fastq'): #count each record (read) and update counter
        before_count += 1
    for record in SeqIO.parse(after_file, 'fastq'): #count each record (read) and update counter
        after_count += 1
    log_file.write(donors_list[i] + ' had ' + str(before_count) + ' read pairs before Bowtie2 filtering and ' + str(after_count) +
                   ' read pairs after.\n')

''' 5 '''
# use bowtie2 output reads to produce one assembly using SPAdes
spades_command = ('spades -k 127 -t 2 --only-assembler ' +
                  '--pe1-1 SRR5660030.1_mapped.1.fastq --pe1-2 SRR5660030.1_mapped.2.fastq ' +
                  '--pe2-1 SRR5660033.1_mapped.1.fastq --pe2-2 SRR5660033.1_mapped.2.fastq ' +
                  '--pe3-1 SRR5660044.1_mapped.1.fastq --pe3-2 SRR5660044.1_mapped.2.fastq ' +
                  '--pe4-1 SRR5660045.1_mapped.1.fastq --pe4-2 SRR5660045.1_mapped.2.fastq ' +
                  '-o spades_assembly/')
os.system(spades_command)

# write the SPAdes command used to the log file
log_file.write('\nSPAdes command used for assembly:\n' + spades_command + '\n')

''' 6 & 7 '''
# calculate the number of contigs with length > 1000; calculate total number of bp in all contigs > 1000bp; write both to log
# also, put all contigs > 1000 in list, will use to find longest one in part 8
contig_file = open('spades_assembly/contigs.fasta', 'r')

contig_list = [] #will hold contigs > 1000bp
num_contigs = 0
total_len = 0
for record in SeqIO.parse(contig_file, 'fasta'): #parse contig_file, count contigs>1000bp, len of all those contigs, and update contig_list
    if len(record.seq) > 1000:
        num_contigs += 1
        total_len += len(record.seq)
        contig_list.append(str(record.seq))

#write info to log
log_file.write('\nThere are ' + str(num_contigs) + ' contigs > 1000 bp in the assembly.\n')
log_file.write('\nThere are ' + str(total_len) + ' bp in the assembly.\n\n')

''' 8 '''
# retrieve longest contig from spades assembly and write to a file for use in blast
long_contig_file = open('longest_contig.fasta', 'w')
contig_list.sort(key=len, reverse=True) #sort list by size in descending order (largest = first)
long_contig_file.write('>longest_contig_spades\n' + contig_list[0])
long_contig_file.close()

# get sequences from the Betaherpesvirinae subfamily using biopython (Entrez esearch and efetch)
# first, use esearch to search for the virus tax with the refseq filter (change retmax bcs default is 20 which is too low)
handle = Entrez.esearch(db='nucleotide', term='txid10357[Organism:exp]', retmax=100000)
# read the handle w Entrez and extract list of ids that match the search term
record = Entrez.read(handle)
ids_list = record.get('IdList')
# next, use efetch and list of ids to get records from ncbi, parse with seqIO
handle = Entrez.efetch(db='nucleotide', id=ids_list, rettype='fasta')
records = list(SeqIO.parse(handle, 'fasta'))
# write the records to a file in fasta format
betaherpes_seqs = open('betaherpesvirinae_seqs.fasta', 'w')
for record in records:
    betaherpes_seqs.write(record.format('fasta'))
betaherpes_seqs.close()

# make local blast db from sequences of betaherpesvirinae subfamily retrieved into a file above
os.system('makeblastdb -in betaherpesvirinae_seqs.fasta -out HCMV_blastdb -title HCMV_blastdb -dbtype nucl')

# blastn the longest contig from spades assembly against the created db (oufmt 6 is tab)
os.system('blastn -query longest_contig.fasta -db HCMV_blastdb -out HCMV_blastn_results.txt ' +
          '-outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"')

# write the top 10 hits (first 10 hits for now) to the log file (\t)
b_results = open('HCMV_blastn_results.txt', 'r')

log_file.write('sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n') #write header, sep by tab
for i in range(10): #iterate 10 times to read 10 lines from blast output aka 10 hits
    log_file.write(str(b_results.readline())) #write the line to the logfile. already tab sep due to outfmt blast option
log_file.close()
