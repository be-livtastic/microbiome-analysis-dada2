# Install SRA toolkit first, then:
# Create a directory for the project
mkdir -p ~/microbiome_practice/raw_data
cd ~/microbiome_practice/raw_data

# Download the SRA run info to see all samples
# Go to: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA762524
# Click on "SRA Experiments" to see the list of samples and their SRA run accessions (SRR numbers)


sudo apt install sra-toolkit # Install the SRA toolkit

# Download a few samples first to test pipeline 
prefetch SRR15862371 
fastq-dump --split-files SRR15862371 #

prefetch SRR15862375
fastq-dump --split-files SRR15862375

prefetch SRR15862387
fastq-dump --split-files SRR15862387

prefetch SRR15862389
fastq-dump --split-files SRR15862389

# For paired-end data, this creates:
# SRRxxxxxxx_1.fastq (forward reads)
# SRRxxxxxxx_2.fastq (reverse reads)
#NCBI SRA merges the paired-end reads into a single file, so we use --split-files to separate them back into forward and reverse reads.
