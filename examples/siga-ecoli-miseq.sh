#! /bin/bash -x

#
# Example assembly of 150bp E. coli reads
#

# Download the E. coli FASTQ files from Illumina's website
wget ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz
wget ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz

IN1=MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz
IN2=MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz

#
# Parameters
#

# Program paths
if [ -z ${siga_main} ]; then
    siga_main="${HOME}/siga/src/siga"
fi

# The number of threads to use
CPU=8

# Correction k-mer 
CORRECTION_K=41

# The minimum overlap to use when computing the graph.
# The final assembly can be performed with this overlap or greater
MIN_OVERLAP=85

# The overlap value to use for the final assembly
ASSEMBLE_OVERLAP=111

# Branch trim length
TRIM_LENGTH=150

#
# Dependency checks
#

# Check the required programs are installed and executable
prog_list="${siga_main}"
for prog in $prog_list; do
    hash $prog 2>/dev/null || { echo "Error $prog not found. Please place $prog on your PATH or update the *_BIN variables in this script"; exit 1; }
done 

# Check the files are found
file_list="$IN1 $IN2"
for input in $file_list; do
    if [ ! -f $input ]; then
        echo "Error input file $input not found"; exit 1;
    fi
done

#
# Preprocessing
#

# Preprocess the data to remove ambiguous basecalls
${siga_main} preprocess --pe-mode=1 -o reads.pp.fastq $IN1 $IN2

#
# Error Correction
#

# Build the index that will be used for error correction
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
${siga_main} index -t $CPU --no-reverse reads.pp.fastq

# Perform k-mer based error correction.
# The k-mer cutoff parameter is learned automatically.
${siga_main} correct -k $CORRECTION_K --learn -t $CPU -o reads.ec.fastq reads.pp.fastq

#
# Primary (contig) assembly
#

# Index the corrected data.
${siga_main} index -t $CPU reads.ec.fastq

# Compute the structure of the string graph
${siga_main} overlap -m $MIN_OVERLAP -t $CPU reads.ec.filter.pass.fa

# Perform the contig assembly
${siga_main} assemble -m $ASSEMBLE_OVERLAP --min-branch-length $TRIM_LENGTH -o primary reads.ec.filter.pass.asqg.gz
