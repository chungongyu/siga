#!/bin/sh
#
CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

########################
# env
########################
#export siga_main="${HOME}/siga/src/siga"
#export siga_log4cxx="${HOME}/siga/src/log4cxx.properties"
if [ -z ${siga_main} ]; then
    siga_main="${HOME}/siga/src/siga"
fi
if [ -z ${siga_log4cxx} ]; then
    siga_log4cxx="${HOME}/siga/src/log4cxx.properties"
fi

coverage=50
sigma=150
insert_size=1000
read_len=150

help() {
    echo "usage: `basename $0` -p <dirname> -r <datadir> -x <insert_size> -c <coverage> -d <sigma> -l <read_length>"
    exit $1
}

while getopts 'p:r:x:c:d:l:h' OPT; do
    case $OPT in
        p) dirname="${OPTARG}";;
        r) datadir="${OPTARG}";;
        x) insert_size="${OPTARG}";;
        c) coverage="${OPTARG}";;
        d) sigma="${OPTARG}";;
        l) read_len="${OPTARG}";;
        h) help 0;;
        ?) help 1;;
    esac
done

if [ -z ${dirname} ]; then
    help 1
fi
if [ -z ${datadir} ]; then
    help 1
fi

echo "-----------------------------"
echo "dirname: ${dirname}"
echo "datadir: ${datadir}"
echo "-----------------------------"

R1="${datadir}/read_${read_len}_${coverage}_${insert_size}_${sigma}_paired_R1"
R2="${datadir}/read_${read_len}_${coverage}_${insert_size}_${sigma}_paired_R2"
f="${dirname}/read_${read_len}_${coverage}_${insert_size}_${sigma}_paired_R"

mkdir -p ${dirname}

########################
#
# Assemble paired read
#
########################

${siga_main} preprocess -c ${siga_log4cxx} --pe-mode=1 --pe-orientation=ff --no-primer-check -o ${f}.fastq ${R1}.fasta ${R2}.fasta
${siga_main} index -c ${siga_log4cxx} -t 1 -p ${f} "${f}.fastq"
${siga_main} overlap -c ${siga_log4cxx} -t 16 --batch-size=10000 -m 100 --no-opposite-strand -p ${f} "${f}.fastq"
#${siga_main} overlap -c ${siga_log4cxx} -t 16 --batch-size=10000 -m 100 -p ${f} "${f}.fastq"
${siga_main} assemble -c ${siga_log4cxx} -t 16 --batch-size=10000 -m 100 --pe-mode=1 --max-distance=100 -p ${f} ${f}.asqg.gz

########################
#
# Assemble simple contigs
#
########################
${siga_main} index -c ${siga_log4cxx} -t 1 -p "${f}-contigs" "${f}-contigs.fa"
${siga_main} rmdup -c ${siga_log4cxx} -t 8 -p "${f}-contigs" "${f}-contigs.fa"
${siga_main} index -c ${siga_log4cxx} -t 1 -p "${f}-contigs.rmdup" "${f}-contigs.rmdup.fa"
${siga_main} overlap -c ${siga_log4cxx} -t 8 -m 10 --no-opposite-strand -p "${f}-contigs.rmdup" "${f}-contigs.rmdup.fa"
#${siga_main} overlap -c ${siga_log4cxx} -t 8 -m 10 -p "${f}-contigs.rmdup" "${f}-contigs.rmdup.fa"
${siga_main} assemble -c ${siga_log4cxx} --pe-mode=0 -m 100 -p "${f}-contigs.rmdup" "${f}-contigs.rmdup.asqg.gz"

########################
#
# Evaluate
#
########################
cat ${f}-contigs.rmdup-contigs.fa | python ${CWD}/contigs_mapping.py 300 ${datadir}/new_ref.fa fasta ${dirname}/unmatched_contigs > ${dirname}/siga-contigs_${insert_size}.stats

########################
#
# Draw graph
#
########################
zcat ${f}-contigs.rmdup-graph.asqg.gz | awk -f ${CWD}/graphviz.awk | dot -Tjpg -o ${f}-contigs.rmdup.jpg

cat ${dirname}/siga-contigs_${insert_size}.stats
