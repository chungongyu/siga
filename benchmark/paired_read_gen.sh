#!/bin/sh
#

COVERAGE=50
SIGMA=100
INSERT_SIZE=1000
READ_LENGTH=150

help() {
    echo "usage: `basename $0` -p <dirname> -x <insert_size> -c <coverage> -d <sigma> -l <read_length>"
    exit $1
}

while getopts 'p:x:c:d:l:h' OPT; do
    case $OPT in
        p) DIRNAME="$OPTARG";;
        x) INSERT_SIZE="$OPTARG";;
        c) COVERAGE="$OPTARG";;
        d) SIGMA="$OPTARG";;
        l) READ_LENGTH="$OPTARG";;
        h) help 0;;
        ?) help 1;;
    esac
done

if [ -z $DIRNAME ]; then
    help 1
fi

echo ${DIRNAME}
mkdir -p ${DIRNAME}

python paired_read_gen.py ${DIRNAME}/new_ref.fa ${READ_LENGTH} ${COVERAGE} ${INSERT_SIZE} ${SIGMA} | python fasta.py make_file ${DIRNAME}/read_${READ_LENGTH}_${COVERAGE}_${INSERT_SIZE}_${SIGMA}_paired_R1.fasta ${DIRNAME}/read_${READ_LENGTH}_${COVERAGE}_${INSERT_SIZE}_${SIGMA}_paired_R2.fasta
