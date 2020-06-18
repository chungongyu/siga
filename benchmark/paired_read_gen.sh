#!/bin/sh
#
CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

coverage=50
sigma=100
insert_size=1000
read_len=150

help() {
    echo "usage: `basename $0` -p <dirname> -x <insert_size> -c <coverage> -d <sigma> -l <read_length>"
    exit $1
}

while getopts 'p:x:c:d:l:h' OPT; do
    case ${OPT} in
        p) dirname="${OPTARG}";;
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

echo "-----------------------------"
echo "dirname: ${dirname}"
echo "-----------------------------"
mkdir -p ${dirname}

python ${CWD}/paired_read_gen.py ${dirname}/new_ref.fa ${read_len} ${coverage} ${insert_size} ${sigma} | python ${CWD}/fasta.py make_file ${dirname}/read_${read_len}_${coverage}_${insert_size}_${sigma}_paired_R1.fasta ${dirname}/read_${read_len}_${coverage}_${insert_size}_${sigma}_paired_R2.fasta
