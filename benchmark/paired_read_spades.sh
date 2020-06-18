#!/bin/sh
#
CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

########################
# env
########################
#
#export spades_main="python ${HOME}/SPAdes-3.6.2/spades.py"
#
if [ -z ${spades_main} ]; then
    spades_main="python ${HOME}/SPAdes-3.6.2/spades.py"
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
        p) dirname="$OPTARG";;
        r) datadir="$OPTARG";;
        x) insert_size="$OPTARG";;
        c) coverage="$OPTARG";;
        d) sigma="$OPTARG";;
        l) read_len="$OPTARG";;
        h) help 0;;
        ?) help 1;;
    esac
done

if [ -z $dirname ]; then
    help 1
fi
if [ -z $datadir ]; then
    help 1
fi

echo "-----------------------------"
echo "dirname: ${dirname}"
echo "datadir: ${datadir}"
echo "-----------------------------"

R1="${datadir}/read_${read_len}_${coverage}_${insert_size}_${sigma}_paired_R1"
R2="${datadir}/read_${read_len}_${coverage}_${insert_size}_${sigma}_paired_R2"

########################
#
# Assemble paired read
#
########################
${spades_main} --careful --only-assembler -k 21,27,33,55,77,99,127 --pe1-ff -o ${dirname} -1 ${R1}.fasta -2 ${R2}.fasta

########################
#
# Evaluate
#
########################
cat ${dirname}/contigs.fasta | python ${CWD}/contigs_mapping.py 300 ${datadir}/new_ref.fa fasta ${dirname}/unmatched_contigs > ${dirname}/spades-contigs_${insert_size}.stats
