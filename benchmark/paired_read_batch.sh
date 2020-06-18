#!/bin/sh
#
CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

help() {
    echo "usage: `basename $0` -e <tool> -p <dirname> -r <datadir> -x <insert_size_list> -c <coverage_list> -d <sigma_list> -l <read_length_list>"
    exit $1
}

coverage_list=50
sigma_list=100
insert_size_list=1000
read_len_list=150

while getopts 'e:p:r:x:c:d:l:h' OPT; do
    case $OPT in
        e) tool="${OPTARG}";;
        p) dirname="${OPTARG}";;
        r) datadir="${OPTARG}";;
        x) insert_size_list="${OPTARG}";;
        c) coverage_list="${OPTARG}";;
        d) sigma_list="${OPTARG}";;
        l) read_len_list="${OPTARG}";;
        h) help 0;;
        ?) help 1;;
    esac
done

if [ -z ${tool} ]; then
    help 1
fi
if [ -z ${dirname} ]; then
    help 1
fi
if [ -z ${datadir} ]; then
    help 1
fi

for p in $(ls ${datadir}); do
    mkdir -p ${datadir}/${p}
    for insert_size in "${insert_size_list}"; do
        for coverage in "${coverage_list}"; do
            for sigma in "${sigma_list}"; do
                for read_len in "${read_len_list}"; do
                    sh ${tool} -p ${dirname}/${p} -r ${datadir}/${p} -c ${coverage} -x ${insert_size} -d ${sigma} -l ${read_len}
                done
            done
        done
    done
done
