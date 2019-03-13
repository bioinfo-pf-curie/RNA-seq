#!/bin/bash

##
## This script writes the expected list of file at the beginning of the pipeline
## This list can then be checked at the end of the process
## 

function usage {
    echo -e "usage : makeOutputList.sh -c CONFIG -p PREFIX -o ODIR"
    echo -e "Use option -h|--help for more information"
}

while getopts "c:p:o:hv" OPT
do
    case $OPT in
        c) CONF=$OPTARG;;
	p) PREFIX=$OPTARG;;
        o) ODIR=$OPTARG;;
        h) usage ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
             exit 1
             ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done


if [[ -z $CONF  || -z ${PREFIX} || -z $ODIR ]]; then
    usage
    exit
fi

## Load config file                                                                                                                                                                                        
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/utils.inc.sh
read_config $CONF

## Output
ofile=${ODIR}/logs/output_list.tsv
prefix_dir=$(basename $ODIR)
if [ -e ${ofile} ]; then
    /bin/rm -f ${ofile}
fi

## Mapping
if [ ! -z ${BOWTIE_RRNA_INDEX} ]; then
    prefix_bam=${PREFIX}_norRNA
else
    prefix_bam=${PREFIX}
fi
echo -e "${prefix_dir}/mapping/${prefix_bam}_markdup.bam" >> ${ofile}

## Preseq
echo -e "${prefix_dir}/preseq/${prefix_bam}.preseq" >> ${ofile}

## Stats
echo -e "${prefix_dir}/stats/${PREFIX}.stats" >> ${ofile}

## Counts
echo -e "${prefix_dir}/counts/${prefix_bam}_counts.csv" >> ${ofile}

## Strand_check
if [ -z ${STRANDED} ]; then
    echo -e "${prefix_dir}/strand_check/${PREFIX}.rseqc" >> ${ofile}
fi
