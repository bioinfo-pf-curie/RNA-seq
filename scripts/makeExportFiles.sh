#!/bin/bash

##
## This script is an internal script that is only used on the Institut Curie data system
##

print_export_persample()
{
    local ofile=$1
    kdi_bioapp="raw transcriptomics"

    if [ ! -z ${BOWTIE_RRNA_INDEX} ]; then
	## BAM files
	echo -e "backup/{SAMPLE}/mapping/{SAMPLE}_norRNA.bam\tanalysis_results/{SAMPLE}/mapping/{SAMPLE}.bam\tbam\tprocessed\t$kdi_bioapp" > ${ofile}
	echo -e "backup/{SAMPLE}/mapping/{SAMPLE}_norRNA.bam.bai\tanalysis_results/{SAMPLE}/mapping/{SAMPLE}.bam.bai\tbam\tprocessed\t$kdi_bioapp" >> ${ofile}
	
	## Counts
	echo -e "backup/{SAMPLE}/counts/{SAMPLE}_norRNA_counts.csv\tanalysis_results/{SAMPLE}/counts/{SAMPLE}_noRNA_counts.csv\ttxt\tprocessed\t$kdi_bioapp" >> ${ofile}

   else
	## BAM files
	echo -e "backup/{SAMPLE}/mapping/{SAMPLE}.bam\tanalysis_results/{SAMPLE}/mapping/{SAMPLE}.bam\tbam\tprocessed\t$kdi_bioapp" > ${ofile}
        echo -e "backup/{SAMPLE}/mapping/{SAMPLE}.bam.bai\tanalysis_results/{SAMPLE}/mapping/{SAMPLE}.bam.bai\tbam\tprocessed\t$kdi_bioapp" >> ${ofile}

	## Counts
	echo -e "backup/{SAMPLE}/counts/{SAMPLE}_counts.csv\tanalysis_results/{SAMPLE}/counts/{SAMPLE}_counts.csv\ttxt\tprocessed\t$kdi_bioapp" >> ${ofile}
    fi

    ## Strand check
    echo -e "backup/{SAMPLE}/strand_check/{SAMPLE}.rseqc\tanalysis_results/{SAMPLE}/strand_check/{SAMPLE}.rseqc\ttxt\tprocessed\t$kdi_bioapp" >> ${ofile}       
}

print_export_perrun()
{
    local ofile=$1
    local odir=$2
    
    outfiles=(report/report.html report/tablecounts_raw.csv report/tablecounts_tpm.csv)
    for item in ${outfiles[*]}
    do
	f=$(basename ${item})
	if [ -e ${odir}/${item} ]; then
	    echo -e "backup/${item}\treport/${f}\treport\tprocessed\tQuality Control" >> ${ofile}
	fi
    done
}

function usage {
    echo -e "usage : makeExportFile.sh -c CONFIG -o ODIR"
    echo -e "Use option -h|--help for more information"
}

while getopts "c:o:hv" OPT
do
    case $OPT in
        c) CONF=$OPTARG;;
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


if [[ -z $CONF  || -z $ODIR ]]; then
    usage
    exit
fi

## Load config file                                                                                                                                                                                        
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/utils.inc.sh
read_config $CONF

echo -e "Preparing KDI export files ..."
echo

mkdir -p ${ODIR}/export

## Merge list of output files
if [ -e ${ODIR}/export/output_list.tsv ]; then
    /bin/rm -f ${ODIR}/export/output_list.tsv
fi

find ${ODIR} -name "output_list.tsv" -not -path "*/export/*"  | while read olist; do
    cat $olist >> ${ODIR}/export/output_list.tsv
done

## Add report file
echo "report/report.html" >> ${ODIR}/export/output_list.tsv

## Export File for KDI
print_export_persample ${ODIR}/export/exportToKDI.conf
print_export_perrun ${ODIR}/export/exportToKDI-RunFiles.conf ${ODIR}

