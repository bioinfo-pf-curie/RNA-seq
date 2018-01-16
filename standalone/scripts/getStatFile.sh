## getStatFile.sh
##
## Copyright (c) 2017 Institut Curie                               
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

## Template of stat file for NGS data analysis
## sample_ID;bioID;#cluster;%good_qual_bases_R1;%good_qual_bases_R2;#reads;
## #reads_on_rRNA;#reads_not_on_rRNA;unique_hits;unique_hits_on_trs;intronic_hits
## #properly_paired;non_dup;mean_depth_coverage

VERSION=0.0.1

function usage {
    echo -e "usage : ./getStatFile.sh -f FORWARD -b BAM -c CONFIG [-r REVERSE] [-x R_RNA_BAM] [-g GTF] [-s SAMPLE_ID]"
    echo -e "Use option -h|--help for more information"
}

function version {
    echo -e "$SOFT version $VERSION"
    exit
}


function help {
    usage;
    echo
    echo "getStatFile.sh $VERSION"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -f INPUT: input forward fastq file"
    echo "   -b BAM: bam file of aligned reads"
    echo "   -c CONFIG: configuration file for RNA processing"
    echo "   [-r INPUT]: input reverse fastq file"
    echo "   [-x R_RNA_BAM]: bam file of rRNA mapping"
    echo "   [-g GTF]: transcriptome  gtf file"
    echo "   [-s SAMPLE_ID]: biosample ID"
    echo "   [-h]: help"
    echo "   [-v]: version"
    exit;
}

while getopts "f:b:c:r:x:g:s:hv" OPT
do
    case $OPT in
        f) FORWARD=$OPTARG;;
	b) BAM=$OPTARG;;
	c) CONF=$OPTARG;;
	r) REVERSE=$OPTARG;;
        x) R_RNA_BAM=$OPTARG;;
        g) GTF=$OPTARG;;
	s) SAMPLE_ID=$OPTARG;;
        v) version ;;
        h) help ;;
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


if [[ -z $CONF  || -z $FORWARD || -z $BAM ]]; then
    usage
    exit
fi

## Load config file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/utils.inc.sh
read_config $CONF

## sampleID
if [ ! -z ${REVERSE} ]; then
    ## Get common part of R1/R2 + remove norRNA suffix + remove .R(1/2) suffix
    fastqID=$(basename $(printf "%s\n%s\n" "${FORWARD}" "${REVERSE}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' | sed -e 's/[\._]*$//'))
else
    fastqID=$(basename ${FORWARD} | sed -e 's/.fastq\(.gz\)//' | sed -e 's/[\_.]*R*[12]*$//')
fi

## #cluster
if [[ ${FORWARD} =~ "gz" ]];then
    nb_cluster=$(($(zcat ${FORWARD} | wc -l)/4))
elif [[ ${FORWARD} =~ "fastq" ]];then
    nb_cluster=$(($(wc -l < ${FORWARD})/4))
else
    die "ERROR : Wrong file type in input for fastq file; file: ${FORWARD}"
fi
## #reads
if [ ! -z ${REVERSE} ]; then
    if [[ ${REVERSE} =~ "gz" ]];then
        nb_reverse=$(($(zcat ${REVERSE} | wc -l)/4))
	elif [[ ${REVERSE} =~ "fastq" ]];then
		nb_reverse=$(($(wc -l < ${REVERSE})/4))
	else
		die "ERROR : Wrong file type in input for fastq file; file: ${REVERSE}"  
    fi 
else
    nb_reverse=0
fi

nb_reads=$((${nb_cluster} + ${nb_reverse}))

## rRNA
if [ ! -z $R_RNA_BAM ]; then
    nb_rrna=$(${SAMTOOLS_PATH}/samtools view -F 4 -c $R_RNA_BAM)
    nb_no_rrna=$(($nb_reads - $nb_rrna))
else
    nb_rrna=NA
    nb_no_rrna=$nb_reads
fi

#################
##
## Mapping Stats
##
#################
UNIQUE_BAM=$(echo $BAM | sed -e 's/.bam$/_unique.bam/')

if [ ${MAPPING_TOOL} == "TOPHAT2" ]; then 
    
    #https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
    aligned=$(${SAMTOOLS_PATH}/samtools view -f 0x40 -c ${BAM})
  
    ${SAMTOOLS_PATH}/samtools view -bq50  ${BAM} > ${UNIQUE_BAM}
    ubam=$(${SAMTOOLS_PATH}/samtools view -c ${UNIQUE_BAM})

    mbam=$(${SAMTOOLS_PATH}/samtools view -F 0x100  ${BAM} | awk '($5!=50){print $1}' | sort -u | wc -l)
  
elif [ ${MAPPING_TOOL} == "STAR" ]; then

    aligned=$(${SAMTOOLS_PATH}/samtools view -F 0x100 -c ${BAM})
  
    ## STAR
    #Proper Paired aligns:
    #N12=$(${SAMTOOLS_PATH}/samtools view -q255 -c -f 0x2 $BAM)
    #End 1 aligns:
    #N1=$(${SAMTOOLS_PATH}/samtools view -q255 -c -F 0x2 -f 0x40 $BAM)
    #End 2 aligns:
    #N2=$(${SAMTOOLS_PATH}/samtools view -q255 -c -F 0x2 -f 0x80 $BAM)
    #The total number of uniquely aligned read is then:
    #N12/2+N1+N2 which should exactly agree with the number in the Log.final.out
    #ubam=$(($N12/2 + $N1 + $N2))

    ${SAMTOOLS_PATH}/samtools view -bq255  ${BAM} > ${UNIQUE_BAM}
    ubam=$( ${SAMTOOLS_PATH}/samtools view -c ${UNIQUE_BAM})

    ## Multiple Hits
    mbam=$(${SAMTOOLS_PATH}/samtools view -F 0x100  ${BAM} | awk '($5!=255){print $1}' | sort -u | wc -l)

    if [ ! -z ${REVERSE} ]; then
	mbam=$(( $mbam * 2 ))
    fi
fi

## duplicates
MDUP_BAM=$(echo $BAM | sed -e 's/.bam$/_markdup.bam/')
ndup=NA
if [ -e ${MDUP_BAM} ]; then
    ndup=$(${SAMTOOLS_PATH}/samtools view -f 1024 ${MDUP_BAM} | wc -l)
fi


##############
##
## Annotations
##
##############
if  [ ! -z $GTF ]; then
    EXON_BED=$(echo $GTF | sed -e 's/.gtf$/_exon.bed/')
    INTRON_BED=$(echo $GTF | sed -e 's/.gtf$/_intron.bed/')
    INTER_BED=$(echo $GTF | sed -e 's/.gtf$/_inter.bed/')
    
    ## hits_on_trs
    ## overlap the exons (at least 1-base)
    if [ -e ${EXON_BED} ]; then
	nb_exon=$(${BEDTOOLS_PATH}/intersectBed -a ${BAM} -b ${EXON_BED} | ${SAMTOOLS_PATH}/samtools view -c -)
    fi
    ## overlap the introns but not the exons as they are already counted before
    if [ -e ${INTRON_BED} ]; then
	nb_intron=$(${BEDTOOLS_PATH}/intersectBed -a ${BAM} -b ${INTRON_BED} | ${BEDTOOLS_PATH}/intersectBed -a stdin -b ${EXON_BED} -v | ${SAMTOOLS_PATH}/samtools view -c -)
    fi
    ## are fully intergenic
    if [ -e ${INTER_BED} ]; then
	nb_inter=$(${BEDTOOLS_PATH}/intersectBed -a ${BAM} -b ${INTER_BED} -f 1 | ${SAMTOOLS_PATH}/samtools view -c -)
    fi
else
    nb_exon=NA
    nb_intron=NA
    nb_inter=NA
fi

## verbose
echo -e "Sample_identifier,Biological_identifier,Number_of_cluster,Number_of_reads,Number_of_rRNA_reads,Number_of_non_rRNA_reads,Number_of_aligned_reads,Number_of_uniquely_aligned_reads,Number_of_multiple_aligned_reads,Number_of_duplicates,Number_of_exonic_reads,Number_of intronic_reads,Number_of_intergenics_reads"
echo -e ${fastqID},${SAMPLE_ID},${nb_cluster},${nb_reads},${nb_rrna},${nb_no_rrna},${aligned},${ubam},${mbam},${ndup},${nb_exon},${nb_intron},${nb_inter}



