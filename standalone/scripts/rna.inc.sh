set -o pipefail  # trace ERR through pipes
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value

## rna.inc.sh
##
## Copyright (c) 2017 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

## $1 = input file(s)
## $2 = output path
## $3 = log dir
fastqc_func()
{
    echo -e "Running FastQC ..."
    
    local out=$2/fastqc
    local log=$3/fastqc.log
    mkdir -p ${out}

    local cmd="${FASTQC_PATH}/fastqc -t 4 -j ${JAVA_PATH}/java -o ${out} $1 2> $log"
    exec_cmd ${cmd}
}


## $1 = input files
## $2 = output dir
## $3 = log dir
rRNA_mapping_func()
{
    echo -e "Running rRNA mapping ..."

    if [[ -z ${BOWTIE_RRNA_IDX} ]]; then
	die "rRNA indexes file not set. Exit"
    fi
    
    local log=$3/bowtie_rRNA.log
    local out=$2/rRNA_mapping
    mkdir -p ${out}

    #need to figure out --sam output ${bowtie_sam}
    bowtie_sam=${out}/$(basename $1 | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*/.sam/')
    bowtie_bam=${out}/$(basename $1 | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*/.bam/')

    inputs=($1)
    if [[ ${#inputs[@]} -eq 1 ]]; then
        if [[ ${inputs[0]} =~ \.gz ]]; then
            cmd_in="<(gzip -cd ${inputs[0]})"
        elif [[ ${inputs[0]} =~ \.fastq ]]; then
            cmd_in="${inputs[0]}"
        else
            die "ERROR : Wrong file type in input for Bowtie1; file: ${inputs[0]}"
        fi
        local unmapped=$(basename ${inputs[0]} | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*/_norRNA.fastq/')
        cmd_un="--un ${out}/${unmapped}"
    
    elif [[ ${#inputs[@]} -eq 2 ]]; then
    	if [[ ${inputs[0]} =~ \.gz ]]; then
            if [[ ${inputs[1]} =~ \.gz ]];then
                cmd_in="-1 <(gzip -cd ${inputs[0]}) -2 <(gzip -cd ${inputs[1]})"
            elif [[  ${inputs[1]} =~ \.fastq ]];then
                cmd_in="-1 <(gzip -cd ${inputs[0]}) -2 ${inputs[1]}"
            else
                die "ERROR : Wrong file type in input for Bowtie1; file: ${inputs[1]}"
            fi
        elif [[ ${inputs[0]} =~ \.fastq ]];then
            if [[ ${inputs[1]} =~ \.gz ]];then
                cmd_in="-1 ${inputs[0]} -2 <(gzip -cd ${inputs[1]})"
            elif [[ ${inputs[1]} =~ \.fastq ]];then
                cmd_in="-1 ${inputs[0]} -2 ${inputs[1]}"
            else
                die "ERROR : Wrong file type in input for Bowtie1; file: ${inputs[1]}"
            fi
        else
            die "ERROR : Wrong file type in input for Bowtie1; file: ${inputs[0]}"
        fi
        
	local unmapped=$(basename ${inputs[0]} | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*/_norRNA.fastq/')
        cmd_un="--un ${out}/${unmapped}"
	#local fq1_out="${out}/$(basename ${unmapped} | sed -e 's/.fastq/_1.fastq/')"
	#local fq2_out="${out}/$(basename ${unmapped} | sed -e 's/.fastq/_2.fastq/')"
    else
	die "Bowtie1 -  found more than two input files !"
    fi

    local cmd="${BOWTIE_PATH}/bowtie ${BOWTIE_OPTS} -p 4 ${cmd_un} --sam ${BOWTIE_RRNA_IDX} ${cmd_in} ${bowtie_sam}"
    exec_cmd ${cmd} 2> $log

    local cmd="${SAMTOOLS_PATH}/samtools view -bS ${bowtie_sam} > ${bowtie_bam}"
    exec_cmd ${cmd} 2>> $log

    ## zip output file
    for fastq in $(echo ${out}/*.fastq)
    do
        local cmd="gzip -f $fastq"
        exec_cmd ${cmd} 2>> $log
    done
    
    ## clean
    rm -f ${bowtie_sam}
}

## $1 = input files
## $2 = output dir
tophat2_func()
{

    echo -e "Running Tophat ..."

    if [[ $STRANDED == "reverse" ]]; then
	local stranded_opt="--library-type fr-firststrand"
    elif [[ $STRANDED == "yes" ]]; then
	local stranded_opt="--library-type fr-secondstrand"
    elif [[ $STRANDED == "no" ]]; then
	local stranded_opt="--library-type fr-unstranded"
    else
	die "Unknown STRANDED parameter ; must be reverse/yes/no"
    fi

    local log=$3/tophat2.log
    local out=$2/mapping
    mkdir -p ${out}
    inputs=($1)
    local out_mapped=$(basename ${inputs[0]} | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*/.bam/') 

    local cmd="${TOPHAT2_PATH}/tophat2 ${TOPHAT2_OPTS} --GTF ${TRANSCRIPTS_GTF} ${stranded_opt} -o ${out} ${TOPHAT2_IDX_PATH} $1"
    exec_cmd ${cmd} 2> $log 

    cmd="mv ${out}/accepted_hits.bam ${out}/${out_mapped}"
    exec_cmd ${cmd} 2>> $log
}


## $1 = input file(s)
## $2 = output dir
star_func()
{
    echo -e "Running STAR mapping ..."
    local log=$3/star.log
    local cmd="${STAR_PATH}/STAR --genomeDir ${STAR_IDX_PATH} --sjdbGTFfile ${TRANSCRIPTS_GTF} --sjdbOverhang 151"

    ## input type
    inputs=($1)
    if [[ ${#inputs[@]} -eq 1 ]]; then
	if [[ ${inputs[0]} =~ \.gz ]]; then
            cmd_in="--readFilesIn <(gzip -cd ${inputs[0]})"
        else
            cmd_in="--readFilesIn $1"
        fi
    elif [[ ${#inputs[@]} -eq 2 ]]; then
	if [[ ${inputs[0]} =~ \.gz ]]; then
		cmd_in="--readFilesIn <(gzip -cd ${inputs[0]}) <(gzip -cd ${inputs[1]})"
	else
		cmd_in="--readFilesIn $1 $2"
	fi
    fi

    ## sample id
    if [[ ! -z ${SAMPLE_ID} ]]; then
	cmd="$cmd --outSAMattrRGline ID:$SAMPLE_ID SM:$3 LB:Illumina PL:Illumina"
    fi

    ## estimate counts during mapping
    if [[ ${COUNT_TOOL} == "STAR" ]]; then
	cmd="$cmd --quantMode GeneCounts"
    fi

    local out=$2/mapping
    mkdir -p ${out}
    local out_mapped=$(basename ${inputs[0]} | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*//')

    cmd="$cmd $cmd_in --outFileNamePrefix ${out}/${out_mapped} ${STAR_OPTS}"
    exec_cmd ${cmd} 2> $log

    cmd="mv ${out}/${out_mapped}Aligned.sortedByCoord.out.bam ${out}/${out_mapped}.bam"
    exec_cmd ${cmd} 2>> $log
}

## $1 = input file
## $2 = output dir
htseq_func()
{
    echo -e "Running HTSeq ..."
    ## For stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature.
    ## For stranded=yes and single-end reads, the read has to be mapped to the same strand as the feature.
    ## For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For stranded=reverse, these rules are reversed.
    
    if [[ $STRANDED == "reverse" ]]; then
        local stranded_opt="-s reverse"

    elif [[ $STRANDED == "yes" ]]; then
        local stranded_opt="-s yes"

    elif [[ $STRANDED == "no" ]]; then
	local stranded_opt="-s no"

    else
        die "Unknown STRANDED parameter ; must be reverse/yes/no"	
    fi

    local log=$3/htseq.log
    local out=$2/counts
    mkdir -p ${out}
    local out_count=$(basename $1 | sed -e 's/.bam$/_htseq.csv/')

    cmd="${HTSEQ_PATH}/htseq-count ${HTSEQ_OPTS} $stranded_opt $1 ${TRANSCRIPTS_GTF} > ${out}/${out_count}"
    exec_cmd ${cmd}  2> $log
}

## $1 = input file
## $2 = output dir
featurecounts_func()
{
    echo -e "Running FeatureCounts ..."
    ## -s <int>
    ## Perform strand-specific read counting. Acceptable values:
    ## 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
    ## 0 by default.


    if [[ $STRANDED == "reverse" ]]; then
        local stranded_opt="-s 2"

    elif [[ $STRANDED == "yes" ]]; then
        local stranded_opt="-s 1"

    elif [[ $STRANDED == "no" ]]; then
        local stranded_opt="-s 0"

    else
	die "Unknown STRANDED parameter ; must be reverse/yes/no"
    fi

    local log=$3/featurecounts.log
    local out=$2/counts
    mkdir -p ${out}
    local out_count=$(basename $1 | sed -e 's/.bam$/_featurecounts.csv/')

    cmd="${FEATURECOUNTS_PATH}/featureCounts -a ${TRANSCRIPTS_GTF} -o ${out}/${out_count} ${FEATURECOUNTS_OPTS} $stranded_opt $1"
    exec_cmd ${cmd} 2> $log
}

starcounts_func()
{
    ## In theory, the counts have been run together with the mapping
    ## So we just want to move them into the counts folder
    local log=$3/star.log
    local out=$2/counts
    mkdir -p ${out}
    local out_count=$(echo $1 | sed -e 's/.bam$/ReadsPerGene.out.tab/')

    cmd="mv ${out_count} ${out}/"
    exec_cmd ${cmd} 2> $log

}


## $1 = input file
## $2 = output dir
## $3 log dir
rseqc_func()
{
    echo -e "Running RSeqc ..."
    ## check library strandedness; default sampling is 200000 reads with mapq>30

    local log=$3/rseqc.log
    local out=$2/strand_check
    mkdir -p ${out}
    outfile=$(basename $1 ".bam").rseqc

    if [[ -z ${GENE_BED} ]]; then
        die "Check that transcript reference bed file exists in CONFIG. Exit"
    fi

    cmd="${PYTHON_PATH}/python ${RSEQC_PATH}/infer_experiment.py -i $1 -r ${GENE_BED} > ${out}/${outfile}"
    exec_cmd ${cmd} 2> $log
}

## $1 = input file
## $2 =  log dir
markdup(){
    echo -e "Mark duplicates ..."

    local log=$2/markdup.log
    local prefix=$(echo $1 | sed -e 's/.bam$//')  
    cmd="java -jar ${PICARD_PATH} MarkDuplicates I=$1 O=${prefix}_markdup.bam REMOVE_DUPLICATES=false M=${prefix}_METRIC"
    exec_cmd ${cmd}  2> $log
}

mapping_stat(){
    echo -e "Running mapping stats ..."

    local log=$6/getStatFile.log
    local output=$5/stats
    mkdir -p $output
    inputs=($1)

    if [[ ${#inputs[@]} -eq 1 ]]; then
	cmd_input="-f ${inputs[0]}"

    elif [[ ${#inputs[@]} -eq 2 ]]; then
	cmd_input="-f ${inputs[0]} -r ${inputs[1]}"

    fi

    outfile=$(basename ${inputs[0]} | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*/.stats/')

    if [ ! -z ${SAMPLE_ID} ]; then
	cmd="bash ${SCRIPTS_PATH}/getStatFile.sh $cmd_input -c $2 -b $3 -x $4 -g $TRANSCRIPTS_GTF -s ${SAMPLE_ID} > $output/$outfile"
    else
	cmd="bash ${SCRIPTS_PATH}/getStatFile.sh $cmd_input -c $2 -b $3 -x $4 -g $TRANSCRIPTS_GTF > $output/$outfile"
    fi
    exec_cmd ${cmd}  2> $log
}



