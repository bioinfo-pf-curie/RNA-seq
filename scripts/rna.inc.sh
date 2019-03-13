## rna.inc.sh
## Functions of RNA-seq pipeline
##
## Copyright (c) 2017-2018 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details
##
##
##
##
## GLOBAL VARIABLE ARE WRITTEN IN UPPER CASE
## local variable are written in lower case

check_env()
{
    if [ -z ${NB_PROC} ]; then
        export NB_PROC=4
    fi
}


## $1 = input file(s)
## $2 = output path
## $3 = log dir
fastqc_func()
{
    check_env
    local out=$2/fastqc
    local log=$3/fastqc.log
    mkdir -p ${out}
    
    echo -e "Running FastQC ..."
    echo -e "Logs: ${log}"
    echo

    local cmd="${FASTQC_PATH}/fastqc -t ${NB_PROC} -j ${JAVA_PATH}/java -o ${out} $1 2> $log"
    exec_cmd ${cmd} > ${log} 2>&1
}


## $1 = input files
## $2 = output dir
## $3 = log dir
rRNA_mapping_func()
{
    check_env
    if [[ -z ${BOWTIE_RRNA_INDEX} ]]; then
	die "rRNA indexes file not set. Exit"
    fi
    
    local log=$3/bowtie_rRNA.log
    local out=$2/rRNA_mapping
    mkdir -p ${out}

    ## Logs
    echo -e "Running rRNA mapping ..."
    echo -e "Logs: $log"
    echo

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
	die "Bowtie1 -  found more than two input files. Exit"
    fi

    local cmd="${BOWTIE_PATH}/bowtie ${BOWTIE_OPTS} -p ${NB_PROC} ${cmd_un} --sam ${BOWTIE_RRNA_INDEX} ${cmd_in} ${bowtie_sam}"
    exec_cmd ${cmd} > $log 2>&1

    local cmd="${SAMTOOLS_PATH}/samtools view -@ ${NB_PROC} -bS ${bowtie_sam} > ${bowtie_bam}"
    exec_cmd ${cmd} >> $log 2>&1

    ## zip output file
    for fastq in $(echo ${out}/*.fastq)
    do
        local cmd="gzip -f $fastq"
        exec_cmd ${cmd} >> $log 2>&1
    done
    
    ## clean
    cmd="/bin/rm -f ${bowtie_sam}"
    exec_cmd ${cmd} >> $log 2>&1
}

## $1 = input files
## $2 = output dir
tophat2_func()
{
    check_env
    local log=$3/tophat2.log
    local out=$2/mapping
    mkdir -p ${out}
    inputs=($1)
 
    echo -e "Running Tophat ..."
    echo -e "Logs: $log"
    echo

    if [ -z ${TRANSCRIPTS_GTF} ]; then
        die "TRANSCRIPT_GTF not defined - cannot run Tophat2 !"
    fi

    ## Check if bowtie2 is in the path
    which bowtie2 > /dev/null
    if [ $? != "0" ]; then
        die "Tophat2 launcher - Bowtie2 is not the $PATH !"
    fi

    local cmd="${TOPHAT2_PATH}/tophat2 -p ${NB_PROC}"
    if [[ $STRANDED == "reverse" ]]; then
	local stranded_opt="--library-type fr-firststrand"
    elif [[ $STRANDED == "yes" ]]; then
	local stranded_opt="--library-type fr-secondstrand"
    elif [[ $STRANDED == "no" ]]; then
	local stranded_opt="--library-type fr-unstranded"
    else
	die "Unknown STRANDED parameter ; must be reverse/yes/no"
    fi

    ## sample id
    if [[ ! -z ${SAMPLE_ID} ]]; then
        cmd="$cmd --rg-id ${SAMPLE_ID} --rg-sample ${SAMPLE_ID} --rg-library Illumina --rg-platform Illumina --rg-platform-unit ${SAMPLE_ID}"
    fi

    local out_prefix=${out}/$(get_fastq_prefix ${inputs[0]})

    local cmd="$cmd ${TOPHAT2_OPTS} --GTF ${TRANSCRIPTS_GTF} ${stranded_opt} -o ${out} ${TOPHAT2_INDEX} $1"
    exec_cmd ${cmd} > $log 2>&1

    cmd="mv ${out}/accepted_hits.bam ${out_prefix}.bam"
    exec_cmd ${cmd} >> $log 2>&1

    ## Clean log folder
    if [ -d ${out}/logs ]; then
        cmd="/bin/rm -rf ${out}/logs"
        exec_cmd ${cmd} >> $log 2>&1
    fi
}


## STAR
star_func()
{
    check_env

    local log=$3/mapping.log
    local cmd="${STAR_PATH}/STAR --genomeDir ${STAR_INDEX} --outSAMtype BAM SortedByCoordinate --runThreadN ${NB_PROC} --runMode alignReads"

    ## logs
    echo -e "Running STAR mapping ..."
    echo -e "Logs: $log"
    if [ -z ${TRANSCRIPTS_GTF+x} ]; then
        echo -e "Warning: GTF file not defined in STAR mapping"
    else
        cmd="$cmd --sjdbGTFfile ${TRANSCRIPTS_GTF}"
    fi
    echo

    ## input type
    inputs=($1)
    if [[ ${#inputs[@]} -eq 1 ]]; then
        if [[ ${inputs[0]} =~ \.gz ]]; then
            cmd_in="--readFilesIn <(gzip -cd ${inputs[0]})"
            #cmd_in="--readFilesIn ${inputs[0]} --readFilesCommand zcat"
        else
            cmd_in="--readFilesIn ${inputs[0]}"
	    fi
    elif [[ ${#inputs[@]} -eq 2 ]]; then
        if [[ ${inputs[0]} =~ \.gz ]]; then
            cmd_in="--readFilesIn <(gzip -cd ${inputs[0]}) <(gzip -cd ${inputs[1]})"
            #cmd_in="--readFilesIn ${inputs[0]} ${inputs[1]} --readFilesCommand zcat"
        else
            cmd_in="--readFilesIn ${inputs[0]} ${inputs[1]}"
        fi
    fi

    ## sample id
    if [[ ! -z ${SAMPLE_ID} ]]; then
        cmd="$cmd --outSAMattrRGline ID:${SAMPLE_ID} SM:${SAMPLE_ID} LB:Illumina PL:Illumina"
    fi

    ## estimate counts during mapping     
    if [[ ${COUNT_TOOL} == "STAR" ]]; then 
	cmd="$cmd --quantMode GeneCounts"     
    fi

    local out=$2/mapping
    mkdir -p ${out}
    local out_prefix=$(get_fastq_prefix ${inputs[0]})

    cmd="$cmd $cmd_in --outFileNamePrefix ${out}/ --outTmpDir ${out}/tmp ${STAR_OPTS}"
    exec_cmd ${cmd} > $log 2>&1

    cmd="mv ${out}/Aligned.sortedByCoord.out.bam ${out}/${out_prefix}.bam"
    exec_cmd ${cmd} >> $log 2>&1
}


## $1 = input file(s)
## $2 = output dir
## $3 = log file
## $4 = work on the first N reads
bowtie2_rseqc_func()
{
    check_env
    local log=$3
    echo -e "Running fast bowtie2 mapping on $4 first reads ..." >> ${log} 2>&1

    bowtie2_opts="-p ${NB_PROC} --fast --end-to-end --reorder -u $4"
    if [[ -z ${BOWTIE2_INDEX} ]]; then
        die "Bowtie2 indexes not detected. Exit"
    fi

    ## input type
    inputs=($1)
    if [[ ${#inputs[@]} -eq 1 ]]; then
            cmd_in="-U ${inputs[0]}"
    elif [[ ${#inputs[@]} -eq 2 ]]; then
            cmd_in="-1 ${inputs[0]} -2 ${inputs[1]}"
    fi

    ## Warning - do not store the bowtie logs to avoid multiQC reporting
    cmd="bowtie2 ${bowtie2_opts} -x ${BOWTIE2_INDEX} $cmd_in > $2/subsample.bam"
    exec_cmd ${cmd} >> $log 2>&1
}


## $1 = input BAM file
## $2 = log directory
bamindex_func()
{
    check_env
    local log=$2/index.log
    echo -e "Indexing BAM file ...."
    echo -e "Logs: ${log}"
    echo

    local cmd="${SAMTOOLS_PATH}/samtools index $1"
    exec_cmd ${cmd} > $log 2>&1
}


## $1 = input file
## $2 = output dir
## $3 = log dir
htseq_func()
{
    check_env
    local log=$3/htseq.log
    local out=$2/counts
    mkdir -p ${out}
    local out_count=$(basename $1 | sed -e 's/.bam$/_counts.csv/')
    
    echo -e "Running HTSeq ..."
    echo -e "Logs: ${log}"
    echo
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

    cmd="${HTSEQ_PATH}/htseq-count ${HTSEQ_OPTS} $stranded_opt $1 ${TRANSCRIPTS_GTF} > ${out}/${out_count}"
    exec_cmd ${cmd} > $log 2>&1
}


## $1 = input file
## $2 = output dir
## $3 = log dir
featurecounts_func()
{    
    check_env
    local log=$3/featurecounts.log
    local out=$2/counts
    mkdir -p ${out}
    local out_count=$(basename $1 | sed -e 's/.bam$/_counts.csv/')

    echo -e "Running FeatureCounts ..."
    echo -e "Logs: ${log}"
    echo
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

    cmd="${FEATURECOUNTS_PATH}/featureCounts -T ${NB_PROC} -a ${TRANSCRIPTS_GTF} -o ${out}/${out_count} ${FEATURECOUNTS_OPTS} $stranded_opt $1"
    exec_cmd ${cmd} > $log 2>&1
}


## $1 = input file
## $2 = output dir
## $3 = log dir
starcounts_func()
{
    ## In theory, the counts have been run together with the mapping
    ## So we just want to move them into the counts folder
    local log=$3/star.log
    local out=$2/counts
    mkdir -p ${out}

    local out_count=$(basename $1 | sed -e 's/.bam$/_counts.csv/')

    echo -e "Getting counts from the STAR mapper ..."
    echo
    #local out_count=$(echo $1 | sed -e 's/.bam$/ReadsPerGene.out.tab/')
    local out_star=$(dirname $1)/ReadsPerGene.out.tab

    cmd="mv ${out_star} ${out}/${out_count}"
    exec_cmd ${cmd} > $log 2>&1
}


## $1 = input file
## $2 = output dir
## $3 log dir
rseqc_func()
{
    ## check library strandedness; default sampling is 200000 reads with mapq>30
    local log=$3/rseqc.log
    local out=$2/strand_check
    local n=200000

    echo -e "Running RSeQC ..."
    echo -e "Logs: ${log}"
    echo

    mkdir -p ${out}
    mkdir -p ${out}/tmp
    if [[ -z ${GENE_BED} ]]; then
        die "Check that transcript reference bed file exists in CONFIG. Exit"
    fi
    if [[ -z ${BOWTIE2_INDEX} ]]; then
        die "Bowtie2 indexes required for RSeQC. Exit"
    fi

    ## Infer experiments from bam files
    if [[ $1 =~ ".bam$" ]]; then
        echo -e "Input bam file detected ..." > ${log} 2>&1
        outfile=$(get_fastq_prefix ${1}).rseqc
        input=$1
    ## Infer experiments from fastq files
    else
        echo -e "Input fastq files detected ..." > ${log} 2>&1
        inputs=($1)
        outfile=$(get_fastq_prefix ${inputs[0]}).rseqc

        ## Single end
        if [[ ${#inputs[@]} -eq 1 ]]; then
            bowtie2_rseqc_func ${inputs[0]} ${out}/tmp ${log} ${n}

        ## paired-end
        elif [[ ${#inputs[@]} -eq 2 ]]; then
            bowtie2_rseqc_func "${inputs[0]} ${inputs[1]}" ${out}/tmp ${log} ${n}
        fi
        input=${out}/tmp/subsample.bam
    fi

    cmd="python ${RSEQC_PATH}/infer_experiment.py -i $input -r ${GENE_BED} > ${out}/${outfile}"
    exec_cmd ${cmd} >> $log 2>&1

    ## clean mapping folder
    cmd="/bin/rm -rf ${out}/tmp/subsample.bam"
    exec_cmd ${cmd} >> $log 2>&1
}


## $1 reseqc output
parse_rseqc_output()
{
    rseqout=$1
    ret="undetermined"

    if [ ! -e $rseqout ]; then
        die "File $rseqout was not find ! Exit"
    fi

    ##PE
    if [[ $(grep -c "PairEnd" $rseqout) -ne 0 ]]; then
	nb_fail=$(grep "failed" $rseqout | awk -F": " '{print $2}')
	nb_fs=$(echo "$nb_fail > 0.5" | bc)
	if [ $nb_fs -eq 1 ]; then ret="undetermined"; fi

	nb_fr=$(grep "1++" $rseqout | awk -F": " '{print $2}') ## fr-secondstrand = yes
	nb_rf=$(grep "2++" $rseqout | awk -F": " '{print $2}') ## fr-firststrand = reverse

	if [[ ! -z $nb_fr && ! -z $nb_rf ]]; then
	    nb_yes=$(echo "$nb_fr - $nb_rf > 0.5" | bc)
	    nb_rev=$(echo "$nb_fr - $nb_rf < -0.5" | bc)
	    
	    if [ $nb_rev -eq 1 ];then
		ret="reverse"
	    elif [ $nb_yes -eq 1 ];then
		ret="yes"
	    else
		ret="no"
	    fi
	fi
    else
    ##SE
	nb_fail=$(grep "failed" $rseqout | awk -F": " '{print $2}')
        nb_fs=$(echo "$nb_fail > 0.5" | bc)
        if [ $nb_fs -eq 1 ]; then ret="undetermined"; fi

        nb_ss=$(grep "++" $rseqout | awk -F": " '{print $2}') ## fr-secondstrand = yes
        nb_ds=$(grep "+-" $rseqout | awk -F": " '{print $2}') ## fr-firststrand = reverse

        if [[ ! -z $nb_ss && ! -z $nb_ds ]]; then
	    nb_yes=$(echo "$nb_ss - $nb_ds > 0.5" | bc)
            nb_rev=$(echo "$nb_ss - $nb_ds < -0.5" | bc)
	
            if [ $nb_rev -eq 1 ]; then
		ret="reverse"
            elif [ $nb_yes -eq 1 ];then
		ret="yes"
            else
		ret="no"
            fi
	fi
    fi
    
    echo $ret
}


## $1 = input file
## $2 = output dir
## $3 log dir
preseq_func()
{
    check_env
    local log=$3/preseq.log
    local out=$2/preseq
    mkdir -p ${out}
    outfile=${out}/$(basename $1 ".bam").preseq
    
    echo -e "Running saturation curve with preseq ..."
    echo -e "Logs: ${log}"
    echo
    
    cmd="${PRESEQ_PATH}/preseq c_curve -o ${outfile} -v -B $1"
    exec_cmd ${cmd} > $log 2>&1
}


## $1 = input file
## $2 = log dir
## $3 = if == 1, output file erase input file
markdup()
{
    check_env
    local log=$2/markdup.log
    local prefix=$(echo $1 | sed -e 's/.bam$//')
 
    echo -e "Mark duplicates ..."
    echo -e "Logs: ${log}"
    echo

    ## VALIDATION_STRINGENCY=LENIENT - That way, it will report all the errors it sees to STDIN, but it will finish the job anyway.
    cmd="${JAVA_PATH}/java -jar ${PICARD_BIN} MarkDuplicates I=$1 O=${prefix}_markdup.bam M=${prefix}_METRIC VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false PG=null"
    exec_cmd ${cmd} > $log 2>&1
}


## $1 = input fastq file(s)
## $2 = configuration file 
## $3 = input BAM file
## $4 = output dir
## $5 = log directory
## $6 = rRNA BAM file (optional)
mapping_stat(){
    local log=$5/getStatFile.log
    local output=$4/stats
    mkdir -p $output
    inputs=($1)
    
    echo -e "Generate StatFile ..."
    echo -e "Logs: ${log}"
    echo
    
    if [[ ${#inputs[@]} -eq 1 ]]; then
	cmd_input="-f ${inputs[0]}"

    elif [[ ${#inputs[@]} -eq 2 ]]; then
	cmd_input="-f ${inputs[0]} -r ${inputs[1]}"

    fi

    opts="-b $3 -c $2"
    outfile=$(basename ${inputs[0]} | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*/.stats/')
    if [[ ! -z $ANNOT_DIR && ! -z $TRANSCRIPTS_GTF ]]; then
	prefix=$(basename ${TRANSCRIPTS_GTF} | sed -e 's/.gtf//')
	opts="${opts} -g ${ANNOT_DIR}/${prefix}";
    fi
    if [ ! -z $6 ]; then opts="${opts} -x $6"; fi
    if [ ! -z ${SAMPLE_ID} ]; then opts="$opts -s ${SAMPLE_ID}"; fi 
    
    cmd="bash ${SCRIPTS_PATH}/getStatFile.sh $cmd_input $opts > $output/$outfile"
    exec_cmd ${cmd} > $log 2>&1
}



