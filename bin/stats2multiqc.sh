#!/bin/bash

splan=$1
aligner=$2
is_pe=$3


## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_id,Sample_name,Number_of_reads,Number_of_rRNA,Percent_of_rRNA,Strandness,Number_of_aligned_reads,Percent_of_aligned_reads,Number_of_uniquely_aligned_reads,Percent_uniquely_aligned_reads,Number_of_multiple_aligned_reads,Percent_multiple_aligned,Number_of_duplicates,Percent_duplicates" > mq.stats

for sample in $all_samples
do
    ##id
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    ##n_reads
    if [ -d "rrna" ]; then
	n_reads=$(grep "reads processed" rrna/${sample}.log | cut -d: -f 2 | sed -e 's/ //')
    elif [ $aligner == "star" ]; then
	n_reads=$(grep "Number of input reads" alignment/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
    elif [[ $aligner == "hisat2" && $is_pe == "1" ]]; then
	n_reads=$(grep "Total pairs" alignment/${sample}.hisat2_summary.txt | cut -d: -f2 | sed -e 's/ //g')
    elif [[ $aligner == "hisat2" && $is_pe == "0" ]]; then
	n_reads=$(grep "Total reads" alignment/${sample}.hisat2_summary.txt | cut -d: -f2 | sed -e 's/ //g')
    fi

    ##n_rRNA
    if [ -d "rrna" ]; then
	n_rrna=$(grep "reads with at least one" rrna/${sample}.log | cut -d: -f 2 | tr "(|)" " " | cut -d" " -f2)
	p_rrna=$(grep "reads with at least one" rrna/${sample}.log | cut -d: -f 2 | tr "(|)" " " | cut -d" " -f4 | sed -e 's/%//')
    else
	n_rrna='NA'
	p_rrna='NA'
    fi

    if [ $is_pe == "1" ]; then
	if [ $aligner == "star" ]; then
	    n_unique=$(grep "Uniquely mapped reads number" alignment/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	    n_multi=$(grep "Number of reads mapped to multiple loci" alignment/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	    n_mapped=$(($n_unique + $n_multi))
	elif [ $aligner == "hisat2" ]; then
	    n_unique=$(grep " 1 time" alignment/${sample}.hisat2_summary.txt | cut -d: -f 2 | sed -e 's/ //g' | awk -F"(" 'BEGIN{s=0}{s=s+$1}END{print s}')
	    n_multi=$(grep ">1 time" alignment/${sample}.hisat2_summary.txt | cut -d: -f 2 | sed -e 's/ //g' | awk -F"(" 'BEGIN{s=0}{s=s+$1}END{print s}')
	    n_mapped=$(($n_unique + $n_multi))
	else
	    echo -e "Aligner not yet supported"
	    exit 1
	    fi
    else
	if [ $aligner == "star" ]; then
            n_unique=$(grep "Uniquely mapped reads number" alignment/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
            n_multi=$(grep "Number of reads mapped to multiple loci" alignment/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
            n_mapped=$(($n_unique + $n_multi))
	elif [ $aligner == "hisat2" ]; then
	    n_unique=$(grep " 1 time" alignment/${sample}.hisat2_summary.txt | cut -d: -f 2 | sed -e 's/ //g' | awk -F"(" 'BEGIN{s=0}{s=s+$1}END{print s}')
	    n_multi=$(grep ">1 time" alignment/${sample}.hisat2_summary.txt | cut -d: -f 2 | sed -e 's/ //g' | awk -F"(" 'BEGIN{s=0}{s=s+$1}END{print s}')
	    n_mapped=$(($n_unique + $n_multi))
	else
            echo -e "Aligner not yet supported"
            exit 1
            fi
    fi

    ##n_dup
    n_dup="NA"
    p_dup="NA"
    if [ -d "picard" ]; then
	if ls picard/${sample}Aligned*markDups_metrics.txt  1>/dev/null 2>&1; then 
	    n_dup_pair=$(grep -a2 "## METRICS" picard/${sample}Aligned*markDups_metrics.txt | tail -1 | awk -F"\t" '{print $7}')
	    n_dup_single=$(grep -a2 "## METRICS" picard/${sample}Aligned*markDups_metrics.txt | tail -1 | awk -F"\t" '{print $6}')
	    n_dup_optical=$(grep -a2 "## METRICS" picard/${sample}Aligned*markDups_metrics.txt | tail -1 | awk -F"\t" '{print $8}')
	    n_dup=$(( $n_dup_pair * 2 + $n_dup_single + $n_dup_optical ))
	    p_dup=$(echo "${n_dup} ${n_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	fi
    fi

    ## Strandness
    strandness=$(cat strandness/${sample}_strandness.txt)

    ## Calculate percentage
    p_mapped=$(echo "${n_mapped} ${n_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    if [ $n_unique != 'NA' ]; then  
	p_unique=$(echo "${n_unique} ${n_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    else
	p_unique='NA'
    fi
    if [ $n_multi != 'NA' ]; then  
	p_multi=$(echo "${n_multi} ${n_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ') 
    else
	p_multi='NA'
    fi


    echo -e ${sample},${sname},${n_reads},${n_rrna},${p_rrna},${strandness},${n_mapped},${p_mapped},${n_unique},${p_unique},${n_multi},${p_multi},${n_dup},${p_dup} >> mq.stats

done
