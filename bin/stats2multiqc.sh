#!/bin/bash

splan=$1
aligner=$2
is_pe=$3


## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_id,Sample_name,Number_of_frag,Number_of_rRNA,Percent_of_rRNA,Strandness,Number_of_aligned,Percent_of_aligned,Number_of_uniquely_aligned,Percent_uniquely_aligned,Number_of_multiple_aligned,Percent_multiple_aligned,Number_of_duplicates,Percent_duplicates,Percent_saturation" > mq.stats

for sample in $all_samples
do
    ##id
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    ##n_frag
    n_frag='NA'
    if [ -d "rrna" ]; then
	n_frag=$(grep "reads processed" rrna/${sample}.log | cut -d: -f 2 | sed -e 's/ //')
    elif [ $aligner == "star" ]; then
	n_frag=$(grep "Number of input reads" alignment/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
    elif [[ $aligner == "hisat2" && $is_pe == "1" ]]; then
	n_frag=$(grep "Total pairs" alignment/${sample}.hisat2_summary.txt | cut -d: -f2 | sed -e 's/ //g')
    elif [[ $aligner == "hisat2" && $is_pe == "0" ]]; then
	n_frag=$(grep "Total reads" alignment/${sample}.hisat2_summary.txt | cut -d: -f2 | sed -e 's/ //g')
    elif [[ $aligner == "salmon" ]]; then
	n_frag=$(grep "num_processed" counts/${sample}/aux_info/meta_info.json | cut -f2 -d: | sed -e "s/ \|,//g")
    fi

    ##n_rRNA
    if [ -d "rrna" ]; then
	n_rrna=$(grep "reads with at least one" rrna/${sample}.log | cut -d: -f 2 | tr "(|)" " " | cut -d" " -f2)
	p_rrna=$(grep "reads with at least one" rrna/${sample}.log | cut -d: -f 2 | tr "(|)" " " | cut -d" " -f4 | sed -e 's/%//')
    else
	n_rrna='NA'
	p_rrna='NA'
    fi

    ## All values are in pairs for paired-end data /reads for single-end data
    if [ $is_pe == "1" ]; then
	if [ $aligner == "star" ]; then
	    n_unique=$(grep "Uniquely mapped reads number" alignment/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	    n_multi=$(grep "Number of reads mapped to multiple loci" alignment/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	    n_mapped=$(($n_unique + $n_multi))
	elif [ $aligner == "hisat2" ]; then
	    n_unique=$(grep " 1 time" alignment/${sample}.hisat2_summary.txt | cut -d: -f 2 | sed -e 's/ //g' | awk -F"(" 'BEGIN{s=0}{s=s+$1}END{print s}')
	    n_multi=$(grep ">1 time" alignment/${sample}.hisat2_summary.txt | cut -d: -f 2 | sed -e 's/ //g' | awk -F"(" 'BEGIN{s=0}{s=s+$1}END{print s}')
	    n_mapped=$(($n_unique + $n_multi))
	elif [ $aligner == "salmon" ]; then
	    n_mapped=$(grep "num_mapped" counts/SRR1106775_1/aux_info/meta_info.json | cut -f2 -d: | sed -e "s/ \|,//g")
	else
	    n_unique='NA'
	    n_multi='NA'
	    n_mapped='NA'
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
	elif [ $aligner == "salmon" ]; then
	    n_mapped=$(grep "num_mapped" counts/SRR1106775_1/aux_info/meta_info.json | cut -f2 -d: | sed -e "s/ \|,//g")
	else
	    n_unique='NA'
	    n_multi='NA'
	    n_mapped='NA'
        fi
    fi

    ##n_dup
    n_dup="NA"
    p_dup="NA"
    if [ -d "picard" ]; then
	if ls picard/${sample}*markDups_metrics.txt  1>/dev/null 2>&1; then 
	    n_dup_pair=$(grep -a2 "## METRICS" picard/${sample}*markDups_metrics.txt | tail -1 | awk -F"\t" '{print $7}')
	    n_dup_single=$(grep -a2 "## METRICS" picard/${sample}*markDups_metrics.txt | tail -1 | awk -F"\t" '{print $6}')
	    n_dup_optical=$(grep -a2 "## METRICS" picard/${sample}*markDups_metrics.txt | tail -1 | awk -F"\t" '{print $8}')
	    ## Picard values are always in reads
	    n_dup=$(( $n_dup_pair * 2 + $n_dup_single + $n_dup_optical ))
	    if [ $is_pe == "1" ]; then
		n_dup=$(( $n_dup / 2 ))
	    fi
	    p_dup=$(echo "${n_dup} ${n_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	fi
    fi

    ## Strandness
    strandness=$(cat strandness/${sample}_strandness.txt | awk -F, '{print $2}')

    ## Calculate percentage
    p_mapped='NA'
    p_unique='NA'
    p_multi='NA'
    if [[ $n_mapped != 'NA' && $n_frag != 'NA' ]]; then
	p_mapped=$(echo "${n_mapped} ${n_frag}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    fi
    if [ $n_unique != 'NA' ]; then  
	p_unique=$(echo "${n_unique} ${n_frag}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    fi
    if [ $n_multi != 'NA' ]; then  
	p_multi=$(echo "${n_multi} ${n_frag}" | awk ' { printf "%.*f",2,$1*100/$2 } ') 
    fi


    if [ -e preseq/${sample}_extrap_ccurve.txt ]; then
	p_sat=$(awk -v nbreads=${n_mapped} 'BEGIN{p=0} NR>1 && $1>nbreads && nbreads>p{val=$2} {p=$1; x=$2} END{printf "%.*f",2,val/x*100}' preseq/${sample}_extrap_ccurve.txt)
    else
	p_sat='NA'
    fi

    echo -e ${sample},${sname},${n_frag},${n_rrna},${p_rrna},${strandness},${n_mapped},${p_mapped},${n_unique},${p_unique},${n_multi},${p_multi},${n_dup},${p_dup},${p_sat} >> mq.stats

done
