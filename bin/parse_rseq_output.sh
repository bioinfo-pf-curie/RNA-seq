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
rseq_output_file=$1

## $1 reseqc output
parse_rseqc_output()
{
    rseqout=$1
    #ret="undetermined"

    if [[ -z $rseqout || ! -e $rseqout ]]; then
        echo "ERROR - input was not find ! Exit"
	exit 1
    fi

    nb_fail=$(grep "failed" $rseqout | awk -F": " '{print $2}')
    if (( $(echo "$nb_fail 0.5" | awk '{print ($1 > $2)}') )); then
	ret="undetermined";
	echo "$ret"
    fi

    if [[ $(grep -c "PairEnd" $rseqout) -ne 0 ]]; then
	nb_ss=$(grep "1++" $rseqout | awk -F": " '{print $2}') ## fr-secondstrand = yes = forward
	nb_ds=$(grep "2++" $rseqout | awk -F": " '{print $2}') ## fr-firststrand = reverse
    else
	nb_ss=$(grep "++" $rseqout | awk -F": " '{print $2}') ## fr-secondstrand = yes = forward
	nb_ds=$(grep "+-" $rseqout | awk -F": " '{print $2}') ## fr-firststrand = reverse
    fi
	

    if [[ ! -z $nb_ss && ! -z $nb_ds ]]; then
	nb_diff=$(echo "$nb_ss $nb_ds" | awk '{print ($1 - $2)}')
		
	if (( $(echo "$nb_diff 0.5" | awk '{print ($1 > $2)}') )); then
	    ret="forward"
	elif (( $(echo "$nb_diff -0.5" | awk '{print ($1 < $2)}') )); then
	    ret="reverse"
	else
	    ret="no"
	fi
    fi

    if [ ! -z ${ret} ]; then
	echo "$ret"
    else
	echo "ERROR - unknown data type !"
	exit 1
    fi
}

parse_rseqc_output $rseq_output_file
