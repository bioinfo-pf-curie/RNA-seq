#!/bin/bash
## Nicolas Servant
## Run the pipeline for a set of sample using PBS -t option

cd $PBS_O_WORKDIR

## run:
## qsub -v SAMPLE_PLAN='./SAMPLE_PLAN',ODIR='/data/tmp/test_multi',CONFIG='CONFIG_star' run_qsub.sh

## Set PATHS
SCRIPTS_PATH=`dirname $0`
ABS_SCRIPTS_PATH=`cd "$SCRIPTS_PATH"; pwd`
BIN_PATH="$ABS_SCRIPTS_PATH/../bin/"
VERSION=$(bash ${BIN_PATH}/RNApip -v)


if [[ -z ${SAMPLE_PLAN} || -z ${ODIR} || -z ${CONFIG} ]];then
    echo -e "Error - Missing parameter(s). Stop."
    exit -1
fi

id=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $1}' ${SAMPLE_PLAN})
bioid=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $2}' ${SAMPLE_PLAN})
forward=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $3}' ${SAMPLE_PLAN})
reverse=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $4}' ${SAMPLE_PLAN})

mkdir -p ${ODIR}/${id}


echo -e "--------------------"
echo -e "Running RNA pipeline cluster mode (v${VERSION})"
when=$(date +%Y-%m-%d)
echo -e "Date: ${when}"
where=$(hostname)
echo -e "Host: ${where}"
echo
echo -e "Id: ${id}"
echo -e "Sample_id: ${bioid}"
echo -e "Forward: ${forward}"
echo -e "Reverse:  ${reverse}"
echo -e "Ouput: ${ODIR}"
echo -e "Config: ${CONFIG}"
echo -e "Log: ${ODIR}/${id}/rnapip.log"
echo -e "--------------------"
echo

## Run RNA pipeline per sample
if [ ! -z "$reverse" ]; then
    bash ${BIN_PATH}/RNApip -f ${forward} -r ${reverse} -o ${ODIR}/${id} -c ${CONFIG} -s ${bioid} > ${ODIR}/${id}/rnapip.log 
else
    bash ${BIN_PATH}/RNApip -f ${forward} -o ${ODIR}/${id} -c ${CONFIG} -s ${bioid} > ${ODIR}/${id}/rnapip.log
fi
