#!/bin/bash

usage="$(basename "$0") [-h] [-d] [-p] [-t] -- deploy git pipeline

where:
    -h show this help text
    -d specify dir_out installation
    -p specify git pipeline (ssh://git@gitlab.curie.fr:2222/xxxxx/yyyyyyyy.git)
    -t specify tag to fetch (tag or commit id)

Main use case :
 deploy_git_pipeline.sh -d /bioinfo/local/curie/ngs-data-analysis/centos -p ssh://git@gitlab.curie.fr:2222/data-analysis/RNA-seq.git -t 3561888d"

while getopts ":h:d:p:t:" options; do
    case "${options}" in
        d)
            DIROUT=${OPTARG}
            ;;
        p)
            URLPROJET=${OPTARG}
            ;;
        t)
            TAG=${OPTARG}
            ;;
        *)
            echo "$usage"
            exit 0
            ;;
    esac
done
shift $((OPTIND-1))

# Preparation commande 
cmd_git="mkdir -p ${DIROUT} && cd ${DIROUT} && echo '### git init started ###' && git init -q &&  echo '### git remote started ###' && git remote add origin ${URLPROJET} &&  echo '### git fetch started ###' && git fetch &&  echo '### git checkout started ###' && git checkout -q -f ${TAG} &&  echo '### git clean started ###' && git clean -q -d -x -f && rm -Rf .git && echo '### local users_bioinfo configuration ###' && mv ${DIROUT}/conf/users_bioinfo.config ${DIROUT}/conf/local.config && echo ${TAG} > version"

# Lancement commande de deploiement 
eval ${cmd_git}
#echo ${cmd_git}
