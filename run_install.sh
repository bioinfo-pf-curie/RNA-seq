#!/bin/bash
#set -eu
#
#  Copyright (c) 2018, U900, Institut Curie
#  Copyright (c) 2018, Philippe La Rosa (Equipe HPC)
#
# run_install.sh : This file is part of NGS=>RNASEQ
#
# Objectifs :
# Generateur d'architecture et d'environnement (espace d'exec) et lancement
# du pipeline via le WorkFlow Manager NextFlow. 
# run_main.sh effectue la creation des repertoires et fichiers 
# qui sont indispensables au lancement du pipeline via le script nextflow, 
# qui est lui-meme chargé de faire tourner un run SAFIR_ALL sur le cluster 
# de calcul sous Centos, pour les differentes taches du pipeline : integration 
# preliminaires, analyse NGS, integration KDI.
#
# Contexte d'usages :
# Dans le cadre d'un deploiement dans l'espace cible d'installation canonique
# c.a.d theoriquement installé dans : /bioinfo/pipelines/safir_all/${env}/bin 
# run_main.sh sera lancé à partir d'un des noued de soummission, soit  :
# par l'equipe data management ou via un automate (RIMS, jenkins, ...) et sous
# le user applicatif safir_all, ceci afin de cadrer/valider avec les contraintes
# liées aux droits sur le projet/data.
# Cependant run_main.sh peut aussi etre lance par n'importe quel user ayant bien sur
# les droits d'acces aux datai (sequencage, ...) et ceci afin de permettre 
# d'effectuer des tests de validations (jenkins, ...).
#
# Parametres Entree :
# -p : le nom du projet (exemple : SAFIR_ALL)
# -r : le nom du run (exemple : L153)
# -e : l'environnement (dev/valid/prod)
# -i : le chemin complet du repertoire resultat du sequencage, reads par exemple : 
# /data/transfert/Illumina/LHASSA_MISEQ/170320_M04391_0153_000000000-B4LBP
# [-o] : la racine chemin du repertoire ou seront sauvegarde l'ensemble des fichiers
# et repertoires lors du run (defaut : /data/tmp/NGS_RUN_TEMP)
# -q : nom de la file de job sur le cluster de calcul
# [-h] : usage
#
# Sorties :
# Soumission au cluster via qsub, plus :
# affichage console (stdout) exemple pour le run V295 du projet SAFIR_ALL :
# run_main.sh : create_nxf_script for projet SAFIR_ALL and run V295 :
# /data/tmp/NGS_RUN_TEMP/SAFIR_ALL_V295_1529940592001/SAFIR_ALL.nf
# make_run_safir_all_create_nxf_work_dir for projet SAFIR_ALL and run V295
# WORK_DIR = /data/tmp/NGS_RUN_TEMP/SAFIR_ALL_V295_1529940592001 
# REPORTING_DIR = /data/tmp/NGS_RUN_TEMP/SAFIR_ALL_V295_1529940592001/results
# make_run_safir_all_create_nxf_configs for projet SAFIR_ALL and run V295 :
# CONFIGS_PIPELINE_PATH : /data/tmp/NGS_RUN_TEMP/SAFIR_ALL_V295_1529940592001/configs/pipeline/V295.conf
# CONFIGS_NXF_PATH : /data/tmp/NGS_RUN_TEMP/SAFIR_ALL_V295_1529940592001/configs/nextflow/SAFIR_ALL-nf.config
#

SAFIR_ALL_DEBUG=${SAFIR_ALL_DEBUG:=0}; [[ "$SAFIR_ALL_DEBUG" == 'x' ]] && set -x

if [[ $TERM && $TERM != 'dumb' ]]
  then
    if command -v tput &>/dev/null
      then
        GREEN=$(tput setaf 2; tput bold)
        YELLOW=$(tput setaf 3)
        RED=$(tput setaf 1)
        NORMAL=$(tput sgr0)
    fi
fi

function echo_red() {
    >&2 echo -e "$RED$*$NORMAL"
}

function echo_green() {
    echo -e "$GREEN$*$NORMAL"
}

function echo_yellow() {
    >&2 echo -e "$YELLOW$*$NORMAL"
}

function die() {
  echo_red "$*"
  exit 1
}

### usage ###
function usage (){
    echo -e "\nUsage: $0"
    echo -e "\n [Options]"    
    echo -e "\t-p : nom du projet (exemple : SAFIR_ALL)"     
    echo -e "\t-e : environnement used to launch the operational test (either dev, valid or prod)" 
    echo -e "\t-r : nom du run (exemple : Z112)"  
    echo -e "\t-i : chemin complet du repertoire resultat de l'étape de pre-proc (exemple : /data/tmp_app/services/bcl2fastq/Z112-SAFIR_ALL"         
    echo -e "\t-o : racine chemin du repertoire ou seront sauvegarde l'ensemble des fichiers pour toutes les etapes : d'integration initiale et d'analyses"            
    echo -e "\t-q : nom de la file de job sur le cluster de calcul"        
    echo -e "\t-h : usage"          
    echo -e "\n\n [Example]: \n\t# $0 -p SAFIR_ALL -e dev -r Z112 -i /data/tmp_app/services/bcl2fastq/Z112-SAFIR_ALL -o /data/tmp/NGS_RUN_TEMP -q diag"
    exit 1
}

((!$#)) && echo "Il n'y a pas d'arguments!!" && usage

if [[ ($# < 11) || ($# > 12) ]]
then
    echo "Nombre d'arguments incorrect ($#) !!"
    usage
fi

while getopts "p:e:r:i:o:q:v" optionName; do
case "$optionName" in

p) project_name="$OPTARG";;
e) env="$OPTARG";;
r) run="$OPTARG";;
i) preproc_dir="$OPTARG";;
o) racine_work_dir="$OPTARG";;
q) queue="$OPTARG";;
h) usage;;
*) usage;;
esac
done

function safir_all_create_env() {
   echo " 1) create_env"
   # path par defaut du repertoire racine d'exec du pipeline
   RACINE_WORK_DIR_DEF="/data/tmp/NGS_RUN_TEMP"
   [[ ! $racine_work_dir ]] && racine_work_dir=${RACINE_WORK_DIR_DEF}
   WORK_DIR=${WORK_DIR:="${racine_work_dir}/${project_name}_${run}_${name_end}"}
   LOGNAME=${LOGNAME:="inconnu"}
   [[ ! $RACINE_PIPELINES_DIR ]] && RACINE_PIPELINES_DIR="/bioinfo/pipelines/rnaseq/${env}"
   NXF_DIR=./
   NXF_NAME=${NXF_NAME:="main.nf"}
   #NXF_BIN_DIR=${NXF_BIN_DIR:="/bioinfo/local/build/Centos/nextflow/nextflow-0.31.0.4885"}
   NXF_BIN_DIR=${NXF_BIN_DIR:="nextflow"}
   CONFIGS_DIR=${CONFIGS_DIR:="configs"}
   CONFIGS_NXF_DIR=${CONFIGS_NXF_DIR:="nextflow"}
   SCRIPTS_RUN_NXF_RACINE_NAME=${SCRIPTS_RUN_NXF_RACINE_NAME:="run"}
   LOCAL_SCRIPTS_PATH=${LOCAL_SCRIPTS_PATH:="${WORK_DIR}/${run}.sh"}
   BIN_DIR=${BIN_DIR:="bin"}
   DOCS_DIR=${DOCS_DIR:="docs"}
   ASSETS_DIR=${ASSETS_DIR:="assets"}
   LOCAL_LOG_DIR=${LOCAL_LOG_DIR:="LOG"}
   PIPELINES_PATH=${PIPELINES_PATH:="/bioinfo/pipelines"}
}


function safir_all_create_nxf_work_dir() {
   echo " 2) create_nxf_work_dir and configurations files"
   # creation architecture
   mkdir -p ${WORK_DIR}
   # copie des elements dans l'espace de travail
   cp -r ${RACINE_PIPELINES_DIR}/${CONFIGS_DIR} ${WORK_DIR}
   cp ${RACINE_PIPELINES_DIR}/nextflow.config ${WORK_DIR}
   # creation script de lancement cd nextflow
   cat <<COMM > ${LOCAL_SCRIPTS_PATH}
genome=\$1
cd  ${WORK_DIR}
# lancement du pipeline sur le cluster 

echo "cd \$(pwd);export NXF_OPTS='-Xms1g -Xmx4g';nextflow run main.nf -resume -profile curieu900 -c nextflow.config --genome \${genome} --aligner star --reads '${preproc_dir}' --queue '${queue}' --featureCounts" | qsub -q ${queue} -N nfcore-rnaseq
COMM

}

function safir_all_install_nxf_script() {
   echo " 3) install_nxf_script"
   [[ -e ${RACINE_PIPELINES_DIR}/${NXF_NAME} ]] || die "${RACINE_PIPELINES_DIR}/${NXF_NAME} n'est pas accessible"
   ln -s ${RACINE_PIPELINES_DIR}/${NXF_NAME} ${WORK_DIR}
   ln -s ${RACINE_PIPELINES_DIR}/${BIN_DIR} ${WORK_DIR}
   ln -s ${RACINE_PIPELINES_DIR}/${DOCS_DIR} ${WORK_DIR}
   ln -s ${RACINE_PIPELINES_DIR}/${ASSETS_DIR} ${WORK_DIR}
}

function nxf_pipeline() {
   echo " 4) nxf_pipeline_nfcore-rnaseq"
   
   echo_green "WORK_DIR = ${WORK_DIR}"
   echo_green "LOCAL_SCRIPTS_PATH = ${LOCAL_SCRIPTS_PATH}"
   echo_green "RACINE_PIPELINES_DIR = ${RACINE_PIPELINES_DIR}"
   echo_green "CONFIG_NXF_PATH = ${WORK_DIR}/${CONFIGS_DIR}"

   d=$(date)
}


function nxf_date() {
    local ts=$(date +%s%3N); [[ $ts == *3N ]] && date +%s000 || echo $ts
}

function get_date() {
    local ts=$(date +%d-%m-%y) && echo $ts
}
name_end=$(nxf_date)
date_run=$(get_date)

##
echo_yellow "#### $0 for projet ${project_name} env:${env} and run:${run} by ${LOGNAME}"
#
# creation des variables d'environnement nécessaires au lancement des outils du pipeline (fasqc, ...)
# creation ligne de comm pour la creation du fichier de configuration du pipeline creation du fichier 
# de configuration specifique au run safir_all ${WORK_DIR}/${CONFIGS_DIR}/${CONFIGS_NXF_DIR}
safir_all_create_env

#
# creation du repertoire d'exec : doit contenir l'ensemble des composants necessaire au lancement du pipeline
# creation repertoires resultats, configs, script, creation script bash de lancement 
# de nextflow pointant sur les parametres.
safir_all_create_nxf_work_dir

# ajout lien vers le script nxf du pipeline dans WORK_DIR
safir_all_install_nxf_script

#
# lancement de nexflow pour le pipeline
nxf_pipeline


