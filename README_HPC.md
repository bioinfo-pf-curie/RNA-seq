*Copyright (c) 2018, U900, HPC Institut Curie*  

**nfcore/rnaseq (instance curieU900)**
====================

## Objectif
Tests et validation d'une implémentation du pipeline **nfcore/rnaseq** en interne à la plate-forme et avec nos propres data et fichiers de configurations.

## Introduction
Afin de faciliter l'exploitation, le lancement est constitué en deux étapes; pour la première c'est le script "run\_install.sh" qui à été créé et implémenté à la racine du projet. 

Il effectue les tâches ci-dessous : 
 - Generateur d'architecture et d'environnement d'execution du pipeline (espace d'exec) : 

Il effectue la creation des répertoires et fichiers qui sont indispensables au lancement du pipeline via le script nextflow : **main.nf**, qui est lui-meme chargé de faire tourner un run sur le cluster de calcul, pour les differentes tâches du pipeline.
 - Géneration dans l'espace d'exec d'un script de lancement du pipeline (deuxieme étape)

## Options Entrees
```
-p : nom du projet d'analyse NGS  
-e : environnement utilisé (dev, valid ou prod)  
-r : nom du run  
-i : le chemin complet des données reads
-o : le chemin complet racine de l'espace d'exec
-q : nom de la file de job sur le cluster de calcul (exemple : diag, batch)      
-h : usage (affichage de l'usage)  
```
## Sorties(STDOUT)
affichage console (stdout) exemple pour le lancement d'un run d'un projet :
```
#### run_install.sh for projet RNASEQ env:dev and run:G261 by plarosa
 1) create_env
 2) create_nxf_work_dir and configurations files
 3) install_nxf_script
 4) nxf_pipeline_nfcore-rnaseq
WORK_DIR = /data/tmp/NGS_RUN_TEMP/RNASEQ_G261_1552391265888
LOCAL_SCRIPTS_PATH = /data/tmp/NGS_RUN_TEMP/RNASEQ_G261_1552391265888/G261.sh
RACINE_PIPELINES_DIR = /bioinfo/pipelines/sandbox/dev/nfcore/rnaseq
CONFIG_NXF_PATH = /data/tmp/NGS_RUN_TEMP/RNASEQ_G261_1552391265888/conf

```

## Lancements premiere étape 
### Installation espace d'exec

1. Ligne de commande exemple 

```
bash run_install.sh  -p RNASEQ -e dev -r G261 -i "/bioinfo/local/curie/ngs-data-analysis/testdataset_pipelines/RNA-seq/Gendrel2012/*{1,2}.fastq.gz" -o /data/tmp/NGS_RUN_TEMP -q diag

```    
                      
## Lancements du RUN
### mode nominal

1. se connecter à **calcsub**
2. lancer **qsub -I** (pour se retrouver sur un des noeud du cluster)
2. se positionner dans l'espace d'exec (WORK_DIR voir Sorties(STDOUT))  
3. Lancer le run avec la version du génome souhaité; à noter que le HG19 et HG38 sont implémentés dans la configuration 

Ligne de commande exemple 
```
bash G261.sh HG19

```
Script nextflow lancé 
```
nextflow run main.nf -resume -profile curieu900 -c nextflow.config --genome HG19 --reads '/bioinfo/local/curie/ngs-data-analysis/testdataset_pipelines/RNA-seq/Gendrel2012/*{1,2}.fastq.gz' --queue diag --email Philippe.La-Rosa@curie.fr

```    

A noter qui est possible d'ajouter son email dans le script si l'on souhaite être notifié du résultat par email 

dans ce cas il faut juste ajouter par exemple : --email Philippe.La-Rosa@curie.fr
Exemple 
```
genome=$1
cd  /data/tmp/NGS_RUN_TEMP/RNASEQ_G261_1552391265888
# lancement du pipeline sur le cluster

echo "cd $(pwd);export NXF_OPTS='-Xms1g -Xmx4g';nextflow run main.nf -resume -profile curieu900 -c nextflow.config --genome ${genome}  --reads '/bioinfo/local/curie/ngs-data-analysis/testdataset_pipelines/RNA-seq/Gendrel2012/*{1,2}.fastq.gz' --queue 'diag' --email 'Philippe.La-Rosa@curie.fr'" | qsub -q diag -N nfcore-rnaseq

```
