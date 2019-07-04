#!/bin/bash
env=$1
dir_conda_modules=/data/modules/pipelines/rnaseq/${env}
env_conda_name=rnaseq-2.0
### installation conda
bash /bioinfo/local/build/Centos/miniconda/Miniconda3-4.5.12-Linux-x86_64.sh -b -p ${dir_conda_modules}/conda

export PATH=${dir_conda_modules}/conda/bin:$PATH
echo $PATH
## Installation des tools via conda dans /data/modules/pipelines/rnaseq/$env 
conda env create -p ${dir_conda_modules}/conda/envs/${env_conda_name} -f environment.yml 2>&1 | tee ${dir_conda_modules}/install_env_conda.log 
