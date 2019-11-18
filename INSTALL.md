# INSTALL 

## Prerequisite

Check that `cmake` with  `version 3` is installed on you computer.
In some distributions (such as CentOS) it is available with the command `cmake3`

## Installation steps

* Let's assume you have cloned the git repository in the folder `${HOME}/RNA-seq`

* Create the folder `${HOME}/build`

* `cd ${HOME}/build`

### Display the options available for installation

The different options are displayed using the following command in the **Cache values** section:

* `cmake -N -LH ../RNA-seq/`

The options for the **A**nalysis **P**ipeline start with the prefix **AP**


### Set the option for the installation

Default options have to be replaced otherwise, the pipeline will not work.

* `cmake ../RNA-seq -DCMAKE_INSTALL_PREFIX=${HOME}/install -DAP_SINGULARITY_IMAGE_PATH=/path/to/images -DAP_ANNOTATION_PATH=/data/annotations/pipelines`

Other options can be provided.

### Installation

* `make; make install`

The pipeline will be available in ${HOME}/install. 
The config files with the defined options will be available in `${HOME}/install/conf`.

## For developpers

To avoid multiple installations for testing while developping, the developpers can copy the content of `${HOME}/install/conf` to `${HOME}/RNA-seq/conf` and run the pipeline directly from `${HOME}/RNA-seq/conf`.




