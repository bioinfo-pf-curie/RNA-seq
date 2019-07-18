# Installation

This help page was originally cloned from the [nf-core](https://nf-co.re/) project.
To start using the nf-core/mypipeline pipeline, follow the steps below:

1. [Install Nextflow](#1-install-nextflow)
2. [Install the pipeline](#2-install-the-pipeline)
3. [Pipeline configuration](#3-pipeline-configuration)
    * [Cluster environment](#31-cluster-environment)
    * [Software dependencies](#32-software-dependencies)
    * [Other configurations](#33-other-configurations)
4. [Reference genomes](#4-reference-genomes)

## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

See [nextflow.io](https://www.nextflow.io/) for further instructions on how to install and configure Nextflow.

## 2) Install the pipeline

The pipeline does not require any specific installation. Just download/clone and transfer the pipeline files manually:

```bash
wget https://github.com/bioinfo-pf-curie/rnaseq/archive/master.zip
mkdir -p ~/ic-pipelines/rnaseq/
unzip master.zip -d ~/ic-pipelines/rnaseq/
cd ~/my_data/
nextflow run ~/ic-pipelines/rnaseq/main.nf --help
```

If you would like to make changes to the pipeline, it's best to make a fork on GitHub and then clone the files. Once cloned you can run the pipeline directly as above.


## 3) Pipeline configuration
By default, the pipeline loads a basic server configuration [`conf/base.config`](../conf/base.config)
This uses a number of sensible defaults for process requirements and is suitable for running on a simple (if powerful!) local server.

Be warned of two important points about this default configuration:

1. The default profile uses the `local` executor
    * All jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.
    * See the [nextflow docs](https://www.nextflow.io/docs/latest/executor.html) for information about running with other hardware backends. Most job scheduler systems are natively supported.
2. Nextflow will expect all software to be installed and available on the `PATH`
    * It's expected to use an additional config profile for singularity or conda support. See below.

#### 3.1) Cluster environment

In order to run the pipeline on a computational cluster, we define the `cluster` profile.
In this case, instead of running locally, each process is submitted to the cluster through your workflow management system.

```bash
## Send the jobs to a cluster
nextflow run main.nf -profile cluster
```

Note that in this case, the master job will still be run locally (whereas all processes will be submitted to the cluster).
It could therefore be useful to also submit the master job to the cluster

```bash
## PBS/Torque submission
echo "nextflow run main.nf -profile cluster" | qsub -l "mem=1gb,nodes=1:ppn=1"
```

#### 3.2) Software dependencies

If you do not want to locally install all software dependencies, then `conda` or `singularity` are your best friends.
The process is very similar: running the pipeline with the option `-profile singularity` or `-profile conda` tells Nextflow to enable singularity (resp. conda) for this run. 
If you're not able to use Singularity, you can instead use conda to manage the software requirements.
This is slower and less reproducible than the above, but is still better than having to install all requirements yourself!

In the current version of the pipeline, we offer several options:

```bash
nextflow run main.nf -profile curie
```

Will run the pipeline using a global conda environment available from our local institutional cluster (see conf/curie.config).

```bash
nextflow run main.nf -profile singularity
```

Will run the pipeline using singularity images for each processes. These images have to be defined in the conf/singularity.config.
Note that so far, there is no way to automatically fetch these images from a web server.


```bash
nextflow run main.nf -profile conda
```

Will build a new conda environment from the .yaml file before running the pipeline.
The pipeline ships with a conda environment file and nextflow has built-in support for this.
To use it first ensure that you have conda installed (we recommend [miniconda](https://conda.io/miniconda.html))

#### 3.3) Other configurations

See [`docs/configuration/adding_your_own.md`](configuration/adding_your_own.md)

## 4) Reference genomes

See [`docs/configuration/reference_genomes.md`](configuration/reference_genomes.md)
