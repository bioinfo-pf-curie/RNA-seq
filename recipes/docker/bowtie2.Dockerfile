FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/RNA-seq.git"
LABEL gitCommit="bc28c568b4184a54d8245a025513ae53ee38d0dd / devel"

# real path from baseDir: bowtie2.yml
ADD bowtie2.yml /opt/bowtie2.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/bowtie2.yml \
&& echo "source activate bowtie2-2.4.2" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/bowtie2-2.4.2/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
