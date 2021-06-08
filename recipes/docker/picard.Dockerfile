FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/RNA-seq.git"
LABEL gitCommit="bc28c568b4184a54d8245a025513ae53ee38d0dd / devel"

# real path from baseDir: picard.yml
ADD picard.yml /opt/picard.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/picard.yml \
&& echo "source activate picard-2.25.3" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/picard-2.25.3/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
