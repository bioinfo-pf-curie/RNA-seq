FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/RNA-seq.git"
LABEL gitCommit="bc28c568b4184a54d8245a025513ae53ee38d0dd / devel"

# real path from baseDir: deeptools.yml
ADD deeptools.yml /opt/deeptools.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/deeptools.yml \
&& echo "source activate deeptools-3.5.1" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/deeptools-3.5.1/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
