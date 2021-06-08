FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/RNA-seq.git"
LABEL gitCommit="bc28c568b4184a54d8245a025513ae53ee38d0dd / devel"

# real path from baseDir: htseq.yml
ADD htseq.yml /opt/htseq.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/htseq.yml \
&& echo "source activate HTSeqCounts-0.13.5" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/HTSeqCounts-0.13.5/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
