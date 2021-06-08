FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/RNA-seq.git"
LABEL gitCommit="bc28c568b4184a54d8245a025513ae53ee38d0dd / devel"

# real path from baseDir: star.yml
ADD star.yml /opt/star.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/star.yml \
&& echo "source activate star-2.7.6a" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/star-2.7.6a/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
