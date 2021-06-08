FROM centos:7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/RNA-seq.git"
LABEL gitCommit="bc28c568b4184a54d8245a025513ae53ee38d0dd / devel"

RUN yum install -y which \
&& yum clean all

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
