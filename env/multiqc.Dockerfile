FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

RUN conda install bioconda::multiqc
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  
