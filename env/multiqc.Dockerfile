FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

RUN conda install python=3.7
RUN conda install -c bioconda -c conda-forge multiqc
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  
