FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>


RUN conda install bioconda::fastp=0.20.0
RUN conda install bioconda::multiqc
RUN conda install r::r=3.6.0
RUN conda install bioconda::bwa=0.7.17
RUN conda install bioconda::samtools=1.9
RUN conda install bioconda::picard=2.21.3
RUN conda install r-gsheet
RUN Rscript -e "install.packages('roperators',dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('fuzzyjoin', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
# RUN Rscript -e "install.packages('devtools', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('tidyverse', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  
