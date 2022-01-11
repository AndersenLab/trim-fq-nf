FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

RUN Rscript -e "install.packages('roperators',dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('fuzzyjoin', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('tidyverse', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"

RUN apt-get --allow-releaseinfo-change update && \
	apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*
