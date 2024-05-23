FROM mambaorg/micromamba:1.5.0
LABEL Author: Mike Sauria <mike.sauria@jhu.edu>

COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml .
RUN \
    micromamba install -n base -f conda.yml -y \
	&& micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN Rscript -e "install.packages('roperators', dependencies=TRUE, repos='http://cran.us.r-project.org')"

USER root

RUN apt-get --allow-releaseinfo-change update && \
	apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
