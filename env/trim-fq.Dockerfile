FROM continuumio/miniconda3
RUN apt-get update && apt-get install -y procps && \
    apt-get clean
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels r
RUN conda create -n trim.fq \
                        bioconda::fastp=0.20.0 \
                        bioconda::multiqc=1.8 \
                        r=3.6.0 \
    && conda clean -a
ENV PATH /opt/conda/envs/trim-fq/bin:$PATH
# Use libhts.so from conda
ENV LD_LIBRARY_PATH /opt/conda/envs/trim-fq/lib
RUN conda env export --name trim-fq > trim-fq.yml

LABEL Name="trim-fq" Author="Katie Evans"