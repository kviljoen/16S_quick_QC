FROM debian:jessie

ENV CONDA_INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
#Exports conda path
ENV PATH $PATH:/opt/conda/bin/

RUN printf "deb http://archive.debian.org/debian/ jessie main\ndeb-src http://archive.debian.org/debian/ jessie main\ndeb http://security.debian.org jessie/updates main\ndeb-src http://security.debian.org jessie/updates main" > /etc/apt/sources.list

RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    wget \
    ca-certificates \
    bash \
    procps \
    wget \
    curl \
    gzip \
    perl && \
    wget --quiet https://repo.continuum.io/miniconda/${CONDA_INSTALLER} && \
    /bin/bash /${CONDA_INSTALLER} -b -p /opt/conda && \
    /opt/conda/bin/conda install --yes conda && \
    conda install conda-build && \
    #Update conda and uses it to install software used by YAMP
    #that is required to use YAMP on AWS Batch && \
    conda install -c bioconda -y bbmap fastqc multiqc &&\
    conda install -c conda-forge -y awscli &&\
    conda clean --yes --tarballs --packages --source-cache
   
################## UCT Hex specific ###########################
RUN mkdir -p /researchdata/fhgfs
RUN mkdir -p /scratch/DB/bio/YAMP
