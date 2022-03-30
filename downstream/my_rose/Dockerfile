FROM ubuntu:20.04
LABEL maintainer "Tatsuya Hattori"

USER root
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && apt-get install -y \
    bash \
    curl \
    g++ \
    gcc \
    gfortran \
    git \
    grep \
    less \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    libreadline-dev \
    make \
    python2 \
    r-base \
    tar \
    unzip \
    wget \
    xorg-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python2 /usr/bin/python

# bcftools-1.2
WORKDIR /usr/local/src
RUN wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2 \
    && tar xvf bcftools-1.2.tar.bz2 \
    && cd bcftools-1.2 \
    && make \
    && make install \
    && cd htslib-1.2.1 \
    && make \
    && make install

# samtools-1.13
WORKDIR /usr/local/src
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar jxvf samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    ./configure && \
    make && \
    make install 


# add user
ARG username=hattori
ARG wkdir=/home/${username}
RUN useradd -m ${username}
COPY entrypoint.sh /
ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash"]


# ROSE
WORKDIR ${wkdir}
RUN git clone https://bitbucket.org/young_computation/rose.git && \
    wget http://younglab.wi.mit.edu/ROSE/ROSE_DATA.zip && \
    unzip ROSE_DATA.zip && \
    mv ROSE_DATA/* rose/

RUN chown -cR ${username}:${username} ${wkdir}

# for MRR detection
WORKDIR /usr/local/src
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz \
    && tar -zxvf bedtools-2.29.1.tar.gz \
    && cd bedtools2 \
    && make

RUN add-apt-repository universe \
    && apt-get update && apt-get install -y \
    python-setuptools \
    python3-setuptools \
    python3-pip

RUN python3 -m pip install numpy==1.8.0 scipy==0.13.3 pysam==0.9.1.4


WORKDIR ${wkdir}
RUN git clone https://bitbucket.org/YichaoCai/rose_strict_share.git
