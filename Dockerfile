FROM ubuntu:18.04
RUN apt-get update && apt-get -y upgrade
RUN apt-get -y -qq install wget unzip make gcc g++ zlib1g-dev libbz2-dev liblzma-dev

# Download and install boost
RUN wget -q https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.bz2 && \
    tar xjf boost_1_81_0.tar.bz2 && \
    cd ./boost_1_81_0 && \
    ./bootstrap.sh --with-libraries=filesystem,iostreams,system,program_options && \
    ./b2 -j 8 install

# Download and install zstd
RUN wget -q https://github.com/facebook/zstd/releases/download/v1.5.4/zstd-1.5.4.tar.gz && \
    tar xvf zstd-1.5.4.tar.gz && \
    cd zstd-1.5.4 && \
    make -j 8

# Download and install htslib
RUN wget -q https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
    tar xjf htslib-1.17.tar.bz2 && \
    cd htslib-1.17 && \
    ./configure && \
    make -j 8 && \
    make install

# Set env vars
ENV CPATH=/zstd-1.5.4/lib:/htslib-1.17
ENV LIBRARY_PATH=/zstd-1.5.4/lib:/usr/local/lib
ENV LD_LIBRARY_PATH=/zstd-1.5.4/lib:/usr/local/lib

# Download and install Gamma-SMC
RUN wget -q https://github.com/regevs/gamma_smc/releases/download/v0.1-alpha/gamma_smc-v0.1-alpha.zip && \
    unzip gamma_smc-v0.1-alpha.zip && \
    cd gamma_smc-main && \
    make
    
