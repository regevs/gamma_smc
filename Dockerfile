FROM ubuntu:22.04

# Install packages and set the g++/gcc defaults
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y -qq install wget unzip make gcc-12 g++-12 bzip2 zlib1g-dev libbz2-dev liblzma-dev && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100 && \
    update-alternatives --install /usr/bin/cc gcc /usr/bin/gcc-12 100 && \
    apt-get -y clean all

# Download and install boost
RUN wget -q https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.bz2 && \
    tar xjf boost_1_81_0.tar.bz2 && \
    cd ./boost_1_81_0 && \
    ./bootstrap.sh --with-libraries=filesystem,iostreams,system,program_options --with-toolset=gcc && \
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

# Copy and build Gamma-SMC, use a general x86-64 architecture
COPY . /home/
RUN cd /home/ &&\
    make clean && make MARCH=x86-64 && \
    cp /home/bin/gamma_smc /usr/local/bin/gamma_smc

# Define entry point
ENTRYPOINT ["gamma_smc"]

# How to build image:
# $ docker build -t rs42/gamma_smc .
# How to convert to singularity:
# $ rm rs42_gamma_smc.sif && bsub -o ~/temp/out.txt -e ~/temp/err.txt -R "select[mem>=4000] rusage[mem=4000]" -M 4000 'remote-docker-singularity.py rs42/gamma_smc'
# How to push image:
# $ docker tag rs42/gamma_smc docker.io/regevsch/gamma_smc:v0.2-alpha    
# $ docker push docker.io/regevsch/gamma_smc:v0.2-alpha    
