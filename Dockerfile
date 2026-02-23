# syntax=docker/dockerfile:1
FROM ubuntu:22.04

LABEL maintainer="scrumpis"
ENV DEBIAN_FRONTEND=noninteractive

# Core dependencies
RUN apt-get update && apt-get install -y \
    wget bzip2 ca-certificates git curl unzip \
    python3 python3-venv python3-pip \
    build-essential make gcc g++ zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-gnutls-dev \
    r-base r-cran-ggplot2 \
    minimap2 mummer trf cd-hit ncbi-blast+ gnuplot \
    samtools bedtools \
    librsvg2-bin \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install R package 'rideogram'
RUN R -e "install.packages('rideogram', repos='http://cran.us.r-project.org')"

# Install Miniconda for tidk
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    /opt/conda/bin/conda update -n base -c defaults conda -y

ENV PATH="/opt/conda/bin:${PATH}"

# Install tidk via conda
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -y tidk && \
    conda clean --all -y

# Install quarTeT
WORKDIR /opt
RUN git clone https://github.com/aaranyue/quarTeT.git
WORKDIR /opt/quarTeT
RUN pip3 install .

# Make quartet accessible globally
ENV PATH="/opt/quarTeT:/opt/conda/bin:$PATH"

# Install Circos 0.69-9
WORKDIR /opt

RUN apt-get update && apt-get install -y \
    curl ca-certificates libgd-dev cpanminus && \
    rm -rf /var/lib/apt/lists/*

# Securely download Circos 0.69-9
RUN wget --no-check-certificate https://circos.ca/distribution/circos-0.69-9.tgz && \
    ls -lh circos-0.69-9.tgz && \
    tar -xzf circos-0.69-9.tgz && \
    mv circos-0.69-9 circos && \
    rm circos-0.69-9.tgz

# Ensure /bin/env exists for Perl scripts
RUN [ -f /bin/env ] || ln -s /usr/bin/env /bin/env

ENV PATH="/opt/circos/bin:${PATH}"

# Install required Perl modules
RUN cpanm --notest \
    Config::General     \
    Font::TTF::Font     \
    GD                  \
    Math::Bezier        \
    Math::Round         \
    Math::VecStat       \
    Params::Validate    \
    Regexp::Common      \
    Set::IntSpan        \
    Statistics::Basic   \
    SVG                 \
    Text::Format

# Validate the installation
RUN circos -modules && circos -version

# Set default behavior
CMD ["-h"]
(base) john7932@ANR-YQ96QLT6HC iwgc-circos % more Dockerfile
# syntax=docker/dockerfile:1
FROM ubuntu:22.04

LABEL maintainer="scrumpis"
ENV DEBIAN_FRONTEND=noninteractive

# Core dependencies
RUN apt-get update && apt-get install -y \
    wget bzip2 ca-certificates git curl unzip \
    python3 python3-venv python3-pip \
    build-essential make gcc g++ zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-gnutls-dev \
    r-base r-cran-ggplot2 \
    minimap2 mummer trf cd-hit ncbi-blast+ gnuplot \
    samtools bedtools \
    librsvg2-bin \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install R package 'rideogram'
RUN R -e "install.packages('rideogram', repos='http://cran.us.r-project.org')"

# Install Miniconda for tidk
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    /opt/conda/bin/conda update -n base -c defaults conda -y

ENV PATH="/opt/conda/bin:${PATH}"

# Install tidk via conda
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -y tidk && \
    conda clean --all -y

# Install quarTeT
WORKDIR /opt
RUN git clone https://github.com/aaranyue/quarTeT.git
WORKDIR /opt/quarTeT
RUN pip3 install .

# Make quartet accessible globally
ENV PATH="/opt/quarTeT:/opt/conda/bin:$PATH"

# Install Circos 0.69-9
WORKDIR /opt

RUN apt-get update && apt-get install -y \
    curl ca-certificates libgd-dev cpanminus && \
    rm -rf /var/lib/apt/lists/*

# Securely download Circos 0.69-9
RUN wget --no-check-certificate https://circos.ca/distribution/circos-0.69-9.tgz && \
    ls -lh circos-0.69-9.tgz && \
    tar -xzf circos-0.69-9.tgz && \
    mv circos-0.69-9 circos && \
    rm circos-0.69-9.tgz

# Ensure /bin/env exists for Perl scripts
RUN [ -f /bin/env ] || ln -s /usr/bin/env /bin/env

ENV PATH="/opt/circos/bin:${PATH}"

# Install required Perl modules
RUN cpanm --notest \
    Config::General     \
    Font::TTF::Font     \
    GD                  \
    Math::Bezier        \
    Math::Round         \
    Math::VecStat       \
    Params::Validate    \
    Regexp::Common      \
    Set::IntSpan        \
    Statistics::Basic   \
    SVG                 \
    Text::Format

# Validate the installation
RUN circos -modules && circos -version

# Set default behavior
CMD ["-h"]
