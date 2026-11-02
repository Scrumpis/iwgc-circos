# syntax=docker/dockerfile:1
# Dependency image: run the versioned checkout through a bind mount at /data.
# amd64 only (the pinned Miniconda installer is x86_64).
FROM ubuntu:22.04 AS core
LABEL maintainer="scrumpis"
ENV DEBIAN_FRONTEND=noninteractive

# Signed Ubuntu packages supply Circos and Perl dependencies without an insecure
# direct tarball download or unversioned CPAN installations.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates wget bzip2 git curl \
    python3 gawk samtools bedtools circos liblist-moreutils-perl minimap2 librsvg2-bin fontconfig \
    && rm -rf /var/lib/apt/lists/*

# Make the exact licensed fallback available to SVG rasterizers as well as Circos.
COPY fonts/LiberationSans-Bold.ttf /usr/local/share/fonts/iwgc/LiberationSans-Bold.ttf
COPY fonts/LICENSE.LiberationSans /usr/share/doc/iwgc-circos-font/LICENSE.LiberationSans
RUN fc-cache -f /usr/local/share/fonts/iwgc

WORKDIR /data
CMD ["bash"]

FROM core AS full
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-pip python3-setuptools mummer trf cd-hit ncbi-blast+ \
    gnuplot r-base r-base-dev r-cran-ggplot2 libxml2-dev librsvg2-dev \
    && rm -rf /var/lib/apt/lists/*

# Optional quarTeT plotting dependency, retained for compatibility. te --noplot
# is the path used by iwgc-circos. Stop on install errors rather than R exit 0.
# r-base-dev supplies make/compilers, pkg-config and PNG/JPEG/ICU headers;
# XML and rsvg additionally need their development libraries above.
RUN Rscript -e 'install.packages("RIdeogram", repos="https://cloud.r-project.org"); stopifnot(requireNamespace("RIdeogram", quietly=TRUE))'

ARG MINICONDA=Miniconda3-py312_25.1.1-2-Linux-x86_64.sh
ARG MINICONDA_SHA256=4766d85b5f7d235ce250e998ebb5a8a8210cbd4f2b0fea4d2177b3ed9ea87884
RUN wget -q "https://repo.anaconda.com/miniconda/${MINICONDA}" -O /tmp/miniconda.sh \
    && echo "${MINICONDA_SHA256}  /tmp/miniconda.sh" | sha256sum -c - \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh
# Use only these channels; do not silently accept Anaconda default-channel terms.
RUN /opt/conda/bin/conda install -y --override-channels -c conda-forge -c bioconda tidk=0.2.0 \
    && /opt/conda/bin/conda clean --all -y
ENV PATH="/opt/conda/bin:${PATH}"

ARG QUARTET_REF=91b197f3992debafcaa409e82601374495b7e40e
RUN git clone https://github.com/aaranyue/quarTeT.git /opt/quarTeT \
    && git -C /opt/quarTeT checkout --detach "$QUARTET_REF" \
    && /usr/bin/pip3 install --no-deps /opt/quarTeT
ENV PATH="/opt/quarTeT:${PATH}"

RUN samtools --version && bedtools --version && circos -modules && circos -version \
    && tidk --version && quartet.py -h && rsvg-convert --version
WORKDIR /data
CMD ["bash"]
