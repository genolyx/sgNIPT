###############################################################################
#  Single Gene NIPT Pipeline - Docker Image
###############################################################################
#  Panel: Twist UCL_SingleGeneNIPT_TE-96276661_hg38
#
#  Includes:
#    - BWA-MEM2, samtools, bcftools, freebayes, mosdepth, fastp
#    - Nextflow runtime
#    - Python 3 with numpy, scipy, pysam, python-pptx
#    - LibreOffice headless (PPTX→PDF conversion)
#    - MultiQC, SnpSift (optional annotation)
#
#  Build:
#    docker build -t sgnipt:latest -f docker/Dockerfile .
###############################################################################

FROM ubuntu:22.04 AS base

LABEL maintainer="Single Gene NIPT Team"
LABEL description="Single Gene NIPT Pipeline for cfDNA-based prenatal screening"
LABEL panel="Twist UCL_SingleGeneNIPT_TE-96276661_hg38"

# Prevent interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

# ── System packages ──────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl ca-certificates gnupg2 \
    build-essential cmake pkg-config \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    libncurses5-dev libdeflate-dev \
    python3 python3-pip python3-dev \
    openjdk-17-jre-headless \
    libreoffice-impress libreoffice-common \
    fonts-noto-cjk fonts-noto-core \
    pigz tabix gzip unzip \
    procps less \
    && rm -rf /var/lib/apt/lists/*

# ── BWA-MEM2 ────────────────────────────────────────────────────────────────
ARG BWA_MEM2_VERSION=2.2.1
RUN wget -qO- https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VERSION}/bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2 \
    | tar -xjf - -C /opt/ \
    && ln -s /opt/bwa-mem2-${BWA_MEM2_VERSION}_x64-linux/bwa-mem2 /usr/local/bin/bwa-mem2

# ── samtools / bcftools / htslib ─────────────────────────────────────────────
ARG SAMTOOLS_VERSION=1.20
RUN wget -qO- https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    | tar -xjf - -C /tmp/ \
    && cd /tmp/samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) && make install \
    && rm -rf /tmp/samtools-${SAMTOOLS_VERSION}

ARG BCFTOOLS_VERSION=1.20
RUN wget -qO- https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    | tar -xjf - -C /tmp/ \
    && cd /tmp/bcftools-${BCFTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) && make install \
    && rm -rf /tmp/bcftools-${BCFTOOLS_VERSION}

# ── fastp ────────────────────────────────────────────────────────────────────
ARG FASTP_VERSION=0.23.4
RUN wget -q https://github.com/OpenGene/fastp/releases/download/v${FASTP_VERSION}/fastp \
    -O /usr/local/bin/fastp \
    && chmod +x /usr/local/bin/fastp

# ── freebayes ────────────────────────────────────────────────────────────────
RUN wget -q https://github.com/freebayes/freebayes/releases/download/v1.3.7/freebayes-1.3.7-linux-amd64-static.gz \
    -O /tmp/freebayes.gz \
    && gunzip /tmp/freebayes.gz \
    && mv /tmp/freebayes /usr/local/bin/freebayes \
    && chmod +x /usr/local/bin/freebayes

# ── mosdepth ─────────────────────────────────────────────────────────────────
ARG MOSDEPTH_VERSION=0.3.8
RUN wget -q https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth \
    -O /usr/local/bin/mosdepth \
    && chmod +x /usr/local/bin/mosdepth

# ── Nextflow ─────────────────────────────────────────────────────────────────
ARG NEXTFLOW_VERSION=24.04.4
RUN wget -q https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VERSION}/nextflow-${NEXTFLOW_VERSION}-all \
    -O /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow

# ── MultiQC ──────────────────────────────────────────────────────────────────
RUN pip3 install --no-cache-dir multiqc==1.22

# ── Python dependencies ─────────────────────────────────────────────────────
RUN pip3 install --no-cache-dir \
    numpy \
    scipy \
    pandas \
    pysam \
    python-pptx \
    pypdf \
    matplotlib \
    jinja2

# ── SnpSift (optional annotation) ───────────────────────────────────────────
ARG SNPSIFT_VERSION=5.2
RUN mkdir -p /opt/snpsift \
    && wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
       -O /tmp/snpeff.zip \
    && unzip -q /tmp/snpeff.zip -d /opt/ \
    && rm /tmp/snpeff.zip \
    && echo '#!/bin/bash\njava -jar /opt/snpEff/SnpSift.jar "$@"' > /usr/local/bin/SnpSift \
    && chmod +x /usr/local/bin/SnpSift

# ── Pipeline code ───────────────────────────────────────────────────────────
COPY bin/ /Work/SgNIPT/pipeline/
RUN chmod +x /Work/SgNIPT/pipeline/scripts/*.py

# ── Entrypoint ──────────────────────────────────────────────────────────────
COPY docker/entrypoint.sh /Work/SgNIPT/entrypoint.sh
RUN chmod +x /Work/SgNIPT/entrypoint.sh

# ── Create mount point directories ──────────────────────────────────────────
RUN mkdir -p \
    /Work/SgNIPT/fastq \
    /Work/SgNIPT/data \
    /Work/SgNIPT/config \
    /Work/SgNIPT/analysis \
    /Work/SgNIPT/output \
    /Work/SgNIPT/log

WORKDIR /Work/SgNIPT

ENTRYPOINT ["/Work/SgNIPT/entrypoint.sh"]
