###############################################################################
#  Single Gene NIPT Pipeline - Docker Image
###############################################################################
#  Nextflow + Python analysis pipeline for cfDNA single-gene disorder screening
#
#  Panel: Twist UCL_SingleGeneNIPT_TE-96276661_hg38
#    - 5,138 target regions (~1.07 Mb), 183 genes, 22 chromosomes
#    - No Y-chromosome targets → ChrX-based supplementary FF
#    - 29 zero-probe regions flagged in QC
#
#  Base image: mambaorg/micromamba (conda/mamba for bioinformatics tools)
#  All bioinformatics tools installed via environment.yml (bioconda channel)
#
#  Build:
#    docker build -t sgnipt:latest -f docker/Dockerfile .
###############################################################################

FROM mambaorg/micromamba:1.5.8

USER root

LABEL maintainer="Single Gene NIPT Team"
LABEL description="Single Gene NIPT Pipeline for cfDNA-based prenatal screening"
LABEL panel="Twist UCL_SingleGeneNIPT_TE-96276661_hg38"

# ── System dependencies ─────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    procps \
    openjdk-17-jre-headless \
    ca-certificates \
    # LibreOffice headless for PPTX→PDF report conversion
    libreoffice-impress libreoffice-common \
    fonts-noto-cjk fonts-noto-core \
    && rm -rf /var/lib/apt/lists/*

# ── Install Nextflow ────────────────────────────────────────────────────────
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow && \
    nextflow -version

# ── Create work directories (volumes mounted at runtime) ────────────────────
RUN mkdir -p /Work/SgNIPT/fastq /Work/SgNIPT/data /Work/SgNIPT/config \
    /Work/SgNIPT/analysis /Work/SgNIPT/output /Work/SgNIPT/log

# ── Conda environment with pipeline tools ───────────────────────────────────
COPY docker/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes && \
    rm /tmp/environment.yml

# ── Pipeline code ───────────────────────────────────────────────────────────
#  Copies: main.nf, nextflow.config, modules/*.nf, scripts/*.py
COPY bin/ /Work/SgNIPT/pipeline/
RUN chmod +x /Work/SgNIPT/pipeline/scripts/*.py

# ── Entrypoint ──────────────────────────────────────────────────────────────
COPY docker/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# ── Environment ─────────────────────────────────────────────────────────────
ENV TZ=Asia/Seoul
ENV PYTHONUNBUFFERED=1
ENV NXF_HOME=/tmp/.nextflow

WORKDIR /Work/SgNIPT/pipeline

ENTRYPOINT ["/entrypoint.sh"]
