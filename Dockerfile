# syntax=docker/dockerfile:1
FROM rocker/r-ver:4.5.0@sha256:f97b4944221bdc309583dfbc37dc0bf2654cfbd1f330cded047aeec4ae5df7e3

LABEL org.opencontainers.image.title="receptor-signaling" \
      org.opencontainers.image.description="RDS single-cell receptor association with negative-binomial mixed models" \
      org.opencontainers.image.source="https://github.com/kjgaulton/receptor_signaling"

ENV CRAN_SNAPSHOT="https://packagemanager.posit.co/cran/2025-06-01" \
    BIOCONDUCTOR_VERSION="3.21" \
    OMP_NUM_THREADS="1" OPENBLAS_NUM_THREADS="1" MKL_NUM_THREADS="1" \
    HOME="/tmp" TMPDIR="/tmp"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       ca-certificates cmake pkg-config libnlopt-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
       libglpk-dev libgmp3-dev libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

COPY container/install.R /tmp/install.R
RUN Rscript --vanilla /tmp/install.R && rm /tmp/install.R

WORKDIR /opt/receptor-signaling
COPY R/ R/
COPY scripts/run_nb.R scripts/run_nb.R
COPY tests/test_nb.R tests/test_nb.R
COPY container/entrypoint.sh /usr/local/bin/receptor-signaling
COPY container/config.txt container/config.txt
COPY README.md CONTAINER.md ./
RUN chmod 755 /usr/local/bin/receptor-signaling \
    && mkdir -p /work && chown 10001:10001 /work

USER 10001:10001
RUN receptor-signaling --self-test
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/receptor-signaling"]
CMD ["--help"]
