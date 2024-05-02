FROM mambaorg/micromamba:1.5-jammy

LABEL \
  author="Jacob Munro" \
  description="Container for fastlin" \
  maintainer="Bahlo Lab"

# install os deps
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
      procps \
  && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# install env with micromamba
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# copy barcodes file
ADD https://raw.githubusercontent.com/rderelle/barcodes-fastlin/8ed3d83f8f54414a7930cfe211ecc43fc0f153af/MTBC_barcodes.tsv /fastlin/MTBC_barcodes.tsv
RUN chmod a+r /fastlin/MTBC_barcodes.tsv

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC