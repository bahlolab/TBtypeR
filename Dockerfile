FROM mambaorg/micromamba:1.5-jammy

LABEL \
  author="Jacob Munro" \
  description="Container for TBtypeR" \
  maintainer="Bahlo Lab"

# install os deps
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends procps \
  && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# install env with micromamba
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Install TBtypeR R package
COPY . /tmp/TBtypeR
RUN /opt/conda/bin/R --slave --vanilla -e \
    "devtools::install(pkg = '/tmp/TBtypeR', force = TRUE, upgrade = 'never')"

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC
