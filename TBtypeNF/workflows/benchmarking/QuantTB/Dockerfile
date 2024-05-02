FROM mambaorg/micromamba:1.5-jammy

LABEL \
  author="Jacob Munro" \
  description="Container for quanttb" \
  maintainer="Bahlo Lab"

# install os deps
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
      procps zip unzip \
    && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# install env with micromamba
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

RUN mkdir -p /tmp/quanttb/scripts/ && cd /tmp \
    && unzip -p /opt/conda/lib/python2.7/site-packages/quanttb-1.1-py2.7.egg quanttb/scripts/classify.py \
    | sed 's:java -jar:java -Xmx7G -jar:' | sed 's:-t 4:-t 2:' | sed 's:--threads 4:--threads 2:' > quanttb/scripts/classify.py \
    && zip /opt/conda/lib/python2.7/site-packages/quanttb-1.1-py2.7.egg quanttb/scripts/classify.py

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC