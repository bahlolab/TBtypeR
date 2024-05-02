FROM mambaorg/micromamba:1.5-jammy

LABEL \
  author="Jacob Munro" \
  description="Container for SplitStrains" \
  maintainer="Bahlo Lab"

# install os deps
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
      git \
      procps \
  && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# install env with micromamba
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# install SplitStrains commit 6850e65
RUN git clone https://github.com/WGS-TB/SplitStrains.git /SplitStrains \
    && cd /SplitStrains \
    && git checkout f936cb2 \
    && echo '#!/usr/bin/env python3' > /opt/conda/bin/splitStrains \
    && cat /SplitStrains/splitStrains.py >> /opt/conda/bin/splitStrains \
    && chmod +x /opt/conda/bin/splitStrains

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC