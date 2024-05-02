FROM mambaorg/micromamba:1.5-jammy

LABEL \
  author="Jacob Munro" \
  description="Container for MixInfect" \
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

# install mixinfect script
ADD https://raw.githubusercontent.com/bensobkowiak/MixInfect/79a73156efa4a1d218a8d06e32c135e048c0387e/MixInfect_FreeBayes.R /tmp/MixInfect_FreeBayes.R
RUN echo "#!/usr/bin/env Rscript" > /opt/conda/bin/MixInfectFB.R \
    && sed '85d;515d;520d;521d;522d;523d;524d;' /tmp/MixInfect_FreeBayes.R >> /opt/conda/bin/MixInfectFB.R \
    && echo "args <- commandArgs(trailingOnly=TRUE)" >> /opt/conda/bin/MixInfectFB.R \
    && echo "MixInfect(args[1], output=args[2], Qual=as.integer(args[3]), LowCov=as.integer(args[4]))" >> /opt/conda/bin/MixInfectFB.R \
    && chmod +x /opt/conda/bin/MixInfectFB.R

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC