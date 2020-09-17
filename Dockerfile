FROM nfcore/base:1.10.2
LABEL authors="Peter Kruczkiewicz" \
      description="Docker image containing all software requirements for the peterk87/nf-illmap pipeline"

RUN apt-get update \
 && apt-get install -y procps curl \
 && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-illmap-1.0.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-illmap-1.0.1 > nf-illmap-1.0.1.yml
