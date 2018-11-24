FROM nfcore/base
LABEL authors="Anthony Underwood" \
      description="Docker image containing all requirements for nf-core/bacterialmappingphylogeny pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-bacterialmappingphylogeny-1.0dev/bin:$PATH
