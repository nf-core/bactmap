FROM nfcore/base
LABEL authors="Anthony Underwood" \
      description="Docker image containing all requirements for nf-core/bactmap pipeline"

RUN apt update; apt install -y gcc bc procps

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

COPY bin/pipeline_python_scripts/filtered_bcf_to_fasta.py /usr/local/bin/ 
RUN chmod 755 /usr/local/bin/filtered_bcf_to_fasta.py

COPY bin/pipeline_python_scripts/calculate_fraction_of_non_GATC_bases.py /usr/local/bin/ 
RUN chmod 755 /usr/local/bin/calculate_fraction_of_non_GATC_bases.py

RUN wget https://download.asperasoft.com/download/sw/connect/3.8.1/ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz -O aspera-connect.tar.gz; \
    mkdir /aspera; tar xvfz aspera-connect.tar.gz -C /aspera ; \
    /aspera/ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh; \
    mv /root/.aspera /.aspera; \
    echo "[aspera]\nASPERA_BIN  = /.aspera/connect/bin/ascp\nASPERA_PRIVATE_KEY = /.aspera/connect/etc/asperaweb_id_dsa.openssh\nASPERA_OPTIONS =\nASPERA_SPEED = 100M" > /aspera.ini

RUN cd /root; git clone https://github.com/enasequence/enaBrowserTools.git; mv enaBrowserTools /enaBrowserTools

ENV PATH /opt/conda/envs/nf-core-bactmap-1.0dev/bin:$PATH
ENV PATH /enaBrowserTools/python3:$PATH