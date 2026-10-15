# Reproducible image: CRISPResso2 2.3.4 + this wrapper.
#   docker build -t crispr-indel .
#   docker run --rm -v "$PWD":/data -w /data crispr-indel \
#       --locus /opt/crispr-indel/config/PROCR_g1.yaml --control WT.fastq.gz KO.fastq.gz -o results
FROM mambaorg/micromamba:1.5.8
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/crispr-indel
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-deps /opt/crispr-indel
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "crispr-indel"]
