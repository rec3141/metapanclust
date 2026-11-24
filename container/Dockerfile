# metapanclust tools image. Build context: a directory holding spacedust/bin/{spacedust,foldseek} (static builds).
# docker build -t metapanclust-tools:0.2 . ; SIF via apptainer build ... docker-archive://tools.tar
FROM mambaorg/micromamba:2-ubuntu24.04
USER root
RUN apt-get update && apt-get install -y --no-install-recommends libatomic1 libgomp1 ca-certificates && rm -rf /var/lib/apt/lists/*
RUN micromamba install -y -n base -c conda-forge -c bioconda python=3.12 ppanggolin mmseqs2 pyrodigal pandas numpy scipy scikit-learn && micromamba clean -a -y
COPY spacedust/bin/spacedust spacedust/bin/foldseek /usr/local/bin/
ENV PATH=/opt/conda/bin:$PATH
# 0.2: modernised binny (rec3141/binny, branch modern) for the CAMI binning comparison
RUN micromamba install -y -n base -c conda-forge -c bioconda opentsne hdbscan networkx matplotlib-base pyarrow joblib pyhmmer git \
    && /opt/conda/bin/pip install --no-deps "binny @ git+https://github.com/rec3141/binny.git@5ef17c1" \
    && /opt/conda/bin/binny --help >/dev/null && micromamba clean -a -y
