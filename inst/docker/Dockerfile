# =============================================================================
# scONCO Docker image
#
# Provides a reproducible R 4.4 + Seurat v5 + Bioconductor 3.20 environment
# pre-loaded with scONCO, its v1.1 dependencies, and the comparator tools
# used in benchmarks/tools/. Includes Python 3 + scanpy for cross-language
# interop (e.g., M3 trajectory module in v2.0).
#
# Build:
#   docker build -t kgdlab/sconco:1.1 -f inst/docker/Dockerfile .
# Run (interactive R):
#   docker run --rm -it kgdlab/sconco:1.1 R
# Run a benchmark:
#   docker run --rm -v $(pwd):/work kgdlab/sconco:1.1 \
#     Rscript benchmarks/run_all.R --datasets GSE163973 --tools scMRMA,scONCO --dry-run
# =============================================================================

# --- Base ---------------------------------------------------------------
FROM rocker/r-ver:4.4.1

LABEL maintainer="Chia-Jung Chang <charlene.cj.chang@gmail.com>"
LABEL description="scONCO: pan-cancer single-cell annotation framework"
LABEL version="1.1.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV RENV_PATHS_CACHE=/opt/renv-cache

# --- System dependencies ------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libgsl-dev libcairo2-dev libxt-dev \
    libhdf5-dev libgdal-dev libproj-dev libudunits2-dev \
    libglpk-dev libgmp-dev \
    git curl wget pandoc \
    python3 python3-pip python3-venv \
    && rm -rf /var/lib/apt/lists/*

# --- R packages: BiocManager + CRAN core --------------------------------
RUN R -e "install.packages(c('BiocManager','remotes','renv','devtools','testthat'), \
            repos='https://cloud.r-project.org', Ncpus = 4)"

# --- Bioconductor + Seurat v5 ------------------------------------------
RUN R -e "BiocManager::install(version = '3.20', ask = FALSE, update = FALSE)"
RUN R -e "install.packages(c( \
            'Seurat','SeuratObject','dplyr','tibble','tidyr','stringr', \
            'ggplot2','Matrix','pheatmap','pROC','easyPubMed','optparse'), \
            repos='https://cloud.r-project.org', Ncpus = 4)"
RUN R -e "BiocManager::install(c( \
            'SingleR','celldex','ComplexHeatmap','SummarizedExperiment', \
            'SingleCellExperiment','CelliD','SCINA'), \
            ask = FALSE, update = FALSE)"

# --- GitHub packages ----------------------------------------------------
RUN R -e "remotes::install_github('JiaLiVUMC/scMRMA', upgrade = 'never')"
RUN R -e "install.packages('scCATCH', repos='https://cloud.r-project.org')"
# GPTCelltype optional (needs OPENAI_API_KEY at runtime):
RUN R -e "remotes::install_github('Winnie09/GPTCelltype', upgrade = 'never')"

# --- Python deps (for M3 trajectory module + cross-language interop) ----
RUN pip3 install --no-cache-dir --break-system-packages \
        scanpy scvelo cellrank anndata h5py numpy pandas matplotlib

# --- Workdir + install scONCO itself -----------------------------------
WORKDIR /opt/scONCO
COPY DESCRIPTION NAMESPACE LICENSE README.md /opt/scONCO/
COPY R          /opt/scONCO/R/
COPY inst       /opt/scONCO/inst/
COPY tests      /opt/scONCO/tests/
COPY vignettes  /opt/scONCO/vignettes/

RUN R -e "remotes::install_local('/opt/scONCO', upgrade='never')"

# --- CLI link -----------------------------------------------------------
RUN ln -sf $(Rscript -e "cat(system.file('cli','sconco.R',package='scONCO'))") \
           /usr/local/bin/sconco && \
    chmod +x /usr/local/bin/sconco

# --- Final ---------------------------------------------------------------
WORKDIR /work
CMD ["R", "--no-save"]
