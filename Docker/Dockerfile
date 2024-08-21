FROM bioconductor/bioconductor_docker:RELEASE_3_17

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root


RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends \
    cmake \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN Rscript -e 'install.packages(c("optparse", "remotes"))'
RUN Rscript -e 'install.packages(c("foreach","doParallel"))'
RUN Rscript -e 'install.packages(c("tidyverse"))'
RUN Rscript -e 'devtools::install_version("spatstat", version = "1.64-1")'
RUN Rscript -e 'install.packages("ranger")'
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'remotes::install_github("cran/Matrix", ref = "1.7-0")'
RUN Rscript -e 'BiocManager::install("glmGamPoi")'
RUN Rscript -e 'install.packages("ggpubr")'
RUN Rscript -e 'install.packages("hdf5r")'
RUN Rscript -e 'install.packages("pheatmap")'
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'install.packages("remotes")'
RUN Rscript -e 'install.packages("future.apply")'
RUN Rscript -e 'install.packages("progressr")'
RUN Rscript -e 'install.packages("ranger")'

# Install seurat v5 packages
RUN Rscript -e 'remotes::install_github("satijalab/seurat", "seurat5")'
RUN Rscript -e 'remotes::install_github("satijalab/azimuth", "master")'
RUN Rscript -e 'remotes::install_github("satijalab/seurat-data", "master")'
#RUN Rscript -e 'remotes::install_version("SeuratObject", "5.0.2")'
RUN Rscript -e 'remotes::install_github("satijalab/seurat-wrappers", "master")'
RUN Rscript -e 'remotes::install_github("stuart-lab/signac", "master")'
RUN Rscript -e 'remotes::install_github("bnprks/BPCells")'


WORKDIR /opt

# Run container as non-root
RUN groupadd -g 10001 scClustOpt && \
   useradd -u 10000 -g scClustOpt scClustOpt 

USER scClustOpt:scClustOpt

CMD ["/bin/bash"]

# Copy dockerfile into the image
COPY Dockerfile /Dockerfile
