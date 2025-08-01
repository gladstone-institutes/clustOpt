FROM bioconductor/bioconductor_docker:RELEASE_3_19

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

COPY . /opt/clustOpt
WORKDIR /opt/clustOpt

# Remove renv files that cause build issues
RUN rm -rf renv/ renv.lock .Rprofile


# Install GitHub dependencies first
RUN R -e "remotes::install_github('bnprks/BPCells/r')"

# Install the local package with all dependencies
RUN R -e "devtools::install('.', dependencies = TRUE, upgrade = 'never')"

# Install strongly suggested package
RUN R -e  "BiocManager::install('glmGamPoi')"

# Install extra packages for Rscripts
RUN R -e 'install.packages(c("optparse", "readr"), repos = "https://cloud.r-project.org")'

# Default command
CMD ["/bin/bash"]

# Copy the Dockerfile into the image for reference
COPY Dockerfile /Dockerfile