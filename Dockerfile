FROM --platform=linux/amd64 bioconductor/bioconductor_docker:RELEASE_3_19

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

COPY . /opt/clustOpt
WORKDIR /opt/clustOpt

# renv is only used for local development of the package
RUN R -e 'renv::deactivate()'

# Install pak for faster package installation 
RUN R -e 'install.packages("pak", repos = "https://cloud.r-project.org")'

# Install using using pak (update the DESCRIPTION for new builds)
RUN R -e 'pak::pkg_install(pkg = ".", dependencies = TRUE)'

# Install optparse to make passing args to Rscripts easier
RUN R -e 'pak::pkg_install(pkg = "optparse")'

# Default command
CMD ["/bin/bash"]

# Copy the Dockerfile into the image for reference
COPY Dockerfile /Dockerfile