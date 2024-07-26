FROM bioconductor/bioconductor_docker:RELEASE_3_18

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root


RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends \
    cmake \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Use renv to install packages
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

WORKDIR /opt/clustOpt
RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json
COPY renv.lock renv.lock

RUN R -e "renv::restore()"

# Copy the package
COPY . /opt/clustOpt

# Install the package
RUN R -e "devtools::install()"

# Run container as non-root
RUN groupadd -g 10001 notroot && \
   useradd -u 10000 -g notroot notroot 

USER notroot:notroot

CMD ["/bin/bash"]

# Copy dockerfile into the image
COPY Dockerfile /Dockerfile