FROM rocker/shiny:latest

LABEL maintainer="Kruttika Dabke <kdabke@emila.org>"

# system libraries of general use
RUN apt-get update && apt-get install --no-install-recommends -y \
    pandoc \   
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl3 \
    libxml2 \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# system library dependency for the example app
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

# basic shiny functionality
RUN R -q -e "install.packages(c('shiny', 'rmarkdown', 'BiocManager'))"


RUN R -e "BiocManager::install('ComplexHeatmap', ask = FALSE)"

RUN R -e "install.packages(c(\
  'reticulate', \
  'png', \
  'grid', \
  'gridExtra', \
  'colorspace', \
  'readxl', \
  'Matrix', \
  'tibble', \
  'dplyr', \
  'scales', \
  'stringr', \
  'ggplot2', \ 
  'vegan', \
  'base64enc', \
  'shinyjs', \
  'shinyWidgets'\
), repos='https://cloud.r-project.org')"



ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p ${CONDA_DIR} && \
    rm ~/miniconda.sh

# Install Miniconda using reticulate
RUN R -e "reticulate::install_miniconda(update = TRUE)"

# Add Conda to PATH
ENV PATH="${CONDA_DIR}/bin:${PATH}"

# Setup environment for R to find conda
#RUN echo "export PATH=${CONDA_DIR}/bin:$PATH" >> /etc/R/Renviron.site

# Copy your environment.yml file
COPY my-rdkit-env2.yml /tmp/my-rdkit-env2.yml

# Create conda environment from the file
RUN conda env create -f /tmp/my-rdkit-env2.yml


# copy the app to the image
RUN mkdir /root/chempass
COPY chempass /root/chempass

# Copy your custom Rprofile.site
COPY Rprofile.site /usr/local/lib/R/etc/

# Append conda path to Rprofile.site after copying
RUN echo 'Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", "/opt/conda/bin"))' >> /usr/local/lib/R/etc/Rprofile.site

EXPOSE 3838

# uncomment this line to run on bare shiny instead of shinyproxy
CMD ["R", "-q", "-e", "shiny::runApp('/root/chempass')"]
