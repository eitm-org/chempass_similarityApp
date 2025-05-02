FROM rocker/shiny:latest

LABEL maintainer="Kruttika Dabke <kdabke@emila.org>"

# system libraries of general use
RUN apt-get update && apt-get install --no-install-recommends -y \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.1 \
    libxml2 \
    && rm -rf /var/lib/apt/lists/*

# system library dependency for the example app
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

# basic shiny functionality
RUN R -q -e "install.packages(c('shiny', 'rmarkdown'))"

# install dependencies of the example app
RUN R -q -e "install.packages('Rmpfr')"

# copy the app to the image
RUN mkdir /root/chempass
COPY chempass /root/chempass

COPY Rprofile.site /usr/local/lib/R/etc/

EXPOSE 3838

# uncomment this line to run on bare shiny instead of shinyproxy
# CMD ["R", "-q", "-e", "shiny::runApp('/root/chempass')"]
