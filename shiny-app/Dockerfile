
# File: Dockerfile
FROM rocker/shiny-verse:latest


# FROM rocker/rstudio:4.2.0
# ENV SHINY_SERVER_VERSION 1.5.18.979
# RUN /rocker_scripts/install_shiny_server.sh


RUN apt-get update && apt-get install -y --no-install-recommends \
        bedtools \
        cmake \
        curl \
	lbzip2 \
        libbz2-dev \
        libicu-dev \
        liblzma-dev \
        libmagick++-dev \
        libpq-dev \
        libproj-dev \
        libssl-dev \
        libv8-dev \
        libxml2-dev \
        tcl-dev \
	libharfbuzz-dev \
	libfribidi-dev \
        tk-dev \
  && install2.r --error --deps TRUE \
    shinyWidgets


RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
ENV RENV_VERSION 0.16.0
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY  app-eas/renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library

RUN --mount=type=cache,target=/var/cache/renv R -e "renv::restore()"
# RUN R -e "renv::restore()"

# remove demo apps
RUN rm -rf /srv/shiny-server/*

# copy in the app
COPY app/ /srv/shiny-server/app/
COPY app-eas/ /srv/shiny-server/app-eas/
CMD ["/usr/bin/shiny-server"]
