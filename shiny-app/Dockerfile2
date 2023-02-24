#syntax=docker/dockerfile:1.2

FROM rocker/rstudio:4.2.0

ENV RENV_VERSION 0.16.0
# ENV SHINY_SERVER_VERSION 1.5.18.979
RUN /rocker_scripts/install_shiny_server.sh

RUN apt-get update && apt-get install -y --no-install-recommends \
        bedtools \
        cmake \
        curl \
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
    && rm -rf /var/lib/apt/lists/*


RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY app/renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library

# mounting cache can speed up if running locally. If you comment this in, 
# be sure to comment out the other RUN R command.
RUN --mount=type=cache,target=/var/cache/renv R -e "renv::restore()"


# RUN R -e "renv::restore()"

# remove demo apps
RUN rm -rf /srv/shiny-server/*

# copy in the app
COPY app/ /srv/shiny-server/

# copy in container files # TODO: move up
COPY container_files/ /



CMD ["/usr/bin/shiny-server"]
