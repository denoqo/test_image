# File: Dockerfile
FROM rocker/shiny-verse:latest

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    lbzip2 \
  && install2.r --error --deps TRUE \
    shinyWidgets


RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY app/renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library

RUN --mount=type=cache,target=/var/cache/renv R -e "renv::restore()"


# remove demo apps
RUN rm -rf /srv/shiny-server/*

# copy in the app
COPY app/ /srv/shiny-server/
