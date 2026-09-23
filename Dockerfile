# syntax=docker/dockerfile:1

#FROM docker.io/library/ubuntu:noble
FROM bioconductor/shiny:3.23


COPY scripts/install_dependency.R /scripts/install_dependency.R
RUN Rscript /scripts/install_dependency.R

COPY Shiny_AMMOA_local /Shiny_AMMOA_local
EXPOSE 3838
CMD ["Rscript", "-e", "shiny::runApp('/Shiny_AMMOA_local', host='0.0.0.0', port=3838)"]
