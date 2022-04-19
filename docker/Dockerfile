
# docker build --build-arg version=0.0.2.6 -t test.grab .

FROM rocker/tidyverse:4.0.0

LABEL maintainer="wenjianb@pku.edu.cn"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3-pip \
    libzstd-dev

RUN R -e "devtools::install_github('GeneticAnalysisinBiobanks/GRAB', ref='dev1')"

RUN echo "Complete installing GRAB package\n"

RUN echo ${version}

