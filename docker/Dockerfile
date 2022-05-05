
# cd /home/wenjianb/Docker/GRAB
# docker build --build-arg version=0.0.2.8 -t grab_0.0.2.8 .
# docker tag grab_0.0.2.8 wenjianb/grab:0.0.2.8
# docker push wenjianb/grab:0.0.2.8

FROM rocker/tidyverse:4.0.0

LABEL maintainer="wenjianb@pku.edu.cn"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3-pip \
    libzstd-dev

RUN R -e "devtools::install_github('GeneticAnalysisinBiobanks/GRAB', ref='dev1')"

RUN echo "Complete installing GRAB package\n"

CMD echo "Example: docker run -v /home/wenjianb/Docker/GRAB/:/grab/ grab:0.0.2.8 Rscript /grab/code.R"

