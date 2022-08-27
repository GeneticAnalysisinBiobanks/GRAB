
# cd /home/wenjianb/Docker/GRAB
# docker build --build-arg version=0.0.3.3 -t wenjianb/grab:0.0.3.3 .
# docker push wenjianb/grab:0.0.3.3
# docker save wenjianb/grab:0.0.3.3 > grab_0.0.3.3.tar
# gzip grab_0.0.3.3.tar

FROM rocker/tidyverse:4.0.0

LABEL maintainer="wenjianb@pku.edu.cn"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3-pip \
    libzstd-dev

WORKDIR /grab
ADD code.R .
ADD GRAB.Region.R .
ADD GRAB.ReadGeno.R .

RUN R -e "devtools::install_github('GeneticAnalysisinBiobanks/GRAB', ref='main')"

RUN echo "Complete installing GRAB package\n"

CMD echo "Example: docker run -v /home/wenjianb/Docker/GRAB/:/grab/ wenjianb/grab:0.0.3.3 Rscript /grab/code.R"

