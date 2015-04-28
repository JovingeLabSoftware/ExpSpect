# Dockerfile to configure dev environment for ExpSpect project
# Branching off of here: https://github.com/rocker-org/rocker/tree/master/rstudio
FROM rocker/rstudio

MAINTAINER Andrew Borgman "ndrwsbrgmn@gmail.com"

# install and configure s3cmd
RUN wget -O- -q http://s3tools.org/repo/deb-all/stable/s3tools.key | sudo apt-key add -
RUN wget -O/etc/apt/sources.list.d/s3tools.list http://s3tools.org/repo/deb-all/stable/s3tools.list
RUN sudo apt-get update && sudo apt-get install s3cmd
COPY .s3cmd /root/.s3cmd

# configure our R environment
RUN R -e "install.packages(c('R6', 'RCurl', 'RJSONIO', 'devtools'), repos='http://cran.rstudio.com/')"
RUN R -e 'source("http://bioconductor.org/biocLite.R");biocLite();biocLite("rhdf5");biocLite("hgu133plus2.db");'
RUN R -e "devtools::install_github('JovingeLabSoftware/ExpSpect')"

# grab our data from ec2
RUN s3cmd get s3://data.lincscloud.org/l1000/level4/zspc_n1328098x22268.gctx /data/zspc_n1328098x22268.gctx
RUN s3cmd get  s3://data.lincscloud.org/l1000/metadata/inst.info /data/zspc_n1328098x22268.gctx

