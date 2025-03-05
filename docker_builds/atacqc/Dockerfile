FROM linearparadox/basegcp
RUN apt-get update
RUN apt install --yes r-base libssl-dev libcurl4-openssl-dev libxml2-dev libgsl-dev
RUN R -e "install.packages('BiocManager', dependencies=T, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('ATACseqQC', ask=FALSE)"
WORKDIR /fastq
COPY atacqc.sh /fastq/atacqc.sh
RUN chmod +x /fastq/atacqc.sh
COPY atacqc.R /fastq/atacqc.R
ENTRYPOINT ["bash", "/fastq/atacqc.sh"]
