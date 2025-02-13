FROM ubuntu
WORKDIR tools
RUN apt-get update
RUN apt-get install -y build-essential wget curl unzip
COPY bowtie2-2.5.4-source.zip .
RUN unzip bowtie2-2.5.4-source.zip
RUN cd /tools/bowtie2-2.5.4 && make static-libs && make STATIC_BUILD=1
