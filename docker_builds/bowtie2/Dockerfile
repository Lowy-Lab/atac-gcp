FROM linearparadox/basegcp
LABEL project="appendix-atac"
WORKDIR /tools
RUN wget -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-aarch64.zip && \
    unzip bowtie2.zip && rm bowtie2.zip && \
    echo 'export PATH=$PATH:/tools/bowtie2-2.5.4-linux-aarch64' >> ~/.bashrc && \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xvf samtools-1.21.tar.bz2 && rm samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && ./configure --disable-bz2 --prefix /tools/samtools && make && make install && \
    echo 'export PATH=$PATH:/tools/samtools/bin' >> ~/.bashrc
WORKDIR /fastq
COPY align.sh /fastq/align.sh
RUN chmod 777 align.sh
ENV PATH="/tools/bowtie2-2.5.4-linux-aarch64:/tools/samtools/bin:${PATH}"
ENTRYPOINT ["bash", "-i", "/fastq/align.sh"]