FROM linearparadox/basegcp:latest
WORKDIR /tools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xvf samtools-1.21.tar.bz2 && rm samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && ./configure --disable-bz2 --prefix /tools/samtools && make && make install && \
    echo 'export PATH=$PATH:/tools/samtools/bin' >> ~/.bashrc
WORKDIR /fastq
COPY filter.sh /fastq/filter.sh
RUN chmod 777 filter.sh
ENV PATH="/tools/samtools/bin:${PATH}"
COPY filter_mt.py /fastq/
ENTRYPOINT ["bash", "/fastq/filter.sh"]