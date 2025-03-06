FROM linearparadox/samtools
WORKDIR /tools
RUN wget https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar && \
    apt-get update && apt-get install -y default-jdk python3-venv python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    python3 -m venv /tools/pandas && \
    /tools/pandas/bin/pip install pandas gcsfs google-cloud-storage
COPY get_min_complex.py /fastq
COPY find_subsample.py /fastq
COPY picard.sh /fastq
WORKDIR /fastq
ENTRYPOINT ["bash", "/fastq/picard.sh"]