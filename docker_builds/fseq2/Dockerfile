FROM linearparadox/samtools
WORKDIR /tools
RUN apt-get update && \
    apt-get install -y software-properties-common build-essential bedtools && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3.8 python3.8-dev python3.8-venv python3.8-distutils pkg-config libhdf5-dev
RUN python3.8 -m venv fseq2 && \
    fseq2/bin/pip3 install --upgrade pip wheel setuptools && \
    fseq2/bin/pip3 install "numpy==1.19.5" && \
    fseq2/bin/pip3 install "fseq2"

ENV PATH="/tools/fseq2/bin:${PATH}"
WORKDIR /fastq
COPY call_peaks.sh /fastq
ENTRYPOINT ["bash", "/fastq/call_peaks.sh"]