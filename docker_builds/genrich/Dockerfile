FROM linearparadox/samtools
WORKDIR /tools
RUN apt-get update && apt-get install -y git build-essential python3-venv python3-pip zlib1g-dev && \
    git clone https://github.com/jsh58/Genrich.git
WORKDIR /fastq
COPY call_peaks.sh /fastq
ENTRYPOINT ["bash", "/fastq/call_peaks.sh"]