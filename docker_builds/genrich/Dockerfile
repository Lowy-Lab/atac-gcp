FROM linearparadox/samtools
WORKDIR /tools
RUN apt-get update && apt-get install -y git build-essential python3-venv python3-pip && \
    git clone https://github.com/jsh58/Genrich.git && \
    make Genrich/Genrich && \
    python3 -m venv /tools/macs3 && \
    /tools/macs3/bin/pip install --upgrade macs3
WORKDIR /fastq
COPY callpeaks.sh /fastq
ENTRYPOINT ["./callpeaks.sh"]