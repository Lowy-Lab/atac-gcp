FROM linearparadox/basegcp
LABEL project="appendix-atac"
WORKDIR /tools
COPY pipx_global.sh /tools/pipx_global.sh
RUN wget -O fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc.zip && rm fastqc.zip && \
    chmod +x FastQC/fastqc && \
    echo 'export PATH=$PATH:/tools/FastQC' >> ~/.bashrc && \
    apt-get install -y python3-dev default-jre && \
    pipx ensurepath && \
    curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    echo "export PATH=$PATH:/tools/TrimGalore-0.6.10" >> ~/.bashrc
RUN chmod +x /tools/pipx_global.sh && /tools/pipx_global.sh "cutadapt"
WORKDIR /fastq
COPY trimgalore.sh /fastq/trimgalore.sh
RUN chmod +x /fastq/trimgalore.sh
ENTRYPOINT ["bash", "-i", "/fastq/trimgalore.sh"]