FROM ubuntu
WORKDIR /tools
RUN apt-get update
RUN apt-get install -y build-essential wget curl unzip apt-transport-https ca-certificates gnupg default-jre pipx pigz && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && apt-get install -y google-cloud-cli && \  
    wget -O fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc.zip && rm fastqc.zip && \
    chmod +x FastQC/fastqc && \
    echo 'export PATH=$PATH:/tools/FastQC' >> ~/.bashrc && \
    pipx ensurepath && \
    pipx install cutadapt && \
    pipx ensurepath
WORKDIR /fastq
COPY trimgalore.sh /fastq/trimgalore.sh
RUN chmod +x /fastq/trimgalore.sh
ENTRYPOINT ["/fastq/trimgalore.sh"]



