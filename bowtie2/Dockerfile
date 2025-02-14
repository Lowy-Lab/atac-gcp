FROM ubuntu
WORKDIR /tools
RUN apt-get update
RUN apt-get install -y build-essential wget curl unzip apt-transport-https ca-certificates gnupg curl && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && apt-get install -y google-cloud-cli && \
    wget -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-aarch64.zip && \
    unzip bowtie2.zip && rm bowtie2.zip && \
    echo 'export PATH=$PATH:/tools/bowtie2-2.5.4-linux-aarch64' >> ~/.bashrc