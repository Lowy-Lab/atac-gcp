FROM ubuntu
LABEL infrastructure="gcp"
WORKDIR /tools
RUN apt-get update
RUN apt-get install -y build-essential wget curl unzip apt-transport-https ca-certificates gnupg curl libncurses5-dev zlib1g-dev liblzma-dev pigz python3 pipx python3.12-venv dos2unix && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && apt-get install -y google-cloud-cli                      