docker build -t linearparadox/basegcp:latest --no-cache base-gcp/
docker push linearparadox/basegcp:latest
docker build -t linearparadox/trimgalore:atac --no-cache trimgalore/
docker push linearparadox/trimgalore:atac
docker build -t linearparadox/bowtie2:atac --no-cache bowtie2/
docker push linearparadox/bowtie2:atac