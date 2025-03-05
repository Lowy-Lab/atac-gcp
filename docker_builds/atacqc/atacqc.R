suppressMessages(library(ATACseqQC))
bam <- "filtered_sorted.bam"
index <- "filtered_sorted.bam.bai"
bam_dups <- readsDupFreq(bam, index)
complexity<-estimateLibComplexity(bam_dups, times=100, interpolate.sample.sizes=seq(0.1, 1, by=0.01))
complexity <- complexity[complexity$relative.size <= 1,]
write.csv(complexity, file="library_complexity.csv", row.names=FALSE)
png("fragment_size_distribution.png", width=10, height=7, unit="in", res=300)
fragSize <- fragSizeDist(bam, "frag-size-distribution")
dev.off()