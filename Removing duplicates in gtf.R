transcripts <- read.delim("transcripts.gtf")

head(transcripts)

JustTranscripts <- subset(transcripts, transcripts$transcript == "transcript")
deduppedTranscripts <- subset(JustTranscripts, !duplicated(JustTranscripts[,9]))
            
write.table(deduppedTranscripts, file = "deduppedTranscripts.gtf", sep="\t", col.names=F, row.names=F,quote=F)
head(deduppedTranscripts)
