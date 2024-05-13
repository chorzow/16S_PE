library(dada2); packageVersion("dada2")

args <- commandArgs(trailingOnly = TRUE)

print("Start 16S_script.R")
directory <- args[1] ### The first input after start pipeline
path <- paste0(directory, "/trimmed")
# list.files(path)

fnFs <- sort(list.files(path, pattern = ".1P.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = ".2P.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "(.[12]P.fastq)"), function(x) x[1])
print(sample.names)

filtFs <- file.path(directory, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(directory, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

print("Do output of samples")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE) # TODO: On Windows set multithread=FALSE

head(out)

rdataFile <- "./metagenomes.RData"
save.image(rdataFile)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)
save.image(rdataFile)
plotErrors(errR, nominalQ = TRUE)
save.image(rdataFile)
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
dadaFs[[1]]
save.image(rdataFile)

# merge and get sequence table
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                      maxMismatch = 20,
                      verbose = TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)

#Filtering seqtab table

print("######################")
print("Reads table dimension:")
dim(seqtab)
thresh <- args[2] ### The second input after start pipeline
seqtab = seqtab[,colSums(seqtab)>thresh]
print("")
print("######################")
print(paste0("Reads table dimension after filtration: ", thresh))
dim(seqtab)

table(nchar(getSequences(seqtab)))
save.image(rdataFile)
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
save.image(rdataFile)
sum(seqtab.nochim) / sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- as.data.frame(cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
                rowSums(seqtab.nochim)
                ))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedF", "merged", "nonchim")
rownames(track) <- sample.names
track$passed_frac <- track$nonchim / track$input 
head(track)
write.table(track, file = file.path(directory, "track_stats_DADA2.tsv"))
save.image(rdataFile)

print("Save information:frequency ASV from seqtab.nochim:")
#save information:frequency ASV from seqtab.nochim
seqtab.nochim.export <- t(seqtab.nochim)
rownames(seqtab.nochim.export) <- paste0("ASV_", seq(nrow(seqtab.nochim.export)))
head(seqtab.nochim.export)
write.table(seqtab.nochim.export, 
            file = file.path(directory, "seqtabnochim_n.tsv"), 
            col.names = NA, sep = '\t'
            )

#save information: make fasta file with ASV sequences
library("ape")
fastaMx <- as.matrix(colnames(seqtab.nochim))
rownames(fastaMx) <- rownames(seqtab.nochim.export)
write.dna(fastaMx, file = file.path(directory, "seqtabnochim.fasta"), 
          format = "fasta", nbcol = -1, colw = 700)
          
print("######################")
print("DADA2_part_1.R: completed.")
