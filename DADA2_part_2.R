library("ape")
library(dada2); packageVersion("dada2")

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  path <- args[1]
}

directory = main()

#return data to DADA2
seqtab.nochimOTU <- read.table(file = file.path(directory, "results/OTU_table_filtered.tsv"), 
                               header = TRUE, check.names = FALSE, row.names = 2, sep = '\t')

### sorting for decreasing order
seqtab.nochimOTU[1] <- NULL
sum_frequency = data.frame(rowSums(seqtab.nochimOTU))
colnames(sum_frequency) = "Sum_of_frequencies"
seqtab.nochimOTU = cbind(seqtab.nochimOTU, sum_frequency)
seqtab.nochimOTU = seqtab.nochimOTU[order(seqtab.nochimOTU$Sum_of_frequencies, 
                                          decreasing = TRUE, na.last = TRUE),]
seqtab.nochimOTU$Sum_of_frequencies <- NULL

seqtab.nochimOTU <- t(seqtab.nochimOTU)

#taxonomy assignment
print("Start taxonomy assignment...")
taxa <- assignTaxonomy(seqtab.nochimOTU, "/home/lam2/main/silva_database/silva_nr99_v138_train_set.fa.gz", tryRC = T, multithread = T)
taxa <- addSpecies(taxa, "/home/lam2/main/silva_database/silva_species_assignment_v138.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
print('Get a taxonomy table...')
head(taxa.print)
dim(taxa)

print('Save OTU sequences as fasta...')
names(taxa)[names(taxa) == ""] <- "OTU_sequence"
taxa$OTU_name <- paste0("OTU_", seq(nrow(taxa)))
otunames <- paste0(">", taxa$OTU_name)
seqs <- c(rbind(otunames, taxa$OTU_sequence))
write(seqs, file.path(directory, "results/otu.fasta"))


print('Save taxonomy...')
#if required save taxonomy with and without nt
write.table(taxa, sep = "\t", file = file.path(directory, "results/all_OTUs_phylogeny.tsv"), col.names = NA)
write.table(taxa.print, sep = "\t", file = file.path(directory, "results/all_phylogeny.tsv"), col.names = NA)


#save information:frequency OTU from seqtab.nochi
seqtab.nochim.exportOTU <- t(seqtab.nochimOTU)
rownames(seqtab.nochim.exportOTU) <- paste0("OTU_", seq(nrow(seqtab.nochim.exportOTU)))
write.table(seqtab.nochim.exportOTU, sep = "\t", 
            file = file.path(directory, "results/all_OTU_frequency.tsv"), col.names = NA)
