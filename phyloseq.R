### good phyloseq script
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tibble) 
library(RColorBrewer)
library(vegan); packageVersion("vegan")
library(ape)
library(ggpubr)
library(scales)
library(ggradar)

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  path <- args[1]
}

directory = main()

otu_mat = read.csv2(file = paste(directory, '/results/all_OTU_frequency.tsv', sep=''),
                    header=TRUE, check.names = F, sep='\t')
tax_mat = read.csv2(file = paste(directory, '/results/all_phylogeny.tsv', sep=''),
                    header=TRUE, check.names = F, sep='\t')

tax_mat[1] = paste("OTU_", seq(nrow(tax_mat)), sep = '')

sample_names <- colnames(otu_mat)

print(sample_names)

colnames(otu_mat) = c('OTU', tail(sample_names, -1))


OTU_name = tax_mat[1]
colnames(OTU_name) = 'OTU'

rownames(tax_mat) = OTU_name$OTU
tax_mat[1] = NULL
rownames(otu_mat) = OTU_name$OTU
otu_mat[1] = NULL

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)


### get phyloseq object
carbom <- phyloseq(OTU, TAX)
carbom

sample_names(carbom)
rank_names(carbom)


### Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)


### plots figures
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.05) > 0, TRUE)
carbom_abund

total = median(sample_sums(carbom_abund))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_abund = transform_sample_counts(carbom_abund, standf)


fig_1 <- plot_bar(carbom_abund, fill = "Class") + 
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") + 
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 8), axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size = 12))

ggsave(paste(directory, "/results/phyloseq.png", sep=''))

fig_2 <- plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Class", taxa.order = "Class", 
             trans=NULL, low="beige", high="red", na.value="beige") + 
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20),
        axis.text.x = element_text(size = 6, angle = 65, vjust = 1, hjust = 1), 
        axis.title = element_text(size = 20))
ggsave(paste(directory, "/results/heat_map.png", sep=''))


# beta diversity analysis
set.seed(1)

?ordinate # see documentation

carbom.ord <- ordinate(carbom_abund, "NMDS", "bray", k=3)

sample_variables(carbom_abund)

?plot_ordination # see documentation

ord_p <- plot_ordination(carbom_abund, carbom.ord, axes=c(1, 2), type="Sample.common.name", 
                         color="Kits_type") + 
  geom_point(size=3) + scale_shape_manual(values=c(15, 16, 17, 18, 3, 4, 8, 9, 10)) + theme_light() + 
  theme(legend.title = element_text(size=12), legend.text = element_text(size=12),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16))


