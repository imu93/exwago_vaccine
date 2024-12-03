# This script aims to calculate exWAGO TPM normalized values
# Load packages
pacman::p_load(rtracklayer, DGEobj.utils, dplyr)
# Define input variables
args = commandArgs(trailingOnly = T)
gtfFile = args[1]  # gene annotation gtf file 
countFile = args[2] # featureCounts table
exwago_id = args[3] # exWAGO ID

species = sub("\\..*", "", gtfFile) 
outFile = paste0(species, ".exWAGO_TPM_norm.txt")

# Import gtf and extract exons
gtf = import(gtfFile)
genes = gtf[gtf$type == "exon",]
# split by gene ID and count efective gene length
gene_list = split(genes, genes$gene_id)
exonic_lengths = lapply(gene_list, function(exon_group) sum(width(reduce(exon_group))))
# merge in a df
gene_length = do.call(rbind, exonic_lengths) %>% as.data.frame()
colnames(gene_length) = "Length"

# Now read counts  and estimate TPM
counts = read.delim(countFile)
counts = counts %>% as.matrix() 
tmp = convertCounts(counts, geneLength = gene_length$Length, unit = "TPM")
# I could save the whole table but not in this case

exwago_g = tmp[exwago_id,]
names(exwago_g) = NULL
# Write exWAGO tab
write.table(exwago_g %>% as.data.frame(), outFile, sep = "\t", quote = F)