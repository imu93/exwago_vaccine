# This script aims to produce tables with sRNA counts using segmented annotations
pacman::p_load(Rsubread, rtracklayer)
args = commandArgs(trailingOnly = T)
annFile = args[1]
genome =  import(annFile)
threads = 20

species = sub("\\..*", "", annFile)
outFile1 = paste0(species, ".srna_1827.FC.RDS")
outFile2 = paste0(species, ".counts.1827.txt")
# My annotation files has some unstranded categories such as Unknown
# Let's count in both strands just in case
unk_As = genome[grepl("Unk|MITE|Sat|intergenic", genome$class),]
strand(unk_As) = ifelse(strand(unk_As) == "+", "-", "+")
genome = c(genome, unk_As)

# Let's build the SAF fromat
df = data.frame("GeneID"=genome$ID, "Chr"=seqnames(genome), "Start"=start(genome), "End"=end(genome), "Strand"=strand(genome))

# In bams 
bam_list = list.files(pattern = ".*mapped.bam$")

# And count
tb = featureCounts(files = bam_list, annot.ext = df, allowMultiOverlap= TRUE, strandSpecific = 1,
                   fracOverlap= .7, nthreads = threads, useMetaFeatures=TRUE, largestOverlap = TRUE)

# allowMultiOverlap = TRUE means that I will allow feature-feature overlaps
# strandSpecific = 1 means that the count should be guided by strand
# fracOverlap means that at least 70% of each read must align to feature to be assigned
# largestOverlap = TRUE specify that if I have reads aligning to both features, reads will be assigned to the one with the larges overlap
# Bearing in mid that all these bams lack the 'NH' tag all of them will count as unique mappers

# Finally I'll save the table
saveRDS(tb, outFile1)
write.table(tb$counts, outFile2, quote = F, sep = "\t", col.names = T, row.names = T)
