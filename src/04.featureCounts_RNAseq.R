# This script aims produce gene expression tables from bam files
# This script should be run within a the directory where all the bamfiles
# from a single species are located. 

# As input the script also rquere the gtf files of the gene annotaton of each species

# Load packages
pacman::p_load(Rsubread, rtracklayer, dplyr)
# import gtf to extract gene coordinates
args = commandArgs(trailingOnly = T)
gtfFile = args[1] # gtf input file
species = sub("\\..*", "", gtfFile)
outfile1 = paste0(species, ".mRNA.FC.RDS")
outfile2 = paste0(species, ".counts_mRNA.txt")

# read the annotation file
genome =  import(gtfFile)
genome = genome[genome$type == "gene",]

# Create a df for fc
df = data.frame("GeneID"= genome$ID,
                "Chr"=seqnames(genome),
                "Start"=start(genome),
                "End"=end(genome),
                "Strand"=strand(genome))

# deffine input bamfiles
bam_list = list.files(pattern = ".*sorted.bam$") 

# Set reverse estrand for H. bakeri
# The rest of datasets are unstranded
str = ifelse(grepl("heligmosomoides_bakeri|teladorsagia_circumcincta", species), 2, 0)
thr = 40
tb = featureCounts(files = bam_list,
                   annot.ext = df, # annotation
                   allowMultiOverlap= TRUE, 
                   strandSpecific = str,
                   isPairedEnd = T, # all my libraries are PE
                   requireBothEndsMapped = TRUE,  
                   countChimericFragments = FALSE,
                   nthreads = thr,
                   useMetaFeatures=TRUE,
                   fraction = TRUE) # distribute mm reads in fractions

# export fc object
saveRDS(tb, outfile1)
# and the fc table
write.table(tb$counts, outfile2,
            quote = F, sep = "\t", 
            col.names = T, row.names = T)