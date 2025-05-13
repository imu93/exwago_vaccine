# This script aims to produce fn matrices
# Load packages
pacman::p_load(rtracklayer,parallel)
args = commandArgs(trailingOnly = T)
# add suffix for outFile
species=args[1]
outFile = sub("\\..*", "", species)

# Read all Rds files produced with 01.bam2Rds.R
files = list.files(pattern = ".*.M1.Rds$")

# Now I need a function to obtain fn
get_fn_mtx = function(file, ...){
  cat("Processing file:", file, "\n")
  x = readRDS(file)
  lib  = subset(x, seqnames(x) != "MtDNA") # 
  fnuc = substr(lib$seq,1,1)
  lens = split(fnuc, width(lib))
  lis = lapply(lens,table)
  lis = lapply(lis, function(x){x[c("C","T","A","G")]})
  counts = as.data.frame(do.call(rbind, lis))
  return(counts)
}

# Run over all the replicates
fn_mtx = mclapply(files, get_fn_mtx, mc.cores = 20)
names(fn_mtx) = files
# Rave the list
outFile = paste0(outFile,".exwago.fn.all_lengths.Rds")
saveRDS(fn_mtx, outFile)