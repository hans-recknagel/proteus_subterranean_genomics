#script for filtering fasta files with contracted/expanded gene groups, to get the longest sequence according to the provided list (closest to human)
#changed to work with all orthogroups

library(Biostrings)
library(ShortRead)

analysis = "eye_loss"
conExp = "contracted"
#conExp = "expanded"


listPath = paste0("./path/to/convergent_sequences/", "species_priority_", analysis, ".txt")
seqPath = paste0("./path/to/convergent_sequences/", analysis)

outPath = paste0("./path/to/convergent_sequences/", analysis, "_out")


spList = read.table(listPath)
spList = spList$V1


seqFilesNames = list.files(seqPath)
seqFiles = list.files(seqPath, full.names = TRUE)

seqs = lapply(seqFiles, readAAStringSet)
names(seqs) = seqFilesNames

counter = 0

filteredSeqs = lapply(1:length(seqs), function(x, write, seqList){
    

  og = seqList[[x]]
  
  ogNames = names(og)
  seqname = gsub(pattern = "(\\.fa)", replacement = "", x = names(seqList)[[x]])
  
  isIn = NULL
  isInLen = length(isIn)
  i = 1
  while (isInLen == 0) {
    counter <<- counter +1
    
    sp = spList[[i]]
    isIn = grep(pattern = sp, x = ogNames)
    isInLen = length(isIn)
    i = i + 1
  }
  
  spFilter = og[isIn]
  spFilterLen = sapply(spFilter, length)
  spFilterLenMax = max(spFilterLen)
  spFilterLenMaxSeq = spFilter[which(spFilterLen == spFilterLenMax)][1]
  
  maxHeader = names(spFilterLenMaxSeq)
  maxHeaderOG = paste(seqname, maxHeader, sep = "-")
  names(spFilterLenMaxSeq) = maxHeaderOG
  
  if (write == TRUE) {
    dir.create(outPath)
    writePath = paste0(outPath, "/", seqname, "_", sp, ".fa")
    writeXStringSet(spFilterLenMaxSeq, writePath)
  }
  
  
  return(spFilterLenMaxSeq)
  
}, write = T, seqList = seqs)







