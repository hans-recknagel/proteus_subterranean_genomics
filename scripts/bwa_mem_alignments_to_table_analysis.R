#bam results from bwa alignments from exp con analyses to table format

#Takes bam files from sam:to_bam.sh and ortholog list tables for input.

library(Rsamtools)
library(dplyr)

bamPath = "./path/to/exp_con_bwa_mem_alignments_analyses/"
ortablesPath = "./path/to/ortholog_location_target_tables/"

bamFiles = list.files(bamPath, pattern = "sorted.*bam", full.names = TRUE)
bamIndex = list.files(bamPath, pattern = ".*csi", full.names = TRUE)
bamIndex = gsub(pattern = "()\\.bai", replacement = "\\1", x = bamIndex)


bamNames = list.files(bamPath, pattern = "*sorted.bam")
bamNames = gsub(pattern = "()\\.bwa.aln.output_sorted.bam", replacement = "\\1", x = bamNames)


bamImport = lapply(1:length(bamFiles), function(x, b, i){
  out = scanBam(b[x])#, index = i[x])
},b = bamFiles, i = bamIndex)

names(bamImport) = bamNames

alignmentTables = lapply(bamImport, function(x){
  
  df = as.data.frame(x)
  df2 = df %>%
    mutate(og = gsub(pattern = ".*-(OG\\d+)-.*", replacement = "\\1", x = qname)) %>%
    mutate(seqLen = nchar(seq)) %>%
    mutate(chr = gsub(pattern = "(^.*\\d+)_sliding.*", replacement = "\\1", x = rname)) %>%
    mutate(chrStart = as.numeric(gsub(pattern = ".*sliding_(\\d+)-\\d+$", replacement = "\\1", x = rname))) %>%
    mutate(chrStart = chrStart - 1) %>%
    mutate(seqStart = chrStart + pos) %>%
    mutate(seqEnd = seqStart + seqLen) %>%
    mutate(query = gsub(pattern = ".*OG\\d+-(.*)$", replacement = "\\1", x = qname)) %>%
    group_by(og, qname) %>%
    arrange(desc(qwidth), .by_group = TRUE) %>%
    slice_max(qwidth)
  
  out = df2 %>%
    select(query, og, chr, seqStart, seqEnd)
  
  out
  
})

orList = list.files(ortablesPath, pattern = ".csv", full.names = TRUE)
orFiles = lapply(orList, read.csv)

alignmentTablesOut = lapply(1:length(alignmentTables), function(x, altables, ortables){
  
  al = altables[[x]]
  or = ortables[[x]]
  
  out = full_join(or, al, by = c("OG_ID" = "og"))
  out
  
}, altables = alignmentTables, ortables = orFiles)

names(alignmentTablesOut) = names(alignmentTables)

sapply(names(alignmentTables), function(x){write.csv(alignmentTablesOut[[x]], file = paste0(bamPath, x, ".csv"), quote = FALSE, row.names = FALSE)})

