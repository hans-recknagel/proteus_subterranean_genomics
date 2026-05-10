#this is a separate chromosome location finding script
#this script uses gff file not bwa mem results to determine the locations

library(ShortRead)
library(dplyr)
library(tidyr)

tablesPath = "./test_data/ortholog_location_target_tables/"
gffPath = "./path/to/gff/file.gff"


ofResultPath = "./test_data/orthofinder_results/"

tablesList = list.files(tablesPath, full.names = TRUE, pattern = ".csv")
tablesListNames = list.files(tablesPath, pattern = ".csv")
tablesListNames = gsub(pattern = "_candidate_gene_selection_analysis.csv", replacement = "", x = tablesListNames)


tables = lapply(tablesList, read.csv)
names(tables) = tablesListNames

tablesOGS = lapply(tables, function(x){
  out = x$OG_ID
  out
  
})

ogs = lapply(1:length(tablesOGS), function(x, ogsIn){
  
  analysisName = names(ogsIn[x])
  ogsPath = paste0(ofResultPath, analysisName, "/Orthogroup_Sequences")
  ogsListFull = list.files(ogsPath, full.names = TRUE)
  ogsListnames = list.files(ogsPath)
  ogsListnames = gsub(pattern = ".fa", replacement = "", x = ogsListnames)
  
  ogsSelection = ogsIn[[x]]
  ogsListnamesSel = ogsListFull[which(ogsListnames %in% ogsSelection)]
  
  outOG = lapply(ogsListnamesSel, readAAStringSet)
  outOgHeaders = lapply(outOG, names)
  names(outOgHeaders) = ogsListnames[which(ogsListnames %in% ogsSelection)]
  
  outOgHeadersPrangNum = lapply(outOgHeaders, function(x){
    outGrep = x[grep("Prang", x)]
    outGrep
  })
  
  outOgHeadersPrang = lapply(outOgHeadersPrangNum, gsub, pattern = "_Prang", replacement = "")
  
  df = data.frame(lapply(outOgHeadersPrang, function(x) {
    x = unlist(x)
    length(x) = max(lengths(outOgHeadersPrang))
    return(x)
  }))
  df2 = pivot_longer(df, cols = 1:ncol(df), names_to = "OG_ID", values_to = "headers")
  df2 = drop_na(df2)
  return(df2)
  
}, ogsIn = tablesOGS)

gff = read.table(gffPath, header = FALSE, sep = "\t", col.names = c("chr", "model", "type", "start", "end", "score", "strand", "phase", "attributes"))

gff2 = gff %>%
  mutate(headers = gsub(pattern = ".*=", replacement = "", x = attributes)) %>%
  mutate(headers = gsub(pattern = ";$", replacement = "", x = headers)) %>%
#  mutate(len = if_else(strand == "+", end - start, start - end))
  mutate(len = end - start) %>%
  filter(type == "mRNA")

ogs2 = lapply(ogs, left_join, gff2, by = "headers")
ogsNA = lapply(ogs2, function(x){
  na = sum(is.na(x$len))
  perc = na/nrow(x)
  perc
})


