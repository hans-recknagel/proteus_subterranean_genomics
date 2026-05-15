#backtranslating OF results

#libraries
library(ShortRead)
library(stringi)
library(parallel)
library(dplyr)
source("./scripts/orthogroups_filter.R")

#Inputs for this function are the otrhogroup table from Orthofinder results, Orthogroup sequence fasta files and an output path.
orthogroups_filterOGs = orthogroups_filter(OGtable = "./path/to/orthogroups_tables/Orthogroups.tsv", OGfastaPath = "./path/to/othhofinder/fasta/results/Orthogroup_Sequences/", OGoutPath = "./path/to/filtered/output/", writeOGs = TRUE)

#Inputs for later are path to CDS files corresponding to the same species used for Orthofinder input, output from the function above and an output path. 
#cds path
cdsPath = "./path/to/cds/files/"
#OG result fastas
ofFolderPath = "./path/to/filtered/output/"
#out fasta Path
finalFastaOut = "./path/to/output/folder/"


#import cds files
cdsFilesFullNames = list.files(cdsPath, full.names = TRUE)
cdsFilesNames = list.files(cdsPath)

cds = lapply(cdsFilesFullNames, readFasta)
names(cds) = cdsFilesNames
cdsSeq = lapply(cds, sread)
cdsHeaders = lapply(cds, ShortRead::id)
cdsHeaders = lapply(cdsHeaders, as.character)



#putting all protein sequences in one dataframe
#first make a data frame
filterOGdf = lapply(orthogroups_filterOGs, lapply, as.data.frame)

#make columns with locus names and orthogroup names
filterOGdf = lapply(filterOGdf, function(x){
  out = lapply(1:length(x), function(y, dfs){
    og = names(dfs[y])
    df = dfs[[y]]
    df$headers = rownames(df)
    df$orthogroup = og
    df
    
  }, dfs = x)
  out
})

#combine dataframes for each number of species
filterOGdf = lapply(filterOGdf, do.call, what = rbind)

#add number of species in DFs
filterOGdf = lapply(1:length(filterOGdf), function(x, dfs){
  nSp = names(dfs[x])
  df = dfs[[x]]
  df$nSpecies = nSp
  df  
}, dfs = filterOGdf)

#combine in one DF
filterOGdf = do.call(rbind, filterOGdf)

#cds data frame
cdsdf = lapply(1:length(cdsSeq), function(x, cdsS, cdsH){
  sp = names(cdsS[x])
  names(cdsS[[x]]) = cdsH[[x]]
  out = as.data.frame(cdsS[[x]])
  out$header = rownames(out)
  out$species = sp
  colnames(out) = c("cdsSeq", "headers", "species")
  out
}, cdsS = cdsSeq, cdsH = cdsHeaders)

#combine cds in one data frame
cdsdf = do.call(rbind, cdsdf)

#joind cds and protein data frames
cdsOGdf = left_join(filterOGdf, cdsdf, by = "headers")

#
cdsOGdf2 = cdsOGdf %>%
  group_by(nSpecies) %>%
  group_split()

cdsOGdf2 = lapply(1:length(cdsOGdf2), function(x, inLst){
  df = inLst[[x]] %>%
    group_by(orthogroup) %>%
    group_split
  
  ogs = unique(inLst[[x]]$orthogroup)
  names(df) = ogs
  df
  
}, inLst = cdsOGdf2)

names(cdsOGdf2) = names(orthogroups_filterOGs)

backtranslatedOGs = lapply(1:length(cdsOGdf2), function(x, lst1){
 #x = 3
  lst1name = names(lst1[x])
  out = lapply(1:length(lst1[[x]]), function(y, lst2, nsp){
    #y = 809
    df = lst2[[y]]
    
    og = unique(df$orthogroup)
    dir.create(finalFastaOut)
    nspFolder = paste0(finalFastaOut, nsp, "/")
    dir.create(nspFolder)

    
    out = DNAStringSet(df$cdsSeq)
    names(out) = df$headers
    
    outfile = paste0(nspFolder, og, ".fa")
    
    #writeFasta(out, outfile)
    
    
  }, lst2 = lst1[[x]], nsp = lst1name)
  
}, lst1 = cdsOGdf2)

