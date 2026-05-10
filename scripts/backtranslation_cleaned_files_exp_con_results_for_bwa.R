#backtranslating OF results

#libraries
library(ShortRead)
library(stringi)
library(parallel)
library(dplyr)
library(tidyr)
source("./orthogroups_filter.R")

#path to list of orthogroups of interest
tablesPath = "./path/to/exp_con_files_for_bwa_analyses/"

cdsPath = "./path/to/cds/files/"
finalFastaOut = "/out/path/"

tables = list.files(tablesPath, full.names = TRUE)
tablesNames = list.files(tablesPath)
tablesNames = gsub(pattern = "\\.txt$", replacement = "", x = tablesNames)

bwaOGs = lapply(tables, read.table, header = TRUE)

bwaOGs2 = lapply(bwaOGs, function(x){
  out = x[,1]
  out
})

names(bwaOGs2) = tablesNames


OGtables = "./test_data/orthogroups_tables/"
OGtablesAllList = list.files(OGtables, full.names = TRUE, pattern = ".*Orthogroups\\.tsv")
OGtablesAllNames = list.files(OGtables, pattern = ".*Orthogroups\\.tsv")

#OGtablesAllList = OGtablesAllList[!grepl(pattern = "phylogeny", x = OGtablesAllList)]
#OGtablesAllNames = OGtablesAllNames[!grepl(pattern = "phylogeny", x = OGtablesAllNames)]

OGTablesAll = lapply(OGtablesAllList, read.table, sep = "\t", header = TRUE)
names(OGTablesAll) = OGtablesAllNames

ogdfList = lapply(1:length(bwaOGs2), function(x, bwadfs, ogdfs){

  bwadf = bwadfs[[x]]
  bwadfName = names(bwadfs)[x]

  if (grepl("cave", bwadfName)) {
    ogdf = ogdfs$cave_species_Orthogroups.tsv
    analysis = "cave_species"
  }else if (grepl("pigm", bwadfName)) {
    ogdf = ogdfs$pigmentation_loss_Orthogroups.tsv
    analysis = "pigmentation_loss"
  }else if (grepl("eye", bwadfName)) {
    ogdf = ogdfs$eye_loss_Orthogroups.tsv
    analysis = "eye_loss"
  }else if (grepl("geno", bwadfName)) {
    ogdf = ogdfs$genome_size_Orthogroups.tsv
    analysis = "genome_size"
  }else if (grepl("long", bwadfName)) {
    ogdf = ogdfs$longevity_Orthogroups.tsv
    analysis = "longevity"
  }

  ogdf2 = ogdf %>%
    filter(Orthogroup %in% bwadf) %>%
    select(Orthogroup, Prang.mark) %>%
    filter(Prang.mark != "")

  if (nrow(ogdf2) == 0) {
    return()

  }else{
    ogdf2 = ogdf2 %>%
      separate_wider_delim(Prang.mark, delim = ",", names_sep = "_", too_few = "align_start") %>%
      pivot_longer(cols = 2:ncol(.), names_to = "Prang", values_to = "header") %>%
      drop_na()
    return(ogdf2)
  }

}, bwadfs = bwaOGs2, ogdfs = OGTablesAll)

names(ogdfList) = tablesNames
tablesNames2 = tablesNames[!sapply(ogdfList,is.null)]
ogdfList = ogdfList[!sapply(ogdfList,is.null)]

orthogroups_filterOGs = lapply(1:length(ogdfList), function(x, ogdfs){
  bwadf = ogdfs[[x]]
  bwadf$header = gsub(pattern = " ", replacement = "", x = bwadf$header)
  bwadfName = names(ogdfs)[x]

  if (grepl("cave", bwadfName)) {
    analysis = "cave_species"
  }else if (grepl("pigm", bwadfName)) {
    analysis = "pigmentation_loss"
  }else if (grepl("eye", bwadfName)) {
    analysis = "eye_loss"
  }else if (grepl("geno", bwadfName)) {
    analysis = "genome_size"
  }else if (grepl("long", bwadfName)) {
    analysis = "longevity"
  }

  filteredOgs = unique(bwadf$Orthogroup)
  fastaPath = paste0("./test_data/orthofinder_results/", analysis, "/Orthogroup_Sequences/")
  fastaList = list.files(fastaPath)
  fastaList = gsub(pattern = "\\.fa$", replacement = "", x = fastaList)

  fastaListFiltered = which(fastaList %in% filteredOgs)

  fastaListPath = list.files(fastaPath, full.names = TRUE)
  fastaListPath = fastaListPath[fastaListFiltered]

  fastas = lapply(fastaListPath, readAAStringSet)
  names(fastas) = filteredOgs

  ids = lapply(fastas, names)

  fastaFilter = lapply(1:length(ids), function(x, idIn, fastaIn, bwas){

    fasta = fastaIn[[x]]
    id = idIn[[x]]
    fastaName = names(fastaIn)[x]

    bwaog = bwas %>%
      filter(Orthogroup == fastaName)

    headers = bwaog$header
    idNums = which(id %in% headers)

    if (length(idNums) != length(headers)) {
      stop("Missing sequences!")
    }

    fastaOut = fasta[idNums]
    fastaOut

  }, idIn = ids, fastaIn = fastas, bwas = bwadf)

  names(fastaFilter) = filteredOgs



  fastaFilter

}, ogdfs = ogdfList)

names(orthogroups_filterOGs) = tablesNames2

#cds path



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
  df$expCon = nSp
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


cdsOGdf2 = cdsOGdf %>%
  group_by(expCon) %>%
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
  lst1name = names(lst1[x])
  out = lapply(1:length(lst1[[x]]), function(y, lst2, nsp){
    df = lst2[[y]]

    og = unique(df$orthogroup)
    dir.create(finalFastaOut)
    nspFolder = paste0(finalFastaOut, nsp, "/")
    dir.create(nspFolder)


    out = DNAStringSet(df$cdsSeq)
    names(out) = df$headers

    outfile = paste0(nspFolder, og, ".fa")

    writeFasta(out, outfile)
    out

  }, lst2 = lst1[[x]], nsp = lst1name)
  out
}, lst1 = cdsOGdf2)


