###hyphy json result to table

library(jsonlite)
library(dplyr)

analysis = "analysis_name"
path = paste0("/path/to/relax/results/", analysis)
#relax

allJsonRelax = list.files(path, recursive = TRUE)
allJsonRelaxFullNames = list.files(path, recursive = TRUE, full.names = TRUE)
allRuns = list.files(path)


counterx = 0
countery = 0
relaxRunDFs = lapply(1:length(allRuns), function(x, allnames, fullnames, runs){
    
  counterx <<- x
  run = runs[[x]]
  runNums = grep(pattern = run, x = fullnames)
  fullnamesx = fullnames[runNums]
  allnamesx = allnames[runNums]
  countery = 0
  
  allRelaxFiles = lapply(fullnamesx, fromJSON)


  
  names(allRelaxFiles) = allnamesx
  allRelaxResults = lapply(1:length(allRelaxFiles), function(y, jsondata){
      
    countery <<- y
    res = jsondata[[y]]
    resName = names(jsondata)[[y]]
    outStats = res$`test results`
    out = as.data.frame(do.call(cbind, outStats))
    hits = gsub(pattern = ".*/(hits_per_\\d{2}_species)/.*", replacement = "\\1", x = resName)
    hits = gsub(pattern = "hits_per_(\\d{2})_species", replacement = "\\1", x = hits)
    
    orthogroup = gsub(pattern = ".*(OG\\d+)\\.fa.*", replacement = "\\1", x = resName)
    
    out$nSp = hits
    out$og = orthogroup
    
    proteusName = names(res$`branch attributes`$`0`)
    proteusGrepl = grepl(pattern = "Prang", x = proteusName)
    proteusIn = TRUE %in% proteusGrepl
    
    out$proteus = proteusIn
    
    return(out)
    
    
  }, jsondata = allRelaxFiles)
  
  allRelaxResultsDF = do.call(rbind, allRelaxResults)
  
  rm(allRelaxFiles)
  
  names(allRelaxResultsDF) = c("LRT", "pvalue", "K", "nSp", "og", "proteus")
  
  allRelaxResultsDF2 = allRelaxResultsDF %>%
    mutate(selection = if_else(pvalue > 0.05, "none", if_else(K < 1, "relaxed", "intensified")))
  
  colnames(allRelaxResultsDF2) = paste0(colnames(allRelaxResultsDF2), "_", x)
  
  return(allRelaxResultsDF2)
  
}, allnames = allJsonRelax, fullnames = allJsonRelaxFullNames, runs = allRuns)

relaxRunDF = do.call(cbind, relaxRunDFs)
selectionCol = relaxRunDF[grep(pattern = "selection", x = colnames(relaxRunDF))]
selectionCol$intensified = rowSums(selectionCol == "intensified")
selectionCol$relaxed = rowSums(selectionCol == "relaxed")
selectionCol$none = rowSums(selectionCol == "none")


selectionCol$selection = apply(selectionCol[,6:8], 1, function(x){
  out = which(x == max(x))
  
  if (length(out) > 1) {
    out = "unknown"
  }else{
    out = names(out)
  }
  
} )

selectionCol = select(selectionCol, c("intensified", "relaxed", "none", "selection"))
relaxRunDF = cbind(relaxRunDF, selectionCol)

means = function(value){
  out = relaxRunDF[grep(pattern = value, x = colnames(relaxRunDF), ignore.case = FALSE)]
  out2 = rowMeans(out)
  return(out2)
  
}

pvalue_mean = means("pvalue")
LRT_mean = means("LRT")
K_mean = means("K")

relaxRunDF$pvalue_mean = pvalue_mean
relaxRunDF$LRT_mean = LRT_mean
relaxRunDF$K_mean = K_mean


outname = paste0("pot/path/hyphy_RELAX_table_fixed_added_", analysis, ".csv")
