

library(jsonlite)
library(dplyr)

analysis = "analysis"
pathToAnalysis = paste0("path/to/absrel/results/", analysis)

allJsonAbsrel = list.files(pathToAnalysis, pattern = ".json", recursive = TRUE)
allJsonAbsrelFullNames = list.files(pathToAnalysis, pattern = ".json", recursive = TRUE, full.names = TRUE)

#uncomment for the correct analysis
#cave_species
#foreground = c("Asmec", "Luden", "Prang", "Sians", "Sirhi", "Trros")
#background = c("Ammex", "Ammel", "Asmes", "Darer", "Gamor", "Hosap", "Leocu", "Mymur", "Nebri", "Plwal", "Sigra", "Tagut", "Thmac", "Xetro", "Zoviv", "Trtib")

#longevity
#foreground = c("Acbae", "Cacar", "Prang", "Elmax", "Hosap", "Lacha", "Myluc", "Hegla", "Ororc", "Rhtyp", "Tecar")
#background = c("Ancor", "Botau", "Capor", "Darer", "Chsab", "Chpic", "Cocri", "Ercal", "Modom", "Ovari", "Ranor", "Sasal", "Scvul", "Tagut")

#eye_loss
#foreground = c("Asmec", "Chasi", "Luden", "Cocri", "Geser", "Hegla", "Miuni", "Nagal", "Prang", "Sians", "Trros")
#background = c("Darer", "Ranor", "Elmax", "Ammex", "Ammel", "Asmes", "Bubuf", "Chlan", "Ectel", "Ereur", "Gamor", "Jajac", "Mymur", "Ornil", "Orcun", "Plwal", "Sigra", "Soara", "Thmac", "Trtib", "Sirhi")

#pigmentation_loss
#foreground = c("Ammex", "Asmec", "Hegla", "Prang", "Luden", "Sirhi", "Sians", "Trros")
#background = c("Ammel", "Asmes", "Darer", "Gamor", "Hosap", "Leocu", "Mymur", "Plwal", "Sigra", "Tagut", "Thmac", "Xetro", "Zoviv", "Trtib", "Chlan", "Ornil", "Scvul")

#genome_size
foreground = c("Plwal", "Ammex", "Prang", "Bubuf", "Nefor", "Prann", "Sctor")
background = c("Asmes", "Asmec", "Darer", "Hosap", "Leocu", "Tagut", "Ornil", "Camil", "Enpus", "Lacha", "Menas", "Pemar", "Pipip", "Xetro", "Zoviv")



allAbsrelFiles = lapply(allJsonAbsrelFullNames, fromJSON)
names(allAbsrelFiles) = allJsonAbsrel

allAbsrelResults = lapply(1:length(allAbsrelFiles), function(x, jsondata){
  
  res = jsondata[[x]]
  resName = names(jsondata)[[x]]
  outStats = res$`branch attributes`$`0`
  speciesOnlyResults = grep(pattern = "Node", x = names(outStats), invert = TRUE)
  outStats = outStats[speciesOnlyResults]
  
  
  
  out = lapply(1:length(outStats), function(y, datain){
    element = datain[[y]]
    pvalue = element$`Corrected P-value`
    return(pvalue)
  }, datain = outStats)
  
  species = gsub(pattern = "^.*_(\\w{5}$)", replacement = "\\1", x = names(outStats))
  
  names(out) = species
  
  hits = gsub(pattern = "/.*$", replacement = "", x = resName)
  orthogroup = gsub(pattern = ".*(OG\\d+)\\.fa.*", replacement = "\\1", x = resName)
  
  out2 = as.data.frame(do.call(cbind, out))
  out2$og = orthogroup
  out2$nSp = hits
  
  return(out2)
}, jsondata = allAbsrelFiles)

allAbsrelResultsDF = do.call(bind_rows, allAbsrelResults)

rm(allAbsrelFiles)

allAbsrelResultsDF = allAbsrelResultsDF[,order(colnames(allAbsrelResultsDF))]



allAbsrelResultsDF2 = allAbsrelResultsDF %>%
  relocate(nSp, .after = last_col()) %>%
  relocate(og, .after = last_col())

allAbsrelResultsDF2f = allAbsrelResultsDF2[, foreground]
allAbsrelResultsDF2b = allAbsrelResultsDF2[, background]

foreSelection = rowSums(allAbsrelResultsDF2f < 0.05, na.rm = TRUE)
backSelection = rowSums(allAbsrelResultsDF2b < 0.05, na.rm = TRUE)

foreAll = rowSums(!is.na(allAbsrelResultsDF2f), na.rm = TRUE)
backAll = rowSums(!is.na(allAbsrelResultsDF2b), na.rm = TRUE)

#foreIns = rowSums(allAbsrelResultsDF2f == 1, na.rm = TRUE)
#backIns = rowSums(allAbsrelResultsDF2b == 1, na.rm = TRUE)

allAbsrelResultsDF2$selectionF = foreSelection
allAbsrelResultsDF2$selectionB = backSelection
allAbsrelResultsDF2$allF = foreAll
allAbsrelResultsDF2$allB = backAll


write.csv(allAbsrelResultsDF2, file = paste0("/out/path/", analysis, "_absrel.csv"), quote = FALSE, row.names = FALSE)


