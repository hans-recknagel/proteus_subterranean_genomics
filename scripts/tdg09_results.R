#script for putting tdg09 results in table form 

library(yaml)


analysis = "analysis_name"

resPath = paste0("/path/to/tdg09/results", analysis, "/tdg09_out")

tdg09files = list.files(resPath, full.names = TRUE)
tdg09filesNames = list.files(resPath)

tdg09filesNames = gsub(pattern = ".tdg09.out", replacement = "", x = tdg09filesNames)

yamlFiles = lapply(tdg09files,  yaml.load_file)
names(yamlFiles) = tdg09filesNames

yamlLen = lapply(yamlFiles, length)
yamlLen = unlist(yamlLen)
yamlLenInts = which(yamlLen %in% 17)

yamlFilesFilter = yamlFiles[yamlLenInts]
 counter = 0

lrt_results = lapply(yamlFilesFilter, function(x){
   counter <<- counter + 1
  out = as.data.frame(matrix(unlist(x$LrtResults), ncol=5, byrow=T))
  return(out)
})

counter = 0
lrt_resultsFinal = lapply(1:length(lrt_results), function(x, yamls){
  counter <<- counter + 1
  df = yamls[[x]]
  dfName = names(yamls[x])
  
  names(df) = c("site", "deltaLnL", "dof", "lrt", "fdr")
  df$og = dfName
  return(df)
  
}, yamls = lrt_results)


dfOut = do.call(rbind, lrt_resultsFinal)

outFile = paste0("./path/to/output/tdg09_outdf_", analysis, ".csv")

write.csv(dfOut, file = outFile, quote = FALSE, row.names = FALSE)
