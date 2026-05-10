library(parallel)
library(stringr)

ogTablePath = "/path/to/orthogroups_tables/Orthogroups.tsv"

ogTablesPerSpecies = read.table(ogTablePath, sep = "\t", header = TRUE)

nonNcbiColsAll = c("Ammex.mark", "Darer.mark", "Nefor.mark", "Plwal.mark", "Prang.mark", "Acbae.mark", "Menas.mark", "Pipip.mark", "Scvul.mark")

colnamesAll = which(colnames(ogTablesPerSpecies) %in% nonNcbiColsAll)
ogTablesPerSpecies2 = ogTablesPerSpecies[,-colnamesAll]

cl = makeCluster(4)
clusterEvalQ(cl = cl, library(stringr))
ogTablesPerSpecies3 = parApply(cl = cl, ogTablesPerSpecies2, c(1,2), str_replace_all, pattern = "_[:upper:]{1}[:lower:]{4}", replacement = "")
stopCluster(cl)


write.table(ogTablesPerSpecies3, "./path/to/orthogroups_tables/orthogroups_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
