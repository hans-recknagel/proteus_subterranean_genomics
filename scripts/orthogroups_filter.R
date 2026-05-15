#Functions in this script are called by other scripts, it is not run by itself.
#getting OG fiter in function

#orthogroups_filter

orthogroups_filter = function(OGtable, OGfastaPath, OGoutPath, writeOGs = TRUE){
  
  
  #libraries
  library(stringi)
  library(dplyr)
  library(ShortRead)
  
  #reading orthogroup table
  orthogroups = read.table(file = OGtable, sep = "\t", header = TRUE)
  
  #number of species in analysis
  nspecies = ncol(orthogroups) - 1
  
  #sequence hits per species in orthogroup
  hitsPerSp = apply(orthogroups[, 2:ncol(orthogroups)], 1, stri_split_fixed, ",", omit_empty = NA)
  
  #getting number of hits per species per orthogroup
  hitsPerSpDF = lapply(hitsPerSp, lapply, function(x){length(x[!is.na(x)])})
  hitsPerSpDF = as.data.frame(do.call(rbind, hitsPerSpDF))
  hitsPerSpDF = cbind(orthogroups$Orthogroup, hitsPerSpDF)
  
  colnames(hitsPerSpDF) = colnames(orthogroups)
  
  #adding columns for number of species with more than one, one and zero hits 
  greater1 = apply(hitsPerSpDF, c(1,2), function(x){x > 1})
  greater1 = apply(greater1, 1, sum, na.rm = TRUE)
  equal1 = apply(hitsPerSpDF, c(1,2), function(x){x == 1})
  equal1 = apply(equal1, 1, sum, na.rm = TRUE)
  equal0 = apply(hitsPerSpDF, c(1,2), function(x){x < 1})
  equal0 = apply(equal0, 1, sum, na.rm = TRUE)
  
  hitsPerSpDF$greater1 = greater1
  hitsPerSpDF$equal1 = equal1
  hitsPerSpDF$equal0 = equal0
  
  #getting 1:1 orthologs, somewhat redundant, they are on OF results anyway
  oneToOne = dplyr::filter(hitsPerSpDF, equal1 == nspecies)
  
  #minimum number of species in orthogroup
  minSpecies = round(nspecies * 0.75)
  
  #star number of species, one less than all, since all would be 1:1 orthogroups
  #nSp = nspecies - 1
  nSp = nspecies
  
  #making tables for 1 - n species down to minimum allowed number, each table separate
  while (nSp >= minSpecies) {
    spTable = filter(hitsPerSpDF, equal1 == nSp)
    #spTable = inner_join(spTable, orthogroups, "Orthogroup")
    assign(paste0("hits_per_", nSp, "_species"), spTable)
    nSp = nSp - 1
  }
  #removing redundant table
  rm(spTable)
  
  #getting all tables with hits per number of species in a list
  oghitsList = mget(ls(pattern = "hits_per_\\d+_species"))
  
  #loop for getting and writing sequences in filtered orthogroups (tits_per_X_species) to disk
  finalPercentageOGs = lapply(1:length(oghitsList), function(y, inList, oghitsListNames){
    
    #getting OGs with species with extra sequences
    ogFIlter = inList[[y]] %>%
      filter(greater1 > 0)
    
    if (nrow(ogFIlter) != 0) {
      #getting filtered orthogroups with one sequence per species
      filteredOGs = apply(ogFIlter, 1, function(x){
        
        #row imported as list, fixing that
        ogSpeciesOver1 = unlist(x)
        #which species in OG have more than one sequence
        ogSpeciesOver1 = which(ogSpeciesOver1[1:(1+nspecies)] > 1)
        
        #column names for species with more than one sequence
        ogSpeciesOver1Names = colnames(ogFIlter)[ogSpeciesOver1]
        #OG name
        ogSpeciesOver1OG = x[[1]]
        
        #getting sequence headers for species with Y 1 sequence per species 
        ogSpeciesOver1Actual = orthogroups %>%
          filter(Orthogroup == ogSpeciesOver1OG) %>%
          select(ogSpeciesOver1Names)
        
        #reading orthogroups
        ogPath = paste0(OGfastaPath, ogSpeciesOver1OG, ".fa")
        og = readAAStringSet(ogPath)
        
        #rettong OG headers
        ogNames = names(og)
        #getting OG headers beginnings
        ogNamesStart = gsub(pattern = " .*$", replacement = "", x = ogNames)
        
        #getting multiple sequence per species headers in a vector
        dupls = c(ogSpeciesOver1Actual)
        dupls = unlist(lapply(dupls, stri_split_fixed, ","))
        dupls = gsub(pattern = " ", replacement = "", x = dupls)
        
        #which are not multiple sequence per vector
        uniqSeqs = which(!(ogNamesStart %in% dupls))
        
        #getting final sequences and filtered OGs
        ogOut = og[uniqSeqs]
        ogOut
      })
      
      #naming new list
      names(filteredOGs) = ogFIlter$Orthogroup
      
      #getting one sequence and no sequence species
      ogFilter2 = inList[[y]] %>%
        filter(greater1 == 0)
      
      #reading in the rest of the OGs
      zeroOGs = lapply(ogFilter2$Orthogroup, function(x){
        
        ogPath = paste0(OGfastaPath, x, ".fa")
        og = readAAStringSet(ogPath)
        og
        
      })
      
      #combining all OGs in one list
      names(zeroOGs) = ogFilter2$Orthogroup
      outOGs = append(filteredOGs, zeroOGs)
      
    }else{
      
      #getting one sequence and no sequence species
      ogFilter2 = inList[[y]] %>%
        filter(greater1 == 0)
      
      #reading in the rest of the OGs
      zeroOGs = lapply(ogFilter2$Orthogroup, function(x){
        
        ogPath = paste0(OGfastaPath, x, ".fa")
        og = readAAStringSet(ogPath)
        og
        
      })
      
      #combining all OGs in one list
      names(zeroOGs) = ogFilter2$Orthogroup
      outOGs = zeroOGs
      
    }
    
    
    
    #writing OGs
    if (writeOGs == TRUE) {
      lapply(1:length(outOGs), function(x, fastaList){
        
        #og name
        ogListName = oghitsListNames[[y]]
        
        #new directory
        dirName = paste0(OGoutPath, ogListName)
        dir.create(dirName)
        
        #writing fata file in new directory
        outName = names(fastaList)[[x]]
        outFile = paste0(OGoutPath, ogListName, "/", outName, "_filtered.fa")
        
        writeXStringSet(fastaList[[x]], outFile)
        
      }, fastaList = outOGs)
    }
    
    #returning OGs
    outOGs
    
  }, inList = oghitsList, oghitsListNames = names(oghitsList))
  
  names(finalPercentageOGs) = names(oghitsList)
  return(finalPercentageOGs)
  
}

orthogroups_filter_percentage = function(OGtable, OGfastaPath, OGoutPath){
  
  
  #libraries
  library(stringi)
  library(dplyr)
  library(ShortRead)
  
  #reading orthogroup table
  orthogroups = read.table(file = OGtable, sep = "\t", header = TRUE)
  
  #number of species in analysis
  nspecies = ncol(orthogroups) - 1
  
  #sequence hits per species in orthogroup
  hitsPerSp = apply(orthogroups[, 2:ncol(orthogroups)], 1, stri_split_fixed, ",", omit_empty = NA)
  
  #getting number of hits per species per orthogroup
  hitsPerSpDF = lapply(hitsPerSp, lapply, function(x){length(x[!is.na(x)])})
  hitsPerSpDF = as.data.frame(do.call(rbind, hitsPerSpDF))
  hitsPerSpDF = cbind(orthogroups$Orthogroup, hitsPerSpDF)
  
  colnames(hitsPerSpDF) = colnames(orthogroups)
  
  #adding columns for number of species with more than one, one and zero hits 
  greater1 = apply(hitsPerSpDF, c(1,2), function(x){x > 1})
  greater1 = apply(greater1, 1, sum, na.rm = TRUE)
  equal1 = apply(hitsPerSpDF, c(1,2), function(x){x == 1})
  equal1 = apply(equal1, 1, sum, na.rm = TRUE)
  equal0 = apply(hitsPerSpDF, c(1,2), function(x){x < 1})
  equal0 = apply(equal0, 1, sum, na.rm = TRUE)
  
  hitsPerSpDF$greater1 = greater1
  hitsPerSpDF$equal1 = equal1
  hitsPerSpDF$equal0 = equal0
  
  #getting 1:1 orthologs, somewhat redundant, they are on OF results anyway
  oneToOne = dplyr::filter(hitsPerSpDF, equal1 == nspecies)
  
  #minimum number of species in orthogroup
  minSpecies = round(nspecies * 0.75)
  
  #star number of species, one less than all, since all would be 1:1 orthogroups
  #nSp = nspecies - 1
  nSp = nspecies
  
  #making tables for 1 - n species down to minimum allowed number, each table separate
  while (nSp >= minSpecies) {
    spTable = filter(hitsPerSpDF, equal1 == nSp)
    #spTable = inner_join(spTable, orthogroups, "Orthogroup")
    assign(paste0("hits_per_", nSp, "_species"), spTable)
    nSp = nSp - 1
  }
  #removing redundant table
  rm(spTable)
  
  hitsPerSpDFFiltered = hitsPerSpDF %>%
    filter(equal1 >= minSpecies)
  
  percMissingTotal = hitsPerSpDF %>%
    filter(equal1 >= minSpecies) %>%
    group_by(equal1) %>%
    summarise(n1 = n()) %>%
    mutate(diff = max(equal1) - equal1) %>%
    mutate(multiSp = diff * n1) %>%
    summarise(perc = sum(multiSp) / (max(equal1) * sum(n1))) %>%
    pull()
  
  percMissingGr1 = hitsPerSpDF %>%
    filter(equal1 >= minSpecies) %>%
    group_by(greater1) %>%
    summarise(n1 = n()) %>%
    mutate(multiSp = greater1 * n1) %>%
    summarise(perc = sum(multiSp) / (max(hitsPerSpDFFiltered$equal1) * nrow(hitsPerSpDFFiltered))) %>%
    pull()
  
  percMissingGrEq0= hitsPerSpDF %>%
    filter(equal1 >= minSpecies) %>%
    group_by(equal0) %>%
    summarise(n1 = n()) %>%
    mutate(multiSp = equal0 * n1) %>%
    summarise(perc = sum(multiSp) / (max(hitsPerSpDFFiltered$equal1) * nrow(hitsPerSpDFFiltered))) %>%
    pull()  
  
  out = data.frame("percMissingTotal" = percMissingTotal, "percMissingGr1" = percMissingGr1, "percMissingGrEq0" = percMissingGrEq0)
  
  return(out)
  
}

species_headers = function(OGtable = "./test_data/orthogroups_tables/cave_species_Orthogroups.tsv", OGfastaPath = "./test_data/orthofinder_results_old_headers/cave_species_fixed_headers_results/cave_species_fixed_headers_results"){
  
  #libraries
  library(tidyr)
  library(stringi)
  library(dplyr)
  library(ShortRead)
  
  #reading orthogroup table
  orthogroups = read.table(file = OGtable, sep = "\t", header = TRUE)
  #added for right order of species
  orthogroupsSortCols = c(1, order(colnames(orthogroups)[2:ncol(orthogroups)]) + 1)
  orthogroups = orthogroups[,orthogroupsSortCols]
  
  #number of species in analysis
  nspecies = ncol(orthogroups) - 1
  
  #sequence hits per species in orthogroup
  hitsPerSp = apply(orthogroups[, 2:ncol(orthogroups)], 1, stri_split_fixed, ",", omit_empty = NA)
  
  #getting number of hits per species per orthogroup
  hitsPerSpDF = lapply(hitsPerSp, lapply, function(x){length(x[!is.na(x)])})
  hitsPerSpDF = as.data.frame(do.call(rbind, hitsPerSpDF))
  hitsPerSpDF = cbind(orthogroups$Orthogroup, hitsPerSpDF)
  
  colnames(hitsPerSpDF) = colnames(orthogroups)
  
  #adding columns for number of species with more than one, one and zero hits 
  greater1 = apply(hitsPerSpDF, c(1,2), function(x){x > 1})
  greater1 = apply(greater1, 1, sum, na.rm = TRUE)
  equal1 = apply(hitsPerSpDF, c(1,2), function(x){x == 1})
  equal1 = apply(equal1, 1, sum, na.rm = TRUE)
  equal0 = apply(hitsPerSpDF, c(1,2), function(x){x < 1})
  equal0 = apply(equal0, 1, sum, na.rm = TRUE)
  
  hitsPerSpDF$greater1 = greater1
  hitsPerSpDF$equal1 = equal1
  hitsPerSpDF$equal0 = equal0
  
  #getting 1:1 orthologs, somewhat redundant, they are on OF results anyway
  oneToOne = dplyr::filter(hitsPerSpDF, equal1 == nspecies)
  
  #minimum number of species in orthogroup
  minSpecies = round(nspecies * 0.75)
  
  #star number of species, one less than all, since all would be 1:1 orthogroups
  #nSp = nspecies - 1
  nSp = nspecies
  
  #making tables for 1 - n species down to minimum allowed number, each table separate
  while (nSp >= minSpecies) {
    spTable = filter(hitsPerSpDF, equal1 == nSp)
    #spTable = inner_join(spTable, orthogroups, "Orthogroup")
    assign(paste0("hits_per_", nSp, "_species"), spTable)
    nSp = nSp - 1
  }
  #removing redundant table
  rm(spTable)
  
  #getting all tables with hits per number of species in a list
  oghitsList = mget(ls(pattern = "hits_per_\\d+_species"))
  
  #loop for getting and writing sequences in filtered orthogroups (tits_per_X_species) to disk
  finalPercentageOGs = lapply(1:length(oghitsList), function(y, inList, oghitsListNames){
    
    #getting OGs with species with extra sequences
    
    ogFIlter = inList[[y]] %>%
      filter(greater1 > 0)
    
    #print(names(inList)[[y]])
    
    if (nrow(ogFIlter) != 0) {
      
      #getting filtered orthogroups with one sequence per species
      filteredOGs = apply(ogFIlter, 1, function(x){
        
        #row imported as list, fixing that
        ogSpeciesOver1 = unlist(x)
        #which species in OG have more than one sequence
        ogSpeciesOver1 = which(ogSpeciesOver1[1:(1+nspecies)] > 1)
        
        #column names for species with more than one sequence
        ogSpeciesOver1Names = colnames(ogFIlter)[ogSpeciesOver1]
        #OG name
        ogSpeciesOver1OG = x[[1]]
        
        #getting sequence headers for species with Y 1 sequence per species 
        ogSpeciesOver1Actual = orthogroups %>%
          filter(Orthogroup == ogSpeciesOver1OG) %>%
          select(ogSpeciesOver1Names)
        
        #reading orthogroups
        ogPath = paste0(OGfastaPath, ogSpeciesOver1OG, ".fa")
        #print(ogPath)
        
        og = readAAStringSet(ogPath)
        
        #getting OG headers
        ogNames = names(og)
        #getting OG headers beginnings
        ogNamesStart = gsub(pattern = " .*$", replacement = "", x = ogNames)
        
        #getting multiple sequence per species headers in a vector
        dupls = c(ogSpeciesOver1Actual)
        dupls = unlist(lapply(dupls, stri_split_fixed, ","))
        dupls = gsub(pattern = " ", replacement = "", x = dupls)
        
        #which are not multiple sequence per vector
        uniqSeqs = which(!(ogNamesStart %in% dupls))
        
        #getting final sequences and filtered OGs
        ogOut = og[uniqSeqs]
        ogOut
      })
      
      #naming new list
      names(filteredOGs) = ogFIlter$Orthogroup
      
      #getting one sequence and no sequence species
      ogFilter2 = inList[[y]] %>%
        filter(greater1 == 0)
      
      #reading in the rest of the OGs
      zeroOGs = lapply(ogFilter2$Orthogroup, function(x){
        
        ogPath = paste0(OGfastaPath, x, ".fa")
        og = readAAStringSet(ogPath)
        og
        
      })
      
      #combining all OGs in one list
      names(zeroOGs) = ogFilter2$Orthogroup
      outOGs = append(filteredOGs, zeroOGs)
      
      #returning OGs
      outOGs
    }else{
      #getting one sequence and no sequence species
      ogFilter2 = inList[[y]] %>%
        filter(greater1 == 0)
      
      #reading in the rest of the OGs
      zeroOGs = lapply(ogFilter2$Orthogroup, function(x){
        
        ogPath = paste0(OGfastaPath, x, ".fa")
        og = readAAStringSet(ogPath)
        og
        
      })
      
      #combining all OGs in one list
      names(zeroOGs) = ogFilter2$Orthogroup
      outOGs = zeroOGs
      
      #returning OGs
      outOGs
    }
    
    
    
  }, inList = oghitsList, oghitsListNames = names(oghitsList))
  
  #making a long table for all hits tables
  logTablesHits = lapply(oghitsList, function(x){
    
    out = x %>%
      select(-c("greater1", "equal1", "equal0")) %>%
      pivot_longer(cols = !Orthogroup) %>%
      filter(value == 1)
    out
    
  })
  
  #naming fasta list
  names(finalPercentageOGs) = names(oghitsList)
  
  #making data frames from fastas
  finalPercentageOGdfs1 = lapply(finalPercentageOGs, lapply, as.data.frame)
  
  #adding OG info and headers to fasta data frames
  finalPercentageOGdfs2 = lapply(1:length(finalPercentageOGdfs1), function(y, list1){
    
    insideDf = lapply(1:length(list1[[y]]), function(x, list2){
      
      ogName = names(list2)[[x]]
      df = list2[[x]]
      df$headers = row.names(df)
      df$og = ogName
      
      df
      
    }, list2 = list1[[y]])
    insideDf
  }, list1 = finalPercentageOGdfs1 )
  
  names(finalPercentageOGdfs2) = names(oghitsList)
  
  #combining into lists
  finalPercentageOGdfsF = lapply(finalPercentageOGdfs2, do.call, what = rbind)
  
  #combining
  spHeaderList = lapply(1:length(finalPercentageOGdfsF), function(x, headerDFs, speciesDFs){
    
    headDF = headerDFs[[x]]
    spDF = speciesDFs[[x]]
    speciesN = names(headerDFs)[[x]]
    
    out = cbind(headDF, spDF)
    out = out[,c("headers", "name", "og")]
    out$Nspecies = speciesN
    out
  }, headerDFs = finalPercentageOGdfsF, speciesDFs = logTablesHits)
  
  #putting all in one table
  spHeaderDF = do.call(rbind, spHeaderList)
  
  return(spHeaderDF)
}
