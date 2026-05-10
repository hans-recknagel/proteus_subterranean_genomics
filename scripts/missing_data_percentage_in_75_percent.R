
source("./scripts/orthogroups_filter.R")

cave_species = orthogroups_filter_percentage(OGtable = "./path/to/orthogroups_tables/Orthogroups.tsv", OGfastaPath = "./path/to/Orthogroup_Sequences/", OGoutPath = "./path/to/filtered/orthogroups/")
eye_loss = orthogroups_filter_percentage(OGtable = "./path/to/orthogroups_tables/Orthogroups.tsv", OGfastaPath = "./path/to/Orthogroup_Sequences/", OGoutPath = "./path/to/filtered/orthogroups/")
longevity = orthogroups_filter_percentage(OGtable = "./path/to/orthogroups_tables/Orthogroups.tsv", OGfastaPath = "./path/to/Orthogroup_Sequences/", OGoutPath = "./path/to/filtered/orthogroups/")
pigmentation_loss = orthogroups_filter_percentage(OGtable = "./path/to/orthogroups_tables/Orthogroups.tsv", OGfastaPath = "./path/to/Orthogroup_Sequences/", OGoutPath = "./path/to/filtered/orthogroups/")


listPerc = list("cave_species" = cave_species, "eye_loss" = eye_loss, "longevity" = longevity, "pigmentation_loss" = pigmentation_loss)
dfPerc = do.call(rbind, listPerc)
dfPerc$analysis = rownames(dfPerc)

write.csv(dfPerc, file = "species_perc_missing_data.csv", row.names = FALSE, quote = FALSE)
