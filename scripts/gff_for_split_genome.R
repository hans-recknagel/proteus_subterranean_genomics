#script for adapting proteus gff annotation file to split genome

#Inputs are the gff file and a list of headers from the split genome fasta file. Output path can be set at the end of the script.

library(dplyr)
library(tidyr)

gffPath = "./path/to/ggf/file.gff"
splitHeadersPath = "./path/to/genome/fasta/headers/list/fasta.headers.txt"

gffIn = read.table(gffPath, header = FALSE, sep = "\t", col.names = c("chr", "model", "type", "start", "end", "score", "strand", "phase", "attributes"))
splitHeaders = read.table(splitHeadersPath, header = FALSE, col.names = "headersSplit")

splitHeadersSc = splitHeaders %>%
  mutate(chrSc = gsub(pattern = "\\d+.*", replacement = "", x = headersSplit)) %>%
  mutate(chrSc = gsub(pattern = ">", replacement = "", x = chrSc)) %>%
  filter(chrSc == "scaffold") %>%
  mutate(headersSplit = gsub(pattern = ">", replacement = "", x = headersSplit)) %>%
  mutate(chr = gsub(pattern = "_sliding.*", replacement = "", x = headersSplit))
  

gffOut = gffIn %>%
  #starting 
  mutate(chrStart = as.integer(start / 1E9) * 1E9 + 1) %>%
  mutate(chrStart = format(chrStart, scientific = FALSE)) %>%
  mutate(chrEnd = as.integer(end / 1E9) * 1E9 + 1E9) %>%
  mutate(chrEnd = format(chrEnd, scientific = FALSE)) %>%
  mutate(chrStart2 = as.integer(start / 1E9) * 1E9) %>%
  mutate(chrStart2 = format(chrStart2, scientific = FALSE)) %>%
  mutate(chrEnd2 = as.integer(end / 1E9) * 1E9) %>%
  mutate(chrEnd2 = format(chrEnd2, scientific = FALSE)) %>%
  mutate(sameChr = chrStart2 == chrEnd2) %>%
  mutate(chrStart2 = if_else(sameChr == FALSE, as.numeric(chrEnd2) + 1, as.numeric(chrStart2))) %>%
  mutate(chrEnd2 = if_else(sameChr == FALSE, as.numeric(chrEnd2) + 1E9, as.numeric(chrEnd2))) %>%
  mutate(chrEnd = if_else(sameChr == FALSE, as.numeric(chrEnd2) - 1E9, as.numeric(chrEnd))) %>%
  mutate(rn = row_number()) %>%
  mutate(count = if_else(sameChr == FALSE, 2, 1)) %>%
  uncount(count) %>%
  group_by(rn) %>%
  mutate(chrStart3 = 1:n()) %>%
  mutate(chrStart = if_else(chrStart3 == 2, as.numeric(chrStart) + 1E9, as.numeric(chrStart))) %>%
  mutate(chrEnd = if_else(chrStart3 == 2, as.numeric(chrEnd2), as.numeric(chrEnd))) %>%
  mutate(chrStart = format(chrStart, scientific = FALSE)) %>%
  mutate(chrEnd = format(chrEnd, scientific = FALSE)) %>%
  mutate(attributes = if_else(sameChr == FALSE, gsub(pattern = "evm.model", replacement = "split_evm.model", x = attributes), attributes)) %>%
  mutate(chrEnd = if_else(chr == "chr1" & chrEnd == "5000000000", "4107498685", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr2" & chrEnd == "4000000000", "3543809019", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr3" & chrEnd == "4000000000", "3217785387", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr4" & chrEnd == "4000000000", "3068009127", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr5" & chrEnd == "3000000000", "2645194432", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr6" & chrEnd == "3000000000", "2025613667", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr7" & chrEnd == "2000000000", "1981864174", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr8" & chrEnd == "2000000000", "1427536823", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr9" & chrEnd == "2000000000", "1412067188", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr10" & chrEnd == "2000000000", "1288210443", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr11" & chrEnd == "2000000000", "1221434848", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr12" & chrEnd == "2000000000", "1206315881", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr13" & chrEnd == "2000000000", "1068532332", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr14" & chrEnd == "1000000000", "978015207", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr15" & chrEnd == "1000000000", "935486047", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr16" & chrEnd == "1000000000", "821545095", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr17" & chrEnd == "1000000000", "819422292", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr18" & chrEnd == "1000000000", "802104798", chrEnd)) %>%
  mutate(chrEnd = if_else(chr == "chr19" & chrEnd == "1000000000", "799243929", chrEnd))

#several steps to increase efficiency
gffOut2 = gffOut %>%
  ungroup %>%
  select(-c(chrStart2, chrEnd2, rn)) %>%
  mutate(chrSc = gsub(pattern = "\\d+", replacement = "", x = chr)) %>%
  mutate(chrSc = gsub(pattern = "\\d+", replacement = "", x = chr)) %>%
  left_join(splitHeadersSc, by = "chr") %>%
  select(-chrSc.y) %>%
  mutate(chrStart = gsub(pattern = " ", replacement = "", x = chrStart)) %>%
  mutate(chrEnd = gsub(pattern = " ", replacement = "", x = chrEnd)) %>%
  mutate(newChr = if_else(chrSc.x == "chr", paste0(chr, "_sliding_", chrStart, "-", chrEnd), headersSplit)) %>%
  mutate(newStart = if_else(as.numeric(chrStart) > 1, as.numeric(start) - as.numeric(chrStart) + 1, as.numeric(start))) %>%
  mutate(newEnd = if_else(as.numeric(chrStart) > 1, as.numeric(end) - as.numeric(chrStart) + 1, as.numeric(end))) %>%
  mutate(newStart = if_else(chrStart3 == 2, 1, as.numeric(newStart))) %>%
  select(-headersSplit) %>%
  mutate(newEnd = if_else(sameChr == FALSE & chrStart3 == 1, 1E9, newEnd))

gffOut3 = gffOut2 %>%
  select(newChr, model, type, newStart, newEnd, score, strand, phase, attributes)
    
write.table(gffOut3, "path/to/output/split_file.gff", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)





















    
