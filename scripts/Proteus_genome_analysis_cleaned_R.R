########################### PROTEUS SUBTERRANEAN ADAPTATION IN GENOMES ###########################

library(ape)
library(biomaRt)
library(caper)
library(chromPlot)
library(circlize)
library(clusterProfiler)
library(data.table) 
library(dplyr)
library(enrichplot)
library(geiger)
library(GenomicRanges)
library(ggimage)
library(ggstar)
library(ggtree)
library(ggtreeExtra)
library(gridExtra)
library(httr)
library(imager)
library(jsonlite)
library(lattice)
library(limma)
library(lintr)
library(magick)
library(nlme)
library(patchwork)
library(phylolm)
library(phytools)
library(readr)
library(Rsamtools)
library(rstatix)
library(SuperExactTest)
library(TDbook)
library(tidyr)
library(tidyverse)
library(topGO)
library(venn) 
library(VennDiagram)
library(zoo)

######### GENE FAMILY EXPANSION & CONTRACTION ######### 
### prepare ultrametric trees, downloaded from timetree

### CAVE SPECIES ### 
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/01_cave_species")
cav.spec <- phytools::read.newick("cave_species_list_mod.tre")
new.tree <- force.ultrametric(cav.spec, method=c("extend"))
plot(new.tree)
plot(cav.spec)
is.ultrametric(new.tree)
write.tree(new.tree, file="cave_species_list_mod_corr_R.tre")
cav.spec <- phytools::read.newick("cave_species_list_mod_corr_newick.txt")
cav.spec <- phytools::read.newick("cave_species_list_mod_corr_newick_asty_mod.txt")
plot(cav.spec)
is.ultrametric(cav.spec)
new.tree <- force.ultrametric(cav.spec, method=c("extend"))
write.tree(new.tree, file="cave_species_list_mod_corr_newick_asty_mod.tre")
cav.spec <- phytools::read.newick("cave_species_list_mod_corr_newick_asty_mod.tre")
plot(cav.spec)
is.ultrametric(cav.spec)

### PIGMENTATION LOSS ### 
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/02_pigmentation_loss")
pigm.loss <- phytools::read.newick("pigmentation_species_list_mod.tre")
is.ultrametric(pigm.loss)
new.tree <- force.ultrametric(pigm.loss, method=c("extend"))
plot(new.tree)
plot(pigm.loss)
is.ultrametric(new.tree)
write.tree(new.tree, file="pigmentation_species_list_mod_corr_R.tre")

### EYE LOSS ### 
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/03_eye_loss")
eye.loss <- phytools::read.newick("eye_loss_species_list.tre")
is.ultrametric(eye.loss)

### LONGEVITY ### 
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/04_longevity")
longevity <- phytools::read.newick("longevity_species_list_abbrev.tre")
is.ultrametric(longevity)
plot(longevity)
new.tree <- force.ultrametric(longevity, method=c("extend"))
plot(new.tree)
is.ultrametric(new.tree)
write.tree(new.tree, file="longevity_species_list_abbrev_corr_R.tre")

### test different lambdas for sharks, fishes, tetrapods, salamanders, mammals
## choose model with lowest likelihood value
#### test if number of expansions & contractions differ between groups (e.g. pigmented vs. non-pigmented)
## cave species: 355209, 353926, 349461; lambda 3 best
## pigmentation:378324, 376598, 370388;  lambda 3
## eye loss: 407743, 404792, 403945, 398884; lambda 4
## longevity: 255284, 254601, 251561, 253678; lambda 3

#### EXPANSION VS CONTRACTION DIFFERENCES (PER SE) #### 
### CAVE SPECIES ###
# remove hashtag (#) before Taxon_ID before loading table
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/01_cave_species/03_filtered_50_lambda3")
data <- read.table("Base_clade_results.txt", header=TRUE)
data$Taxon_ID <- sub("^<+", "", data$Taxon_ID)
data$Taxon_ID <- gsub("<.*", "", data$Taxon_ID)
data$Taxon_ID <- sub(">", "", data$Taxon_ID)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/01_cave_species/")
FvsB <- read.csv("F_vs_B_cave_species.csv", header=TRUE)
data <- merge(data,FvsB,by="Taxon_ID")
data$DvsI<-data$Decrease/data$Increase
t.test(Increase~F.node, var.equal = FALSE, data=data)
t.test(Decrease~F.node, var.equal = FALSE, data=data)
t.test(DvsI~F.node, var.equal = FALSE, data=data) #*0.02349
data <- data[order(data$F.node, decreasing = TRUE),]
#write.csv(data, "summary_tree_cave_species_leaf_nodes_overview.csv", row.names=F)
# extract mean
data.sum.cave.leaf <- data %>% group_by(F.node) %>% 
  summarise_at(vars("Increase", "Decrease", "DvsI"), mean)
### create boxplot
# grouped boxplot of absolute difference
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/01_cave_species/")
FvsB <- read.csv("summary_tree_cave_species_leaf_nodes_overview_for_R_no_outliers.csv", header=TRUE)
str(FvsB)
p1 <- ggplot(FvsB, aes(x=cave_or_not, y=Number, fill=gene_family_fate)) + 
  geom_boxplot(coef=2) +
  scale_fill_manual(values=c("#56B4E9", "#F5B452")) + 
  scale_y_continuous(name = "Number of gene families", limits = c(0, 4000)) +
  theme_classic()
# difference in ratio of expanded/contracted
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/01_cave_species/")
FvsB <- read.csv("summary_tree_cave_species_leaf_nodes_overview.csv", header=TRUE)
str(FvsB)
p2 <- ggplot(FvsB, aes(x=F.node, y=IvsD, fill=F.node)) + 
  geom_boxplot(coef=2) +
  scale_y_continuous(name = "Rate of expanded to contracted gene families", limits = c(0, 3)) +
  scale_fill_manual(values=c("#8C8C8C", "#ECD5E0")) + 
  theme_classic()

#### FIGURE 2b,c ####
grid.arrange(p1,p2, ncol=2)

### PIGMENTATION LOSS ###
## pigmentation:378324, 376598, 370388;  lambda 3
# remove hashtag (#) before Taxon_ID before loading table
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/02_pigmentation_loss/04_filtered_50_lambda3")
data <- read.table("Base_clade_results.txt", header=TRUE)
str(data)
data$Taxon_ID <- sub("^<+", "", data$Taxon_ID)
data$Taxon_ID <- gsub("<.*", "", data$Taxon_ID)
data$Taxon_ID <- sub(">", "", data$Taxon_ID)
str(data)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/02_pigmentation_loss/")
FvsB <- read.csv("F_vs_B_pigmentation_loss.csv", header=TRUE)
str(FvsB)
data <- merge(data,FvsB,by="Taxon_ID")
str(data)
data$DvsI<-data$Decrease/data$Increase
str(data)
t.test(Increase~F.node, var.equal = FALSE, data=data)
t.test(Decrease~F.node, var.equal = FALSE, data=data)
t.test(DvsI~F.node, var.equal = FALSE, data=data) #*0.0121
data <- data[order(data$F.node, decreasing = TRUE),]
#write.csv(data, "summary_tree_pigmentation_loss_leaf_nodes_overview.csv", row.names=F)
# extract mean
data.sum.pig.leaf <- data %>% group_by(F.node) %>% 
  summarise_at(vars("Increase", "Decrease", "DvsI"), mean)

### EYE LOSS ###
## eye loss: 407743, 404792, 403945, 398884; lambda 4
# remove hashtag (#) before Taxon_ID before loading table
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/03_eye_loss/05_filtered_50_lambda4")
data <- read.table("Base_clade_results.txt", header=TRUE)
str(data)
data$Taxon_ID <- sub("^<+", "", data$Taxon_ID)
data$Taxon_ID <- gsub("<.*", "", data$Taxon_ID)
data$Taxon_ID <- sub(">", "", data$Taxon_ID)
str(data)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/03_eye_loss/")
FvsB <- read.csv("F_vs_B_eye_loss.csv", header=TRUE)
str(FvsB)
data <- merge(data,FvsB,by="Taxon_ID")
str(data)
data$DvsI<-data$Decrease/data$Increase
str(data)
t.test(Increase~F.node, var.equal = FALSE, data=data)
t.test(Decrease~F.node, var.equal = FALSE, data=data)
t.test(DvsI~F.node, var.equal = FALSE, data=data) 
data <- data[order(data$F.node, decreasing = TRUE),]
#write.csv(data, "summary_tree_eye_loss_leaf_nodes_overview.csv", row.names=F)
# extract mean
data.sum.eye.leaf <- data %>% group_by(F.node) %>% 
  summarise_at(vars("Increase", "Decrease", "DvsI"), mean)

### LONGEVITY ###
## longevity: 255284, 254601, 251561, 253678; lambda 3
# remove hashtag (#) before Taxon_ID before loading table
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/04_longevity/04_filtered_50_lambda3")
data <- read.table("Base_clade_results.txt", header=TRUE)
str(data)
data$Taxon_ID <- sub("^<+", "", data$Taxon_ID)
data$Taxon_ID <- gsub("<.*", "", data$Taxon_ID)
data$Taxon_ID <- sub(">", "", data$Taxon_ID)
str(data)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/04_longevity/")
FvsB <- read.csv("F_vs_B_longevity.csv", header=TRUE)
str(FvsB)
data <- merge(data,FvsB,by="Taxon_ID")
str(data)
data$DvsI<-data$Decrease/data$Increase
str(data)
t.test(Increase~F.node, var.equal = FALSE, data=data)
t.test(Decrease~F.node, var.equal = FALSE, data=data)
t.test(DvsI~F.node, var.equal = FALSE, data=data) 
data <- data[order(data$F.node, decreasing = TRUE),]
#write.csv(data, "summary_tree_longevity_leaf_nodes_overview.csv", row.names=F)
# extract mean
data.sum.long.leaf <- data %>% group_by(F.node) %>% 
  summarise_at(vars("Increase", "Decrease", "DvsI"), mean)

### TABLE FOR ALL 
## summarize all mean values in one table
summary.all.traits.leaf.nodes <- bind_rows(data.sum.cave.leaf, data.sum.pig.leaf, data.sum.eye.leaf, data.sum.long.leaf)

#### EXPANSION VS CONTRACTION DIFFERENCES IN FAST EVOLVING FAMILIES #### 
### combine results and put into table averages for all values and analyses + significance

### CAVE SPECIES: FAST EVOLVING GENE FAMILIES ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/01_cave_species/03_filtered_50_lambda3")
# Probability of rapid change
cavspec.prob.L3 <- read.table("Base_branch_probabilities.tab", header=TRUE)
names(cavspec.prob.L3) <- gsub("X\\.|\\..*", '', names(cavspec.prob.L3))
cavespec.prob.L3.bin <- cavspec.prob.L3 %>% mutate(across(where(is.numeric), ~if_else(. > 0.05, "NO", "YES")))
# Gene family count file (positive = exp, negative = cont)
cavspec.count.L3 <- read.table("Base_change.tab", header=TRUE)
names(cavspec.count.L3) <- gsub("X\\.|\\..*", '', names(cavspec.count.L3))
cavspec.count.L3[,-1] <- as.data.frame(apply(cavspec.count.L3[,-1], c(1,2), function(x) ifelse(x>0,"EXP", ifelse(x<0,"CON","0"))))
cavespec.count.L3.bin <- cavspec.count.L3
# Combine fast evolving and exp/con matrices
fast.evol <- cavespec.prob.L3.bin[order(cavespec.prob.L3.bin$FamilyID),]
expa.contr <- cavespec.count.L3.bin[order(cavespec.count.L3.bin$FamilyID),]
OGs <- data.frame(FamilyID = expa.contr[,1])
new.tab <- merge(fast.evol, OGs, by="FamilyID", all.y=TRUE)
new.tab[is.na(new.tab)] <- "NO"
comb <- new.tab
for (i in 1:nrow(new.tab)) {for (j in 1:ncol(new.tab)) {if (new.tab[i,j] == "YES") comb[i,j] <- expa.contr[i,j]}}
# Prepare expanded table
expanded.table <- comb
expanded.table[expanded.table == "CON"] <- "NO"
exp.table.leaves <- expanded.table[,1:23]
# Prepare contracted table
comb <- new.tab
for (i in 1:nrow(new.tab)) {for (j in 1:ncol(new.tab)) {if (new.tab[i,j] == "YES") comb[i,j] <- expa.contr[i,j]}}
comb[comb == "EXP"] <- "NO"
contracted.table <- comb
contr.table.leaves <- contracted.table[,1:23]

### FOR EXPANDED OGs ####
data.trans <- pivot_longer(exp.table.leaves, cols = 2:ncol(exp.table.leaves), names_to = "species", values_to = "expan")
str(data.trans)
dim(exp.table.leaves)
dim(data.trans)
new.table.exp <- data.trans %>%
  mutate(type  = if_else(species %in% c("Luden", "Asmec", "Trros", "Sirhi","Sians","Prang"), "cave", "surface"))
new.table.exp
str(new.table.exp)
new.table.exp.chi <- new.table.exp[-c(2)]
# convert new.table.con.chi to data.table
final.table.exp.chi <- as.data.frame(table(new.table.exp.chi))
str(final.table.exp.chi)
### across all families
overall <- new.table.exp[,3:4]
str(overall)
overall <- overall[, c("type", "expan")]
new.table.exp = table(overall)
print(new.table.exp)
str(new.table.exp)
chisq.test(new.table.exp)
prop.table(table(overall$expan, overall$type), margin=2)*100
### create a Figure for this
# Make a stacked barplot--> it will be in %!
# use percentages and base R barplot
my.table <- prop.table(table(overall$expan, overall$type), margin=2)*100
col <- c("#F5B452", "#D3D3D3")
cave.surf <- barplot(my.table, col=col , border="white", legend = TRUE, xlab="group")
str(my.table)
# Data Visualization Of Contingency Table With ggplot2 (Stacked Bar Graph):
df <- as.data.frame(t(my.table))
## panel d
p1 <- ggplot(data = df, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat = "identity") + 
  labs(x = "\n Answer", y = "Percentage \n", 
       title = "Gene family expansion \n",
       fill = "Gene family evolution") +
  scale_fill_manual(values=c("#F5B452", "#D3D3D3")) +
  theme_classic()
p1
### for each OG individually
expanded.table.2 <- final.table.exp.chi %>%
  group_by(FamilyID) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(type, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$expan
    return(M)}))
expanded.table.3 <- expanded.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
expanded.table.3
exp.ratios <- pivot_wider(final.table.exp.chi, names_from = type, values_from = Freq)
str(exp.ratios)
# calculate EXP/NO ratio in cave species
cave <- exp.ratios %>% 
  group_by(FamilyID, expan) %>% 
  summarise(cave = sum(cave)) %>% 
  mutate(cave.ratio = if_else(expan == "EXP", cave/lead(cave), 1))
cave <- cave[cave$expan!= "NO",]
cave <- cave[,-c(2:3)]
# calculate EXP/NO ratio in surface species
surface <- exp.ratios %>% 
  group_by(FamilyID, expan) %>% 
  summarise(surface = sum(surface)) %>% 
  mutate(surface.ratio = if_else(expan == "EXP", surface/lead(surface), 1))
surface <- surface[surface$expan!= "NO",]
surface <- surface[,-c(2:3)]
expanded.OGs.significance.ratios <- merge(expanded.table.3, cave, by="FamilyID")
expanded.OGs.significance.ratios.2 <- merge(expanded.OGs.significance.ratios, surface, by="FamilyID")
str(expanded.OGs.significance.ratios.2)
#write.table(expanded.OGs.significance.ratios.2,"fisher_expanded_OGs_ratios_cave_surface.csv", sep=",", quote=FALSE, row.names=FALSE)
### check if sharing in cave is higher than in surface species
## load data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/07_main_results")
shared.OGs.expanded.cave.species <- read.csv("fisher_expanded_OGs_ratios_cave_surface.csv")
# Convert significance to logical
shared.OGs.expanded.cave.species$significant <- shared.OGs.expanded.cave.species$significant == "yes"
# Create column to indicate direction (TRUE = more in cave, FALSE = more in surface)
shared.OGs.expanded.cave.species$cave_more <- shared.OGs.expanded.cave.species$DIFF > 0
# Contingency table: cave vs. surface by significance
cont_table <- table(DirectionIsCave = shared.OGs.expanded.cave.species$cave_more,
                    IsSignificant = shared.OGs.expanded.cave.species$significant)
print(cont_table)
## statistical test
chisq.test(cont_table)
## cave species have significantly more shared expanded genes than surface species: X-squared = 674.45, df = 1, p-value < 2.2e-1

#### FOR CONTRACTED OGs ####
data.trans = pivot_longer(contr.table.leaves, cols = 2:ncol(contr.table.leaves), names_to = "species", values_to = "contr")
str(data.trans)
new.table.con <- data.trans %>%
  mutate(type  = if_else(species %in% c("Luden", "Asmec", "Trros", "Sirhi","Sians","Prang"), "cave", "surface"))
new.table.con
new.table.con.chi <- new.table.con[-c(2)]
# convert new.table.con.chi to data.table
final.table.con.chi <- as.data.frame(table(new.table.con.chi))
str(final.table.con.chi)
### across all families
overall <- new.table.con[,3:4]
str(overall)
overall <- overall[, c("type", "contr")]
new.table.con = table(overall)
print(new.table.con)
str(new.table.con)
chisq.test(new.table.con)
prop.table(table(overall$contr, overall$type), margin=2)*100
### make the Figure
# use percentages and base R barplot
my.table <- prop.table(table(overall$contr, overall$type), margin=2)*100
cave.surf <- barplot(my.table, col=col , border="white", legend = TRUE, xlab="group")
str(my.table)
# Data Visualization Of Contingency Table With ggplot2 (Stacked Bar Graph):
df <- as.data.frame(t(my.table))
## panel e
p2 <- ggplot(data = df, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat = "identity") + 
  labs(x = "\n Answer", y = "Percentage \n", title = "Gene family contraction \n", fill = "Gene family evolution") +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  theme_classic()

#### FIGURE 2 d,e ####
### Put contraction and expansion together
grid.arrange(p1,p2,ncol=2)

#### CHECK CONVERGENCE OF INDIVIDUAL OGS #### 
contracted.table.2 <- final.table.con.chi %>%
  group_by(FamilyID) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(type, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$contr
    return(M)}))
contracted.table.3 <- contracted.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
contracted.table.3
contr.ratios <- pivot_wider(final.table.con.chi, names_from = type,
                            values_from = Freq)
str(contr.ratios)
# calculate CON/NO ratio in cave species
cave <- contr.ratios %>% 
  group_by(FamilyID, contr) %>% 
  summarise(cave = sum(cave)) %>% 
  mutate(cave.ratio = if_else(contr == "CON", cave/lead(cave), 1))
cave <- cave[cave$contr!= "NO",]
cave <- cave[,-c(2:3)]
# calculate CON/NO ratio in surface species
surface <- contr.ratios %>% 
  group_by(FamilyID, contr) %>% 
  summarise(surface = sum(surface)) %>% 
  mutate(surface.ratio = if_else(contr == "CON", surface/lead(surface), 1))
surface <- surface[surface$contr!= "NO",]
surface <- surface[,-c(2:3)]
contracted.OGs.significance.ratios <- merge(contracted.table.3, cave, by="FamilyID")
contracted.OGs.significance.ratios.2 <- merge(contracted.OGs.significance.ratios, surface, by="FamilyID")
str(contracted.OGs.significance.ratios.2)
#write.table(contracted.OGs.significance.ratios.2,"fisher_contracted_OGs_ratios_cave_surface.csv", sep=",", quote=FALSE, row.names=FALSE)
## load data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/07_main_results")
shared.OGs.contracted.cave.species <- read.csv("fisher_contracted_OGs_ratios_cave_surface.csv")
# Convert significance to logical
shared.OGs.contracted.cave.species$significant <- shared.OGs.contracted.cave.species$significant == "yes"
# Create column to indicate direction (TRUE = more in cave, FALSE = more in surface)
shared.OGs.contracted.cave.species$cave_more <- shared.OGs.contracted.cave.species$DIFF > 0
# Contingency table: cave vs. surface by significance
cont_table <- table(DirectionIsCave = shared.OGs.contracted.cave.species$cave_more,
                    IsSignificant = shared.OGs.contracted.cave.species$significant)
print(cont_table)
## statistical test
chisq.test(cont_table)
## cave species have significantly more shared contracted genes than surface species: X-squared = 344.09, df = 1, p-value < 2.2e-16

### make a Figure for OGs that are shared 
# Create the data
df <- data.frame(
  environment = c("cave", "cave", "surface", "surface"),
  genefam = c("expanded", "contracted", "expanded", "contracted"),
  values = c(56, 35, 0, 0))
# Set the desired order of environments: surface first, then cave
df$environment <- factor(df$environment, levels = c("surface", "cave"))

#### FIGURE 2g ####
# Make the grouped barplot
ggplot(df, aes(x = environment, y = values, fill = genefam)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("expanded" = "#F5B452", "contracted" = "#56B4E9")) +
  labs(title = "Gene family expansions and contractions by environment", x = "Environment", y = "Number of gene families", fill = "Gene family type") +
  theme_classic() +
  geom_text(aes(label = values), position = position_dodge(width = 0.7), vjust = -0.3, size = 3)

### PIGMENTATION LOSS ### 
## pigmentation:378324, 376598, 370388;  lambda 3
### Probability of rapid change (expansion OR contraction)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/02_pigmentation_loss/04_filtered_50_lambda3")
pigm.prob.L3 <- read.table("Base_branch_probabilities.tab", header=TRUE)
str(pigm.prob.L3)
## remove unwanted  characters
names(pigm.prob.L3) <- gsub("X\\.",'', names(pigm.prob.L3))
str(pigm.prob.L3)
names(pigm.prob.L3) <- gsub("\\..*","", names(pigm.prob.L3))
str(pigm.prob.L3)
pigm.prob.L3
pigm.prob.L3.bin <- pigm.prob.L3 %>% 
  mutate(across(where(is.numeric), function(x) if_else(x > 0.05, "NO", "YES")))
### Gene family count file (positive = exp, negative = cont)
pigm.count.L3 <- read.table("Base_change.tab", header=TRUE)
str(pigm.count.L3)
## remove unwanted  characters
names(pigm.count.L3) <- gsub("X\\.",'', names(pigm.count.L3))
str(pigm.count.L3)
names(pigm.count.L3) <- gsub("\\..*","", names(pigm.count.L3))
str(pigm.count.L3)
pigm.count.L3[,-1] = apply(pigm.count.L3[,-1], c(1,2), function(x) {ifelse(any(x > 0), "EXP", ifelse(any(x < 0), "CON", "0"))})
pigm.count.L3[,-1]  = as.data.frame(pigm.count.L3[,-1])
pigm.count.L3.bin <- pigm.count.L3
#### check which families have significantly evolved faster and expanded or contracted 
fast.evol <- pigm.prob.L3.bin
expa.contr <- pigm.count.L3.bin
# order by FamilyID
fast.evol <- fast.evol[order(fast.evol$FamilyID),]
expa.contr <- expa.contr[order(expa.contr$FamilyID),]
OGs<- data.frame(FamilyID= expa.contr[,1])
str(expa.contr)
str(fast.evol)
new.tab <- merge(fast.evol,OGs, by="FamilyID", all.y=TRUE)
new.tab[order(new.tab$FamilyID),]
new.tab[is.na(new.tab)] <- "NO"
## combine matrices and take value in expa.contr on condition that fast.evol has entry YES
comb = new.tab
for (i in 1:nrow(new.tab)) {
  for (j in 1:ncol(new.tab)) {
    if (new.tab[i,j] == "YES") {
      comb[i,j] = expa.contr[i,j]}}}
comb
## write file with all expanded gene families
comb[comb == "CON"] <- "NO"
expanded.table <- comb
str(expanded.table)
# only extract columns with leaf nodes 
exp.table.leaves <- expanded.table[,1:26]
str(exp.table.leaves)
#write.table(comb, file="expanded_fast_evolving_OGs.csv", sep=",", row.names=F)
## write file with all contracted gene families
comb = new.tab
for (i in 1:nrow(new.tab)) {
  for (j in 1:ncol(new.tab)) {
    if (new.tab[i,j] == "YES") {
      comb[i,j] = expa.contr[i,j]}}}
comb[comb == "EXP"] <- "NO"
contracted.table <- comb
write.table(comb, file="contracted_fast_evolving_OGs.csv", sep=",", row.names=F)
str(contracted.table)
# only extract columns with leaf nodes 
contr.table.leaves <- contracted.table[,1:26]
str(contr.table.leaves)

### FOR EXPANDED OGs ####
data.trans <- pivot_longer(exp.table.leaves, cols = 2:ncol(exp.table.leaves), names_to = "species", values_to = "expan")
str(data.trans)
dim(exp.table.leaves)
dim(data.trans)
new.table.exp <- data.trans %>%
  mutate(type  = if_else(species %in% c("Hegla", "Luden", "Asmec", "Trros", "Sirhi", "Ammex", "Sians", "Prang"), "pigm_loss", "pigm"))
new.table.exp
str(new.table.exp)
new.table.exp.chi <- new.table.exp[-c(2)]
# convert new.table.con.chi to data.table
final.table.exp.chi <- as.data.frame(table(new.table.exp.chi))
str(final.table.exp.chi)
### across all families
overall <- new.table.exp[,3:4]
str(overall)
overall <- overall[, c("type", "expan")]
new.table.exp = table(overall)
print(new.table.exp)
str(new.table.exp)
chisq.test(new.table.exp)
prop.table(table(overall$expan, overall$type), margin=2)*100
### for each OG individually
expanded.table.2 <- final.table.exp.chi %>%
  group_by(FamilyID) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(type, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$expan
    return(M)}))
expanded.table.3 <- expanded.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
expanded.table.3
exp.ratios <- pivot_wider(final.table.exp.chi, names_from = type,
                            values_from = Freq)
str(exp.ratios)
# calculate EXP/NO ratio in pigmentation loss species
pigm_loss <- exp.ratios %>% 
  group_by(FamilyID, expan) %>% 
  summarise(pigm_loss = sum(pigm_loss)) %>% 
  mutate(pigm_loss.ratio = if_else(expan == "EXP", pigm_loss/lead(pigm_loss), 1))
pigm_loss <- pigm_loss[pigm_loss$expan!= "NO",]
pigm_loss <- pigm_loss[,-c(2:3)]
# calculate EXP/NO ratio in pigmented species
pigm <- exp.ratios %>% 
  group_by(FamilyID, expan) %>% 
  summarise(pigm = sum(pigm)) %>% 
  mutate(pigm.ratio = if_else(expan == "EXP", pigm/lead(pigm), 1))
pigm <- pigm[pigm$expan!= "NO",]
pigm <- pigm[,-c(2:3)]
expanded.OGs.significance.ratios <- merge(expanded.table.3, pigm_loss, by="FamilyID")
expanded.OGs.significance.ratios.2 <- merge(expanded.OGs.significance.ratios, pigm, by="FamilyID")
str(expanded.OGs.significance.ratios.2)
#write.table(expanded.OGs.significance.ratios.2,"fisher_expanded_OGs_ratios_pigmentation_loss.csv", sep=",", quote=FALSE, row.names=FALSE)

#### FOR CONTRACTED OGs ####
data.trans = pivot_longer(contr.table.leaves, cols = 2:ncol(contr.table.leaves), names_to = "species", values_to = "contr")
str(data.trans)
new.table.con <- data.trans %>%
  mutate(type  = if_else(species %in% c("Hegla", "Luden", "Asmec", "Trros", "Sirhi", "Ammex", "Sians", "Prang"), "pigm_loss", "pigm"))
new.table.con
new.table.con.chi <- new.table.con[-c(2)]
# convert new.table.con.chi to data.table
final.table.con.chi <- as.data.frame(table(new.table.con.chi))
str(final.table.con.chi)
### across all families
overall <- new.table.con[,3:4]
str(overall)
overall <- overall[, c("type", "contr")]
new.table.con = table(overall)
print(new.table.con)
str(new.table.con)
chisq.test(new.table.con)
prop.table(table(overall$contr, overall$type), margin=2)*100
### both expanded and contracted OGs occur more frequently in cave species
### For each OG individually
contracted.table.2 <- final.table.con.chi %>%
  group_by(FamilyID) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(type, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$contr
    return(M)}))
contracted.table.3 <- contracted.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
contracted.table.3
contr.ratios <- pivot_wider(final.table.con.chi, names_from = type,
                            values_from = Freq)
str(contr.ratios)
# calculate CON/NO ratio in pigmentation loss species
pigm_loss <- contr.ratios %>% 
  group_by(FamilyID, contr) %>% 
  summarise(pigm_loss = sum(pigm_loss)) %>% 
  mutate(pigm_loss.ratio = if_else(contr == "CON", pigm_loss/lead(pigm_loss), 1))
pigm_loss <- pigm_loss[pigm_loss$contr!= "NO",]
pigm_loss <- pigm_loss[,-c(2:3)]
# calculate CON/NO ratio in pigmented species
pigm <- contr.ratios %>% 
  group_by(FamilyID, contr) %>% 
  summarise(pigm = sum(pigm)) %>% 
  mutate(pigm.ratio = if_else(contr == "CON", pigm/lead(pigm), 1))
pigm <- pigm[pigm$contr!= "NO",]
pigm <- pigm[,-c(2:3)]
contracted.OGs.significance.ratios <- merge(contracted.table.3, pigm_loss, by="FamilyID")
contracted.OGs.significance.ratios.2 <- merge(contracted.OGs.significance.ratios, pigm, by="FamilyID")
str(contracted.OGs.significance.ratios.2)
#write.table(contracted.OGs.significance.ratios.2,"fisher_contracted_OGs_ratios_pigmentation_loss_check.csv", sep=",", quote=FALSE, row.names=FALSE)

### EYE LOSS ###
## eye loss: 407743, 404792, 403945, 398884; lambda 4
### Probability of rapid change (expansion OR contraction)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/03_eye_loss/05_filtered_50_lambda4")
eyeloss.prob.L3 <- read.table("Base_branch_probabilities.tab", header=TRUE)
str(eyeloss.prob.L3)
## remove unwanted  characters
names(eyeloss.prob.L3) <- gsub("X\\.",'', names(eyeloss.prob.L3))
str(eyeloss.prob.L3)
names(eyeloss.prob.L3) <- gsub("\\..*","", names(eyeloss.prob.L3))
str(eyeloss.prob.L3)
eyeloss.prob.L3
eyeloss.prob.L3.bin <- eyeloss.prob.L3 %>% 
  mutate(across(where(is.numeric), function(x) if_else(x > 0.05, "NO", "YES")))
### Gene family count file (positive = exp, negative = cont)
eyeloss.count.L3 <- read.table("Base_change.tab", header=TRUE)
str(eyeloss.count.L3)
## remove unwanted  characters
names(eyeloss.count.L3) <- gsub("X\\.",'', names(eyeloss.count.L3))
str(eyeloss.count.L3)
names(eyeloss.count.L3) <- gsub("\\..*","", names(eyeloss.count.L3))
str(eyeloss.count.L3)
eyeloss.count.L3[,-1] = apply(eyeloss.count.L3[,-1], c(1,2), function(x) {ifelse(any(x > 0), "EXP", ifelse(any(x < 0), "CON", "0"))})
eyeloss.count.L3[,-1]  = as.data.frame(eyeloss.count.L3[,-1])
eyeloss.count.L3.bin <- eyeloss.count.L3
#### check which families have significantly evolved faster and expanded or contracted 
fast.evol <- eyeloss.prob.L3.bin
expa.contr <- eyeloss.count.L3.bin
# order by FamilyID
fast.evol <- fast.evol[order(fast.evol$FamilyID),]
expa.contr <- expa.contr[order(expa.contr$FamilyID),]
OGs<- data.frame(FamilyID= expa.contr[,1])
str(expa.contr)
str(fast.evol)
new.tab <- merge(fast.evol,OGs, by="FamilyID", all.y=TRUE)
new.tab[order(new.tab$FamilyID),]
new.tab[is.na(new.tab)] <- "NO"
## combine matrices and take value in expa.contr on condition that fast.evol has entry YES
comb = new.tab
for (i in 1:nrow(new.tab)) {
  for (j in 1:ncol(new.tab)) {
    if (new.tab[i,j] == "YES") {
      comb[i,j] = expa.contr[i,j]}}}
comb
## write file with all expanded gene families
comb[comb == "CON"] <- "NO"
expanded.table <- comb
str(expanded.table)
# only extract columns with leaf nodes 
exp.table.leaves <- expanded.table[,1:33]
str(exp.table.leaves)
#write.table(comb, file="expanded_fast_evolving_OGs.csv", sep=",", row.names=F)
## write file with all contracted gene families
comb = new.tab
for (i in 1:nrow(new.tab)) {
  for (j in 1:ncol(new.tab)) {
    if (new.tab[i,j] == "YES") {
      comb[i,j] = expa.contr[i,j]}}}
comb[comb == "EXP"] <- "NO"
contracted.table <- comb
#write.table(comb, file="contracted_fast_evolving_OGs.csv", sep=",", row.names=F)
str(contracted.table)
# only extract columns with leaf nodes 
contr.table.leaves <- contracted.table[,1:33]
str(contr.table.leaves)

### FOR EXPANDED OGs ####
data.trans <- pivot_longer(exp.table.leaves, cols = 2:ncol(exp.table.leaves), names_to = "species", values_to = "expan")
str(data.trans)
dim(contr.table.leaves)
dim(data.trans)
new.table.exp <- data.trans %>%
  mutate(type  = if_else(species %in% c("Luden",  "Asmec",  "Sians",  "Miuni",  "Geser",  "Prang",  "Cocri",  "Nagal",  "Hegla",  "Chasi",  "Trros"), "eye_loss", "eye"))
new.table.exp
str(new.table.exp)
new.table.exp.chi <- new.table.exp[-c(2)]
# convert new.table.con.chi to data.table
final.table.exp.chi <- as.data.frame(table(new.table.exp.chi))
str(final.table.exp.chi)
#write.table(final.table.exp.chi,"contingency_tests_table_expanded.csv", sep=",")
### across all families
overall <- new.table.exp[,3:4]
str(overall)
overall <- overall[, c("type", "expan")]
new.table.exp = table(overall)
print(new.table.exp)
str(new.table.exp)
chisq.test(new.table.exp)
prop.table(table(overall$expan, overall$type), margin=2)*100
### for each OG individually
expanded.table.2 <- final.table.exp.chi %>%
  group_by(FamilyID) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(type, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$expan
    return(M)}))
expanded.table.3 <- expanded.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
expanded.table.3
exp.ratios <- pivot_wider(final.table.exp.chi, names_from = type,
                            values_from = Freq)
str(exp.ratios)
# calculate EXP/NO ratio in eye loss species
eye_loss <- exp.ratios %>% 
  group_by(FamilyID, expan) %>% 
  summarise(eye_loss = sum(eye_loss)) %>% 
  mutate(eye_loss.ratio = if_else(expan == "EXP", eye_loss/lead(eye_loss), 1))
eye_loss <- eye_loss[eye_loss$expan!= "NO",]
eye_loss <- eye_loss[,-c(2:3)]
# calculate EXP/NO ratio in species with eyes
eye <- exp.ratios %>% 
  group_by(FamilyID, expan) %>% 
  summarise(eye = sum(eye)) %>% 
  mutate(eye.ratio = if_else(expan == "EXP", eye/lead(eye), 1))
eye <- eye[eye$expan!= "NO",]
eye <- eye[,-c(2:3)]
expanded.OGs.significance.ratios <- merge(expanded.table.3, eye_loss, by="FamilyID")
expanded.OGs.significance.ratios.2 <- merge(expanded.OGs.significance.ratios, eye, by="FamilyID")
str(expanded.OGs.significance.ratios.2)
#write.table(expanded.OGs.significance.ratios.2,"fisher_expanded_OGs_ratios_eye_loss.csv", sep=",", quote=FALSE, row.names=FALSE)

#### FOR CONTRACTED OGs ####
data.trans = pivot_longer(contr.table.leaves, cols = 2:ncol(contr.table.leaves), names_to = "species", values_to = "contr")
str(data.trans)
new.table.con <- data.trans %>%
  mutate(type  = if_else(species %in% c("Luden", "Asmec", "Sians", "Miuni", "Geser", "Prang", "Cocri", "Nagal", "Hegla", "Chasi", "Trros"), "eye_loss", "eye"))
new.table.con
new.table.con.chi <- new.table.con[-c(2)]
# convert new.table.con.chi to data.table
final.table.con.chi <- as.data.frame(table(new.table.con.chi))
str(final.table.con.chi)
### across all families
overall <- new.table.con[,3:4]
str(overall)
overall <- overall[, c("type", "contr")]
new.table.con = table(overall)
print(new.table.con)
str(new.table.con)
chisq.test(new.table.con)
prop.table(table(overall$contr, overall$type), margin=2)*100
### both expanded and contracted OGs occur more frequently in cave species
### For each OG individually
contracted.table.2 <- final.table.con.chi %>%
  group_by(FamilyID) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(type, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$contr
    return(M)}))
contracted.table.3 <- contracted.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
contracted.table.3
contr.ratios <- pivot_wider(final.table.con.chi, names_from = type,
                            values_from = Freq)
str(contr.ratios)
# calculate CON/NO ratio in eye loss species
eye_loss <- contr.ratios %>% 
  group_by(FamilyID, contr) %>% 
  summarise(eye_loss = sum(eye_loss)) %>% 
  mutate(eye_loss.ratio = if_else(contr == "CON", eye_loss/lead(eye_loss), 1))
eye_loss <- eye_loss[eye_loss$contr!= "NO",]
eye_loss <- eye_loss[,-c(2:3)]
# calculate CON/NO ratio in species with eyes
eye <- contr.ratios %>% 
  group_by(FamilyID, contr) %>% 
  summarise(eye = sum(eye)) %>% 
  mutate(eye.ratio = if_else(contr == "CON", eye/lead(eye), 1))
eye <- eye[eye$contr!= "NO",]
eye <- eye[,-c(2:3)]
contracted.OGs.significance.ratios <- merge(contracted.table.3, eye_loss, by="FamilyID")
contracted.OGs.significance.ratios.2 <- merge(contracted.OGs.significance.ratios, eye, by="FamilyID")
str(contracted.OGs.significance.ratios.2)
#write.table(contracted.OGs.significance.ratios.2,"fisher_contracted_OGs_ratios_eye_loss.csv", sep=",", quote=FALSE, row.names=FALSE)

### LONGEVITY
## longevity: 255284, 254601, 251561, 253678; lambda 3
### Probability of rapid change (expansion OR contraction)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/04_longevity/04_filtered_50_lambda3")
longev.prob.L3 <- read.table("Base_branch_probabilities.tab", header=TRUE)
str(longev.prob.L3)
## remove unwanted  characters
names(longev.prob.L3) <- gsub("X\\.",'', names(longev.prob.L3))
str(longev.prob.L3)
names(longev.prob.L3) <- gsub("\\..*","", names(longev.prob.L3))
str(longev.prob.L3)
longev.prob.L3
longev.prob.L3.bin <- longev.prob.L3 %>% 
  mutate(across(where(is.numeric), function(x) if_else(x > 0.05, "NO", "YES")))
### Gene family count file (positive = exp, negative = cont)
longev.count.L3 <- read.table("Base_change.tab", header=TRUE)
str(longev.count.L3)
## remove unwanted  characters
names(longev.count.L3) <- gsub("X\\.",'', names(longev.count.L3))
str(longev.count.L3)
names(longev.count.L3) <- gsub("\\..*","", names(longev.count.L3))
str(longev.count.L3)
longev.count.L3[,-1] = apply(longev.count.L3[,-1], c(1,2), function(x) {ifelse(any(x > 0), "EXP", ifelse(any(x < 0), "CON", "0"))})
longev.count.L3[,-1]  = as.data.frame(longev.count.L3[,-1])
longev.count.L3.bin <- longev.count.L3
#### check which families have significantly evolved faster and expanded or contracted 
fast.evol <- longev.prob.L3.bin
expa.contr <- longev.count.L3.bin
# order by FamilyID
fast.evol <- fast.evol[order(fast.evol$FamilyID),]
expa.contr <- expa.contr[order(expa.contr$FamilyID),]
OGs<- data.frame(FamilyID= expa.contr[,1])
str(expa.contr)
str(fast.evol)
new.tab <- merge(fast.evol,OGs, by="FamilyID", all.y=TRUE)
new.tab[order(new.tab$FamilyID),]
new.tab[is.na(new.tab)] <- "NO"
## combine matrices and take value in expa.contr on condition that fast.evol has entry YES
comb = new.tab
for (i in 1:nrow(new.tab)) {
  for (j in 1:ncol(new.tab)) {
    if (new.tab[i,j] == "YES") {
      comb[i,j] = expa.contr[i,j]}}}
comb
## write file with all expanded gene families
comb[comb == "CON"] <- "NO"
expanded.table <- comb
str(expanded.table)
# only extract columns with leaf nodes 
exp.table.leaves <- expanded.table[,1:26]
str(exp.table.leaves)
#write.table(comb, file="expanded_fast_evolving_OGs.csv", sep=",", row.names=F)
## write file with all contracted gene families
comb = new.tab
for (i in 1:nrow(new.tab)) {
  for (j in 1:ncol(new.tab)) {
    if (new.tab[i,j] == "YES") {
      comb[i,j] = expa.contr[i,j]}}}
comb[comb == "EXP"] <- "NO"
contracted.table <- comb
#write.table(comb, file="contracted_fast_evolving_OGs.csv", sep=",", row.names=F)
str(contracted.table)
# only extract columns with leaf nodes 
contr.table.leaves <- contracted.table[,1:26]
str(contr.table.leaves)

### FOR EXPANDED OGs ####
data.trans <- pivot_longer(exp.table.leaves, cols = 2:ncol(exp.table.leaves), names_to = "species", values_to = "expan")
str(data.trans)
dim(contr.table.leaves)
dim(data.trans)
new.table.exp <- data.trans %>%
  mutate(type  = if_else(species %in% c("Rhtyp", "Cacar", "Acbae", "Prang", "Ororc", "Myluc", "Hegla", "Hosap", "Elmax", "Tecar", "Lacha"), "long.life", "short.life"))
new.table.exp
str(new.table.exp)
new.table.exp.chi <- new.table.exp[-c(2)]
# convert new.table.con.chi to data.table
final.table.exp.chi <- as.data.frame(table(new.table.exp.chi))
str(final.table.exp.chi)
### across all families
overall <- new.table.exp[,3:4]
str(overall)
overall <- overall[, c("type", "expan")]
new.table.exp = table(overall)
print(new.table.exp)
str(new.table.exp)
chisq.test(new.table.exp)
my.table <- prop.table(table(overall$expan, overall$type), margin=2)*100
# Data Visualization Of Contingency Table With ggplot2 (Stacked Bar Graph):
df <- as.data.frame(t(my.table))

#### FIGURE 4b ####
p2 <- ggplot(data = df, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat = "identity") + 
  labs(x = "\n Answer", y = "Percentage \n", 
       title = "Gene family contraction \n",
       fill = "Gene family evolution") +
  scale_fill_manual(values=c("#F5B452", "#D3D3D3")) +
  theme_classic()
p2

### for each OG individually
expanded.table.2 <- final.table.exp.chi %>%
  group_by(FamilyID) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(type, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$expan
    return(M)}))
expanded.table.3 <- expanded.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
expanded.table.3
exp.ratios <- pivot_wider(final.table.exp.chi, names_from = type,
                            values_from = Freq)
str(exp.ratios)
# calculate EXP/NO ratio in long-lived species
long.life <- exp.ratios %>% 
  group_by(FamilyID, expan) %>% 
  summarise(long.life = sum(long.life)) %>% 
  mutate(long.life.ratio = if_else(expan == "EXP", long.life/lead(long.life), 1))
long.life <- long.life[long.life$expan!= "NO",]
long.life <- long.life[,-c(2:3)]
# calculate EXP/NO ratio in short-lived species
short.life <- exp.ratios %>% 
  group_by(FamilyID, expan) %>% 
  summarise(short.life = sum(short.life)) %>% 
  mutate(short.life.ratio = if_else(expan == "EXP", short.life/lead(short.life), 1))
short.life <- short.life[short.life$expan!= "NO",]
short.life <- short.life[,-c(2:3)]
expanded.OGs.significance.ratios <- merge(expanded.table.3, long.life, by="FamilyID")
expanded.OGs.significance.ratios.2 <- merge(expanded.OGs.significance.ratios, short.life, by="FamilyID")
str(expanded.OGs.significance.ratios.2)
#write.table(expanded.OGs.significance.ratios.2,"fisher_expanded_OGs_ratios_longevity.csv", sep=",", quote=FALSE, row.names=FALSE)

#### FOR CONTRACTED OGs ####
data.trans = pivot_longer(contr.table.leaves, cols = 2:ncol(contr.table.leaves), names_to = "species", values_to = "contr")
str(data.trans)
new.table.con <- data.trans %>%
  mutate(type  = if_else(species %in% c("Rhtyp", "Cacar", "Acbae", "Prang", "Ororc", "Myluc", "Hegla", "Hosap", "Elmax", "Tecar", "Lacha"), "long.life", "short.life"))
new.table.con
new.table.con.chi <- new.table.con[-c(2)]
# convert new.table.con.chi to data.table
final.table.con.chi <- as.data.frame(table(new.table.con.chi))
str(final.table.con.chi)
### across all families
overall <- new.table.con[,3:4]
str(overall)
overall <- overall[, c("type", "contr")]
new.table.con = table(overall)
print(new.table.con)
str(new.table.con)
chisq.test(new.table.con)
prop.table(table(overall$contr, overall$type), margin=2)*100
my.table <- prop.table(table(overall$contr, overall$type), margin=2)*100
# Data Visualization Of Contingency Table With ggplot2 (Stacked Bar Graph):
df <- as.data.frame(t(my.table))

#### FIGURE 4c ####
p3 <- ggplot(data = df, aes(x = reorder(Var1, Var2), y = Freq, fill = Var2)) + 
  geom_bar(stat = "identity") + 
  labs(x = "\n Answer", y = "Percentage \n", 
       title = "Gene family expansion \n",
       fill = "Gene family evolution") +
  scale_fill_manual(values=c("#56B4E9", "#D3D3D3")) +
  theme_classic()
p3
### both expanded and contracted OGs occur more frequently in cave species

### For each OG individually
contracted.table.2 <- final.table.con.chi %>%
  group_by(FamilyID) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(type, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$contr
    return(M)}))
contracted.table.3 <- contracted.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
contracted.table.3
contr.ratios <- pivot_wider(final.table.con.chi, names_from = type,
                            values_from = Freq)
str(contr.ratios)
# calculate CON/NO ratio in long-lived species
long.life <- contr.ratios %>% 
  group_by(FamilyID, contr) %>% 
  summarise(long.life = sum(long.life)) %>% 
  mutate(long.life.ratio = if_else(contr == "CON", long.life/lead(long.life), 1))
long.life <- long.life[long.life$contr!= "NO",]
long.life <- long.life[,-c(2:3)]
# calculate CON/NO ratio in short-lived species
short.life <- contr.ratios %>% 
  group_by(FamilyID, contr) %>% 
  summarise(short.life = sum(short.life)) %>% 
  mutate(short.life.ratio = if_else(contr == "CON", short.life/lead(short.life), 1))
short.life <- short.life[short.life$contr!= "NO",]
short.life <- short.life[,-c(2:3)]
contracted.OGs.significance.ratios <- merge(contracted.table.3, long.life, by="FamilyID")
contracted.OGs.significance.ratios.2 <- merge(contracted.OGs.significance.ratios, short.life, by="FamilyID")
str(contracted.OGs.significance.ratios.2)
#write.table(contracted.OGs.significance.ratios.2,"fisher_contracted_OGs_ratios_longevity.csv", sep=",", quote=FALSE, row.names=FALSE)


#### COMPARE IF CAVE SPECIES HAVE RELATIVELY MORE EXPANDED / CONTRCTED OGs THAN EYE LOSS AND PIGMENTATION LOSS SPECIES ####
analysis_table <- data.frame(
  trait = c("cave_vs_surface","cave_vs_surface",
            "eye_loss_vs_eyes","eye_loss_vs_eyes",
            "pigment_loss_vs_pigment","pigment_loss_vs_pigment"),
  status = c("expanded","contracted",
             "expanded","contracted",
             "expanded","contracted"),
  N_total = c(91638, 91638,
              157575, 157575,
              120792, 120792),
  N_event = c(1733, 1882,
              1457, 1645,
              1732, 2356))
## input counts for each analysis
analysis_counts <- data.frame(
  trait = c("cave","eye_loss","pigment_loss"),
  N_focal = c(1733, 1457, 1732),      # expanded genes in focal group
  N_focal_total = c(91638, 157575, 120792),
  N_control = c(2925, 3262, 3095),    # expanded genes in control group
  N_control_total = c(244368, 300825, 256683))
## compare proportions using prop.test (chi-squared for proportions)
analysis_counts$p_value <- NA
analysis_counts$chi2_stat <- NA
for(i in 1:nrow(analysis_counts)){
  counts <- c(analysis_counts$N_focal[i], analysis_counts$N_control[i])
  totals <- c(analysis_counts$N_focal_total[i], analysis_counts$N_control_total[i])
  test <- prop.test(counts, totals)
  analysis_counts$p_value[i] <- test$p.value
  analysis_counts$chi2_stat[i] <- test$statistic}
analysis_counts
## input data (expanded gene counts)
analysis_counts <- data.frame(
  trait = c("cave", "eye_loss", "pigment_loss"),
  N_focal = c(1733, 1457, 1732),
  N_focal_total = c(91638, 157575, 120792),
  N_control = c(2925, 3262, 3095),
  N_control_total = c(244368, 300825, 256683))
## compute proportions and odds ratios
analysis_counts <- analysis_counts %>%
  mutate(
    p_focal = N_focal / N_focal_total,
    p_control = N_control / N_control_total,
    diff_prop = p_focal - p_control,
    odds_ratio = (p_focal / (1 - p_focal)) / (p_control / (1 - p_control)),
    log_or = log(odds_ratio))
analysis_counts
## compare the *effect sizes* (log odds ratios) across analyses
# We’ll test if cave vs surface has a significantly larger log odds ratio
# than eye_loss or pigment_loss.
# Compute standard errors for log odds ratios using prop.test data
get_logor_se <- function(success1, total1, success2, total2) {
  p1 <- success1 / total1
  p2 <- success2 / total2
  1/success1 + 1/(total1 - success1) + 1/success2 + 1/(total2 - success2)}
analysis_counts$se_logor <- sqrt(mapply(get_logor_se,
                                        analysis_counts$N_focal,
                                        analysis_counts$N_focal_total,
                                        analysis_counts$N_control,
                                        analysis_counts$N_control_total))
## Z-tests comparing cave vs the other traits
compare_effects <- function(df, trait1, trait2) {
  diff <- df$log_or[df$trait == trait1] - df$log_or[df$trait == trait2]
  se <- sqrt(df$se_logor[df$trait == trait1]^2 + df$se_logor[df$trait == trait2]^2)
  z <- diff / se
  p <- 2 * (1 - pnorm(abs(z)))
  data.frame(compare = paste(trait1, "vs", trait2),
             logOR_diff = diff, z_value = z, p_value = p)}
comparison_results <- rbind(
  compare_effects(analysis_counts, "cave", "eye_loss"),
  compare_effects(analysis_counts, "cave", "pigment_loss"))
comparison_results


#### DO THE SAME FOR ABSREL GENES UNDER SELECTION ####
## input summary counts for each analysis
analysis_counts_sel <- data.frame(
  trait = c("cave", "pigment_loss", "eye_loss"),
  N_focal = c(2358, 2743, 4265),
  N_focal_total = c(25019, 32157, 53906),
  N_control = c(5743, 5561, 7559),
  N_control_total = c(83539, 81551, 101128))
## compute proportions and odds ratios
analysis_counts_sel <- analysis_counts_sel %>%
  mutate(
    p_focal = N_focal / N_focal_total,
    p_control = N_control / N_control_total,
    odds_ratio = (p_focal / (1 - p_focal)) / (p_control / (1 - p_control)),
    log_or = log(odds_ratio))
## compute standard errors for log-odds ratios
get_logor_se <- function(success1, total1, success2, total2) {
  1/success1 + 1/(total1 - success1) + 1/success2 + 1/(total2 - success2)}
analysis_counts_sel$se_logor <- sqrt(mapply(
  get_logor_se,
  analysis_counts_sel$N_focal,
  analysis_counts_sel$N_focal_total,
  analysis_counts_sel$N_control,
  analysis_counts_sel$N_control_total))
## compare effect sizes between traits (Z-tests)
compare_effects <- function(df, trait1, trait2) {
  diff <- df$log_or[df$trait == trait1] - df$log_or[df$trait == trait2]
  se <- sqrt(df$se_logor[df$trait == trait1]^2 + df$se_logor[df$trait == trait2]^2)
  z <- diff / se
  p <- 2 * (1 - pnorm(abs(z)))
  data.frame(compare = paste(trait1, "vs", trait2), logOR_diff = diff, z_value = z, p_value = p)}
# compare cave vs other regressive traits
comparison_sel <- rbind(compare_effects(analysis_counts_sel, "cave", "pigment_loss"), compare_effects(analysis_counts_sel, "cave", "eye_loss"))
## view results
analysis_counts_sel
comparison_sel

######### SELECTION ANALYSES ######### 
######### ABSREL ######### 
### CAVE SPECIES ### 
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/01_cave_species")
cave.sel.unfil <- read.csv("cave_species_absrel.csv") 
str(cave.sel.unfil)
cave.sel.unfil.1 <- cave.sel.unfil %>%
  mutate(F_nosel = .data[["allF"]] - .data[["selectionF"]])
cave.sel.unfil.2 <- cave.sel.unfil.1 %>%
  mutate(B_nosel = .data[["allB"]] - .data[["selectionB"]])
str(cave.sel.unfil.2)
cave.sel <- cave.sel.unfil.2 %>%
  dplyr::select(og,selectionF,selectionB,F_nosel,B_nosel)
str(cave.sel)
cave.sel <- cave.sel %>% 
  dplyr::rename("OG" = "og", "F_sel" = "selectionF","B_sel" = "selectionB")
str(cave.sel)
#cave.sel <- read.csv("cave_species_absrel_cont.csv")
str(cave.sel)
cave.sel.fisher <- pivot_longer(cave.sel, cols = 2:5, names_to = "name",
             values_to = "Freq", values_drop_na = FALSE)
str(cave.sel.fisher)

### Do F branches have the same proportion of selected genes as B branches?
# Start from raw data
overall <- cave.sel.fisher[, 2:3]
# Summarise counts per category
overall <- overall %>%
  group_by(name) %>%
  summarise(total_counts = sum(Freq), .groups = "drop") %>%
  separate(name, into = c("branch", "selection"), sep = "_")
# Store data for plotting BEFORE turning into xtabs
plot_data <- overall
# Turn into contingency table for statistical test
xtab <- xtabs(total_counts ~ branch + selection, data = overall)
# Chi-squared test
chisq.test(xtab)
# Row-wise percentages (for interpretation)
prop.table(xtab, 1) * 100
# cave species contain significantly more genes that are under selection (x2 = 181, P > 0.0001, 9.4% vs. 6.9%)
### make a Figure of this
# Make a stacked barplot--> it will be in %!
### summarize counts per category
overall <- cave.sel.fisher[, 2:3] %>%
  group_by(name) %>%
  summarise(total_counts = sum(Freq), .groups = "drop") %>%
  separate(name, into = c("branch", "selection"), sep = "_")
# Save for plotting
plot_data <- overall
### create contingency table and run chi-squared test
xtab <- xtabs(total_counts ~ branch + selection, data = overall)
chisq.test(xtab)
# cave species contain significantly more genes that are under selection (x2 = 181, P > 0.0001, 9.4% vs. 6.9%)
### percentages (for interpretation)
overall.perc <- prop.table(xtab, 1) * 100
print(overall.perc)
### base R stacked barplot (%)
# Prepare matrix for barplot
bar_data <- plot_data %>%
  group_by(branch) %>%
  mutate(percent = total_counts / sum(total_counts) * 100) %>%
  ungroup() %>%  # IMPORTANT: ungroup before select()
  dplyr::select(branch, selection, percent) %>%
  pivot_wider(names_from = selection, values_from = percent) %>%
  column_to_rownames("branch") %>%
  as.matrix()

#### FIGURE 2f ####
# Plot using base R
barplot(
  t(bar_data),
  beside = FALSE,
  col = c("#D3D3D3", "#66CDAA"),  # non-selected = gray, selected = orange
  legend = TRUE,
  ylim = c(0, 100),
  ylab = "Percentage of Genes",
  main = "Proportion of Selected vs. Non-selected Genes in F vs. B Branches")

### Do F branches share significantly more genes under selection than B branches?
## Assess for each OG individually
cave.sel.fisher <- cave.sel.fisher %>% 
  separate(name, into = c("branch", "selection"), sep="_")
str(cave.sel.fisher)
cave.sel.fisher$OG <- as.factor(cave.sel.fisher$OG)
cave.sel.fisher$branch <- as.factor(cave.sel.fisher$branch)
cave.sel.fisher$selection <- as.factor(cave.sel.fisher$selection)
str(cave.sel.fisher)
cave.sel.fisher <- as.data.frame(cave.sel.fisher)
str(cave.sel.fisher)
selected.table.2 <- cave.sel.fisher %>%
  group_by(OG) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(branch, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$selection
    return(M)}))
selected.table.3 <- selected.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
selected.table.3
sel.ratios <- pivot_wider(cave.sel.fisher, names_from = branch,
                          values_from = Freq)
colnames(sel.ratios) <- c("OG", "selection","F.branch","B.branch") 
str(sel.ratios)
# calculate proportion of selected F branches 
F.branch.tot <- sel.ratios %>% 
  group_by(OG) %>% 
  summarise(F.branch.tot = sum(F.branch)) 
# calculate proportion of selected B branches 
B.branch.tot <- sel.ratios %>% 
  group_by(OG) %>% 
  summarise(B.branch.tot = sum(B.branch)) 
# merge with data frame
F.branches.merged <- merge(sel.ratios,F.branch.tot, by="OG")  
F.B.branches.merged <- merge(F.branches.merged,B.branch.tot, by="OG")
# calculate if F branches or B branches have higher proportion of selected branches
F.B.sel.ratios <- F.B.branches.merged %>% 
  mutate(F.branch.ratio = F.branch/F.branch.tot) %>% 
  mutate(B.branch.ratio = B.branch/B.branch.tot) %>% 
  filter(selection=="sel") %>%
  mutate(F_vs_B = F.branch.ratio>B.branch.ratio)
F.B.sel.ratios <- F.B.sel.ratios[c("OG","F.branch.ratio","B.branch.ratio","F_vs_B")]
str(F.B.sel.ratios)
sel.OGs.significance.ratios <- merge(selected.table.3, F.B.sel.ratios, by="OG")
str(sel.OGs.significance.ratios)
sel.OGs.significance.ratios <- sel.OGs.significance.ratios %>%
  mutate(significance = pvalue<=0.05)
table <- table(
  CaveHigher = sel.OGs.significance.ratios$F_vs_B,
  Significant = sel.OGs.significance.ratios$significance)
print(table)
chisq.test(table)
## sharing is significantly higher in cave vs. surface species: X-squared = 162.13, df = 1, p-value < 2.2e-16

## combine with eggnog annotations
# read the data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/09_eggnog_annotation")
eggnog <- read.csv("all_OGs_cave_species_MM_9n4wx4so.emapper.annotations.tsv", header = TRUE, sep = "\t")
str(eggnog)
# add a new column with OG, gene symbol, description and GO and KEGG annotations
eggnog <- eggnog %>% 
  separate(query,c("OG", "sequence.ID"), sep = "-", remove = FALSE)
eggnog <- eggnog[,c("OG","Preferred_name","Description","GOs","KEGG_ko")]
str(eggnog)
sel.OGs.annot.significance.ratios <- merge(eggnog, sel.OGs.significance.ratios, by="OG", all.y=TRUE)
str(sel.OGs.annot.significance.ratios)
## write file with shared genes and their annotations
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/01_cave_species")
#write.csv(sel.OGs.annot.significance.ratios,"fisher_selected_OGs_ratios_cave_species.csv", quote=TRUE, row.names=FALSE)
## write file with shared genes significantly selected in F branches
str(sel.OGs.annot.significance.ratios)
cave.species.shared.sel <- sel.OGs.annot.significance.ratios %>%
  filter(pvalue<=0.05) %>%
  filter(F_vs_B=="TRUE")
cave.species.shared.sel.OGs <- cave.species.shared.sel[,c("OG","Preferred_name","Description","pvalue","F_vs_B")]
#write.table(cave.species.shared.sel.OGs,"cave_species_shared_sel_OGs.txt", quote=FALSE, row.names=FALSE)
#write.csv(cave.species.shared.sel.OGs,"cave_species_shared_sel_OGs.csv", quote=FALSE, row.names=FALSE)
### make a Figure for OGs that are shared 
environment <- c("cave","surface")
values <- c(77, 3)
# join to frame
df <- data.frame(environment, values)
print(df)

#### FIGURE 2h ####
ggplot(data = df, aes(x = environment, y = values, fill = environment)) + 
  geom_bar(stat = "identity") + 
  labs(x = "\n environment", y = "Number of shared selected genes \n", 
       title = "Selected genes evolved under convergence \n",
       fill = "Shared selected genes in") +
  scale_fill_manual(values=c("#ECD5E0", "#8C8C8C")) +
  theme_classic()

### PIGMENTATION LOSS ### 
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/02_pigmentation_loss")
cave.sel.unfil <- read.csv("pigmentation_loss_absrel.csv") 
str(cave.sel.unfil)
cave.sel.unfil.1 <- cave.sel.unfil %>%
  mutate(F_nosel = .data[["allF"]] - .data[["selectionF"]])
cave.sel.unfil.2 <- cave.sel.unfil.1 %>%
  mutate(B_nosel = .data[["allB"]] - .data[["selectionB"]])
str(cave.sel.unfil.2)
cave.sel <- cave.sel.unfil.2 %>%
  dplyr::select(og,selectionF,selectionB,F_nosel,B_nosel)
str(cave.sel)
cave.sel <- cave.sel %>% 
  dplyr::rename("OG" = "og", "F_sel" = "selectionF","B_sel" = "selectionB")
str(cave.sel)
cave.sel.fisher <- pivot_longer(cave.sel, cols = 2:5, names_to = "name",
                                values_to = "Freq", values_drop_na = FALSE)
str(cave.sel.fisher)
### Do F branches have the same proportion of selected genes as B branches?
overall <- cave.sel.fisher[,2:3]
str(overall)
overall <- overall %>% group_by(name) %>% summarise(total_counts = sum(Freq)) %>% data.table() 
str(overall)
overall <- overall %>% 
  separate(name, into = c("branch", "selection"), sep="_")
overall <- xtabs(total_counts ~ branch+selection, 
                 data=overall)
str(overall)
chisq.test(overall)
prop.table(overall, 1)*100
# non-pigmented species contain significantly more genes that are under selection (x2 = 99.5, P > 0.0001, 8.5% vs. 6.8%)
### Do F branches share significantly more genes under selection than B branches?
## Assess for each OG individually
cave.sel.fisher <- cave.sel.fisher %>% 
  separate(name, into = c("branch", "selection"), sep="_")
str(cave.sel.fisher)
cave.sel.fisher$OG <- as.factor(cave.sel.fisher$OG)
cave.sel.fisher$branch <- as.factor(cave.sel.fisher$branch)
cave.sel.fisher$selection <- as.factor(cave.sel.fisher$selection)
str(cave.sel.fisher)
cave.sel.fisher <- as.data.frame(cave.sel.fisher)
str(cave.sel.fisher)
selected.table.2 <- cave.sel.fisher %>%
  group_by(OG) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(branch, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$selection
    return(M)}))
selected.table.3 <- selected.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
selected.table.3
sel.ratios <- pivot_wider(cave.sel.fisher, names_from = branch,
                          values_from = Freq)
colnames(sel.ratios) <- c("OG", "selection","F.branch","B.branch") 
str(sel.ratios)
# calculate proportion of selected F branches 
F.branch.tot <- sel.ratios %>% 
  group_by(OG) %>% 
  summarise(F.branch.tot = sum(F.branch)) 
# calculate proportion of selected B branches 
B.branch.tot <- sel.ratios %>% 
  group_by(OG) %>% 
  summarise(B.branch.tot = sum(B.branch)) 
# merge with data frame
F.branches.merged <- merge(sel.ratios,F.branch.tot, by="OG")  
F.B.branches.merged <- merge(F.branches.merged,B.branch.tot, by="OG")
# calculate if F branches or B branches have higher proportion of selected branches
F.B.sel.ratios <- F.B.branches.merged %>% 
  mutate(F.branch.ratio = F.branch/F.branch.tot) %>% 
  mutate(B.branch.ratio = B.branch/B.branch.tot) %>% 
  filter(selection=="sel") %>%
  mutate(F_vs_B = F.branch.ratio>B.branch.ratio)
F.B.sel.ratios <- F.B.sel.ratios[c("OG","F.branch.ratio","B.branch.ratio","F_vs_B")]
str(F.B.sel.ratios)
sel.OGs.significance.ratios <- merge(selected.table.3, F.B.sel.ratios, by="OG")
str(sel.OGs.significance.ratios)
sel.OGs.significance.ratios <- sel.OGs.significance.ratios %>%
  mutate(significance = pvalue<=0.05)
## combine with eggnog annotations
# read the data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/09_eggnog_annotation")
eggnog <- read.csv("all_OGs_pigmentation_loss_MM_fe2v4rjm.emapper.annotations.csv", header = TRUE)
str(eggnog)
# add a new column with OG, gene symbol, description and GO and KEGG annotations
eggnog <- eggnog %>% 
  separate(query,c("OG", "sequence.ID"), sep = "-", remove = FALSE)
eggnog <- eggnog[,c("OG","Preferred_name","Description","GOs","KEGG_ko")]
str(eggnog)
sel.OGs.annot.significance.ratios <- merge(eggnog, sel.OGs.significance.ratios, by="OG", all.y=TRUE)
str(sel.OGs.annot.significance.ratios)
## write file with shared genes and their annotations
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/02_pigmentation_loss")
#write.csv(sel.OGs.annot.significance.ratios,"fisher_selected_OGs_ratios_pigmentation_loss.csv", quote=TRUE, row.names=FALSE)
## write file with shared genes significantly selected in F branches
str(sel.OGs.annot.significance.ratios)
cave.species.shared.sel <- sel.OGs.annot.significance.ratios %>%
  filter(pvalue<=0.05) %>%
  filter(F_vs_B=="TRUE")
cave.species.shared.sel.OGs <- cave.species.shared.sel[,c("OG","Preferred_name","Description","pvalue","F_vs_B")]
#write.table(cave.species.shared.sel.OGs,"pigmentation_loss_shared_sel_OGs.txt", quote=FALSE, row.names=FALSE)
#write.csv(cave.species.shared.sel.OGs,"pigmentation_loss_shared_sel_OGs.csv", quote=FALSE, row.names=FALSE)

### EYE LOSS ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/03_eye_loss")
cave.sel.unfil <- read.csv("eye_loss_absrel.csv") 
str(cave.sel.unfil)
cave.sel.unfil.1 <- cave.sel.unfil %>%
  mutate(F_nosel = .data[["allF"]] - .data[["selectionF"]])
cave.sel.unfil.2 <- cave.sel.unfil.1 %>%
  mutate(B_nosel = .data[["allB"]] - .data[["selectionB"]])
str(cave.sel.unfil.2)
cave.sel <- cave.sel.unfil.2 %>%
  dplyr::select(og,selectionF,selectionB,F_nosel,B_nosel)
str(cave.sel)
cave.sel <- cave.sel %>% 
  dplyr::rename("OG" = "og", "F_sel" = "selectionF","B_sel" = "selectionB")
str(cave.sel)
cave.sel.fisher <- pivot_longer(cave.sel, cols = 2:5, names_to = "name",
                                values_to = "Freq", values_drop_na = FALSE)
str(cave.sel.fisher)
### Do F branches have the same proportion of selected genes as B branches?
overall <- cave.sel.fisher[,2:3]
str(overall)
overall <- overall %>% group_by(name) %>% summarise(total_counts = sum(Freq)) %>% data.table() 
str(overall)
overall <- overall %>% 
  separate(name, into = c("branch", "selection"), sep="_")
overall <- xtabs(total_counts ~ branch+selection, 
                 data=overall)
str(overall)
chisq.test(overall)
prop.table(overall, 1)*100
# cave species contain significantly more genes that are under selection (x2 = 9.5, P > 0.002, 7.9% vs. 7.5%)
### Do F branches share significantly more genes under selection than B branches?
## Assess for each OG individually
cave.sel.fisher <- cave.sel.fisher %>% 
  separate(name, into = c("branch", "selection"), sep="_")
str(cave.sel.fisher)
cave.sel.fisher$OG <- as.factor(cave.sel.fisher$OG)
cave.sel.fisher$branch <- as.factor(cave.sel.fisher$branch)
cave.sel.fisher$selection <- as.factor(cave.sel.fisher$selection)
str(cave.sel.fisher)
cave.sel.fisher <- as.data.frame(cave.sel.fisher)
str(cave.sel.fisher)
selected.table.2 <- cave.sel.fisher %>%
  group_by(OG) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(branch, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$selection
    return(M)}))
selected.table.3 <- selected.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
selected.table.3
sel.ratios <- pivot_wider(cave.sel.fisher, names_from = branch,
                          values_from = Freq)
colnames(sel.ratios) <- c("OG", "selection","F.branch","B.branch") 
str(sel.ratios)
# calculate proportion of selected F branches 
F.branch.tot <- sel.ratios %>% 
  group_by(OG) %>% 
  summarise(F.branch.tot = sum(F.branch)) 
# calculate proportion of selected B branches 
B.branch.tot <- sel.ratios %>% 
  group_by(OG) %>% 
  summarise(B.branch.tot = sum(B.branch)) 
# merge with data frame
F.branches.merged <- merge(sel.ratios,F.branch.tot, by="OG")  
F.B.branches.merged <- merge(F.branches.merged,B.branch.tot, by="OG")
# calculate if F branches or B branches have higher proportion of selected branches
F.B.sel.ratios <- F.B.branches.merged %>% 
  mutate(F.branch.ratio = F.branch/F.branch.tot) %>% 
  mutate(B.branch.ratio = B.branch/B.branch.tot) %>% 
  filter(selection=="sel") %>%
  mutate(F_vs_B = F.branch.ratio>B.branch.ratio)
F.B.sel.ratios <- F.B.sel.ratios[c("OG","F.branch.ratio","B.branch.ratio","F_vs_B")]
str(F.B.sel.ratios)
sel.OGs.significance.ratios <- merge(selected.table.3, F.B.sel.ratios, by="OG")
str(sel.OGs.significance.ratios)
sel.OGs.significance.ratios <- sel.OGs.significance.ratios %>%
  mutate(significance = pvalue<=0.05)
## combine with eggnog annotations
# read the data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/09_eggnog_annotation")
eggnog <- read.csv("all_OGs_eye_loss_MM_fe2v4rjm.emapper.annotations.csv", header = TRUE)
str(eggnog)
# add a new column with OG, gene symbol, description and GO and KEGG annotations
eggnog <- eggnog %>% 
  separate(query,c("OG", "sequence.ID"), sep = "-", remove = FALSE)
eggnog <- eggnog[,c("OG","Preferred_name","Description","GOs","KEGG_ko")]
str(eggnog)
sel.OGs.annot.significance.ratios <- merge(eggnog, sel.OGs.significance.ratios, by="OG", all.y=TRUE)
str(sel.OGs.annot.significance.ratios)
## write file with shared genes and their annotations
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/03_eye_loss")
#write.csv(sel.OGs.annot.significance.ratios,"fisher_selected_OGs_ratios_eye_loss.csv", quote=TRUE, row.names=FALSE)
## write file with shared genes significantly selected in F branches
str(sel.OGs.annot.significance.ratios)
cave.species.shared.sel <- sel.OGs.annot.significance.ratios %>%
  filter(pvalue<=0.05) %>%
  filter(F_vs_B=="TRUE")
cave.species.shared.sel.OGs <- cave.species.shared.sel[,c("OG","Preferred_name","Description","pvalue","F_vs_B")]
#write.table(cave.species.shared.sel.OGs,"eye_loss_shared_sel_OGs.txt", quote=FALSE, row.names=FALSE)
#write.csv(cave.species.shared.sel.OGs,"eye_loss_shared_sel_OGs.csv", quote=FALSE, row.names=FALSE)

### LONGEVITY ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/04_longevity")
cave.sel.unfil <- read.csv("longevity_absrel.csv") 
str(cave.sel.unfil)
cave.sel.unfil.1 <- cave.sel.unfil %>%
  mutate(F_nosel = .data[["allF"]] - .data[["selectionF"]])
cave.sel.unfil.2 <- cave.sel.unfil.1 %>%
  mutate(B_nosel = .data[["allB"]] - .data[["selectionB"]])
str(cave.sel.unfil.2)
cave.sel <- cave.sel.unfil.2 %>%
  dplyr::select(og,selectionF,selectionB,F_nosel,B_nosel)
str(cave.sel)
cave.sel <- cave.sel %>% 
  dplyr::rename("OG" = "og", "F_sel" = "selectionF","B_sel" = "selectionB")
str(cave.sel)
cave.sel.fisher <- pivot_longer(cave.sel, cols = 2:5, names_to = "name",
                                values_to = "Freq", values_drop_na = FALSE)
str(cave.sel.fisher)
### Do F branches have the same proportion of selected genes as B branches?
overall <- cave.sel.fisher[,2:3]
str(overall)
overall <- overall %>% group_by(name) %>% summarise(total_counts = sum(Freq)) %>% data.table() 
str(overall)
overall <- overall %>% 
  separate(name, into = c("branch", "selection"), sep="_")
overall <- xtabs(total_counts ~ branch+selection, 
                 data=overall)
str(overall)
chisq.test(overall)
prop.table(overall, 1)*100
overall.perc <- prop.table(overall, 1)*100
# longevity species contain significantly more genes that are under selection (x2 = 230, P > 0.0001, 9.5% vs. 7.5%)
### make a Figure of this
# Make a stacked barplot--> it will be in %!
# Data Visualization Of Contingency Table With ggplot2 (Stacked Bar Graph):

#### FIGURE 4d ####
df <- as.data.frame(t(overall.perc))
ggplot(data = df, aes(x = branch, y = Freq, fill = selection)) + 
  geom_bar(stat = "identity") + 
  labs(x = "\n Answer", y = "Percentage \n", 
       title = "Genes under selection \n",
       fill = "Positive selection") +
  scale_fill_manual(values=c("#D3D3D3", "#66CDAA")) +
  theme_classic()

### Do F branches share significantly more genes under selection than B branches?
## Assess for each OG individually
cave.sel.fisher <- cave.sel.fisher %>% 
  separate(name, into = c("branch", "selection"), sep="_")
str(cave.sel.fisher)
cave.sel.fisher$OG <- as.factor(cave.sel.fisher$OG)
cave.sel.fisher$branch <- as.factor(cave.sel.fisher$branch)
cave.sel.fisher$selection <- as.factor(cave.sel.fisher$selection)
str(cave.sel.fisher)
cave.sel.fisher <- as.data.frame(cave.sel.fisher)
str(cave.sel.fisher)
selected.table.2 <- cave.sel.fisher %>%
  group_by(OG) %>%
  nest() %>%
  mutate(M = map(data, function(dat){
    dat2 <- dat %>% spread(branch, Freq)
    M <- as.matrix(dat2[, -1])
    row.names(M) <- dat2$selection
    return(M)}))
selected.table.3 <- selected.table.2 %>%
  mutate(pvalue = map_dbl(M, ~fisher.test(.x)$p.value)) %>%
  dplyr::select(-data, -M) %>%
  ungroup()
selected.table.3
sel.ratios <- pivot_wider(cave.sel.fisher, names_from = branch,
                          values_from = Freq)
colnames(sel.ratios) <- c("OG", "selection","F.branch","B.branch") 
str(sel.ratios)
# calculate proportion of selected F branches 
F.branch.tot <- sel.ratios %>% 
  group_by(OG) %>% 
  summarise(F.branch.tot = sum(F.branch)) 
# calculate proportion of selected B branches 
B.branch.tot <- sel.ratios %>% 
  group_by(OG) %>% 
  summarise(B.branch.tot = sum(B.branch)) 
# merge with data frame
F.branches.merged <- merge(sel.ratios,F.branch.tot, by="OG")  
F.B.branches.merged <- merge(F.branches.merged,B.branch.tot, by="OG")
# calculate if F branches or B branches have higher proportion of selected branches
F.B.sel.ratios <- F.B.branches.merged %>% 
  mutate(F.branch.ratio = F.branch/F.branch.tot) %>% 
  mutate(B.branch.ratio = B.branch/B.branch.tot) %>% 
  filter(selection=="sel") %>%
  mutate(F_vs_B = F.branch.ratio>B.branch.ratio)
F.B.sel.ratios <- F.B.sel.ratios[c("OG","F.branch.ratio","B.branch.ratio","F_vs_B")]
str(F.B.sel.ratios)
sel.OGs.significance.ratios <- merge(selected.table.3, F.B.sel.ratios, by="OG")
str(sel.OGs.significance.ratios)
sel.OGs.significance.ratios <- sel.OGs.significance.ratios %>%
  mutate(significance = pvalue<=0.05)
## combine with eggnog annotations
# read the data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/09_eggnog_annotation")
eggnog <- read.csv("all_OGs_longevity_MM_j29p4hp0.emapper.annotations.csv", header = TRUE)
str(eggnog)
# add a new column with OG, gene symbol, description and GO and KEGG annotations
eggnog <- eggnog %>% 
  separate(query,c("OG", "sequence.ID"), sep = "-", remove = FALSE)
eggnog <- eggnog[,c("OG","Preferred_name","Description","GOs","KEGG_ko")]
str(eggnog)
sel.OGs.annot.significance.ratios <- merge(eggnog, sel.OGs.significance.ratios, by="OG", all.y=TRUE)
str(sel.OGs.annot.significance.ratios)
## write file with shared genes and their annotations
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/04_longevity")
#write.csv(sel.OGs.annot.significance.ratios,"fisher_selected_OGs_ratios_longevity.csv", quote=TRUE, row.names=FALSE)
## write file with shared genes significantly selected in F branches
str(sel.OGs.annot.significance.ratios)
cave.species.shared.sel <- sel.OGs.annot.significance.ratios %>%
  filter(pvalue<=0.05) %>%
  filter(F_vs_B=="TRUE")
cave.species.shared.sel.OGs <- cave.species.shared.sel[,c("OG","Preferred_name","Description","pvalue","F_vs_B")]
#write.table(cave.species.shared.sel.OGs,"longevity_shared_sel_OGs.txt", quote=FALSE, row.names=FALSE)
#write.csv(cave.species.shared.sel.OGs,"longevity_shared_sel_OGs.csv", quote=FALSE, row.names=FALSE)


######### PAML POSITIVE SELECTION ######### 
## perform hypergeometric tests as provided in SuperExactTest to test for overlap in selcted genes between species
# input files are lists of selected genes for each species
### CAVE SPECIES ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/02_cave_species")
Asmec <- read.table("Asmec.txt", header=TRUE)
Prang <- read.table("Prang.txt", header=TRUE)
Trros <- read.table("Trros.txt", header=TRUE)
Sians <- read.table("Sians.txt", header=TRUE)
Sirhi <- read.table("Sirhi.txt", header=TRUE)
Luden <- read.table("Luden.txt", header=TRUE)
cave_list <- list(as.character(Asmec$OG_symbol),as.character(Prang$OG_symbol),as.character(Trros$OG_symbol),as.character(Sians$OG_symbol),as.character(Sirhi$OG_symbol),as.character(Luden$OG_symbol))
names(cave_list) <- c("Asmec", "Prang","Trros","Sians","Sirhi","Luden")
str(cave_list) 
(length.gene.sets<-sapply(cave_list,length))
# create list with number of genes (assume 20,000 genes per species)
total=6064
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
# probability density distribution of possible interaction sizes
(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))
common.genes=SuperExactTest::intersect(Asmec[[1]], Prang[[1]], Trros[[1]], Sians[[1]], Sirhi[[1]], Luden[[1]])
(num.observed.overlap=length(common.genes))
fit=MSET(cave_list, n=total, lower.tail=FALSE)
fit$FE
fit$p.value
# Looking for all possible interactions across all species
res_all=supertest(cave_list, n=total)
plot(res_all, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15),keep.empty.intersections=FALSE)
plot(res_all, Layout="landscape", degree=2:6, sort.by="size", margin=c(0.5,5,1,2), keep.empty.intersections=FALSE,show.overlap.size=FALSE,minMinusLog10PValue=4)
#write.csv(summary(res_all)$Table, file="res_all.csv", row.names=FALSE)
summary(res_all)
## extract all significantly shared OGs and all shared by at least 2 species
shared.N <- read.table("cave_species_shared_count.txt", header=TRUE)
str(shared.N)
shared.sign <- read.table("cave_species_significantly_shared.txt", header=TRUE)
str(shared.sign)
shared.sign.N <- merge(shared.N, shared.sign, by="OG")
str(shared.sign.N)
## combine with eggnog annotations
# read the data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/09_eggnog_annotation")
eggnog <- read.csv("all_OGs_cave_species_MM_9n4wx4so.emapper.annotations.tsv", header = TRUE, sep = "\t")
str(eggnog)
# add a new column with OG, gene symbol, description and GO and KEGG annotations
eggnog <- eggnog %>% 
  separate(query,c("OG", "sequence.ID"), sep = "-", remove = FALSE)
eggnog <- eggnog[,c("OG","Preferred_name","Description")]
str(eggnog)
sel.shared.sign.N <- merge(eggnog, shared.sign.N, by="OG", all.y=TRUE)
str(sel.shared.sign.N)
## write file with shared genes and their annotations
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/02_cave_species")
#write.csv(sel.shared.sign.N,"cave_species_sign_shared_genes_PAML.csv", quote=TRUE, row.names=FALSE)

### PIGMENTATION LOSS ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/03_pigmentation_loss")
Ammex <- read.table("Ammex.txt", header=TRUE)
Asmec <- read.table("Asmec.txt", header=TRUE)
Hegla <- read.table("Hegla.txt", header=TRUE)
Luden <- read.table("Luden.txt", header=TRUE)
Prang <- read.table("Prang.txt", header=TRUE)
Sians <- read.table("Sians.txt", header=TRUE)
Sirhi <- read.table("Sirhi.txt", header=TRUE)
Trros <- read.table("Trros.txt", header=TRUE)
cave_list <- list(as.character(Ammex$OG_symbol),as.character(Asmec$OG_symbol),as.character(Hegla$OG_symbol),
                  as.character(Luden$OG_symbol),as.character(Prang$OG_symbol),as.character(Sians$OG_symbol),
                  as.character(Sirhi$OG_symbol),as.character(Trros$OG_symbol))
names(cave_list) <- c("Ammex", "Asmec", "Hegla", "Luden", "Prang", "Sians","Sirhi","Trros")
str(cave_list) 
(length.gene.sets<-sapply(cave_list,length))
# create list with number of genes (assume 20,000 genes per species)
total=5728
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
# probability density distribution of possible interaction sizes
(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))
common.genes=SuperExactTest::intersect(Asmec[[1]], Prang[[1]], Trros[[1]], Sians[[1]], Sirhi[[1]], Luden[[1]])
(num.observed.overlap=length(common.genes))
fit=MSET(cave_list, n=total, lower.tail=FALSE)
fit$FE
fit$p.value
# Looking for all possible interactions across all species
res_all=supertest(cave_list, n=total)
plot(res_all, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15),keep.empty.intersections=FALSE)
plot(res_all, Layout="landscape", degree=3:8, sort.by="size", margin=c(0.5,5,1,2), keep.empty.intersections=FALSE,show.overlap.size=FALSE,minMinusLog10PValue=4)
#write.csv(summary(res_all)$Table, file="res_all.csv", row.names=FALSE)
summary(res_all)
## extract all significantly shared OGs and all shared by at least 2 species
shared.N <- read.table("pigmentation_loss_shared_count.txt", header=TRUE)
str(shared.N)
shared.sign <- read.table("pigmentation_loss_significantly_shared.txt", header=TRUE)
str(shared.sign)
shared.sign.N <- merge(shared.N, shared.sign, by="OG")
str(shared.sign.N)
## combine with eggnog annotations
# read the data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/09_eggnog_annotation")
eggnog <- read.csv("all_OGs_pigmentation_loss_MM_fe2v4rjm.emapper.annotations.csv", header = TRUE)
str(eggnog)
# add a new column with OG, gene symbol, description and GO and KEGG annotations
eggnog <- eggnog %>% 
  separate(query,c("OG", "sequence.ID"), sep = "-", remove = FALSE)
eggnog <- eggnog[,c("OG","Preferred_name","Description")]
str(eggnog)
sel.shared.sign.N <- merge(eggnog, shared.sign.N, by="OG", all.y=TRUE)
str(sel.shared.sign.N)
## write file with shared genes and their annotations
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/03_pigmentation_loss")
#write.csv(sel.shared.sign.N,"pigmentation_loss_sign_shared_genes_PAML.csv", quote=TRUE, row.names=FALSE)

### EYE LOSS ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/04_eye_loss")
Asmec <- read.table("Asmec.txt", header=TRUE)
Chasi <- read.table("Chasi.txt", header=TRUE)
Cocri <- read.table("Cocri.txt", header=TRUE)
Geser <- read.table("Geser.txt", header=TRUE)
Hegla <- read.table("Hegla.txt", header=TRUE)
Luden <- read.table("Luden.txt", header=TRUE)
Miuni <- read.table("Miuni.txt", header=TRUE)
Nagal <- read.table("Nagal.txt", header=TRUE)
Prang <- read.table("Prang.txt", header=TRUE)
Sians <- read.table("Sians.txt", header=TRUE)
Trros <- read.table("Trros.txt", header=TRUE)
cave_list <- list(as.character(Ammex$OG_symbol),as.character(Chasi$OG_symbol),as.character(Cocri$OG_symbol),
                  as.character(Geser$OG_symbol),as.character(Hegla$OG_symbol),as.character(Luden$OG_symbol),
                  as.character(Miuni$OG_symbol),as.character(Nagal$OG_symbol),as.character(Prang$OG_symbol),
                  as.character(Sians$OG_symbol),as.character(Trros$OG_symbol))
names(cave_list) <- c("Ammex", "Chasi", "Cocri", "Geser", "Hegla", "Luden",
                      "Miuni", "Nagal", "Prang", "Sians", "Trros")
str(cave_list) 
(length.gene.sets<-sapply(cave_list,length))
# create list with number of genes (assume 20,000 genes per species)
total=5448
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
# probability density distribution of possible interaction sizes
(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))
common.genes=SuperExactTest::intersect(Asmec[[1]], Chasi[[1]], Cocri[[1]], 
                                       Geser[[1]], Hegla[[1]], Luden[[1]],
                                       Miuni[[1]], Nagal[[1]], Prang[[1]],
                                       Sians[[1]], Trros[[1]])
(num.observed.overlap=length(common.genes))
fit=MSET(cave_list, n=total, lower.tail=FALSE)
fit$FE
fit$p.value
# Looking for all possible interactions across all species
res_all=supertest(cave_list, n=total)
plot(res_all, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15),keep.empty.intersections=FALSE)
plot(res_all, Layout="landscape", degree=6:11, sort.by="size", margin=c(0.5,5,1,2), keep.empty.intersections=FALSE,show.overlap.size=FALSE,minMinusLog10PValue=4)
#write.csv(summary(res_all)$Table, file="res_all.csv", row.names=FALSE)
summary(res_all)
## extract all significantly shared OGs and all shared by at least 2 species
shared.N <- read.table("eye_loss_shared_count.txt", header=TRUE)
str(shared.N)
shared.sign <- read.table("eye_loss_significantly_shared.txt", header=TRUE)
str(shared.sign)
shared.sign.N <- merge(shared.N, shared.sign, by="OG")
str(shared.sign.N)
## combine with eggnog annotations
# read the data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/09_eggnog_annotation")
eggnog <- read.csv("all_OGs_eye_loss_MM_fe2v4rjm.emapper.annotations.csv", header = TRUE)
str(eggnog)
# add a new column with OG, gene symbol, description and GO and KEGG annotations
eggnog <- eggnog %>% 
  separate(query,c("OG", "sequence.ID"), sep = "-", remove = FALSE)
eggnog <- eggnog[,c("OG","Preferred_name","Description")]
str(eggnog)
sel.shared.sign.N <- merge(eggnog, shared.sign.N, by="OG", all.y=TRUE)
str(sel.shared.sign.N)
## write file with shared genes and their annotations
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/04_eye_loss")
#write.csv(sel.shared.sign.N,"eye_loss_sign_shared_genes_PAML.csv", quote=TRUE, row.names=FALSE)

### LONGEVITY ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/05_longevity")
Acbae <- read.table("Acbae.txt", header=TRUE)
Cacar <- read.table("Cacar.txt", header=TRUE)
Elmax <- read.table("Elmax.txt", header=TRUE)
Hegla <- read.table("Hegla.txt", header=TRUE)
Hosap <- read.table("Hosap.txt", header=TRUE)
Lacha <- read.table("Lacha.txt", header=TRUE)
Myluc <- read.table("Myluc.txt", header=TRUE)
Ororc <- read.table("Ororc.txt", header=TRUE)
Prang <- read.table("Prang.txt", header=TRUE)
Rhtyp <- read.table("Rhtyp.txt", header=TRUE)
Tecar <- read.table("Tecar.txt", header=TRUE)
cave_list <- list(as.character(Acbae$OG_symbol),as.character(Cacar$OG_symbol),as.character(Elmax$OG_symbol),as.character(Hegla$OG_symbol),
                  as.character(Hosap$OG_symbol),as.character(Lacha$OG_symbol),as.character(Myluc$OG_symbol),as.character(Ororc$OG_symbol),
                  as.character(Prang$OG_symbol),as.character(Rhtyp$OG_symbol),as.character(Tecar$OG_symbol))
names(cave_list) <- c("Acbae", "Cacar","Elmax","Hegla","Hosap","Lacha","Myluc", "Ororc","Prang","Rhtyp","Tecar")
str(cave_list) 
(length.gene.sets<-sapply(cave_list,length))
# create list with number of genes (assume 20,000 genes per species)
total=8084
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
# probability density distribution of possible interaction sizes
(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))
common.genes=SuperExactTest::intersect(Acbae[[1]], Cacar[[1]], Elmax[[1]], Hegla[[1]], Hosap[[1]], Lacha[[1]],Myluc[[1]], Ororc[[1]], 
                                       Prang[[1]], Rhtyp[[1]], Tecar[[1]])
(num.observed.overlap=length(common.genes))
fit=MSET(cave_list, n=total, lower.tail=FALSE)
fit$FE
fit$p.value
# Looking for all possible interactions across all species
res_all=supertest(cave_list, n=total)
plot(res_all, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15),keep.empty.intersections=FALSE)
plot(res_all, Layout="landscape", degree=6:11, sort.by="size", margin=c(0.5,5,1,2), keep.empty.intersections=FALSE,show.overlap.size=FALSE,minMinusLog10PValue=4)
#write.csv(summary(res_all)$Table, file="longevity_res_all.csv", row.names=FALSE)
summary(res_all)
## extract all significantly shared OGs and all shared by at least 2 species
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/05_longevity")
shared.N <- read.table("longevity_shared_count.txt", header=TRUE)
str(shared.N)
shared.sign <- read.table("longevity_significantly_shared.txt", header=TRUE)
str(shared.sign)
shared.sign.N <- merge(shared.N, shared.sign, by="OG")
str(shared.sign.N)
## combine with eggnog annotations
# read the data
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/06_expansion_contraction/09_eggnog_annotation")
eggnog <- read.csv("all_OGs_longevity_MM_j29p4hp0.emapper.annotations.csv", header = TRUE)
str(eggnog)
# add a new column with OG, gene symbol, description and GO and KEGG annotations
eggnog <- eggnog %>% 
  separate(query,c("OG", "sequence.ID"), sep = "-", remove = FALSE)
eggnog <- eggnog[,c("OG","Preferred_name","Description")]
str(eggnog)
sel.shared.sign.N <- merge(eggnog, shared.sign.N, by="OG", all.y=TRUE)
str(sel.shared.sign.N)
## write file with shared genes and their annotations
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/07_PAML/05_longevity")
#write.csv(sel.shared.sign.N,"longevity_sign_shared_genes_PAML.csv", quote=TRUE, row.names=FALSE)
## write file with shared genes significantly selected in F branches
str(sel.OGs.annot.significance.ratios)
cave.species.shared.sel <- sel.OGs.annot.significance.ratios %>%
  filter(pvalue<=0.05) %>%
  filter(F_vs_B=="TRUE")
cave.species.shared.sel.OGs <- cave.species.shared.sel[,c("OG","Preferred_name","Description","pvalue","F_vs_B")]
#write.table(cave.species.shared.sel.OGs,"longevity_shared_sel_OGs.txt", quote=FALSE, row.names=FALSE)
#write.csv(cave.species.shared.sel.OGs,"longevity_shared_sel_OGs.csv", quote=FALSE, row.names=FALSE)


######### RELAX ######### 
### CAVE SPECIES ###
## read in annotation file
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/01_eggnog_annotation_files")
eggnog.cave <- read.csv("all_OGs_cave_species_MM_9n4wx4so.emapper.annotations.csv", header = TRUE)
str(eggnog.cave)
# add a new column just with OG ID
eggnog.cave <- eggnog.cave %>% 
  separate(query,c("OG_ID", "sequence.ID"), sep = "-", remove = FALSE) %>%
  dplyr::select(OG_ID, Preferred_name,Description,GOs,KEGG_ko)
str(eggnog.cave)
## read in selection file
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/")
relax.cave <- read.table("cave_species_selected_OGs.txt", header=TRUE)
str(relax.cave)
selected.OGs <- merge(relax.cave, eggnog.cave, by="OG_ID", all.x=TRUE)
str(selected.OGs)
intensified <- selected.OGs  %>% filter(selection == "intensified") %>% 
  dplyr::select(OG_ID,Preferred_name,selection) 
relaxed <- selected.OGs  %>% filter(selection == "relaxed") %>% 
  dplyr::select(OG_ID,Preferred_name,selection) 
#write.table(selected.OGs, "cave_species_RELAX_annotated_selected_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(intensified, "cave_species_RELAX_annotated_intensified_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(relaxed, "cave_species_RELAX_annotated_relaxed_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
### check if relaxed selection is more prevalent in cave vs. surface species
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/general_descriptive")
cave.relax <- read.csv("cave_species_gen_des_1-5_combined.csv")
cave.tre <- read.tree("cave_species_list_mod_corr_R.tre")
plot(cave.tre)
str(cave.relax)
# Set rownames to species names (if not already done)
rownames(cave.relax) <- cave.relax$sp
# Check for perfect match in species names
if (!all(cave.tre$tip.label %in% rownames(cave.relax))) {
  stop("Mismatch between tree tip labels and data rownames!")}
# Sort the data to match tree tip order
cave.relax <- cave.relax[cave.tre$tip.label, ]
# Final sanity check
stopifnot(all(rownames(cave.relax) == cave.tre$tip.label))
### Statistical comparison of relaxed_ratio between cave (F) and surface (B) species ###
## simple t-test
# Check factor levels and counts
print(table(cave.relax$fb))
print(levels(as.factor(cave.relax$fb)))
# Summarize relaxed_ratio by group
summary_stats <- cave.relax %>%
  group_by(fb) %>%
  summarise(mean_relaxed = mean(relaxed_ratio), median_relaxed = median(relaxed_ratio), sd_relaxed = sd(relaxed_ratio), n = n())
print(summary_stats)
# Plot relaxed_ratio by group with boxplot + jitter

#### FIGURE 2i ####
ggplot(cave.relax, aes(x = fb, y = relaxed_ratio, color = fb)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  scale_color_manual(values = c("F" = "#E9D4DC", "B" = "#8B8B8B")) +
  labs(title = "Relaxed Ratio by Species Type",
    x = "Species Type (B = surface, F = cave)",
    y = "Relaxed Ratio") +
  theme_classic() +
  theme(legend.position = "none")

# Quick t-test for relaxed_ratio between groups (non-phylogenetic)
t_test_result <- t.test(relaxed_ratio ~ fb, data = cave.relax)
print(t_test_result)
## Phylogenetic ANOVA
phylANOVA_relaxed <- phylANOVA(tree = cave.tre, x = cave.relax$fb, y = cave.relax$relaxed_ratio, nsim = 1000)
print(phylANOVA_relaxed)
## PGLS regression
comp.data <- comparative.data(phy = cave.tre, data = cave.relax, names.col = "sp", vcv = TRUE, na.omit = FALSE)
pgls_model_relaxed <- pgls(relaxed_ratio ~ fb, data = comp.data)
summary(pgls_model_relaxed)
## Phylogenetic linear model (Brownian motion)
phylolm_model_relaxed <- phylolm(formula = relaxed_ratio ~ fb, data = cave.relax, phy = cave.tre, model = "BM")
summary(phylolm_model_relaxed)
## compare overlap between intensified selection and positive selection
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/")
intens <- read.table("cave_species_RELAX_annotated_intensified_genes.txt", header = TRUE)
str(intens)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/01_cave_species")
pos <- read.csv("cave_species_shared_sel_OGs.csv", header = TRUE)
str(pos)
names(pos)[names(pos) == 'OG'] <- 'OG_ID'
intens_pos <- merge(intens,pos, by="OG_ID")
#write.csv(intens_pos, "overlap_intesified_positive_selection_cave_species.csv", row.names=FALSE)


#### GENERATE CHROMOSOME FIGURE OF OLM HIGHLIGHTING GENES UNDER SELECTION IN CAVE SPECIES, GENE AND REPEAT DENSITY ####
# set working directory 
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/12_genome_location/01_proteus")
## load chromosome lengths 
chr_df <- read.table("Proteus_chromosome_lengths.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(chr_df) <- c("chr","length")
chr_df <- chr_df[grepl("^chr", chr_df$chr, ignore.case=TRUE), ]
chr_df$chr <- tolower(chr_df$chr)
chr_df <- chr_df[order(chr_df$length, decreasing=TRUE), ]
chr_df$chr_label <- gsub("chr","",chr_df$chr)  # Only numbers 1-19
### load selection datasets 
pos_sel <- read.csv("cave_species_paml_gen_loc.csv", stringsAsFactors=FALSE)
pos_sel <- pos_sel[grepl("^chr", pos_sel$chr, ignore.case=TRUE), ]
pos_sel$chr <- tolower(pos_sel$chr)
rel_sel <- read.csv("cave_species_relaxed_gen_loc.csv", stringsAsFactors=FALSE)
rel_sel <- rel_sel[grepl("^chr", rel_sel$chr, ignore.case=TRUE), ]
rel_sel$chr <- tolower(rel_sel$chr)
exp_sel <- read.csv("Prot_shared_cave_exp.csv", stringsAsFactors=FALSE)
exp_sel <- exp_sel[grepl("^chr", exp_sel$chr, ignore.case=TRUE), ]
exp_sel$chr <- tolower(exp_sel$chr)
### load gene annotation 
gene_anno <- read.table("gene_annotation_mRNA.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(gene_anno)[1] <- "chr"
gene_anno$chr <- tolower(gene_anno$chr)
colnames(gene_anno)[4:5] <- c("start","end")
### repeat annotation 
repeat_anno <- read.table("DY.HiC.repeat.gff", header=FALSE, stringsAsFactors=FALSE)
colnames(repeat_anno)[1] <- "chr"
repeat_anno$chr <- tolower(repeat_anno$chr)
colnames(repeat_anno)[4:5] <- c("start","end")
### function to calculate averaged density per Mb, averaged over avg_window (e.g. 20 Mb) 
calculate_density_per_mb_avg20 <- function(df, chr_lengths, base_window=1e6, avg_window=20e6) {
  density_list <- list()
  # number of base bins to average
  avg_window_bins <- max(1, round(avg_window / base_window))
  for (chr in chr_lengths$chr) {
    chr_len <- chr_lengths$length[chr_lengths$chr == chr]
    # base 1 Mb bins (ensure final bin ends exactly at chromosome length)
    bins <- seq(0, chr_len, by=base_window)
    if (tail(bins, 1) < chr_len) bins <- c(bins, chr_len)
    pos <- c(df$start[df$chr==chr], df$end[df$chr==chr])
    counts <- hist(pos, breaks=bins, plot=FALSE)$counts
    # convert counts to per-Mb rate (counts per 1 Mb)
    per_mb <- counts / (base_window / 1e6)
    # Average across avg_window (sliding moving average)
    if (length(per_mb) >= avg_window_bins && avg_window_bins > 1) {
      per_mb_avg_ts <- stats::filter(per_mb, rep(1/avg_window_bins, avg_window_bins), sides=2)
      per_mb_avg <- as.numeric(per_mb_avg_ts)           # <-- convert ts -> numeric
      per_mb_avg[is.na(per_mb_avg)] <- 0
    } else {
      per_mb_avg <- per_mb}
    density_list[[chr]] <- data.frame(
      chr = chr,
      start = bins[-length(bins)],
      end = bins[-1],
      density = per_mb_avg)}
  density_df <- do.call(rbind, density_list)
  # Global normalization across all chromosomes (0–1), handle constant case
  min_val <- min(density_df$density, na.rm = TRUE)
  max_val <- max(density_df$density, na.rm = TRUE)
  if (max_val - min_val > 0) {
    density_df$density_norm <- (density_df$density - min_val) / (max_val - min_val)
  } else {
    density_df$density_norm <- 0}
  density_df$density_norm <- pmin(pmax(density_df$density_norm, 0), 1)
  return(density_df)}
### compute densities (averaged over 20 Mb) 
gene_density <- calculate_density_per_mb_avg20(gene_anno, chr_df, base_window=1e6, avg_window=20e6)
repeat_density <- calculate_density_per_mb_avg20(repeat_anno, chr_df, base_window=1e6, avg_window=20e6)
### legend labels 
legend_pos <- unique(pos_sel$analysis_paml)
legend_rel <- unique(rel_sel$analysis_relax_relaxed)
legend_exp <- unique(exp_sel$qstart)
### softer alternating chromosome background 
chr_cols <- rep(c("#F2F2F2", "#E6E6E6"), length.out=nrow(chr_df))
### clear Circos 
circos.clear()
circos.par(start.degree=90, gap.degree=1, cell.padding=c(0,0,0,0))

#### FIGURE 3 ####
### initialize chromosomes 
circos.initialize(factors=chr_df$chr, xlim=cbind(0, chr_df$length))
### Ring 1: Positive selection 
circos.trackPlotRegion(ylim=c(0,1), track.height=0.08, bg.col=chr_cols, bg.border="black",
                       panel.fun=function(x,y){
                         sector <- get.cell.meta.data("sector.index")
                         r <- pos_sel[pos_sel$chr==sector, ]
                         if(nrow(r)>0) circos.segments(x0=r$seqStart, y0=0, x1=r$seqStart, y1=1,
                                                       col="#FF7F50", lwd=0.6)})
### Ring 2: Relaxed selection 
circos.trackPlotRegion(ylim=c(0,1), track.height=0.08, bg.col=chr_cols, bg.border="black",
                       panel.fun=function(x,y){
                         sector <- get.cell.meta.data("sector.index")
                         r <- rel_sel[rel_sel$chr==sector, ]
                         if(nrow(r)>0) circos.segments(x0=r$seqStart, y0=0, x1=r$seqStart, y1=1,
                                                       col="#20B2AA", lwd=0.6)})
### Ring 3: Shared cave expression 
circos.trackPlotRegion(ylim=c(0,1), track.height=0.08, bg.col=chr_cols, bg.border="black",
                       panel.fun=function(x,y){
                         sector <- get.cell.meta.data("sector.index")
                         r <- exp_sel[exp_sel$chr==sector, ]
                         if(nrow(r)>0) circos.segments(x0=r$seqStart, y0=0, x1=r$seqStart, y1=1,
                                                       col="#9370DB", lwd=0.6)})
### Ring 4: Gene density (second innermost) 
track_height <- 0.08
circos.trackPlotRegion(ylim=c(0,1), track.height=track_height, bg.border=NA, bg.col=NA,
                       panel.fun=function(x,y){
                         sector <- get.cell.meta.data("sector.index")
                         r <- gene_density[gene_density$chr==sector, ]
                         if(nrow(r)>0){
                           # pick points every 20 base bins (20 * 1Mb = 20Mb)
                           idx <- seq(1, nrow(r), by=20)
                           x_vals <- (r$start[idx] + r$end[idx]) / 2
                           y_vals <- r$density_norm[idx]
                           circos.lines(x=x_vals, y=y_vals, col="#FFD700", lwd=1.2)
                           circos.polygon(x=c(x_vals, rev(x_vals)), y=c(y_vals, rep(0, length(y_vals))),col=adjustcolor("#FFD700", alpha.f=0.4), border=NA)}
                         if(sector == chr_df$chr[1]){
                           circos.yaxis(side="left", at=seq(0,1,0.2), labels=seq(0,1,0.2),
                                        labels.cex=0.5, sector.index=sector, tick.length=0.01)}})
### Ring 5: Repeat density (innermost) 
circos.trackPlotRegion(ylim=c(0,1), track.height=track_height, bg.border=NA, bg.col=NA,
                       panel.fun=function(x,y){
                         sector <- get.cell.meta.data("sector.index")
                         r <- repeat_density[repeat_density$chr==sector, ]
                         if(nrow(r)>0){
                           idx <- seq(1, nrow(r), by=20)
                           x_vals <- (r$start[idx] + r$end[idx]) / 2
                           y_vals <- r$density_norm[idx]
                           circos.lines(x=x_vals, y=y_vals, col="#87CEEB", lwd=1.2)
                           circos.polygon(x=c(x_vals, rev(x_vals)), y=c(y_vals, rep(0, length(y_vals))), col=adjustcolor("#87CEEB", alpha.f=0.4), border=NA)}
                         if(sector == chr_df$chr[1]){
                           circos.yaxis(side="left", at=seq(0,1,0.2), labels=seq(0,1,0.2),
                                        labels.cex=0.5, sector.index=sector, tick.length=0.01)}})
### chromosome labels outside the outermost ring 
circos.track(ylim=c(0,1), track.height=0.05, bg.border=NA, panel.fun=function(x,y){
  sector <- get.cell.meta.data("sector.index")
  chr_number <- chr_df$chr_label[chr_df$chr == sector]
  circos.text(mean(get.cell.meta.data("xlim")), 1.15, labels=chr_number,
              facing="outside", cex=0.7, adj=c(0.5,0))})
### legend 
legend("topleft",
       legend=c("Positive selection", "Relaxed selection", "Shared cave expression",
                "Gene density", "Repeat density"),
       fill=c("#FF7F50","#20B2AA","#9370DB","#FFD700","#87CEEB"),
       border=NA, bty="n", cex=0.8)

### PIGMENTATION LOSS  ###
## read in annotation file
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/01_eggnog_annotation_files")
eggnog.pigm <- read.csv("all_OGs_pigmentation_loss_MM_fe2v4rjm.emapper.annotations.csv", header = TRUE)
str(eggnog.pigm)
# add a new column just with OG ID
eggnog.pigm <- eggnog.pigm %>% 
  separate(query,c("OG_ID", "sequence.ID"), sep = "-", remove = FALSE) %>%
  dplyr::select(OG_ID, Preferred_name,Description,GOs,KEGG_ko)
str(eggnog.pigm)
## read in selection file
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/")
relax.pigm <- read.table("pigmentation_loss_selected_OGs.txt", header=TRUE)
str(relax.pigm)
selected.OGs <- merge(relax.pigm, eggnog.pigm, by="OG_ID", all.x=TRUE)
str(selected.OGs)
intensified <- selected.OGs  %>% filter(selection == "intensified") %>% 
  dplyr::select(OG_ID,Preferred_name,selection) 
relaxed <- selected.OGs  %>% filter(selection == "relaxed") %>% 
  dplyr::select(OG_ID,Preferred_name,selection) 
#write.table(selected.OGs, "pigmentation_loss_RELAX_annotated_selected_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(intensified, "pigmentation_loss_RELAX_annotated_intensified_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(relaxed, "pigmentation_loss_RELAX_annotated_relaxed_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)

## compare overlap between intensified selection and positive selection
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/")
intens <- read.table("pigmentation_loss_RELAX_annotated_intensified_genes.txt", header = TRUE)
str(intens)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/02_pigmentation_loss")
pos <- read.csv("pigmentation_loss_shared_sel_OGs.csv", header = TRUE)
str(pos)
names(pos)[names(pos) == 'OG'] <- 'OG_ID'
intens_pos <- merge(intens,pos, by="OG_ID")

### EYE LOSS ###
## read in annotation file
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/01_eggnog_annotation_files")
eggnog.eyeloss <- read.csv("all_OGs_eye_loss_MM_fe2v4rjm.emapper.annotations.csv", header = TRUE)
str(eggnog.eyeloss)
# add a new column just with OG ID
eggnog.eyeloss <- eggnog.eyeloss %>% 
  separate(query,c("OG_ID", "sequence.ID"), sep = "-", remove = FALSE) %>%
  dplyr::select(OG_ID, Preferred_name,Description,GOs,KEGG_ko)
str(eggnog.eyeloss)
## read in selection file
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/")
relax.eyeloss <- read.table("eye_loss_selected_OGs.txt", header=TRUE)
str(relax.eyeloss)
selected.OGs <- merge(relax.eyeloss, eggnog.eyeloss, by="OG_ID", all.x=TRUE)
str(selected.OGs)
intensified <- selected.OGs  %>% filter(selection == "intensified") %>% 
  dplyr::select(OG_ID,Preferred_name,selection) 
relaxed <- selected.OGs  %>% filter(selection == "relaxed") %>% 
  dplyr::select(OG_ID,Preferred_name,selection) 
#write.table(selected.OGs, "eye_loss_RELAX_annotated_selected_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(intensified, "eye_loss_RELAX_annotated_intensified_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(relaxed, "eye_loss_RELAX_annotated_relaxed_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
## compare overlap between intensified selection and positive selection
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/")
intens <- read.table("eye_loss_RELAX_annotated_intensified_genes.txt", header = TRUE)
str(intens)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/03_eye_loss")
pos <- read.csv("eye_loss_shared_sel_OGs.csv", header = TRUE)
str(pos)
names(pos)[names(pos) == 'OG'] <- 'OG_ID'
intens_pos <- merge(intens,pos, by="OG_ID")

### LONGEVITY ###
## read in annotation file
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/01_eggnog_annotation_files")
eggnog.longev <- read.csv("all_OGs_longevity_MM_j29p4hp0.emapper.annotations.csv", header = TRUE)
str(eggnog.longev)
# add a new column just with OG ID
eggnog.longev <- eggnog.longev %>% 
  separate(query,c("OG_ID", "sequence.ID"), sep = "-", remove = FALSE) %>%
  dplyr::select(OG_ID, Preferred_name,Description,GOs,KEGG_ko)
str(eggnog.longev)
## read in selection file
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/")
relax.longev <- read.table("longevity_selected_OGs.txt", header=TRUE)
str(relax.longev)
selected.OGs <- merge(relax.longev, eggnog.longev, by="OG_ID", all.x=TRUE)
str(selected.OGs)
intensified <- selected.OGs  %>% filter(selection == "intensified") %>% 
  dplyr::select(OG_ID,Preferred_name,selection) 
relaxed <- selected.OGs  %>% filter(selection == "relaxed") %>% 
  dplyr::select(OG_ID,Preferred_name,selection) 
#write.table(selected.OGs, "longevity_RELAX_annotated_selected_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(intensified, "longevity_RELAX_annotated_intensified_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(relaxed, "longevity_RELAX_annotated_relaxed_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)

## compare overlap between intensified selection and positive selection
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/")
intens <- read.table("longevity_RELAX_annotated_intensified_genes.txt", header = TRUE)
str(intens)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/09_absrel/04_longevity")
pos <- read.csv("longevity_shared_sel_OGs.csv", header = TRUE)
str(pos)
names(pos)[names(pos) == 'OG'] <- 'OG_ID'
intens_pos <- merge(intens,pos, by="OG_ID")


######### GENOMIC LOCATION ANALYSIS ######### 
### CAVE SPECIES ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/12_genome_location")
### Test for co-localisation between relaxed and positive selection, but exclude overlaps
### load gene lists (Example data) and the annotated background genes
cave.genloc <- read.csv("cave_species_candidate_gene_selection_analysis.csv") 
str(cave.genloc)
dim(cave.genloc)
unique.cave.loc <- cave.genloc %>% distinct(OG_ID, .keep_all = TRUE)
dim(unique.cave.loc)
unique.cave.loc <- unique.cave.loc %>% drop_na(query)
str(unique.cave.loc)
#write.csv(unique.cave.loc,"cave_species_selection_genomic_location.csv",quote=FALSE, row.names=FALSE)
## drop scaffolds from selection files
unique.cave.loc <- unique.cave.loc %>%
  filter(!grepl("^scaffold", chr))
str(unique.cave.loc)
## check that all scaffolds are excluded
categories <- unique(unique.cave.loc$chr) 
numberOfCategories <- length(categories) 
### create different dataframes for each selection score
absrel <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_absrel) %>% filter(analysis_absrel != "")
colnames(absrel) <- c("Chrom", "Start", "End", "Name")
paml <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_PAML) %>% filter(analysis_PAML != "")
colnames(paml) <- c("Chrom", "Start", "End", "Name")
intens <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_relax_intens) %>% filter(analysis_relax_intens != "")
colnames(intens) <- c("Chrom", "Start", "End", "Name")
expan <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_cafe_expanded) %>% filter(analysis_cafe_expanded != "")
colnames(expan) <- c("Chrom", "Start", "End", "Name")
contr <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_cafe_contr) %>% filter(analysis_cafe_contr != "")
colnames(contr) <- c("Chrom", "Start", "End", "Name")
relax <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_relax_relaxed) %>% filter(analysis_relax_relaxed != "")
colnames(relax) <- c("Chrom", "Start", "End", "Name")
pigm <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, condition_pigmentation) %>% filter(condition_pigmentation != "")
colnames(pigm) <- c("Chrom", "Start", "End", "Name")
eyeloss <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, condition_eyes) %>% filter(condition_eyes != "")
colnames(eyeloss) <- c("Chrom", "Start", "End", "Name")
## extract regions that contain genes
gene_annot <- read.csv("cave_species_6064_genes_coord.csv") 
str(gene_annot)
## remove NAs
gene_annot <- gene_annot %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))  
dim(paml)
paml <- paml %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))
dim(absrel)
absrel <- absrel %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))
dim(intens)
intens <- intens %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))
dim(relax)
relax <- relax %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))
### create GenomicRanges objects for each file
# Convert the data frames into GRanges objects for efficient interval overlap operations
gene_annot_gr <- GRanges(seqnames = gene_annot$Chrom, ranges = IRanges(start = gene_annot$Start, end = gene_annot$End))
paml_gr <- GRanges(seqnames = paml$Chrom, ranges = IRanges(start = paml$Start, end = paml$End))
absrel_gr <- GRanges(seqnames = absrel$Chrom, ranges = IRanges(start = absrel$Start, end = absrel$End))
intens_gr <- GRanges(seqnames = intens$Chrom, ranges = IRanges(start = intens$Start, end = intens$End))
relax_gr <- GRanges(seqnames = relax$Chrom, ranges = IRanges(start = relax$Start, end = relax$End))
# make a function that calculates a permutation test for two gene lists and checks if they are closer than expected by chance
perform_permutation_test <- function(gene_list1, gene_list2, gene_annot, threshold = 100000, n_permutations = 1000) {
  # Extract common chromosomes
  common_chromosomes <- intersect(as.vector(seqnames(gene_list1)), as.vector(seqnames(gene_list2)))
  common_chromosomes <- unique(as.character(common_chromosomes))
  # Calculate observed pairwise distances
  observed_distances <- c()
  for (chrom in common_chromosomes) {
    # Subset genes on the same chromosome
    list1_chr <- gene_list1[seqnames(gene_list1) == chrom]
    list2_chr <- gene_list2[seqnames(gene_list2) == chrom]
    # Calculate all pairwise distances
    pairwise_distances <- distance(list1_chr, list2_chr)
    observed_distances <- c(observed_distances, as.vector(pairwise_distances))} 
  mean_observed_distance <- mean(observed_distances)
  # Perform permutation test
  random_distances <- numeric(n_permutations)
  for (i in 1:n_permutations) {
    random_all_distances <- c()
    for (chrom in common_chromosomes) {
      # Subset the genome annotation and gene lists by chromosome
      gene_annot_chr <- gene_annot[seqnames(gene_annot) == chrom]
      list1_chr <- gene_list1[seqnames(gene_list1) == chrom]
      list2_chr <- gene_list2[seqnames(gene_list2) == chrom]
      # Ensure that there are enough genes on the chromosome to sample
      if (length(gene_annot_chr) >= length(list2_chr)) {
        # Shuffle gene list 2 by sampling indices from the gene annotation on the same chromosome
        shuffled_indices <- sample(seq_along(gene_annot_chr), length(list2_chr), replace = FALSE)
        shuffled_list2_chr <- gene_annot_chr[shuffled_indices]
      } else {
        # If not enough genes are available, sample with replacement
        shuffled_indices <- sample(seq_along(gene_annot_chr), length(list2_chr), replace = TRUE)
        shuffled_list2_chr <- gene_annot_chr[shuffled_indices]}
      # Calculate pairwise distances between list1_chr and shuffled_list2_chr
      pairwise_random_distances <- distance(list1_chr, shuffled_list2_chr)
      random_all_distances <- c(random_all_distances, as.vector(pairwise_random_distances))}
    # Store the mean of the pairwise distances for this permutation
    random_distances[i] <- mean(random_all_distances)}
  # Calculate p-value
  p_value <- mean(random_distances <= mean_observed_distance)
  # Return results
  return(list(
    mean_observed_distance = mean_observed_distance,
    mean_random_distance = mean(random_distances),
    p_value = p_value))}
## do all pairwise comparisons
gene_list1 <- paml_gr
gene_list2 <- absrel_gr
gene_list3 <- intens_gr  
gene_list4 <- relax_gr
## store the lists in a named list for easy access
gene_lists <- list(
  "PAML" = gene_list1,
  "ABSREL" = gene_list2,
  "INTENS" = gene_list3,
  "RELAX" = gene_list4)
## get all pairwise combinations of the lists
pairwise_combinations <- combn(names(gene_lists), 2, simplify = FALSE)
## run the permutation test for each pair
results <- list()
for (pair in pairwise_combinations) {
  list1_name <- pair[1]
  list2_name <- pair[2]
  # Get the actual GRanges objects for the pair
  list1 <- gene_lists[[list1_name]]
  list2 <- gene_lists[[list2_name]]
  # Perform the permutation test
  result <- perform_permutation_test(list1, list2, gene_annot_gr)
  # Store the result with the pair names
  results[[paste(list1_name, "vs", list2_name)]] <- result}


### Test for clustering of positively selected or relaxed genes
# Calculate pairwise distances using the midpoint of each gene
# Function to calculate pairwise distances by chromosome using midpoints
# Define the main function that will run the analysis
run_clustering_analysis <- function(selection_gr, background_gr, n_permutations = 1000) {
  ### calculate observed distances for the gene list
  observed_distances <- calculate_pairwise_midpoint_distances(selection_gr)
  ### perform permutation test to compare with random expectation
  set.seed(123)  # For reproducibility
  random_distances <- vector("numeric", length = n_permutations)
  for (i in 1:n_permutations) {
    # Randomly sample the same number of genes from the background
    random_genes <- sample(background_gr, length(selection_gr), replace = FALSE)
    # Calculate pairwise distances for the randomly sampled genes
    random_distances[i] <- mean(calculate_pairwise_midpoint_distances(random_genes))}
  ### compare observed vs random
  mean_observed_distance <- mean(observed_distances)
  mean_random_distance <- mean(random_distances)
  ### p-value calculation (how often observed distance is smaller than random)
  p_value <- mean(random_distances >= mean_observed_distance)
  ### return the results as a list
  return(list(
    mean_observed_distance = mean_observed_distance,
    mean_random_distance = mean_random_distance,
    p_value = p_value,
    random_distances = random_distances,
    observed_distances = observed_distances))}
# Function to output the results in a readable format
print_clustering_results <- function(results, analysis_name) {
  cat("\nResults for", analysis_name, "\n")
  cat("Mean observed distance between genes: ", results$mean_observed_distance, "\n")
  cat("Mean random distance between genes: ", results$mean_random_distance, "\n")
  cat("P-value: ", results$p_value, "\n")}
### run the analysis for each selection file and store the results
paml_results <- run_clustering_analysis(paml_gr, gene_annot_gr)
relax_results <- run_clustering_analysis(relax_gr, gene_annot_gr)
intens_results <- run_clustering_analysis(intens_gr, gene_annot_gr)
absrel_results <- run_clustering_analysis(absrel_gr, gene_annot_gr)
### print the results for each analysis
print_clustering_results(paml_results, "PAML")
print_clustering_results(relax_results, "RELAX")
print_clustering_results(intens_results, "INTENS")
print_clustering_results(absrel_results, "ABSREL")
### optionally visualize the distribution of random distances for each analysis
# You can modify the layout of the plots if needed to fit multiple histograms in one output
par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
# PAML plot
hist(paml_results$random_distances, breaks = 30, main = "PAML: Distribution of Random Distances", xlab = "Average Distance (bp)")
abline(v = paml_results$mean_observed_distance, col = "red", lwd = 2, lty = 2)
# RELAX plot
hist(relax_results$random_distances, breaks = 30, main = "RELAX: Distribution of Random Distances", xlab = "Average Distance (bp)")
abline(v = relax_results$mean_observed_distance, col = "red", lwd = 2, lty = 2)
# INTENS plot
hist(intens_results$random_distances, breaks = 30, main = "INTENS: Distribution of Random Distances", xlab = "Average Distance (bp)")
abline(v = intens_results$mean_observed_distance, col = "red", lwd = 2, lty = 2)
# ABSREL plot
hist(absrel_results$random_distances, breaks = 30, main = "ABSREL: Distribution of Random Distances", xlab = "Average Distance (bp)")
abline(v = absrel_results$mean_observed_distance, col = "red", lwd = 2, lty = 2)
# Reset layout
par(mfrow = c(1, 1))


#### MAKE CIRCOS PLOT OF GENOME CO-LOCALIZATION ACROSS SPECIES ####
# merge OG ID candidate genes and localisation on genome
# Read the CSV files
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/12_genome_location/01_proteus")
cave_genes <- read.csv("cave_species_6064_genes.csv", stringsAsFactors = FALSE)
supplementary <- read.csv("Supplementary_Data_X_cave_species_OG_expansion_contraction_selection_genes.csv", stringsAsFactors = FALSE)
# Perform left join based on OG_ID
merged_data <- left_join(supplementary, cave_genes, by = "OG_ID")
# Write the result to a new CSV
#write.csv(merged_data, "merged_cave_species_data.csv", row.names = FALSE)
# Set working directory for genome assembly
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/04_genome_assembly/Prang_Hi-c")
### load chromosome lengths
chrom_sizes <- read.table("Proteus_chromosome_lengths.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("Chrom", "Length")
### keep only chromosomes (exclude scaffolds) and add Start/End columns
chromosomes <- chrom_sizes %>% filter(grepl("^chr", Chrom)) %>% mutate(Start = 1, End = Length) %>% dplyr::select(Chrom, Start, End)
### read gene location files
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/12_genome_location/01_proteus")
positive_genes <- read.csv("cave_species_paml_gen_loc.csv", stringsAsFactors = FALSE)
relaxed_genes  <- read.csv("cave_species_relaxed_gen_loc.csv", stringsAsFactors = FALSE)
### filter genes to only main chromosomes
positive_genes <- positive_genes[positive_genes$chr %in% chromosomes$Chrom, ]
relaxed_genes  <- relaxed_genes[relaxed_genes$chr %in% chromosomes$Chrom, ]
### initialize circos plot
circos.clear()
circos.par(start.degree = 90, gap.after = rep(5, nrow(chromosomes)))
circos.initialize(factors = chromosomes$Chrom, xlim = cbind(chromosomes$Start, chromosomes$End))
### optional: chromosome ideogram
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) { chr = CELL_META$sector.index; circos.text(CELL_META$xcenter, 0.5, chr, facing = "clockwise", niceFacing = TRUE, cex = 0.6) }, bg.col = "grey90", track.height = 0.05)
### inner track: relaxed selection
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.1, bg.border = NA, panel.fun = function(x, y) { chr = CELL_META$sector.index; chr_genes <- relaxed_genes[relaxed_genes$chr == chr, ]; if (nrow(chr_genes) > 0) { circos.rect(chr_genes$seqStart, 0, chr_genes$seqEnd, 1, col = "skyblue", border = NA) } })
### outer track: positive selection
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.1, bg.border = NA, panel.fun = function(x, y) { chr = CELL_META$sector.index; chr_genes <- positive_genes[positive_genes$chr == chr, ]; if (nrow(chr_genes) > 0) { circos.rect(chr_genes$seqStart, 0, chr_genes$seqEnd, 1, col = "salmon", border = NA) } })
### add legend
legend("topleft", legend = c("Relaxed selection", "Positive selection"), fill = c("skyblue", "salmon"), border = NA, bty = "n")


### PIGMENTATION LOSS ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/12_genome_location")
cave.genloc <- read.csv("cave_species_candidate_gene_selection_analysis.csv") 
str(cave.genloc)
dim(cave.genloc)
unique.cave.loc <- cave.genloc %>% distinct(OG_ID, .keep_all = TRUE)
dim(unique.cave.loc)
unique.cave.loc <- unique.cave.loc %>% drop_na(query)
str(unique.cave.loc)
#write.csv(unique.cave.loc,"cave_species_selection_genomic_location.csv",quote=FALSE, row.names=FALSE)
### extract chromosome size
chromosomes <- unique.cave.loc %>%
  group_by(chr) %>%
  summarize(seqStart = min(seqStart), seqEnd = max(seqEnd))
# change to data frame and add annotation
chromosomes <- as.data.frame(chromosomes)
chromosomes$annotation <- "chromosome"
print(chromosomes)
colnames(chromosomes) <- c("Chrom", "Start", "End", "Name")
chromosomes <- chromosomes %>%
  filter(!grepl("^scaffold", Chrom))
print(chromosomes)
str(chromosomes)
chromosome_lengths <- chromosomes %>% dplyr::select(Chrom,End)
colnames(chromosome_lengths) <-c("Chrom", "Length")
### create different dataframes for each selection score
absrel <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_absrel) %>% filter(analysis_absrel != "")
colnames(absrel) <- c("Chrom", "Start", "End", "Name")
paml <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_PAML) %>% filter(analysis_PAML != "")
colnames(paml) <- c("Chrom", "Start", "End", "Name")
intens <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_relax_intens) %>% filter(analysis_relax_intens != "")
colnames(intens) <- c("Chrom", "Start", "End", "Name")
expan <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_cafe_expanded) %>% filter(analysis_cafe_expanded != "")
colnames(expan) <- c("Chrom", "Start", "End", "Name")
contr <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_cafe_contr) %>% filter(analysis_cafe_contr != "")
colnames(contr) <- c("Chrom", "Start", "End", "Name")
relax <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_relax_relaxed) %>% filter(analysis_relax_relaxed != "")
colnames(relax) <- c("Chrom", "Start", "End", "Name")
pigm <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, condition_pigmentation) %>% filter(condition_pigmentation != "")
colnames(pigm) <- c("Chrom", "Start", "End", "Name")
eyeloss <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, condition_eyes) %>% filter(condition_eyes != "")
colnames(eyeloss) <- c("Chrom", "Start", "End", "Name")
## combine and rearrange data
# Merge multiple columns into a single column and drop other column
cave.chroms.sel <- unique.cave.loc %>%
  unite("annotation", analysis_absrel, analysis_relax_intens, analysis_PAML,analysis_cafe_expanded, na.rm = TRUE, sep = "_", remove = TRUE) %>%
  dplyr::select(-analysis_relax_relaxed,-analysis_cafe_contr,-condition_pigmentation,-condition_eyes,-gene_symbol,-qname, -query, -OG_ID)
print(cave.chroms.sel)
# drop NAs and change order of variables
cave.chroms.sel <- cave.chroms.sel %>%
  filter(annotation != "") %>%
  dplyr::select(-annotation,annotation)
print(cave.chroms.sel)
colnames(cave.chroms.sel) <- c("Chrom", "Start", "End", "Name")
### combine both data frames
cave.chroms.sel.comb <- full_join(chromosomes, cave.chroms.sel, by = c("Chrom", "Start", "End", "Name"))
print(cave.chroms.sel.comb)
colnames(cave.chroms.sel.comb) <- c("Chrom", "Start", "End", "Name")
### use chromoplot package for visualization
### start with plotting all genes under selection (absrel, intens, paml)
chromPlot(gaps=chromosomes, annot1 = absrel, figCols=19)
### create different plots for all selection scores
chromPlot(gaps=chromosomes, annot1 = absrel, annot2= intens, annot3= paml, figCols=19)
chromPlot(gaps=chromosomes, annot1 = relax, figCols=19)
### extract regions that contain genes
gene_annot <- read.csv("proteus_annotation_hic_cds.csv") 
str(gene_annot)
gene_annot <- gene_annot %>%
  filter(!grepl("^scaffold", Chrom))
## remove NAs
gene_annot <- gene_annot %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))  
paml <- paml %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))
intens <- intens %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))
dim(intens)
### create GenomicRanges objects for each file
# Convert the data frames into GRanges objects for efficient interval overlap operations
gene_annot_gr <- GRanges(seqnames = gene_annot$Chrom, ranges = IRanges(start = gene_annot$Start, end = gene_annot$End))
paml_gr <- GRanges(seqnames = paml$Chrom, ranges = IRanges(start = paml$Start, end = paml$End))
intens_gr <- GRanges(seqnames = intens$Chrom, ranges = IRanges(start = intens$Start, end = intens$End))
### restrict paml and intens to genic regions (filter by genic regions)
paml_genic <- subsetByOverlaps(paml_gr, gene_annot_gr)
intens_genic <- subsetByOverlaps(intens_gr, gene_annot_gr)
### find the real overlaps between the paml and intens regions within genic regions
real_overlaps <- findOverlaps(paml_genic, intens_genic)
real_overlap_count <- length(unique(queryHits(real_overlaps)))
cat("Number of real overlaps:", real_overlap_count, "\n")
### perform randomization test (Monte Carlo simulation)
set.seed(123)  # For reproducibility
# function to randomize the gene regions within the genic regions
randomize_genes <- function(selection_gr, gene_annot_gr) {
  randomized_gr <- GRanges()
  for (i in seq_along(selection_gr)) {
    # Get the length of the gene region
    gene_length <- width(selection_gr[i])
    # Select a random genic region
    random_gene_region <- sample(gene_annot_gr, 1)
    # Randomly place the gene within the selected genic region
    start_pos <- start(random_gene_region) + sample(0:(width(random_gene_region) - gene_length), 1)
    end_pos <- start_pos + gene_length - 1
    # Create a new randomized gene region
    new_region <- GRanges(seqnames = seqnames(random_gene_region),
                          ranges = IRanges(start = start_pos, end = end_pos))
    randomized_gr <- c(randomized_gr, new_region)}
  return(randomized_gr)}
# number of iterations for randomization test
n_iterations <- 1000
random_overlaps <- numeric(n_iterations)
for (i in 1:n_iterations) {
  # randomize paml genes within genic regions
  paml_randomized <- randomize_genes(paml_genic, gene_annot_gr)
  # calculate overlap with intens genic regions
  random_overlap <- findOverlaps(paml_randomized, intens_genic)
  random_overlaps[i] <- length(unique(queryHits(random_overlap)))}
### compare real overlap with random overlaps
# plot histogram of random overlaps and show where the real overlap lies
hist(random_overlaps, breaks = 10, main = "Randomized Overlap Counts", xlab = "Number of Overlaps")
abline(v = real_overlap_count, col = "red", lwd = 2, lty = 2)
# calculate p-value
p_value <- mean(random_overlaps >= real_overlap_count)
cat("P-value for co-localization:", p_value, "\n")


### EYE LOSS ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/12_genome_location")
cave.genloc <- read.csv("cave_species_candidate_gene_selection_analysis.csv") 
str(cave.genloc)
dim(cave.genloc)
unique.cave.loc <- cave.genloc %>% distinct(OG_ID, .keep_all = TRUE)
dim(unique.cave.loc)
unique.cave.loc <- unique.cave.loc %>% drop_na(query)
str(unique.cave.loc)
#write.csv(unique.cave.loc,"cave_species_selection_genomic_location.csv",quote=FALSE, row.names=FALSE)
### extract chromosome size
chromosomes <- unique.cave.loc %>%
  group_by(chr) %>%
  summarize(
    seqStart = min(seqStart),
    seqEnd = max(seqEnd))
# change to data frame and add annotation
chromosomes <- as.data.frame(chromosomes)
chromosomes$annotation <- "chromosome"
print(chromosomes)
colnames(chromosomes) <- c("Chrom", "Start", "End", "Name")
chromosomes <- chromosomes %>%
  filter(!grepl("^scaffold", Chrom))
print(chromosomes)
str(chromosomes)
chromosome_lengths <- chromosomes %>% dplyr::select(Chrom,End)
colnames(chromosome_lengths) <-c("Chrom", "Length")
### create different dataframes for each selection score
absrel <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_absrel) %>% filter(analysis_absrel != "")
colnames(absrel) <- c("Chrom", "Start", "End", "Name")
paml <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_PAML) %>% filter(analysis_PAML != "")
colnames(paml) <- c("Chrom", "Start", "End", "Name")
intens <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_relax_intens) %>% filter(analysis_relax_intens != "")
colnames(intens) <- c("Chrom", "Start", "End", "Name")
expan <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_cafe_expanded) %>% filter(analysis_cafe_expanded != "")
colnames(expan) <- c("Chrom", "Start", "End", "Name")
contr <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_cafe_contr) %>% filter(analysis_cafe_contr != "")
colnames(contr) <- c("Chrom", "Start", "End", "Name")
relax <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_relax_relaxed) %>% filter(analysis_relax_relaxed != "")
colnames(relax) <- c("Chrom", "Start", "End", "Name")
pigm <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, condition_pigmentation) %>% filter(condition_pigmentation != "")
colnames(pigm) <- c("Chrom", "Start", "End", "Name")
eyeloss <- unique.cave.loc %>% dplyr::select(chr, seqStart, seqEnd, condition_eyes) %>% filter(condition_eyes != "")
colnames(eyeloss) <- c("Chrom", "Start", "End", "Name")

## combine and rearrange data
# Merge multiple columns into a single column and drop other column
cave.chroms.sel <- unique.cave.loc %>%
  unite("annotation", analysis_absrel, analysis_relax_intens, analysis_PAML,analysis_cafe_expanded, na.rm = TRUE, sep = "_", remove = TRUE) %>%
  dplyr::select(-analysis_relax_relaxed,-analysis_cafe_contr,-condition_pigmentation,-condition_eyes,-gene_symbol,-qname, -query, -OG_ID)
print(cave.chroms.sel)
# drop NAs and change order of variables
cave.chroms.sel <- cave.chroms.sel %>%
  filter(annotation != "") %>%
  dplyr::select(-annotation,annotation)
print(cave.chroms.sel)
colnames(cave.chroms.sel) <- c("Chrom", "Start", "End", "Name")
### combine both data frames
cave.chroms.sel.comb <- full_join(chromosomes, cave.chroms.sel, by = c("Chrom", "Start", "End", "Name"))
print(cave.chroms.sel.comb)
colnames(cave.chroms.sel.comb) <- c("Chrom", "Start", "End", "Name")
### use chromoplot package for visualization
### start with plotting all genes under selection (absrel, intens, paml)
chromPlot(gaps=chromosomes, annot1 = absrel, figCols=19)
### create different plots for all selection scores
chromPlot(gaps=chromosomes, annot1 = absrel, annot2= intens, annot3= paml, figCols=19)
chromPlot(gaps=chromosomes, annot1 = relax, figCols=19)

### statistically test for co-localization
# Define a function to check overlap
check_overlap <- function(regions1, regions2) {
  overlaps <- 0
  for (i in 1:nrow(regions1)) {
    for (j in 1:nrow(regions2)) {
      if (regions1$Chrom[i] == regions2$Chrom[j] &&
          regions1$End[i] >= regions2$Start[j] &&
          regions1$Start[i] <= regions2$End[j]) {
        overlaps <- overlaps + 1}}}
  return(overlaps)}

# Calculate observed co-localization
observed_co_localization <- check_overlap(paml, intens)
# Function to shuffle QTL regions within their chromosomes
shuffle_within_chromosomes <- function(qtl_data, chromosome_lengths) {
  shuffled_data <- qtl_data
  for (chr in unique(qtl_data$Chrom)) {
    chr_length <- chromosome_lengths$Length[chromosome_lengths$Chrom == chr]
    chr_indices <- which(qtl_data$Chrom == chr)
    shuffled_starts <- sample(qtl_data$Start[chr_indices])
    shuffled_ends <- shuffled_starts + (qtl_data$End[chr_indices] - qtl_data$Start[chr_indices])
    shuffled_data$Start[chr_indices] <- shuffled_starts
    shuffled_data$End[chr_indices] <- shuffled_ends}
  return(shuffled_data)}
# Perform permutation test
set.seed(123)
n_permutations <- 1000
permuted_co_localizations <- numeric(n_permutations)
for (i in 1:n_permutations) {
  permuted_qtl_trait2 <- shuffle_within_chromosomes(intens, chromosome_lengths)
  permuted_co_localizations[i] <- check_overlap(paml, permuted_qtl_trait2)}

# Calculate p-value
p_value <- mean(permuted_co_localizations >= observed_co_localization)
print(paste("Observed co-localization:", observed_co_localization))
print(paste("P-value:", p_value))

### LONGEVITY ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/12_genome_location")
### Step 1: Load gene lists (Example data) and the annotated background genes
## merge longevity candidate genes from literature with results
long.genloc <- read.csv("longevity_candidate_gene_selection_analysis.csv") 
lit.genes <- read.csv("longevity_lit_genes_OGs.csv", header=TRUE)
str(long.genloc)
str(lit.genes)
dim(long.genloc)
long.genloc <- merge(long.genloc, lit.genes, by = "OG_ID", all.x = TRUE)
str(long.genloc)
dim(long.genloc)
unique.long.loc <- long.genloc %>% distinct(OG_ID, .keep_all = TRUE)
dim(unique.long.loc)
unique.long.loc <- unique.long.loc %>% drop_na(query)
str(unique.long.loc)
#write.csv(unique.long.loc,"cave_species_selection_genomic_location.csv",quote=FALSE, row.names=FALSE)
## drop scaffolds from selection files
unique.long.loc <- unique.long.loc %>%
filter(!grepl("^scaffold", chr))
str(unique.long.loc)
## check that all scaffolds are excluded
categories <- unique(unique.long.loc$chr) 
numberOfCategories <- length(categories) 
### create different dataframes for each selection score
absrel <- unique.long.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_absrel) %>% filter(analysis_absrel != "")
colnames(absrel) <- c("Chrom", "Start", "End", "Name")
paml <- unique.long.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_paml) %>% filter(analysis_paml != "")
colnames(paml) <- c("Chrom", "Start", "End", "Name")
intens <- unique.long.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_relaxed_intensified) %>% filter(analysis_relaxed_intensified != "")
colnames(intens) <- c("Chrom", "Start", "End", "Name")
expan <- unique.long.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_cafe_expanded) %>% filter(analysis_cafe_expanded != "")
colnames(expan) <- c("Chrom", "Start", "End", "Name")
contr <- unique.long.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_cafe_contracted) %>% filter(analysis_cafe_contracted != "")
colnames(contr) <- c("Chrom", "Start", "End", "Name")
relax <- unique.long.loc %>% dplyr::select(chr, seqStart, seqEnd, analysis_relax_relaxed) %>% filter(analysis_relax_relaxed != "")
colnames(relax) <- c("Chrom", "Start", "End", "Name")
long <- unique.long.loc %>% dplyr::select(chr, seqStart, seqEnd, condition) %>% filter(condition != "")
colnames(pigm) <- c("Chrom", "Start", "End", "Name")
## extract regions that contain genes
gene_annot <- read.csv("longevity_5736_genes_coord.csv") 
str(gene_annot)
## remove NAs
gene_annot <- gene_annot %>%
  mutate(Start = as.integer(Start),
         End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))  
dim(paml)
paml <- paml %>%
  mutate(Start = as.integer(Start),
         End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))
dim(absrel)
absrel <- absrel %>%
  mutate(Start = as.integer(Start),
         End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))
dim(intens)
intens <- intens %>%
  mutate(Start = as.integer(Start),
         End = as.integer(End)) %>%
  filter(!is.na(Start) & !is.na(End))

# Step 2: Create GenomicRanges objects for each file
# Convert the data frames into GRanges objects for efficient interval overlap operations
gene_annot_gr <- GRanges(seqnames = gene_annot$Chrom,
                         ranges = IRanges(start = gene_annot$Start, end = gene_annot$End))
paml_gr <- GRanges(seqnames = paml$Chrom,
                   ranges = IRanges(start = paml$Start, end = paml$End))
absrel_gr <- GRanges(seqnames = absrel$Chrom,
                     ranges = IRanges(start = absrel$Start, end = absrel$End))
intens_gr <- GRanges(seqnames = intens$Chrom,
                     ranges = IRanges(start = intens$Start, end = intens$End))

# make a function that calculates a permutation test for two gene lists and checks if they are closer than expected by chance
perform_permutation_test <- function(gene_list1, gene_list2, gene_annot, threshold = 100000, n_permutations = 1000) {
  # Extract common chromosomes
  common_chromosomes <- intersect(
    as.vector(seqnames(gene_list1)),
    as.vector(seqnames(gene_list2)))
  common_chromosomes <- unique(as.character(common_chromosomes))
  # Calculate observed pairwise distances
  observed_distances <- c()
  for (chrom in common_chromosomes) {
    # Subset genes on the same chromosome
    list1_chr <- gene_list1[seqnames(gene_list1) == chrom]
    list2_chr <- gene_list2[seqnames(gene_list2) == chrom]
    # Calculate all pairwise distances
    pairwise_distances <- distance(list1_chr, list2_chr)
    observed_distances <- c(observed_distances, as.vector(pairwise_distances))  # Store distances
  }
  mean_observed_distance <- mean(observed_distances)
  # Perform permutation test
  random_distances <- numeric(n_permutations)
  for (i in 1:n_permutations) {
    random_all_distances <- c()
    for (chrom in common_chromosomes) {
      # Subset the genome annotation and gene lists by chromosome
      gene_annot_chr <- gene_annot[seqnames(gene_annot) == chrom]
      list1_chr <- gene_list1[seqnames(gene_list1) == chrom]
      list2_chr <- gene_list2[seqnames(gene_list2) == chrom]
      # Ensure that there are enough genes on the chromosome to sample
      if (length(gene_annot_chr) >= length(list2_chr)) {
        # Shuffle gene list 2 by sampling indices from the gene annotation on the same chromosome
        shuffled_indices <- sample(seq_along(gene_annot_chr), length(list2_chr), replace = FALSE)
        shuffled_list2_chr <- gene_annot_chr[shuffled_indices]
      } else {
        # If not enough genes are available, sample with replacement
        shuffled_indices <- sample(seq_along(gene_annot_chr), length(list2_chr), replace = TRUE)
        shuffled_list2_chr <- gene_annot_chr[shuffled_indices]}
      # Calculate pairwise distances between list1_chr and shuffled_list2_chr
      pairwise_random_distances <- distance(list1_chr, shuffled_list2_chr)
      random_all_distances <- c(random_all_distances, as.vector(pairwise_random_distances))}
    # Store the mean of the pairwise distances for this permutation
    random_distances[i] <- mean(random_all_distances)}
  # Calculate p-value
  p_value <- mean(random_distances <= mean_observed_distance)
  # Return results
  return(list(
    mean_observed_distance = mean_observed_distance,
    mean_random_distance = mean(random_distances),
    p_value = p_value))}
# Do all pairwise comparisons
gene_list1 <- paml_gr
gene_list2 <- absrel_gr
gene_list3 <- intens_gr  
# Store the lists in a named list for easy access
gene_lists <- list(
  "PAML" = gene_list1,
  "ABSREL" = gene_list2,
  "INTENS" = gene_list3)
# Get all pairwise combinations of the lists
pairwise_combinations <- combn(names(gene_lists), 2, simplify = FALSE)
# Run the permutation test for each pair
results <- list()
for (pair in pairwise_combinations) {
  list1_name <- pair[1]
  list2_name <- pair[2]
  # Get the actual GRanges objects for the pair
  list1 <- gene_lists[[list1_name]]
  list2 <- gene_lists[[list2_name]]
  # Perform the permutation test
  result <- perform_permutation_test(list1, list2, gene_annot_gr)
  # Store the result with the pair names
  results[[paste(list1_name, "vs", list2_name)]] <- result}

### Test for clustering of positively selected or relaxed genes
# Calculate pairwise distances using the midpoint of each gene
# Function to calculate pairwise distances by chromosome using midpoints
# Define the main function that will run the analysis
run_clustering_analysis <- function(selection_gr, background_gr, n_permutations = 1000) {
  # Step 1: Calculate observed distances for the gene list
  observed_distances <- calculate_pairwise_midpoint_distances(selection_gr)
  # Step 2: Perform permutation test to compare with random expectation
  set.seed(123)  # For reproducibility
  random_distances <- vector("numeric", length = n_permutations)
  for (i in 1:n_permutations) {
    # Randomly sample the same number of genes from the background
    random_genes <- sample(background_gr, length(selection_gr), replace = FALSE)
    # Calculate pairwise distances for the randomly sampled genes
    random_distances[i] <- mean(calculate_pairwise_midpoint_distances(random_genes))}
  # Step 3: Compare observed vs random
  mean_observed_distance <- mean(observed_distances)
  mean_random_distance <- mean(random_distances)
  # Step 4: P-value calculation (how often observed distance is smaller than random)
  p_value <- mean(random_distances >= mean_observed_distance)
  # Step 5: Return the results as a list
  return(list(
    mean_observed_distance = mean_observed_distance,
    mean_random_distance = mean_random_distance,
    p_value = p_value,
    random_distances = random_distances,
    observed_distances = observed_distances))}
# Function to output the results in a readable format
print_clustering_results <- function(results, analysis_name) {
  cat("\nResults for", analysis_name, "\n")
  cat("Mean observed distance between genes: ", results$mean_observed_distance, "\n")
  cat("Mean random distance between genes: ", results$mean_random_distance, "\n")
  cat("P-value: ", results$p_value, "\n")}
# Step 6: Run the analysis for each selection file and store the results
paml_results <- run_clustering_analysis(paml_gr, gene_annot_gr)
intens_results <- run_clustering_analysis(intens_gr, gene_annot_gr)
absrel_results <- run_clustering_analysis(absrel_gr, gene_annot_gr)
# Step 7: Print the results for each analysis
print_clustering_results(paml_results, "PAML")
print_clustering_results(intens_results, "INTENS")
print_clustering_results(absrel_results, "ABSREL")
# Step 8: Optionally visualize the distribution of random distances for each analysis
# You can modify the layout of the plots if needed to fit multiple histograms in one output
par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
# PAML plot
hist(paml_results$random_distances, breaks = 30, main = "PAML: Distribution of Random Distances", xlab = "Average Distance (bp)")
abline(v = paml_results$mean_observed_distance, col = "red", lwd = 2, lty = 2)
# INTENS plot
hist(intens_results$random_distances, breaks = 30, main = "INTENS: Distribution of Random Distances", xlab = "Average Distance (bp)")
abline(v = intens_results$mean_observed_distance, col = "red", lwd = 2, lty = 2)
# ABSREL plot
hist(absrel_results$random_distances, breaks = 30, main = "ABSREL: Distribution of Random Distances", xlab = "Average Distance (bp)")
abline(v = absrel_results$mean_observed_distance, col = "red", lwd = 2, lty = 2)
# Reset layout
par(mfrow = c(1, 1))


######### OVERLAP ACROSS ANALYSES #########
######## CHECK CANDIDATE GENES ########
### CAVE SPECIES ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/02_cave_species")
genes.pos <- read.csv("cave_species_absrel_shared.csv", header=TRUE)
genes.intens <- read.csv("cave_species_relax_intensified_OGs.csv", header=TRUE)
genes.paml <- read.csv("cave_species_paml_selected.csv", header=TRUE)
genes.expan <- read.csv("cave_species_cafe_expanded_OGs.csv", header=TRUE)
genes.relax <- read.csv("cave_species_relax_relaxed_OGs.csv", header=TRUE)
genes.contr <- read.csv("cave_species_cafe_contracted_OGs.csv", header=TRUE)
genes.pigment <- read.csv("pigmentation.csv", header=TRUE)
genes.eyes <- read.csv("eyeloss.csv", header=TRUE)
str(genes.pos)
str(genes.intens)
str(genes.paml)
str(genes.expan)
str(genes.relax)
str(genes.contr)
str(genes.pigment)
str(genes.eyes)
merged1 <- merge(genes.pos,genes.intens, by="OG_ID", all=TRUE)
merged2 <- merge(merged1,genes.paml, by="OG_ID", all=TRUE)
merged3 <- merge(merged2,genes.expan, by="OG_ID", all=TRUE)
merged4 <- merge(merged3,genes.relax, by="OG_ID", all=TRUE)
merged5 <- merge(merged4,genes.contr, by="OG_ID", all=TRUE)
str(merged5)
# also annotate with OG
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/01_eggnog_annotation_files")
eggnog.gensize <- read.csv("all_OGs_cave_species_MM_9n4wx4so.emapper.annotations.csv", header = TRUE)
str(eggnog.gensize)
# add a new column just with OG ID
eggnog.gensize <- eggnog.gensize %>% 
  separate(query,c("OG_ID", "sequence.ID"), sep = "-", remove = FALSE) %>%
  dplyr::select(OG_ID, Preferred_name)
str(eggnog.gensize)
colnames(eggnog.gensize)[2] <- "gene_symbol"
str(eggnog.gensize)
merged6 <- merge(merged5,eggnog.gensize, by="OG_ID", all.x=TRUE)
## merge pigment and eye loss genes with OGs from eggnog annotation
merge.pig <- merge(genes.pigment,eggnog.gensize, by="gene_symbol", all.x=TRUE)
merge.egg <- merge(genes.eyes,eggnog.gensize, by="gene_symbol", all.x=TRUE)
str(merge.pig)
str(merge.egg)
merge.pig.eye <- merge(merge.pig,merge.egg, by="OG_ID", all=TRUE)
str(merge.pig.eye)
merge.all <- merge(merged6,merge.pig.eye, by="OG_ID", all=TRUE)
dim(merge.all)
merged7.final <- merge(merge.all,eggnog.gensize, by="OG_ID", all.x=TRUE)
dim(merged7.final)
str(merged7.final)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/02_cave_species")
#write.csv(merged7.final, "cave_species_candidate_gene_selection_analysis.csv", quote=FALSE, row.names=FALSE)
### in excel file remove columns not needed and remove pigmentation and eye loss genes not annotated by any OGs
### analyse file for overlap
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/")
#cave.spec.syn <- read.csv ("cave_species_candidate_gene_selection_analysis.csv", header=TRUE)
### VENN DIAGRAM
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/")
cavespec.syn <- read.csv ("cave_species_candidate_gene_selection_analysis.csv", header=TRUE)
str(cavespec.syn)
# define sets for diagram
absrel <- cavespec.syn %>% 
  filter(analysis_absrel=="absrel_selected") %>% 
  dplyr::select(OG_ID)
intens <- cavespec.syn %>% 
  filter(analysis_relax_intens=="intensified") %>% 
  dplyr::select(OG_ID)
paml <- cavespec.syn %>% 
  filter(analysis_paml=="paml_selected") %>% 
  dplyr::select(OG_ID)
expan <- cavespec.syn %>% 
  filter(analysis_cafe_expanded=="expanded") %>% 
  dplyr::select(OG_ID)
contr <- cavespec.syn %>% 
  filter(analysis_cafe_contr=="contracted") %>% 
  dplyr::select(OG_ID)
relax <- cavespec.syn %>% 
  filter(analysis_relax_relaxed=="relaxed") %>% 
  dplyr::select(OG_ID)
absrel
intens
paml
expan
contr
relax
absrel <- as.list(absrel)
intens <- as.list(intens)
paml <- as.list(paml)
## helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)}
## draw the diagram comparing genes under selection
list3 <- c(absrel,intens,paml)
print(list3)
display_venn(list3,
             category.names = c("absrel" , "intens" , "paml"),
             fill = c("#999999", "#E69F00", "#56B4E9"))

### PIGMENTATION LOSS ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/03_pigmentation_loss")
genes.pos <- read.csv("pigmentation_loss_shared_sel_OGs.csv", header=TRUE)
genes.intens <- read.table("pigmentation_loss_relax_intensified.txt", header=TRUE)
genes.paml <- read.csv("pigmentation_loss_paml_selected.csv", header=TRUE)
genes.expan <- read.csv("pigmentation_loss_cafe_expanded.csv", header=TRUE)
genes.relax <- read.table("pigmentation_loss_relax_relaxed.txt", header=TRUE)
genes.contr <- read.csv("pigmentation_loss_cafe_contracted.csv", header=TRUE)
genes.pigment <- read.csv("pigmentation.csv", header=TRUE)
str(genes.pos)
str(genes.intens)
str(genes.paml)
str(genes.expan)
str(genes.relax)
str(genes.contr)
str(genes.pigment)
merged1 <- merge(genes.pos,genes.intens, by="OG_ID", all=TRUE)
merged2 <- merge(merged1,genes.paml, by="OG_ID", all=TRUE)
merged3 <- merge(merged2,genes.expan, by="OG_ID", all=TRUE)
merged4 <- merge(merged3,genes.relax, by="OG_ID", all=TRUE)
merged5 <- merge(merged4,genes.contr, by="OG_ID", all=TRUE)
str(merged5)
# also annotate with OG
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/01_eggnog_annotation_files")
eggnog.gensize <- read.csv("all_OGs_pigmentation_loss_MM_fe2v4rjm.emapper.annotations.csv", header = TRUE)
str(eggnog.gensize)
# add a new column just with OG ID
eggnog.gensize <- eggnog.gensize %>% 
  separate(query,c("OG_ID", "sequence.ID"), sep = "-", remove = FALSE) %>%
  dplyr::select(OG_ID, Preferred_name)
str(eggnog.gensize)
colnames(eggnog.gensize)[2] <- "gene_symbol"
str(eggnog.gensize)
## merge pigment and eye loss genes with OGs from eggnog annotation
merge.pig <- merge(genes.pigment,eggnog.gensize, by="gene_symbol", all.x=TRUE)
dim(merge.pig)
str(merge.pig)
## merge all
merged6 <- merge(merged5,merge.pig, by="OG_ID", all.x=TRUE)
str(merged6)
## annotate all with eggnog gene symbol
merged7.final <- merge(merged6,eggnog.gensize, by="OG_ID", all.x=TRUE)
dim(merged7.final)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/03_pigmentation_loss")
write.csv(merged7.final, "pigmentation_loss_candidate_gene_selection_analysis.csv", quote=FALSE, row.names=FALSE)
### in excel file remove columns not needed and replace missing gene symbols with NAs
### analyse file for overlap
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/")
pigm.loss.syn <- read.csv ("pigmentation_loss_candidate_gene_selection_analysis.csv", header=TRUE)
### VENN DIAGRAM
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/")
pigloss.syn <- read.csv ("pigmentation_loss_candidate_gene_selection_analysis.csv", header=TRUE)
str(pigloss.syn)
#Define sets for diagram
absrel <- pigloss.syn %>% 
  filter(analysis_absrel=="selected") %>% 
  dplyr::select(OG_ID)
intens <- pigloss.syn %>% 
  filter(analysis_relax_selection=="intensified") %>% 
  dplyr::select(OG_ID)
paml <- pigloss.syn %>% 
  filter(analysis_paml=="paml_selected") %>% 
  dplyr::select(OG_ID)
expan <- pigloss.syn %>% 
  filter(analysis_cafe_expanded=="expanded") %>% 
  dplyr::select(OG_ID)
contr <- pigloss.syn %>% 
  filter(analysis_cafe_contracted=="contracted") %>% 
  dplyr::select(OG_ID)
relax <- pigloss.syn %>% 
  filter(analysis_relax_relaxed=="relaxed") %>% 
  dplyr::select(OG_ID)
absrel
intens
paml
expan
contr
relax
absrel <- as.list(absrel)
intens <- as.list(intens)
paml <- as.list(paml)
## Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)}
#Draw the diagram comparing genes under selection
list3 <- c(absrel,intens,paml)
print(list3)
display_venn(list3,
             category.names = c("absrel" , "intens" , "paml"),
             fill = c("#999999", "#E69F00", "#56B4E9"))

### EYE LOSS ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/04_eye_loss")
genes.pos <- read.csv("eye_loss_absrel_shared.csv", header=TRUE)
genes.intens <- read.csv("eye_loss_relax_intensified.csv", header=TRUE)
genes.paml <- read.csv("eye_loss_paml_selected.csv", header=TRUE)
genes.expan <- read.csv("eye_loss_cafe_expanded.csv", header=TRUE)
genes.relax <- read.table("eye_loss_relax_relaxed.txt", header=TRUE)
genes.contr <- read.csv("eye_loss_cafe_contracted.csv", header=TRUE)
genes.eyes <- read.csv("eyeloss.csv", header=TRUE)
str(genes.pos)
str(genes.intens)
str(genes.paml)
str(genes.expan)
str(genes.relax)
str(genes.contr)
str(genes.eyes)
merged1 <- merge(genes.pos,genes.intens, by="OG_ID", all=TRUE)
merged2 <- merge(merged1,genes.paml, by="OG_ID", all=TRUE)
merged3 <- merge(merged2,genes.expan, by="OG_ID", all=TRUE)
merged4 <- merge(merged3,genes.relax, by="OG_ID", all=TRUE)
merged5 <- merge(merged4,genes.contr, by="OG_ID", all=TRUE)
str(merged5)
# also annotate with OG
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/01_eggnog_annotation_files")
eggnog.gensize <- read.csv("all_OGs_eye_loss_MM_fe2v4rjm.emapper.annotations.csv", header = TRUE)
str(eggnog.gensize)
# add a new column just with OG ID
eggnog.gensize <- eggnog.gensize %>% 
  separate(query,c("OG_ID", "sequence.ID"), sep = "-", remove = FALSE) %>%
  dplyr::select(OG_ID, Preferred_name)
str(eggnog.gensize)
colnames(eggnog.gensize)[2] <- "gene_symbol"
str(eggnog.gensize)
## merge eye loss genes with OGs from eggnog annotation
merge.egg <- merge(genes.eyes,eggnog.gensize, by="gene_symbol", all.x=TRUE)
dim(merge.egg)
str(merge.egg)
## merge all
merged6 <- merge(merged5,merge.egg, by="OG_ID", all.x=TRUE)
str(merged6)
## annotate all with eggnog gene symbol
merged7.final <- merge(merged6,eggnog.gensize, by="OG_ID", all.x=TRUE)
dim(merged7.final)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/04_eye_loss")
write.csv(merged7.final, "eye_loss_candidate_gene_selection_analysis.csv", quote=FALSE, row.names=FALSE)
### in excel file remove columns not needed and replace missing gene symbols with NAs
### analyse file for overlap
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/")
eye.loss.syn <- read.csv ("eye_loss_candidate_gene_selection_analysis.csv", header=TRUE)
### VENN DIAGRAM
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/")
eyeloss.syn <- read.csv ("eye_loss_candidate_gene_selection_analysis.csv", header=TRUE)
str(eyeloss.syn)
# define sets for diagram
absrel <- eyeloss.syn %>% 
  filter(analysis_absrel=="absrel") %>% 
  dplyr::select(OG_ID)
intens <- eyeloss.syn %>% 
  filter(analysis_relax_selection=="intensified") %>% 
  dplyr::select(OG_ID)
paml <- eyeloss.syn %>% 
  filter(analysis_paml=="paml_selected") %>% 
  dplyr::select(OG_ID)
expan <- eyeloss.syn %>% 
  filter(analysis_cafe=="expanded") %>% 
  dplyr::select(OG_ID)
contr <- eyeloss.syn %>% 
  filter(analysis_cafe_contracted=="contracted") %>% 
  dplyr::select(OG_ID)
relax <- eyeloss.syn %>% 
  filter(analysis_relax_relaxed=="relaxed") %>% 
  dplyr::select(OG_ID)
absrel
intens
paml
expan
contr
relax
absrel <- as.list(absrel)
intens <- as.list(intens)
paml <- as.list(paml)
## Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)}
# draw the diagram comparing genes under selection
list3 <- c(absrel,intens,paml)
print(list3)
display_venn(list3,
             category.names = c("absrel" , "intens" , "paml"),
             fill = c("#999999", "#E69F00", "#56B4E9"))

### LONGEVITY ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/05_longevity")
genes.pos <- read.csv("longevity_absrel_shared.csv", header=TRUE)
genes.intens <- read.table("longevity_relax_intensified.txt", header=TRUE)
genes.paml <- read.csv("longevity_paml_selected.csv", header=TRUE)
genes.expan <- read.csv("longevity_cafe_expanded.csv", header=TRUE)
genes.relax <- read.table("longevity_relax_relaxed.txt", header=TRUE)
genes.contr <- read.csv("longevity_cafe_contracted.csv", header=TRUE)
str(genes.pos)
str(genes.intens)
str(genes.paml)
str(genes.expan)
str(genes.relax)
str(genes.contr)
merged1 <- merge(genes.pos,genes.intens, by="OG_ID", all=TRUE)
merged2 <- merge(merged1,genes.paml, by="OG_ID", all=TRUE)
merged3 <- merge(merged2,genes.expan, by="OG_ID", all=TRUE)
merged4 <- merge(merged3,genes.relax, by="OG_ID", all=TRUE)
merged5 <- merge(merged4,genes.contr, by="OG_ID", all=TRUE)
# also annotate with OG
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/08_RELAX/01_eggnog_annotation_files")
eggnog.gensize <- read.csv("all_OGs_longevity_MM_j29p4hp0.emapper.annotations.csv", header = TRUE)
str(eggnog.gensize)
# add a new column just with OG ID
eggnog.gensize <- eggnog.gensize %>% 
  separate(query,c("OG_ID", "sequence.ID"), sep = "-", remove = FALSE) %>%
  dplyr::select(OG_ID, Preferred_name)
str(eggnog.gensize)
colnames(eggnog.gensize)[2] <- "gene_symbol"
str(eggnog.gensize)
merged6 <- merge(merged5,eggnog.gensize, by="OG_ID", all.x=TRUE )
str(merged6)
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/05_longevity")
#write.csv(merged6, "longevity_candidate_gene_selection_analysis.csv", quote=FALSE, row.names=FALSE)
### analyse file for overlap
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/")
longevity.syn <- read.csv ("longevity_candidate_gene_selection_analysis.csv", header=TRUE)
str(longevity.syn)
#Define sets for diagram
absrel <- longevity.syn %>% 
  filter(analysis_absrel=="absrel_selected") %>% 
  dplyr::select(OG_ID)
intens <- longevity.syn %>% 
  filter(analysis_relax_selection=="intensified") %>% 
  dplyr::select(OG_ID)
paml <- longevity.syn %>% 
  filter(analysis_paml=="paml_selected") %>% 
  dplyr::select(OG_ID)
expan <- longevity.syn %>% 
  filter(analysis_cafe_expanded=="expanded") %>% 
  dplyr::select(OG_ID)
contr <- longevity.syn %>% 
  filter(analysis_cafe_contracted=="contracted") %>% 
  dplyr::select(OG_ID)
relax <- longevity.syn %>% 
  filter(analysis_relax_relaxed=="relaxed") %>% 
  dplyr::select(OG_ID)
absrel
intens
paml
expan
contr
relax
absrel <- as.list(absrel)
intens <- as.list(intens)
paml <- as.list(paml)
## Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)}
# draw the diagram comparing genes under selection
list3 <- c(absrel,intens,paml)
print(list3)

#### FIGURE 4e ####
display_venn(list3,
             category.names = c("absrel" , "intens" , "paml"),
             fill = c("#999999", "#E69F00", "#56B4E9"),
             euler.d = FALSE,
             scaled = FALSE)

### collect candidate genes for longevity from published studies
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/05_longevity")
long.table <- read.csv("longevity_literature_candidate_genes.csv", header=TRUE)
str(long.table)
# Summarize the data
summary <- long.table %>%
  group_by(Gene_symbol) %>%
  summarize(
    Reference_Count = n_distinct(Reference),     # Count unique references per gene
    TargetSpecies_Count = n_distinct(species))  # Count unique target species per gene
# View the summary
print(summary)
write.csv(summary, "Longevity_genes_literature.csv", row.names=F, quote=F)
### combine candidate gene list from our study and from literature
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy")
lit.cand.gen <-read.csv("longevity_genes_literature.csv", header=TRUE)
str(lit.cand.gen)
cand.gen <- read.csv("longevity_candidate_gene_selection_analysis.csv",header=TRUE)
str(cand.gen)
result <- merge(cand.gen, lit.cand.gen, by = "gene_symbol", all.x = TRUE)
print(result)
write.csv(result, "longevity_candidate_gene_lit_genes_selection_analysis.csv")


#### HYPERGEOMETRIC TESTS ####
### CAVE SPECIES ###
#### test if genes under relaxed selection in cave species overlap with pigmentation and eyeloss genes
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/02_cave_species")
## define your gene lists
genes_pigment <- read.table("pigmentation_genes.txt")
genes_pigment <- as.list(genes_pigment)
genes_eyes <- read.table("eye_loss_genes.txt")
genes_eyes <- as.list(genes_eyes)
genes_relax <- read.table("cave_species_relax_relaxed.txt") 
genes_relax <- as.list(genes_relax)
genes_paml <- read.table("cave_species_paml_selected.txt")
genes_paml <- as.list(genes_paml)
genes_absrel <- read.table("cave_species_absrel_shared.txt")
genes_absrel <- as.list(genes_absrel)
genes_intens <- read.table("cave_species_relax_intensified.txt")
genes_intens <- as.list(genes_intens)
## define total number of genes
total_genes <- 6064
## calculate the number of overlapping genes
overlap_count_relax_pigment <- length(intersect(genes_relax[[1]], genes_pigment[[1]]))
overlap_count_relax_eyes <- length(intersect(genes_relax[[1]], genes_eyes[[1]]))
overlap_count_relax_paml <- length(intersect(genes_relax[[1]], genes_paml[[1]]))
overlap_count_relax_absrel <- length(intersect(genes_relax[[1]], genes_absrel[[1]]))
overlap_count_paml_pigment <- length(intersect(genes_paml[[1]], genes_pigment[[1]]))
overlap_count_paml_eyes <- length(intersect(genes_paml[[1]], genes_eyes[[1]]))
overlap_count_paml_absrel <- length(intersect(genes_paml[[1]], genes_absrel[[1]]))
overlap_count_paml_intens <- length(intersect(genes_paml[[1]], genes_intens[[1]]))
overlap_count_absrel_pigment <- length(intersect(genes_absrel[[1]], genes_pigment[[1]]))
overlap_count_absrel_eyes <- length(intersect(genes_absrel[[1]], genes_eyes[[1]]))
overlap_count_absrel_intens <- length(intersect(genes_absrel[[1]], genes_intens[[1]]))
overlap_count_intens_pigment <- length(intersect(genes_intens[[1]], genes_pigment[[1]]))
overlap_count_intens_eyes <- length(intersect(genes_intens[[1]], genes_eyes[[1]]))
## calculate the size of each gene list
size_relax <- length(genes_relax[[1]])
size_pigment <- length(genes_pigment[[1]])
size_eyes <- length(genes_eyes[[1]])
size_paml <- length(genes_paml[[1]])
size_absrel <- length(genes_absrel[[1]])
size_intens <- length(genes_intens[[1]])
### helper functions for additional statistical test outputs
# helper function to compute all stats
compute_overlap_stats <- function(A, B, universe_size, label_A, label_B) {
  a <- length(intersect(A, B))                # Overlap
  b <- length(setdiff(B, A))                  # In B not in A
  c <- length(setdiff(A, B))                  # In A not in B
  d <- universe_size - a - b - c              # In neither
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, dimnames = list(c(paste0("In ", label_A), paste0("Not in ", label_A)), 
                                                                       c(paste0("In ", label_B), paste0("Not in ", label_B))))
  # Fisher's exact test
  fisher <- fisher.test(contingency_table)
  # hypergeometric test
  # phyper(q, m, n, k, lower.tail = FALSE)
  # q = a - 1; m = size_A; n = rest of genome; k = size_B
  hyper_pval <- phyper(a - 1, length(A), universe_size - length(A), length(B), lower.tail = FALSE)
  # enrichment
  enrichment <- (a / length(A)) / (length(B) / universe_size)
  list(Label_A = label_A, Label_B = label_B, Overlap = a, Size_A = length(A), Size_B = length(B), Odds_Ratio = fisher$estimate,
    CI_95 = fisher$conf.int, Fisher_P = fisher$p.value, Hyper_P = hyper_pval, Enrichment = enrichment, Contingency_Table = contingency_table)}
## run the tests
stats_relax_pigment <- compute_overlap_stats(genes_relax[[1]], genes_pigment[[1]], total_genes, "Relaxed", "Pigmentation")
stats_relax_eyes <- compute_overlap_stats(genes_relax[[1]], genes_eyes[[1]], total_genes, "Relaxed", "Eye Loss")
## print Results
print("Relaxed vs Pigmentation Genes")
print(stats_relax_pigment)
print("Relaxed vs Eye Loss Genes")
print(stats_relax_eyes)
### perform the hypergeometric test
# phyper(q, m, n, k, lower.tail = FALSE) - calculates the probability of observing q or more successes
## test to check if relaxed genes are likely to be pigment genes
p_value <- phyper(overlap_count_relax_pigment - 1, size_relax, total_genes - size_relax, size_pigment, lower.tail = FALSE)
print(p_value)
## test to check if relaxed genes are likely to be eye genes
p_value <- phyper(overlap_count_relax_eyes - 1, size_relax, total_genes - size_relax, size_eyes, lower.tail = FALSE)
print(p_value)
## test to check if relaxed genes are likely to overlap with genes under selection paml
p_value <- phyper(overlap_count_relax_paml - 1, size_relax, total_genes - size_relax, size_paml, lower.tail = FALSE)
print(p_value)
## test to check if relaxed genes are likely to overlap with genes under selection absrel
p_value <- phyper(overlap_count_relax_absrel - 1, size_relax, total_genes - size_relax, size_absrel, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) genes are likely to be pigment genes
p_value <- phyper(overlap_count_paml_pigment - 1, size_paml, total_genes - size_paml, size_pigment, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to be eye genes
p_value <- phyper(overlap_count_paml_eyes - 1, size_paml, total_genes - size_paml, size_eyes, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to overlap with genes under selection absrel
p_value <- phyper(overlap_count_paml_absrel - 1, size_paml, total_genes - size_paml, size_absrel, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to overlap with genes under selection intens
p_value <- phyper(overlap_count_paml_intens - 1, size_paml, total_genes - size_paml, size_intens, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to be pigmentation genes
p_value <- phyper(overlap_count_absrel_pigment - 1, size_absrel, total_genes - size_absrel, size_pigment, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to be eye genes
p_value <- phyper(overlap_count_absrel_eyes - 1, size_absrel, total_genes - size_absrel, size_eyes, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to overlap with genes under selection intens
p_value <- phyper(overlap_count_absrel_intens - 1, size_absrel, total_genes - size_absrel, size_intens, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (intens) are likely to be pigmentation genes
p_value <- phyper(overlap_count_intens_pigment - 1, size_intens, total_genes - size_intens, size_pigment, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (intens) are likely to be eye genes
p_value <- phyper(overlap_count_intens_eyes - 1, size_intens, total_genes - size_intens, size_eyes, lower.tail = FALSE)
print(p_value)


### PIGMENTATION LOSS ###
#### test if genes under relaxed selection in cave species overlap with pigmentation and eyeloss genes
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/03_pigmentation_loss")
pig_loss_overlap <- read.csv("pigmentation_loss_candidate_gene_selection_analysis.csv", header=TRUE)
str(pig_loss_overlap)
## define your gene lists
absrel <- pig_loss_overlap %>%
  filter(!is.na(analysis_absrel)) %>%  # Exclude rows where 'positive' is NA
  dplyr::select(OG_ID)      # Select only the first column with gene IDs
relax <- pig_loss_overlap %>%
  filter(!is.na(analysis_relax_relaxed)) %>%    
  dplyr::select(OG_ID)      
paml <- pig_loss_overlap %>%
  filter(!is.na(analysis_paml)) %>%   
  dplyr::select(OG_ID)      
intens <- pig_loss_overlap %>%
  filter(!is.na(analysis_relax_selection)) %>%   
  dplyr::select(OG_ID) 
expan <- pig_loss_overlap %>%
  filter(!is.na(analysis_cafe_expanded)) %>%   
  dplyr::select(OG_ID)      
contr <- pig_loss_overlap %>%
  filter(!is.na(analysis_cafe_contracted)) %>%   
  dplyr::select(OG_ID)      
pigment <- pig_loss_overlap %>%
  filter(!is.na(condition)) %>%   
  dplyr::select(OG_ID)      
absrel <- as.list(absrel)
intens <- as.list(intens)
relax <- as.list(relax)
paml <- as.list(paml)
expan <- as.list(expan)
contr <- as.list(contr)
pigment <- as.list(pigment)
## define total number of genes
total_genes <- 5448
## calculate the number of overlapping genes
overlap_count_relax_pigment <- length(intersect(relax[[1]], pigment[[1]]))
overlap_count_relax_paml <- length(intersect(relax[[1]], paml[[1]]))
overlap_count_relax_absrel <- length(intersect(relax[[1]], absrel[[1]]))
overlap_count_paml_pigment <- length(intersect(paml[[1]], pigment[[1]]))
overlap_count_paml_absrel <- length(intersect(paml[[1]], absrel[[1]]))
overlap_count_paml_intens <- length(intersect(paml[[1]], intens[[1]]))
overlap_count_absrel_pigment <- length(intersect(absrel[[1]], pigment[[1]]))
overlap_count_absrel_intens <- length(intersect(absrel[[1]], intens[[1]]))
overlap_count_intens_pigment <- length(intersect(intens[[1]], pigment[[1]]))
## calculate the size of each gene list
size_relax <- length(relax[[1]])
size_pigment <- length(pigment[[1]])
size_paml <- length(paml[[1]])
size_absrel <- length(absrel[[1]])
size_intens <- length(intens[[1]])
### perform the hypergeometric test
# phyper(q, m, n, k, lower.tail = FALSE) - calculates the probability of observing q or more successes
## test to check if relaxed genes are likely to be pigment genes
p_value <- phyper(overlap_count_relax_pigment - 1, size_relax, total_genes - size_relax, size_pigment, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (intens) are likely to be pigmentation genes
p_value <- phyper(overlap_count_intens_pigment - 1, size_intens, total_genes - size_intens, size_pigment, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) genes are likely to be pigment genes
p_value <- phyper(overlap_count_paml_pigment - 1, size_paml, total_genes - size_paml, size_pigment, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to be pigmentation genes
p_value <- phyper(overlap_count_absrel_pigment - 1, size_absrel, total_genes - size_absrel, size_pigment, lower.tail = FALSE)
print(p_value)
## test to check if relaxed genes are likely to overlap with genes under selection paml
p_value <- phyper(overlap_count_relax_paml - 1, size_relax, total_genes - size_relax, size_paml, lower.tail = FALSE)
print(p_value)
## test to check if relaxed genes are likely to overlap with genes under selection absrel
p_value <- phyper(overlap_count_relax_absrel - 1, size_relax, total_genes - size_relax, size_absrel, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to overlap with genes under selection absrel
p_value <- phyper(overlap_count_paml_absrel - 1, size_paml, total_genes - size_paml, size_absrel, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to overlap with genes under selection intens
p_value <- phyper(overlap_count_paml_intens - 1, size_paml, total_genes - size_paml, size_intens, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to overlap with genes under selection intens
p_value <- phyper(overlap_count_absrel_intens - 1, size_absrel, total_genes - size_absrel, size_intens, lower.tail = FALSE)
print(p_value)

### EYE LOSS ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/04_eye_loss")
eye_loss_overlap <- read.csv("eye_loss_candidate_gene_selection_analysis.csv", header=TRUE)
str(eye_loss_overlap)
## define your gene lists
absrel <- eye_loss_overlap %>%
  filter(!is.na(analysis_absrel)) %>%  # Exclude rows where 'positive' is NA
  dplyr::select(OG_ID)      # Select only the first column with gene IDs
relax <- eye_loss_overlap %>%
  filter(!is.na(analysis_relax_relaxed)) %>%    
  dplyr::select(OG_ID)      
paml <- eye_loss_overlap %>%
  filter(!is.na(analysis_paml)) %>%   
  dplyr::select(OG_ID)      
intens <- eye_loss_overlap %>%
  filter(!is.na(analysis_relax_selection)) %>%   
  dplyr::select(OG_ID) 
expan <- eye_loss_overlap %>%
  filter(!is.na(analysis_cafe)) %>%   
  dplyr::select(OG_ID)      
contr <- eye_loss_overlap %>%
  filter(!is.na(analysis_cafe_contracted)) %>%   
  dplyr::select(OG_ID)      
eye <- eye_loss_overlap %>%
  filter(!is.na(condition)) %>%   
  dplyr::select(OG_ID)      
absrel <- as.list(absrel)
intens <- as.list(intens)
relax <- as.list(relax)
paml <- as.list(paml)
expan <- as.list(expan)
contr <- as.list(contr)
eye <- as.list(eye)
## define total number of genes
total_genes <- 5728
## calculate the number of overlapping genes
overlap_count_relax_eye <- length(intersect(relax[[1]], eye[[1]]))
overlap_count_relax_paml <- length(intersect(relax[[1]], paml[[1]]))
overlap_count_relax_absrel <- length(intersect(relax[[1]], absrel[[1]]))
overlap_count_paml_eye <- length(intersect(paml[[1]], eye[[1]]))
overlap_count_paml_absrel <- length(intersect(paml[[1]], absrel[[1]]))
overlap_count_paml_intens <- length(intersect(paml[[1]], intens[[1]]))
overlap_count_absrel_eye <- length(intersect(absrel[[1]], eye[[1]]))
overlap_count_absrel_intens <- length(intersect(absrel[[1]], intens[[1]]))
overlap_count_intens_eye <- length(intersect(intens[[1]], eye[[1]]))
## calculate the size of each gene list
size_relax <- length(relax[[1]])
size_eye <- length(eye[[1]])
size_paml <- length(paml[[1]])
size_absrel <- length(absrel[[1]])
size_intens <- length(intens[[1]])
### perform the hypergeometric test
# phyper(q, m, n, k, lower.tail = FALSE) - calculates the probability of observing q or more successes
## test to check if relaxed genes are likely to be eye genes
p_value <- phyper(overlap_count_relax_eye - 1, size_relax, total_genes - size_relax, size_eye, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (intens) are likely to be eye genes
p_value <- phyper(overlap_count_intens_eye - 1, size_intens, total_genes - size_intens, size_eye, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) genes are likely to be eye genes
p_value <- phyper(overlap_count_paml_eye - 1, size_paml, total_genes - size_paml, size_eye, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to be eye genes
p_value <- phyper(overlap_count_absrel_eye - 1, size_absrel, total_genes - size_absrel, size_eye, lower.tail = FALSE)
print(p_value)
## test to check if relaxed genes are likely to overlap with genes under selection paml
p_value <- phyper(overlap_count_relax_paml - 1, size_relax, total_genes - size_relax, size_paml, lower.tail = FALSE)
print(p_value)
## test to check if relaxed genes are likely to overlap with genes under selection absrel
p_value <- phyper(overlap_count_relax_absrel - 1, size_relax, total_genes - size_relax, size_absrel, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to overlap with genes under selection absrel
p_value <- phyper(overlap_count_paml_absrel - 1, size_paml, total_genes - size_paml, size_absrel, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to overlap with genes under selection intens
p_value <- phyper(overlap_count_paml_intens - 1, size_paml, total_genes - size_paml, size_intens, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to overlap with genes under selection intens
p_value <- phyper(overlap_count_absrel_intens - 1, size_absrel, total_genes - size_absrel, size_intens, lower.tail = FALSE)
print(p_value)


### LONGEVITY ###
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/11_synergy/")
longevity_overlap <- read.csv("longevity_candidate_gene_lit_genes_selection_analysis.csv", header=TRUE)
str(longevity_overlap)
## define your gene lists
absrel <- longevity_overlap %>%
  filter(!is.na(analysis_absrel)) %>%  # Exclude rows where 'positive' is NA
  dplyr::select(OG_ID)      # Select only the first column with gene IDs
relax <- longevity_overlap %>%
  filter(!is.na(analysis_relax_relaxed)) %>%    
  dplyr::select(OG_ID)      
paml <- longevity_overlap %>%
  filter(!is.na(analysis_paml)) %>%   
  dplyr::select(OG_ID)      
intens <- longevity_overlap %>%
  filter(!is.na(analysis_relax_selection)) %>%   
  dplyr::select(OG_ID) 
expan <- longevity_overlap %>%
  filter(!is.na(analysis_cafe_expanded)) %>%   
  dplyr::select(OG_ID)      
contr <- longevity_overlap %>%
  filter(!is.na(analysis_cafe_contracted)) %>%   
  dplyr::select(OG_ID)      
long <- longevity_overlap %>%
  filter(!is.na(condition)) %>%   
  dplyr::select(OG_ID)      
absrel <- as.list(absrel)
intens <- as.list(intens)
relax <- as.list(relax)
paml <- as.list(paml)
expan <- as.list(expan)
contr <- as.list(contr)
long <- as.list(long)
## define total number of genes
total_genes <- 8084
## calculate the number of overlapping genes
overlap_count_relax_long <- length(intersect(relax[[1]], long[[1]]))
overlap_count_relax_paml <- length(intersect(relax[[1]], paml[[1]]))
overlap_count_relax_absrel <- length(intersect(relax[[1]], absrel[[1]]))
overlap_count_paml_long <- length(intersect(paml[[1]], long[[1]]))
overlap_count_paml_absrel <- length(intersect(paml[[1]], absrel[[1]]))
overlap_count_paml_intens <- length(intersect(paml[[1]], intens[[1]]))
overlap_count_absrel_long <- length(intersect(absrel[[1]], long[[1]]))
overlap_count_absrel_intens <- length(intersect(absrel[[1]], intens[[1]]))
overlap_count_intens_long <- length(intersect(intens[[1]], long[[1]]))
# check how many longevity genes are present in the 8084 genes
all.OGs <- read.table("longevity_all_OGs_annotated.txt", header=TRUE)
str(all.OGs)
long.single.copy.OGs <- read.table("longevity_single_copy_OGs.txt", header=TRUE)
str(long.single.copy.OGs)
single.copy.genes.annot <- merge(long.single.copy.OGs, all.OGs, by= "OG", all.x=TRUE)
str(single.copy.genes.annot)
long.cand.lit <- read.csv("longevity_genes_literature.csv", header=TRUE)
str(long.cand.lit)
long.lit.genes.single.copy <- merge(single.copy.genes.annot, long.cand.lit,all= TRUE, by="gene_symbol") 
#write.csv(long.lit.genes.single.copy, "longevity_single_copy_genes_literature_gene_candidates.csv", quote = F, row.names = F)
# overlap of longevity literature genes and is 1293 genes
## calculate the size of each gene list
size_long <- 1293
size_relax <- length(relax[[1]])
size_paml <- length(paml[[1]])
size_absrel <- length(absrel[[1]])
size_intens <- length(intens[[1]])
### perform the hypergeometric test
# phyper(q, m, n, k, lower.tail = FALSE) - calculates the probability of observing q or more successes
## test to check if positively selected (intens) are likely to be longevity genes
p_value <- phyper(overlap_count_intens_long - 1, size_intens, total_genes - size_intens, size_long, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) genes are likely to be longevity genes
p_value <- phyper(overlap_count_paml_long - 1, size_paml, total_genes - size_paml, size_long, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to be longevity genes
p_value <- phyper(overlap_count_absrel_long - 1, size_absrel, total_genes - size_absrel, size_long, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to overlap with genes under selection absrel
p_value <- phyper(overlap_count_paml_absrel - 1, size_paml, total_genes - size_paml, size_absrel, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (paml) are likely to overlap with genes under selection intens
p_value <- phyper(overlap_count_paml_intens - 1, size_paml, total_genes - size_paml, size_intens, lower.tail = FALSE)
print(p_value)
## test to check if positively selected (absrel) are likely to overlap with genes under selection intens
p_value <- phyper(overlap_count_absrel_intens - 1, size_absrel, total_genes - size_absrel, size_intens, lower.tail = FALSE)
print(p_value)
### helper functions for additional statistical test outputs
# helper function to compute all stats
compute_overlap_stats <- function(A, B, universe_size, label_A, label_B) {
  a <- length(intersect(A, B))                # Overlap
  b <- length(setdiff(B, A))                  # In B not in A
  c <- length(setdiff(A, B))                  # In A not in B
  d <- universe_size - a - b - c              # In neither
  contingency_table <- matrix(c(a, b, c, d), nrow = 2,
                              dimnames = list(c(paste0("In ", label_A), paste0("Not in ", label_A)),
                                              c(paste0("In ", label_B), paste0("Not in ", label_B))))
  # Fisher's exact test
  fisher <- fisher.test(contingency_table)
  # Hypergeometric test
  # phyper(q, m, n, k, lower.tail = FALSE)
  # q = a - 1; m = size_A; n = rest of genome; k = size_B
  hyper_pval <- phyper(a - 1, length(A), universe_size - length(A), length(B), lower.tail = FALSE)
  # Enrichment
  enrichment <- (a / length(A)) / (length(B) / universe_size)
  list(Label_A = label_A, Label_B = label_B, Overlap = a, Size_A = length(A), Size_B = length(B), Odds_Ratio = fisher$estimate,
       CI_95 = fisher$conf.int, Fisher_P = fisher$p.value, Hyper_P = hyper_pval, Enrichment = enrichment, Contingency_Table = contingency_table)}
# Run the tests
stats_paml_longlit <- compute_overlap_stats(paml[[1]], long[[1]], total_genes, "Positive selection", "Longevity related")
stats_absrel_longlit <- compute_overlap_stats(absrel[[1]], long[[1]], total_genes, "Convergent positive selection", "Longevity related")
# Print Results
print("Positive vs Longevity Literature Genes")
print(stats_paml_longlit)
print("Convergent positive vs Longeivty Literature Genes")
print(stats_absrel_longlit)

#### MAKE A FIGURE SUMMARISING EXPANDED/CONTRACTED GENE FAMILIES AND GENES UNDER SELECTION #### 
### make phylogenetic radial trees to create Figures ### 
#set working dir
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/13_Figures/05_longevity")
### read in your tree (example using a Newick file)
# tree <- read.tree("longevity_species_list_manip.tre")  # includes some replaced names for taxa not available on phylopic
tree <- read.tree("longevity_species_list_abbrev.tre")
### get images for species from phylopic
### load the data frame for your additional data (e.g., number of expanded/contracted gene families)
# Replace with your own data
info <- read.csv("info.csv")

#### FIGURE 4a ####
p <- ggtree(tree, layout = "circular")
p +
  geom_fruit(data=info, geom=geom_star, mapping=aes(y=id, fill=type),offset= 0.1, size=5, starstroke=0) + 
  geom_fruit(data=info, geom=geom_bar, mapping=aes(x=expansion, y=id, fill= "Expansion"),pwidth =0.3, offset= 1,
    orientation="y", stat="identity", label=info$expansion, axis.params = list(
      axis = 'x', text.size = 2, nbreak = 3, text.angle = -40, vjust = 1, hjust = 0, limits = c(0, 5000)), 
    grid.params = list()) +
  geom_fruit(data=info, geom=geom_bar, mapping=aes(x=contraction, y=id, fill= "Contraction"),pwidth =0.3,
    orientation="y", stat="identity", label=info$contraction, axis.params = list(
      axis = 'x', text.size = 2, nbreak = 2, text.angle = -40, vjust = 1, hjust = 0, limits = c(0, 2000)), 
    grid.params = list()) +  
  geom_fruit(data=info, geom=geom_bar, mapping=aes(x=absrel, y=id, fill= "Absrel"),pwidth =0.3,
    orientation="y", stat="identity", label=info$absrel, axis.params = list(
      axis = 'x', text.size = 2, nbreak = 2, text.angle = -40, vjust = 1, hjust = 0, limits = c(0, 1000)), 
    grid.params = list()) +  
  scale_fill_manual(
    values = c("foreground" = "#CD3333", "background" = "#E0EEEE", "Expansion" = "#E69F00", "Contraction" = "#56B4E9", "Absrel" = "#66CDAA"))


#### CANDIDATE GENE ANALYSIS TDG09 ####
setwd("C:/Users/hans_/Desktop/Projects/Proteus_genome_subterranean_traits/04_analysis/15_convergent_AA/02_tdg09")
### CAVE SPECIES ### 
cand_genes_sel <- read.csv("tdg09_outdf_cave_species.csv", header=TRUE)
str(cand_genes_sel)
names(cand_genes_sel)[names(cand_genes_sel) == "og"] <- "OG_ID"
length(unique(cand_genes_sel$og))
# Sort the data frame by P value (ascending)
cand_genes_sorted <- cand_genes_sel[order(cand_genes_sel$lrt), ]
# Filter for rows with P < 0.05
cand_genes_significant <- subset(cand_genes_sorted, lrt < 0.05)
# load gene symbols 
gene.sym <- read.csv("cave_species_candidate_gene_selection_analysis.csv",header=TRUE)
str(gene.sym)
str(cand_genes_significant)
cand_genes_annotated <- merge(cand_genes_significant, gene.sym, by = "OG_ID", all.x = TRUE)
str(cand_genes_annotated)
# Write the significant results to a new CSV file
#write.csv(cand_genes_annotated, "cave_species_significant_genes.csv", row.names = FALSE)

### PIGMENTATION LOSS ###
cand_genes_sel <- read.csv("tdg09_outdf_pigmentation_loss.csv", header=TRUE)
str(cand_genes_sel)
names(cand_genes_sel)[names(cand_genes_sel) == "og"] <- "OG_ID"
length(unique(cand_genes_sel$og))
# Sort the data frame by P value (ascending)
cand_genes_sorted <- cand_genes_sel[order(cand_genes_sel$lrt), ]
# Filter for rows with P < 0.05
cand_genes_significant <- subset(cand_genes_sorted, lrt < 0.05)
# load gene symbols 
gene.sym <- read.csv("pigmentation_loss_candidate_gene_selection_analysis.csv",header=TRUE)
str(gene.sym)
str(cand_genes_significant)
cand_genes_annotated <- merge(cand_genes_significant, gene.sym, by = "OG_ID", all.x = TRUE)
str(cand_genes_annotated)
# Write the significant results to a new CSV file
#write.csv(cand_genes_annotated, "pigmentation_loss_significant_genes.csv", row.names = FALSE)

### EYE LOSS ###
cand_genes_sel <- read.csv("tdg09_outdf_eye_loss.csv", header=TRUE)
str(cand_genes_sel)
names(cand_genes_sel)[names(cand_genes_sel) == "og"] <- "OG_ID"
length(unique(cand_genes_sel$og))
# Sort the data frame by P value (ascending)
cand_genes_sorted <- cand_genes_sel[order(cand_genes_sel$lrt), ]
# Filter for rows with P < 0.05
cand_genes_significant <- subset(cand_genes_sorted, lrt < 0.05)
# load gene symbols 
gene.sym <- read.csv("eye_loss_candidate_gene_selection_analysis.csv",header=TRUE)
str(gene.sym)
str(cand_genes_significant)
cand_genes_annotated <- merge(cand_genes_significant, gene.sym, by = "OG_ID", all.x = TRUE)
str(cand_genes_annotated)
# Write the significant results to a new CSV file
#write.csv(cand_genes_annotated, "eye_loss_significant_genes.csv", row.names = FALSE)

### LONGEVITY ###
cand_genes_sel <- read.csv("tdg09_outdf_longevity.csv", header=TRUE)
str(cand_genes_sel)
names(cand_genes_sel)[names(cand_genes_sel) == "og"] <- "OG_ID"
length(unique(cand_genes_sel$og))
# Sort the data frame by P value (ascending)
cand_genes_sorted <- cand_genes_sel[order(cand_genes_sel$lrt), ]
# Filter for rows with P < 0.05
cand_genes_significant <- subset(cand_genes_sorted, lrt < 0.05)
# load gene symbols 
gene.sym <- read.csv("longevity_candidate_gene_selection_analysis.csv",header=TRUE)
str(gene.sym)
str(cand_genes_significant)
cand_genes_annotated <- merge(cand_genes_significant, gene.sym, by = "OG_ID", all.x = TRUE)
str(cand_genes_annotated)
# Write the significant results to a new CSV file
#write.csv(cand_genes_annotated, "longevity_significant_genes.csv", row.names = FALSE)


