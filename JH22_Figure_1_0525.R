#####################################################################################
#Figure 1
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in circles Tomato manuscript
#  Revision 05/25
#  srobertsonalberty002@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk 
#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
#required packages #info Phyloseq https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217 
library("phyloseq")
library("ggplot2")
library ("vegan")
library ("PMCMR")
#############################################################

#############################################################
#set working directory-UNIVDUN
setwd("/cluster/db/R_shared/CIRCLES_tomato/")
#############################################################

##################################################################
#Import the pre-processed Phyloseq object
#################################################################

#import the datasets
JH22_250 <-readRDS("JH22_250_processed.rds")
#more info about .rds format https://riptutorial.com/r/example/3650/rds-and-rdata--rda--files
sample_data(JH22_250)

#check reads distribution
hist(sample_sums(JH22_250))
sort(sample_sums(JH22_250))
#5099-73503
###################################################################################################
#Rarefy the dataset for comparison purposes
#Rationale: differences in sequencing depth may result in artificial differences in alpha diversity
#Warning on rarefication: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531 
#Information on data treatment (and rerefication): https://pubmed.ncbi.nlm.nih.gov/28253908/ 
###################################################################################################

#the rarefication step - LOWEST READS ARE 5155 SO RAREFY TO 5000
#JH22_250_rare_5K <- rarefy_even_depth(JH22_250, sample.size = 5000)
#0 OTUS REMOVED
#save the generated object for the reproducibility of the code
#saveRDS(JH22_250_rare_5K, file = "JH22_250_processed_rare_5k_data_check.rds")

#import the datasets
JH22_250_rare <-readRDS("JH22_250_processed_rare_5k_data_check.rds")

###################################################################################################
#Alpha diversity indexes calculation, panel A and B
###################################################################################################

#Index calculations
JH22_250_rare_alpha <-  estimate_richness(JH22_250_rare, measures = c("Observed", "Shannon"))

#inspect the file
JH22_250_rare_alpha 
dim(JH22_250_rare_alpha)

###################################################################################################
#Create new datasets for data visualisation and statistical analysis
###################################################################################################

#Sample information
JH22_map <- as.data.frame(as.matrix(sample_data(JH22_250_rare)))
JH22_map
dim(JH22_map)
colnames(JH22_map)

#Microhabitat
design_Microhabitat <- as.data.frame(JH22_map[, 1])
rownames(design_Microhabitat) <- rownames(JH22_map)
colnames(design_Microhabitat) <- c("Microhabitat")
design_Microhabitat 

#Genotype
design_Genotype <- as.data.frame(JH22_map[, 2])
rownames(design_Genotype) <- rownames(JH22_map)
colnames(design_Genotype) <- c("Genotype")
design_Genotype

#data frame Microhab_Genotype
design_MG <- cbind(design_Microhabitat, design_Genotype)
design_MG

#Observed ASVs [1st column of the rare file]
JH22_Observed <- as.data.frame(JH22_250_rare_alpha[ ,1])
rownames(JH22_Observed) <- rownames(JH22_250_rare_alpha)
colnames(JH22_Observed) <- c("Observed")
JH22_Observed
dim(JH22_Observed)
#129 1

#Shannon [4th column of the rare file]
JH22_Shannon <- as.data.frame(JH22_250_rare_alpha[ ,2])
rownames(JH22_Shannon) <- rownames(JH22_250_rare_alpha)
colnames(JH22_Shannon) <- c("Shannon")
JH22_Shannon
dim(JH22_Shannon)
#129 1

#make sure the order of the samples is the same
rownames(design_MG)
rownames(JH22_Observed)
rownames(JH22_Shannon)

#Combine the dataset field Microhabitat and Observed ASVs
JH22_Observed_MG <- cbind(design_MG, JH22_Observed)
JH22_Observed_MG 

#Combine the dataset field Microhabitat and Shannon
JH22_Shannon_MG <- cbind(design_MG, JH22_Shannon)
JH22_Shannon_MG 

#Order the levels according to a defined order
#Observed
JH22_Observed_MG$Genotype <- ordered(JH22_Observed_MG$Genotype, levels=c("Unplanted","Abundo", "SV5197TP"))
JH22_Observed_MG$Microhabitat <- ordered(JH22_Observed_MG$Microhabitat, levels=c("Bulk", "Rhizosphere", "Root", "Leaf"))

#Shannon
JH22_Shannon_MG$Genotype <- ordered(JH22_Shannon_MG$Genotype, levels=c("Unplanted","Abundo", "SV5197TP"))
JH22_Shannon_MG$Microhabitat <- ordered(JH22_Shannon_MG$Microhabitat, levels=c("Bulk", "Rhizosphere", "Root", "Leaf"))

#Colors (we will be coloring the Genotype)
# colour-blind pallette: https://ristretto.black/are-your-documents-colourblind-friendly/ 

#JH22 
CIRCLES_TOM <- c("#999999", "#009E73", "#F179A7")

#Plotting Observed
p <- ggplot(JH22_Observed_MG, aes(x=Microhabitat, y=Observed, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = CIRCLES_TOM)
p + geom_point(aes(y=Observed, group=Genotype), position = position_dodge(width=0.75))

#Plotting Shannon 
p <- ggplot(JH22_Shannon_MG, aes(x=Microhabitat, y=Shannon, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = CIRCLES_TOM)
p + geom_point(aes(y=Shannon, group=Genotype), position = position_dodge(width=0.75))

###################################################################################################
#Stats: check the microhabitat effect corrected for genotype; same bulk soil included in both sets
###################################################################################################

#subset for Abundo
############

#Observed
JH22_Observed_MG_Abundo <- JH22_Observed_MG[which(JH22_Observed_MG$Genotype !="SV5197TP"), ]

#check the distribution of the data
shapiro.test(JH22_Observed_MG_Abundo$Observed)

#apply a non paramtric test
kruskal.test(Observed ~ Genotype, data = JH22_Observed_MG_Abundo)

#post-hoc
pairwise.wilcox.test(JH22_Observed_MG_Abundo$Observed, JH22_Observed_MG_Abundo$Microhabitat, p.adjust.method = "BH")

#Shannon
JH22_Shannon_MG_Abundo <- JH22_Shannon_MG[which(JH22_Shannon_MG$Genotype !="SV5197TP"), ]

#check the distribution of the data
shapiro.test(JH22_Shannon_MG_Abundo$Shannon)

#apply a non paramtric test
kruskal.test(Shannon ~ Genotype, data = JH22_Shannon_MG_Abundo)

#post-hoc
pairwise.wilcox.test(JH22_Shannon_MG_Abundo$Shannon, JH22_Shannon_MG_Abundo$Microhabitat, p.adjust.method = "BH")

#subset for SV
##############

#Observed
JH22_Observed_MG_SV <- JH22_Observed_MG[which(JH22_Observed_MG$Genotype !="SV"), ]

#check the distribution of the data
shapiro.test(JH22_Observed_MG_SV$Observed)

#apply a non paramtric test
kruskal.test(Observed ~ Genotype, data = JH22_Observed_MG_SV)

#post-hoc
pairwise.wilcox.test(JH22_Observed_MG_SV$Observed, JH22_Observed_MG_SV$Microhabitat, p.adjust.method = "BH")

#Shannon
JH22_Shannon_MG_SV <- JH22_Shannon_MG[which(JH22_Shannon_MG$Genotype !="SV5197TP"), ]

#check the distribution of the data
shapiro.test(JH22_Shannon_MG_SV$Shannon)

#apply a non paramtric test
kruskal.test(Shannon ~ Genotype, data = JH22_Shannon_MG_SV)

#post-hoc
pairwise.wilcox.test(JH22_Shannon_MG_SV$Shannon, JH22_Shannon_MG_SV$Microhabitat, p.adjust.method = "BH")


###################################################################################################
#Beta diversity calculation, panel C
#working on data not rarefied
###################################################################################################

#transform into relative abundance
JH22_250_cpm <- transform_sample_counts(JH22_250,  function(x) 1e+06 * x/sum(x))

#ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
JH22_250.ord <- ordinate(JH22_250_cpm, "NMDS", "bray")

#visualisation
plot_ordination(JH22_250_cpm, JH22_250.ord, type="samples", color="Genotype", shape="Microhabitat")

#re-order the factors
sample_data(JH22_250_cpm)$Genotype <- factor(sample_data(JH22_250_cpm)$Genotype, levels=c("Unplanted","Abundo", "SV5197TP"))
sample_data(JH22_250_cpm)$Microhabitat <- factor(sample_data(JH22_250_cpm)$Microhabitat, levels=c("Bulk","Rhizosphere", "Root", "Leaf"))

#ggplots function to increase effectivness of the visualisation
p = plot_ordination(JH22_250_cpm, JH22_250.ord, type="samples", color="Genotype", shape="Microhabitat") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17, 15, 18))
p = p + scale_colour_manual(values = CIRCLES_TOM)
p

#Permanova
BC_bacteria  <- phyloseq::distance(JH22_250_cpm, "bray")
BC_bacteria 

#here we can use a formula ANOVA-like to identify factors of interest
Stat <- adonis2(BC_bacteria ~ Microhabitat, data= as.data.frame(as.matrix(sample_data(JH22_250_cpm))), permutations = 5000)
Stat

#end 