#####################################################################################
#Figure 2 & Supplementary Figure 3
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in CIRCLES Tomato manuscript
#  Revision 07/25
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
library("phyloseq")
library("DESeq2")
library("UpSetR")
library("ggplot2")
library ("PMCMR")
#############################################################


#############################################################
#set working directory-UNIVDUN
setwd("/cluster/db/R_shared/CIRCLES_tomato/")
#############################################################

##############################################################
#Import the pre-processed Phyloseq object
##############################################################

#import the datasets
JH22_250 <- readRDS("JH22_250_processed.rds")

#################################################################################################################
#Create DESeq objects, one per genotype
#################################################################################################################

#split the phyloseq object
JH22_Abundo <- subset_samples(JH22_250, Genotype != "SV5197TP")
JH22_SV <- subset_samples(JH22_250, Genotype != "Abundo")

#prune ASV with 0 counts upon splitting
JH22_Abundo_2 <- prune_taxa(taxa_sums(JH22_Abundo) > 0, JH22_Abundo) 
JH22_Abundo_2
JH22_SV_2 <- prune_taxa(taxa_sums(JH22_SV) > 0, JH22_SV)
JH22_SV_2

##################################
#Abundo
##################################

#extract count data 
JH22_Abundo_counts <- otu_table(JH22_Abundo_2)
countData = as.data.frame(JH22_Abundo_counts)
colnames(JH22_Abundo_counts)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH22_Abundo_2)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
JH22_Abundo_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Microhabitat)

#execute the differential count analysis with the function DESeq 
JH22_Abundo_test <- DESeq(JH22_Abundo_cds, fitType="local", betaPrior = FALSE) 

#define the ASVs differentially enriched in plant compartments samples
Abundo_root <- results(JH22_Abundo_test, contrast = c("Microhabitat", "Bulk", "Root")) 
Abundo_rhizosphere <- results(JH22_Abundo_test, contrast = c("Microhabitat","Bulk", "Rhizosphere")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05
Abundo_root_FDR005 <- Abundo_root[(rownames(Abundo_root)[which(Abundo_root$padj <0.05)]), ]
Abundo_rhizosphere_FDR005 <- Abundo_rhizosphere[(rownames(Abundo_rhizosphere)[which(Abundo_rhizosphere$padj <0.05)]), ]

#enriched in plant compartments (negative fold change, all plant compartments are at the second term of comparison)
Abundo_root_enriched <-  Abundo_root_FDR005[(rownames(Abundo_root_FDR005)[which(Abundo_root_FDR005$log2FoldChange < 0)]), ]
Abundo_rhizosphere_enriched <-  Abundo_rhizosphere_FDR005[(rownames(Abundo_rhizosphere_FDR005)[which(Abundo_rhizosphere_FDR005$log2FoldChange < 0)]), ]

#stat_info
#write.table(Abundo_root_enriched, file="Ab_root_enriched_ASV_check.txt", sep="\t")
#write.table(Abundo_rhizosphere_enriched, file="Ab_root_rhizosphere_ASV_check.txt", sep="\t")

#Extract the ASVs identifiers, these will be used for a second part of this code
Abundo_root_enriched_ASVs <- rownames(Abundo_root_enriched)
Abundo_rhizosphere_enriched_ASVs <- rownames(Abundo_rhizosphere_enriched)

#How many ASVs?
length(Abundo_root_enriched_ASVs)
length(Abundo_rhizosphere_enriched_ASVs)

#Save the lists of ASVs for parallel analysis
#write(Abundo_root_enriched_ASVs, file = "Abundo_root_enriched_ASVs_data_check.txt")
#write(Abundo_rhizosphere_enriched_ASVs, file = "Abundo_rhizoshere_enriched_ASVs_data_check.txt")

##################################
#SV5197TP
##################################

#extract count data 
JH22_SV_counts <- otu_table(JH22_SV_2)
countData = as.data.frame(JH22_SV_counts)
colnames(JH22_SV_counts)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH22_SV_2)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
JH22_SV_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Microhabitat)

#execute the differential count analysis with the function DESeq 
JH22_SV_test <- DESeq(JH22_SV_cds, fitType="local", betaPrior = FALSE) 

#define the ASVs differentially enriched in plant compartments samples
SV_root <- results(JH22_SV_test, contrast = c("Microhabitat", "Bulk", "Root")) 
SV_rhizosphere <- results(JH22_SV_test, contrast = c("Microhabitat","Bulk", "Rhizosphere")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05
SV_root_FDR005 <- SV_root[(rownames(SV_root)[which(SV_root$padj <0.05)]), ]
SV_rhizosphere_FDR005 <- SV_rhizosphere[(rownames(SV_rhizosphere)[which(SV_rhizosphere$padj <0.05)]), ]

#enriched in plant compartments (negative fold change, all plant compartments are at the second term of comparison)
SV_root_enriched <-  SV_root_FDR005[(rownames(SV_root_FDR005)[which(SV_root_FDR005$log2FoldChange < 0)]), ]
SV_rhizosphere_enriched <-  SV_rhizosphere_FDR005[(rownames(SV_rhizosphere_FDR005)[which(SV_rhizosphere_FDR005$log2FoldChange < 0)]), ]

#stat-info
#write.table(SV_root_enriched, file="SV_root_enriched_ASV_check.txt", sep="\t")
#write.table(SV_rhizosphere_enriched, file="SV_root_rhizosphere_ASV_check.txt", sep="\t")

#Extract the ASVs identifiers, these will be used for a second part of this code
SV_root_enriched_ASVs <- rownames(SV_root_enriched)
SV_rhizosphere_enriched_ASVs <- rownames(SV_rhizosphere_enriched)

#How many ASVs?
length(SV_root_enriched_ASVs)
length(SV_rhizosphere_enriched_ASVs)

#Save the lists of ASVs for parallel analysis
#write(SV_root_enriched_ASVs, file = "SV_root_enriched_ASVs_data_check.txt")
#write(SV_rhizosphere_enriched_ASVs, file = "SV_rhizoshere_enriched_ASVs_data_check.txt")

##########################
#Plotting rhizosphere
#########################

#Abundo
Abundo_rhizosphere_enriched_counts <- as.data.frame(Abundo_rhizosphere_enriched[, 1])
dim(Abundo_rhizosphere_enriched_counts)
rownames(Abundo_rhizosphere_enriched_counts) <- rownames(Abundo_rhizosphere_enriched)
colnames(Abundo_rhizosphere_enriched_counts) <- c("counts_Abundo_rhizosphere")
Abundo_rhizosphere_enriched_counts[Abundo_rhizosphere_enriched_counts > 1] <- 1
dim(Abundo_rhizosphere_enriched_counts)
#SV5197TP
SV_rhizosphere_enriched_counts <- as.data.frame(SV_rhizosphere_enriched[, 1])
dim(SV_rhizosphere_enriched_counts)
rownames(SV_rhizosphere_enriched_counts) <- rownames(SV_rhizosphere_enriched)
colnames(SV_rhizosphere_enriched_counts) <- c("counts_SV_rhizosphere")
SV_rhizosphere_enriched_counts[SV_rhizosphere_enriched_counts > 1] <- 1
dim(SV_rhizosphere_enriched_counts)

#combine the datasets: note they have unequal values
#define a list of unique ASVs
ASV_list <- unique(c(rownames(Abundo_rhizosphere_enriched_counts), rownames(SV_rhizosphere_enriched_counts)))
length(ASV_list)
ASV_list_rhizosphere <- ASV_list
#Prior combining the dataset, we need to account for ASV unevenly distributed (i.e., enriched in one compartment not in others)
#Abundo
Abundo_Rhizosphere_eriched_merging <- as.data.frame(Abundo_rhizosphere_enriched_counts[ASV_list_rhizosphere, ])
colnames(Abundo_Rhizosphere_eriched_merging) <- c("ASV_rhizo_Abundo")
row.names(Abundo_Rhizosphere_eriched_merging) <- as.vector(ASV_list_rhizosphere)
#SV
SV_Rhizosphere_eriched_merging <- as.data.frame(SV_rhizosphere_enriched_counts[ASV_list_rhizosphere, ])
colnames(SV_Rhizosphere_eriched_merging) <- c("ASV_rhizo_SV")
row.names(SV_Rhizosphere_eriched_merging) <- as.vector(ASV_list_rhizosphere)

#Merge the datasets
rhizosphere_ASVs <- cbind(Abundo_Rhizosphere_eriched_merging, SV_Rhizosphere_eriched_merging)
#set NA to 0: as some ASVs won't be present in certain microhabitats and NA is a non numerical character
rhizosphere_ASVs[is.na(rhizosphere_ASVs)] <- 0
dim(rhizosphere_ASVs)

##########################
#Plotting roots
#########################

#Abundo
Abundo_root_enriched_counts <- as.data.frame(Abundo_root_enriched[, 1])
dim(Abundo_root_enriched_counts)
rownames(Abundo_root_enriched_counts) <- rownames(Abundo_root_enriched)
colnames(Abundo_root_enriched_counts) <- c("counts_Abundo_root")
Abundo_root_enriched_counts[Abundo_root_enriched_counts > 1] <- 1
dim(Abundo_root_enriched_counts)
#SV5197TP
SV_root_enriched_counts <- as.data.frame(SV_root_enriched[, 1])
dim(SV_root_enriched_counts)
rownames(SV_root_enriched_counts) <- rownames(SV_root_enriched)
colnames(SV_root_enriched_counts) <- c("counts_SV_root")
SV_root_enriched_counts[SV_root_enriched_counts > 1] <- 1
dim(SV_root_enriched_counts)

#combine the datasets: note they have unequal values
#define a list of unique ASVs
ASV_list <- unique(c(rownames(Abundo_root_enriched_counts), rownames(SV_root_enriched_counts)))
length(ASV_list)
ASV_list_root <- ASV_list
#Prior combining the dataset, we need to account for ASV unevenly distributed (i.e., enriched in one compartment not in others)
#Abundo
Abundo_root_eriched_merging <- as.data.frame(Abundo_root_enriched_counts[ASV_list_root, ])
colnames(Abundo_root_eriched_merging) <- c("ASV_root_Abundo")
row.names(Abundo_root_eriched_merging) <- as.vector(ASV_list_root)
#SV
SV_root_eriched_merging <- as.data.frame(SV_root_enriched_counts[ASV_list_root, ])
colnames(SV_root_eriched_merging) <- c("ASV_root_SV")
row.names(SV_root_eriched_merging) <- as.vector(ASV_list_root)

#Merge the datasets
root_ASVs <- cbind(Abundo_root_eriched_merging, SV_root_eriched_merging)
#set NA to 0: as some ASVs won't be present in certain microhabitats and NA is a non numerical character
root_ASVs[is.na(root_ASVs)] <- 0
dim(root_ASVs)

##########################
#Merge the dataset prior visualisation
#########################
enriched_list <- unique(c(rownames(rhizosphere_ASVs), rownames(root_ASVs)))
length(enriched_list)
#Prior combining the dataset, we need to account for ASV unevenly distributed (i.e., enriched in one compartment not in others)
#rhizosphere
rhizosphere_ASVs_merging <- as.data.frame(rhizosphere_ASVs[enriched_list, ])
row.names(rhizosphere_ASVs_merging) <- as.vector(enriched_list)
dim(rhizosphere_ASVs_merging)
#root
root_ASVs_merging <- as.data.frame(root_ASVs[enriched_list, ])
row.names(root_ASVs_merging) <- as.vector(enriched_list)
dim(root_ASVs_merging)

#Merge the datasets
enriched_ASVs <- cbind(rhizosphere_ASVs_merging, root_ASVs_merging)
#set NA to 0: as some ASVs won't be present in certain microhabitats and NA is a non numerical character
enriched_ASVs[is.na(enriched_ASVs)] <- 0
dim(enriched_ASVs)
colnames(enriched_ASVs)

####################
# Figure 2 graphics
###################

dev.off()
plot_UPSET_belowground <-upset(enriched_ASVs, sets = c("ASV_rhizo_Abundo", "ASV_rhizo_SV", "ASV_root_Abundo", "ASV_root_SV"), sets.bar.color = "#56B4E9",
                               order.by = "freq", sets.x.label = "Enriched ASVs Rhizosphere", mainbar.y.label = "Intersection",)
plot_UPSET_belowground

###################################
# Figure S3 Acinetobacter dominance
##################################
#Enriched ASV
JH22_250_enriched <- prune_taxa(enriched_list, JH22_250)
JH22_250_enriched_plant <- subset_samples(JH22_250_enriched, Microhabitat != "Bulk")
JH22_250_enriched_plant <- subset_samples(JH22_250_enriched_plant, Microhabitat != "Leaf")
#Total ASV
JH22_250_plant <- subset_samples(JH22_250, Microhabitat != "Bulk")
JH22_250_plant <- subset_samples(JH22_250_plant, Microhabitat != "Leaf")

#% of enriched ASV vs. total
sort(sample_sums(JH22_250_enriched_plant)/sample_sums(JH22_250_plant))
#dominant ASV identified

#extract Acinetobacter from the dataset for comparative visualisation across Microhabitats
Acinetobacter_enriched_plant <- prune_taxa("TACAGAGGGTGCAAGCGTTAATCGGATTTACTGGGCGTAAAGCGCGCGTAGGCGGCTAATTAAGTCAAATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGTTAGCTAGAGTGTGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAACACTGACGCTGAGGTGCGAAAGCATGGGGAGCAAACAGG", JH22_250)

#% of enriched Acinetobacter vs. total
Acinetobacter_enriched <- as.data.frame(sample_sums(Acinetobacter_enriched_plant)/sample_sums(JH22_250))
colnames(Acinetobacter_enriched) <- c("Actinetobacter_proportion")

#extract the mapping file
Acinetobacter_map <- as.data.frame(as.matrix(sample_data(JH22_250)))

#merge the dataset
Acinetobacter_plotting <- cbind(Acinetobacter_map, Acinetobacter_enriched)

#Order the factors
Acinetobacter_plotting$Genotype <- ordered(Acinetobacter_plotting$Genotype, levels=c("Unplanted", "Abundo", "SV5197TP"))
Acinetobacter_plotting$Microhabitat <- ordered(Acinetobacter_plotting$Microhabitat, levels=c("Bulk", "Rhizosphere", "Root", "Leaf"))

#colors
CIRCLES_TOM <- c("#999999", "#009E73", "#F179A7")

#Plotting Acinetobacter
dev.off()
p <- ggplot(Acinetobacter_plotting, aes(x=Microhabitat, y=Actinetobacter_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = CIRCLES_TOM)
p + geom_point(aes(y=Actinetobacter_proportion, group=Genotype), position = position_dodge(width=0.75))

#stats
shapiro.test(Acinetobacter_plotting$Actinetobacter_proportion)

#apply a non paramtric test
kruskal.test(Actinetobacter_proportion ~ Microhabitat, data = Acinetobacter_plotting)

#post-hoc
pairwise.wilcox.test(Acinetobacter_plotting$Actinetobacter_proportion, Acinetobacter_plotting$Microhabitat, p.adjust.method = "BH")

