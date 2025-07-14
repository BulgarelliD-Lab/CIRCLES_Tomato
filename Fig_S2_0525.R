#####################################################################################
#TAXONOMY CHARTS
####################################################################################
#############################################################
#
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
library("phyloseq")

#############################################################
#set working directory-UNIVDUN
setwd("/cluster/db/R_shared/CIRCLES_tomato/")
#############################################################

#import the dataset
JH22_calculation  <-readRDS("JH22_250_processed.rds")

#inspect the file
JH22_calculation

#number of reads
sum(sample_sums(JH22_calculation))
max(sample_sums(JH22_calculation))
min(sample_sums(JH22_calculation))
median(sample_sums(JH22_calculation))

###################################################################################################
#arrange data for plotting
#note that here we have two factors
#1) the microhabitat
#2) the genotype
##################################################################################################

#prune bulk soil samples
JH22_calculation_plant <- subset_samples(JH22_calculation, Microhabitat != "Bulk")

#Microhabitat
sample_data(JH22_calculation_plant)$Microhabitat <- factor(sample_data(JH22_calculation_plant)$Microhabitat, levels=c("Rhizosphere", "Root", "Leaf"))
#Genotype
sample_data(JH22_calculation_plant)$Genotype <- factor(sample_data(JH22_calculation_plant)$Genotype, levels=c("Abundo", "SV5197TP"))

#transform in relative abundance for plotting
JH22_cpm <- transform_sample_counts(JH22_calculation_plant,  function(x) 1e+06 * x/sum(x))

#################################################################################################
#subset for most dominant phyla in the plant microbiota
#################################################################################################
JH22_B <-subset_taxa(JH22_calculation_plant, Phylum== "Bacteroidota")
JH22_Actino <-subset_taxa(JH22_calculation_plant, Phylum== "Actinobacteriota")
JH22_P <-subset_taxa(JH22_calculation_plant, Phylum== "Proteobacteria")
JH22_F <-subset_taxa(JH22_calculation_plant, Phylum== "Firmicutes")
JH22_Acido <-subset_taxa(JH22_calculation_plant, Phylum== "Acidobacteriota")
#merge the phyloseq objects
JH22_plotting <- merge_phyloseq(JH22_B, JH22_Actino, JH22_P, JH22_F, JH22_Acido)

#proportion of dominant phyla
mean(sample_sums(JH22_plotting))/mean(sample_sums(JH22_calculation_plant))

#https://joey711.github.io/phyloseq/plot_bar-examples.html
plot_bar(JH22_plotting, "Phylum", facet_grid=Genotype~Microhabitat)


