#####################################################################################
#PREPROCESING
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in circles Tomato manuscript
#  Revision 12/24
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
#############################################################

#############################################################
#set working directory-UNIVDUN
setwd("/cluster/db/R_shared/CIRCLES_tomato/")
#############################################################

##################################################################
#Import the Phyloseq object and familiarise with its structure
#################################################################

#import the datasets
JH22_250 <-readRDS("JH22_dada2_silva_138.1.rds")

#remove blank water samples
JH22_no_water <-subset_samples(JH22_250, (Genotype!="BlankControl"))
JH22_no_water
#173 samples to 167 samples

##################################################################
#Remove ASVs assigned to Chloroplast and Mitochondria but retain NA
#rationale: 16S rRNA primers may amplify plant-derived sequences 
#those may interfere with data analysis (we are after bacteria, not plants)
#Note that we use so-called PCR clamps in our PCR reaction to minimise plant cross amplification
#More info here: https://www.nature.com/articles/nmeth.2634 
#################################################################

#Chloroplast
JH22_no_chlor <-subset_taxa(JH22_no_water, (Order!="Chloroplast") | is.na(Order))
JH22_no_chlor
#Mitochondria
JH22_no_plants <-subset_taxa(JH22_no_chlor, (Family!="Mitochondria") | is.na(Family))
JH22_no_plants

##################################################################
#Prune putative contaminant ASVs
#Rationale: bacteria are everywhere and 16S rRNA amplification is very efficient
#We use a list of bacteria likely representing the contamination we may face in our lab as baseline
#More info here: https://www.frontiersin.org/articles/10.3389/fmicb.2018.01650/full 
##################################################################

#Import the list of contaminat ASVs from JH06 library
Contaminant_ASVs <- read.delim("JH06_contaminant_ASVs_ids.txt", header = FALSE)

#identify the proportion of putative contaminants in the merged object
JH22_contaminants <- intersect(taxa_names(JH22_no_plants),Contaminant_ASVs)
JH22_contaminants
#no known contaminants in the data set
#########################################################################
# Remove ASVs assigned to NA at phylum level
# Rationale: if a sequence cannot be assigned at a very high level such as Phylum is not of use us
#########################################################################

JH22_no_plants_1 <- subset_taxa(JH22_no_plants, Phylum!= "NA")
JH22_no_plants_1

#########################################################################
# inspect reads distribution across objects
# Rationale: if we want to compare like-with-like we should know sample distribution
#########################################################################

sort(sample_sums(JH22_no_plants_1))
#35 samples below 5K, consider whether they are relevant for the investigation

hist(sample_sums(JH22_no_plants_1))
#almost bi-modal distribution

#better visualisation
hist(sample_sums(JH22_no_plants_1), main = paste("Reads distribution, raw microbial data"), xlab = paste("Reads"), ylab = paste ("number of samples"), ylim = c(0,60))

#prune ASV with 0 counts
JH22_Sample_integer <- prune_taxa(taxa_sums(JH22_no_plants_1) > 0, JH22_no_plants_1) 
JH22_Sample_integer
sample_data(JH22_Sample_integer)

################################################################################################################
# Define technical reproducibility
# Rationale: not all PCR products are created equal.. we use technical replicates (i.e., the same DNA preparation)
#sequenced multiple times to define the minimum amount of sequencing reads we will be able to identify across experiments
#and library to library variation
#Further info: https://www.nature.com/articles/nature11336 
#warning: this is a hard-wired piece of code, do not modify
################################################################################################################

#apply threshold (we have 4 microhabitats and two genotypes and ~170 samples, assuming an ASV present in 10% of samples in case is microhabitat & genotype specificity)
JH22_Sample_theshold = filter_taxa(JH22_Sample_integer, function(x) sum(x > 20) > (0.01 *length(x)), TRUE)
JH22_Sample_theshold 
#measure the impact of thresholding
sort(sample_sums(JH22_Sample_theshold))

##ratio filtered reads/total reads
ratio <- sum(sample_sums(JH22_Sample_theshold))/sum(sample_sums(JH22_Sample_integer))*100
ratio
#93.3
#Filter samples with less than 5,000 reads
JH22_Sample_theshold_5K <- subset_samples(JH22_Sample_theshold, sample_sums(JH22_Sample_theshold) > 5000) 
JH22_Sample_theshold_5K

#and another look at sample distribution
sort(sample_sums(JH22_Sample_theshold_5K))

hist(sample_sums(JH22_Sample_theshold_5K))

#better visualisation
hist(sample_sums(JH22_Sample_theshold_5K), main = paste("Reads distribution, post-processing"), xlab = paste("Reads"), ylab = paste ("number of samples"), ylim = c(0,50))

#save the object for the next steps of the analysis
#saveRDS(JH22_Sample_theshold_5K, file = "JH22_250_processed.rds")

#END 