#####################################################################################
#PREPROCESING
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations in PlaBase
#  Revision 12/24
#  d.bulgarelli@dundee.ac.uk 

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
#required packages 
library("phyloseq")

#############################################################
#set working directory: DB cpu
setwd("/cluster/db/R_shared/CN04_plabase/")

#import the constituent components
#counts
CN04_counts <- read.delim("CN04_counts_plabase.txt", row.names =1)
#taxonomies
CN04_tax <- read.delim("CN04_taxonomy_plabase.txt", row.names =1)
#samples
CN04_map <- read.delim("CN04_map_plabase.txt", row.names =1)

#ensuring order of samples is congruent
CN04_map_2 <- CN04_map[colnames(CN04_counts), ]

#convert taxonomies into matrices
CN04_tax <- as.matrix(CN04_tax)

#create the phyloseq object
counts = otu_table(CN04_counts , taxa_are_rows = TRUE)
tax = tax_table(CN04_tax)
samples = sample_data(CN04_map_2)

CN04_data <- phyloseq(counts, tax, samples)
CN04_data

#########################################################################
# inspect reads distribution across objects
# Rationale: if we want to compare like-with-like we should know sample distribution
#########################################################################

#summary
sum (sample_sums(CN04_data))
min(sample_sums(CN04_data))
max(sample_sums(CN04_data))
mean(sample_sums(CN04_data))

#graphical outputs
sort(sample_sums(CN04_data))
hist(sample_sums(CN04_data), main = paste("CN04 annotate genes distribution"), xlab = paste("reads"), ylab = paste ("number of samples"))

##############################################################################################################################################
# Filter low abundance genes *TO BE CHECKED BY DB*
#############################################################################################################################################

CN04_threshold = filter_taxa(CN04_data, function(x) sum(x > 100) > (0.10 *length(x)), TRUE)
CN04_threshold 

##ratio filtered reads/total reads
ratio <- sum(sample_sums(CN04_threshold))/sum(sample_sums(CN04_data))*100
ratio

##############################################################################################################################################
# Generate a new phyloseq object for calculation: the problem is that we do have differences in sequencing depth (i.e., some samples have more some have less)
# A solution could be "down-sampling" everything to the lowest count, a process called rarefying. 
# However, rarefying is a "contentious" matter in microbe science
# https://doi.org/10.1371/journal.pcbi.1003531 
#############################################################################################################################################

#remove samples with less than 400000
CN04_400K = prune_samples(sample_sums(CN04_threshold)>=400000, CN04_threshold)

CN04_400K_rare <- rarefy_even_depth(CN04_400K, 400000)
#ignore the warning, the object will be saved
CN04_400K_rare

#collapse at level 3 for congruence of metabolic pathways
CN04_400K_rare_L3 <- tax_glom(CN04_400K_rare, taxrank = "Level3")

#collapse at level 5 for congruence of metabolic pathways
CN04_400K_rare_L5 <- tax_glom(CN04_400K_rare, taxrank = "Level5")

#save the RDS object for the reproducibility of the code
#saveRDS(CN04_400K_rare , file = "CN04_400K_calculation.rds")
#saveRDS(CN04_400K_rare_L3 , file = "CN04_400K_L3_calculation.rds")
#saveRDS(CN04_400K_rare_L5 , file = "CN04_400K_L5_calculation.rds")

#quick test betadiversity
CN04_400K_rare.cap <- ordinate(CN04_400K_rare, "CAP", "bray", ~ Microhabitat)
#note the formula to specify what factor(s) look for

plot_ordination(CN04_400K_rare, CN04_400K_rare.cap, color="Genotype", shape = "Microhabitat")
#note on the proportion of variance explained

#Rhizosphere effect
BC_micrhoabitat <- phyloseq::distance(CN04_400K_rare, "bray")
BC_micrhoabitat

#here we can use a formula ANOVA-like to identify factors of interest
Stat_early <- adonis2(BC_micrhoabitat ~ Microhabitat, data= as.data.frame(as.matrix(sample_data(CN04_400K_rare))), permutations = 5000)
Stat_early
#Genotype effect ~67%

#Run a DESeq analysis to identify the top differentially accumulated genes

