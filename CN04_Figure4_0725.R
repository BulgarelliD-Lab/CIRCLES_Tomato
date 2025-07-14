#####################################################################################
#Figure metagenomics X and supplementary Figure X
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations in PlaBase
#  Revision 12/24
#  amosca001@dundee.ac.uk
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
library("vegan")
library("ggplot2")
library ("DESeq2")
library ("PMCMR")

#############################################################
#set working directory: DB cpu
setwd("/cluster/db/R_shared/CN04_plabase/")
#set working directory: Alexandros
#TBD

#import the datasets
CN04 <-readRDS("CN04_400K_calculation.rds")
CN04
tax_table(CN04)


###################################################################################################
#arrange data for plotting
##################################################################################################

#Genotype
sample_data(CN04)$Genotype <- factor(sample_data(CN04)$Genotype, levels=c("Unplanted","Abundo", "SV5197TP"))

#color coding
DB_cols <- c("#999999", "#009E73", "#F179A7")

###################################################################################################
#Panels A & B functional richness & evenness
###################################################################################################
#extract data for statistical analysis
CN04_L5_alpha <-  estimate_richness(CN04, measures = c("Observed", "Shannon"))

#Genotype
design_Genotype <- as.data.frame(as.matrix(sample_data(CN04)))

#generate a new data frame richness genotype
design_L5_genotype <- cbind(design_Genotype, CN04_L5_alpha)
design_L5_genotype

#re-order the factors
design_L5_genotype$Genotype <- factor(design_L5_genotype$Genotype, levels=c("Unplanted","Abundo", "SV5197TP"))

#graphical output richness
p <- ggplot(design_L5_genotype, aes(x=Genotype, y=Observed, fill=Genotype)) + geom_boxplot()
p + geom_jitter(size=5,shape=21, position=position_jitter(0.1))+ scale_fill_manual(values = DB_cols)

#graphical output evenness
p <- ggplot(design_L5_genotype, aes(x=Genotype, y=Shannon, fill=Genotype)) + geom_boxplot()
p + geom_jitter(size=5,shape=21, position=position_jitter(0.1))+ scale_fill_manual(values = DB_cols)

#statistical analysis
#check the distribution of the data]
shapiro.test(design_L5_genotype$Observed)
shapiro.test(design_L5_genotype$Shannon)

#non parametric test
kruskal.test(Observed ~ Genotype, data = design_L5_genotype)
kruskal.test(Shannon ~ Genotype, data = design_L5_genotype)

#post-hoc
pairwise.wilcox.test(design_L5_genotype$Observed, design_L5_genotype$Genotype, p.adjust.method = "BH")
pairwise.wilcox.test(design_L5_genotype$Shannon, design_L5_genotype$Genotype, p.adjust.method = "BH")

###################################################################################################
#Panel C
###################################################################################################

CN04.ord <- ordinate(CN04, "NMDS", "bray")
plot_ordination(CN04, CN04.ord, type="samples", color="Genotype", shape="Microhabitat")

#ggplots function to increase effectivness of the visualisation
p = plot_ordination(CN04, CN04.ord, type="samples", color="Genotype", shape="Microhabitat") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("NMDS data Functional annotation, Bray distance")

#2D stress = 0.05
#Permanova
#this calculation should use the same distance used to build the graphical output
BC_functional  <- phyloseq::distance(CN04, "bray")
BC_functional

#here we can use a formula ANOVA-like to identify factors of interest
Stat_functional <- adonis2(BC_functional ~ Microhabitat, data= as.data.frame(as.matrix(sample_data(CN04))), permutations = 5000)
Stat_functional

###################################################################################################
#Panel D functional enrichment 
###################################################################################################

#############################################################
# DESeq Calculation to identify function sustaining the diversification
#prepare the material: countData the numbers; colDaata the mapping file
#############################################################

countData = as.data.frame(otu_table(CN04))
class(countData)
colnames(countData)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(CN04)))
class(colData)
rownames(colData)

#construct a DESeq dataset combining count data and sample information
#A DESeqDataSet object must have an associated design formula  The formula should be a tilde (???) followed by the variables of interest. In this case the column "Description" in the desing file depicts the variable of interest
CN04_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData, design= ~ Genotype)

#execute the differential count analysis with the function DESeq 
CN04_cds_test <- DESeq(CN04_cds, fitType="local", betaPrior=FALSE)

#Measure the Rhizosphere effect corrected for genotype
#Abundo
Bulk_Abundo  <- results(CN04_cds_test , contrast = c("Genotype",  "Unplanted", "Abundo"))
#inspect the results file
colnames(Bulk_Abundo)

# extract  genera whose adjusted p.value in a given comparison is below 0.01. 
Bulk_Abundo_FDR_001 <- Bulk_Abundo[(rownames(Bulk_Abundo)[which(Bulk_Abundo$padj <0.01)]), ]
dim(Bulk_Abundo_FDR_001)

#identify genera enriched in Abundo: first term of comparison, negative fold change
Abundo_enriched <-  Bulk_Abundo_FDR_001[(rownames(Bulk_Abundo_FDR_001)[which(Bulk_Abundo_FDR_001$log2FoldChange < 0)]), ]
dim(Abundo_enriched)

#write.table(Abundo_enriched, file = "Abundo_enriched_data_check.txt")

#SV
Bulk_SV  <- results(CN04_cds_test , contrast = c("Genotype",  "Unplanted", "SV5197TP"))
#inspect the results file
colnames(Bulk_SV)

# extract  genera whose adjusted p.value in a given comparison is below 0.01. 
Bulk_SV_FDR_001 <- Bulk_SV[(rownames(Bulk_SV)[which(Bulk_SV$padj <0.01)]), ]
dim(Bulk_SV_FDR_001)

#identify genera enriched in SV: first term of comparison, negative fold change
SV_enriched <-  Bulk_SV_FDR_001[(rownames(Bulk_SV_FDR_001)[which(Bulk_SV_FDR_001$log2FoldChange < 0)]), ]
dim(SV_enriched)

#write.table(SV_enriched, file = "SV_enriched_data_check.txt")

#################################
#Bar plot representation
################################

#identify functions commonly enriched in the two rhizosphere
Rhizo_function <- as.vector(intersect(rownames(Abundo_enriched), rownames(SV_enriched)))
Rhizo_function

#identify functions a) rhizosphere enriched and b) differentially regulated among genotypes
Genotype_function <- as.vector(intersect(Rhizo_function, rownames(SV_Abundo_FDR_001)))
length(Genotype_function)

#identify functions not part of the "core" rhizo enriched
Genotype_function_2 <- as.vector(setdiff(rownames(SV_Abundo_FDR_001), Rhizo_function))
length(Genotype_function_2)

#create a phyloseq object with only rhizosphere enriched function
CN04_rhizo_enriched <- prune_taxa(Rhizo_function, CN04)

#agglomerate data at level 4
CN04_rhizo_enriched_L4 <- tax_glom(CN04_rhizo_enriched, taxrank = "Level4")
CN04_rhizo_enriched_L4

#top 15 function
CN04_rhizo_enriched_top50 <- prune_taxa(names(sort(taxa_sums(CN04_rhizo_enriched_L4),TRUE)[1:15]), CN04_rhizo_enriched_L4)
CN04_rhizo_enriched_top50

#inspect the proportion of the top50 overall
sum(sample_sums(CN04_rhizo_enriched_top50))/sum(sample_sums(CN04_rhizo_enriched_L4))

#plot bar plots
p <- plot_bar(CN04_rhizo_enriched_top50, "Genotype", fill="Genotype", facet_grid=~Level4)
p + scale_fill_manual(values = DB_cols)


