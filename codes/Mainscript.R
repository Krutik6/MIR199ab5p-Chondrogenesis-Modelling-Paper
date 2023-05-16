#prep
BiocManager::install("topGO")
#load libraries
libs<-c("DESeq2", "tximport", "readr", "biomaRt", "miRNAtap", "miRNAtap.db",
        "org.Hs.eg.db", "topGO", "ggplot2", "biomaRt", "RColorBrewer", 
        "clusterProfiler", "dplyr", "ggrepel", "ggbiplot")
lapply(libs, require, character.only = TRUE)
#setwd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load functions.R
source("../codes/Functions.R")
###############################################################################
# File path and which txi to import
load("../TPM/MSC_mir199_txiScaledTPM.RData")
txi <- txiScaledTPM
#DESeq2 setup
SampleTable <- read.table('../TPM/MSC_mir199_sampleTable.txt', row.names = 1, 
                          header = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = SampleTable,
                                   design = ~ Condition)
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]
#Check for difference
par(mar=c(8,5,2,2))
#boxplot(counts(ddsTxi))
#boxplot(counts(dds))
#Set conditions
dds$Condition <- factor(dds$Condition, levels = c("hpC_Day_0",
                                                  "hpC_Day_1",
                                                  "hp199a_Day_0",
                                                  "hp199a_Day_1",
                                                   "hp199b_Day_0",
                                                  "hp199b_Day_1"))
#DESeq2
dds <- DESeq(dds)
resultsNames(dds)
#get DE
hpc0_199a0 <- getDE(dds = dds, Cont1 = "hp199a_Day_0", Cont2 = "hpC_Day_0")
hpc0_199b0 <- getDE(dds = dds, Cont1 = "hp199b_Day_0", Cont2 = "hpC_Day_0")
hpc1_199a1 <- getDE(dds = dds, Cont1 = "hp199a_Day_1", Cont2 = "hpC_Day_1")
hpc1_199b1 <- getDE(dds = dds, Cont1 = "hp199b_Day_1", Cont2 = "hpC_Day_1")
#how many are significantly expressed
hpc0_199a0_sig <- hpc0_199a0 %>% filter(padj < 0.05)
hpc0_199b0_sig <- hpc0_199b0 %>% filter(padj < 0.05)
hpc1_199a1_sig <- hpc1_199a1 %>% filter(padj < 0.05)
hpc1_199b1_sig <- hpc1_199b1 %>% filter(padj < 0.05)
#check for specific genes
hpc0_199a0_sig$names <- rownames(hpc0_199a0_sig)
hpc0_199b0_sig$names <- rownames(hpc0_199b0_sig)
hpc1_199a1_sig$names <- rownames(hpc1_199a1_sig)
hpc1_199b1_sig$names <- rownames(hpc1_199b1_sig)
#for day 0 which are the same
Id0 <- intersect(rownames(hpc0_199a0_sig), rownames(hpc0_199b0_sig))
Id1 <- intersect(rownames(hpc1_199a1_sig), rownames(hpc1_199b1_sig))
#data for enrichr
hpc1_199a1_sig100 <- hpc1_199a1_sig %>% 
  arrange(padj)
hpc1_199a1_sig100 <- hpc1_199a1_sig100[1:100,]
write.csv(hpc1_199a1_sig100, "../enrichr/199a_100.csv")

hpc1_199b1_sig100 <- hpc1_199b1_sig %>% 
  arrange(padj)
hpc1_199b1_sig100 <- hpc1_199b1_sig100[1:100,]
write.csv(hpc1_199b1_sig100, "../enrichr/199b_100.csv")
#which were up and down
hpc0_199a0_sig_pos <- hpc0_199a0_sig %>% filter(log2FoldChange >= 0)
hpc0_199a0_sig_neg <- hpc0_199a0_sig %>% filter(log2FoldChange <= 0)
hpc0_199b0_sig_pos <- hpc0_199b0_sig %>% filter(log2FoldChange >= 0)
hpc0_199b0_sig_neg <- hpc0_199b0_sig %>% filter(log2FoldChange <= 0)
hpc1_199a1_sig_pos <- hpc1_199a1_sig %>% filter(log2FoldChange >= 0)
hpc1_199a1_sig_neg <- hpc1_199a1_sig %>% filter(log2FoldChange <= 0)
hpc1_199b1_sig_pos <- hpc1_199b1_sig %>% filter(log2FoldChange >= 0)
hpc1_199b1_sig_neg <- hpc1_199b1_sig %>% filter(log2FoldChange <= 0)
#contrast control D1 against control D0
hpc1_hpc0 <- getDE(dds = dds, Cont1 = "hpC_Day_1", Cont2 = "hpC_Day_0")
hpc1_hpc0_sig <- hpc1_hpc0 %>% filter(padj < 0.05)
hpc1_hpc0_sig$name <- rownames(hpc1_hpc0_sig)
hpc1_hpc0_sig_pos <- hpc1_hpc0_sig %>% filter(log2FoldChange >= 0)
hpc1_hpc0_sig_neg <- hpc1_hpc0_sig %>% filter(log2FoldChange <= 0)
hpc1_hpc0_sig100 <- hpc1_hpc0_sig %>% 
  arrange(padj)
hpc1_hpc0_sig100 <- hpc1_hpc0_sig100[1:100,]
write.csv(hpc1_hpc0_sig100, "../enrichr/HPC_100.csv")
################################################################################
#Get miRNA targets for 199a and 199b
miR199b5p_Targets <- getTargets(microRNA = "miR-199b-5p")
miR199a5p_Targets <- getTargets(microRNA = "miR-199a-5p")
################################################################################
# Which genes in both 199a and 199b do we have confidence in? 
Day1Targets <- confidentTargets(DEres1 = hpc1_199a1, DEres2 = hpc1_199b1, 
                                targets1 = miR199a5p_Targets,
                                targets2 = miR199b5p_Targets)
Day0Targets <- confidentTargets(DEres1 = hpc0_199a0, DEres2 = hpc0_199b0, 
                                targets1 = miR199b5p_Targets,
                                targets2 = miR199a5p_Targets)
################################################################################
# And which ones are found in D1, D0 and in both inhibitions?
GeneNames <- listOfGenes(D0 = Day0Targets, D1 = Day1Targets)
################################################################################
MSC_199a_targets <- getLogfc(D0 = hpc0_199a0, D1 = hpc1_199a1,
                             Targets = GeneNames)
MSC_199b_targets <- getLogfc(D0 = hpc0_199b0, D1 = hpc1_199b1,
                             Targets = GeneNames)
################################################################################
# write.table(x = MSC_miR199a_targets, file = 'MSC_199a_targets.txt', quote = FALSE)
# write.table(x = MSC_miR199b_targets, file = 'MSC_199b_targets.txt', quote = FALSE)
################################################################################
# Scatterplots/ dotplots
data <- PrepForScatt()
Scatt(L = "a", titleString = "miR-199a-5p KD", fileString = "../scatter/199a.jpeg")
Scatt(L = "b", titleString = "miR-199b-5p KD", fileString = "../scatter/199b.jpeg")
###############################################################################
# #functional analysis
a0 <- lessthan(X = hpc0_199a0, val = 0.1)
a1 <- lessthan(X = hpc1_199a1, val = 0.05)
b0 <- lessthan(X = hpc0_199b0, val = 0.05)
b1 <- lessthan(X = hpc1_199b1, val = 0.05)
###############
#GSE_GO
a0_go <- getGo(listedGenes = a0)
b0_go <- getGo(listedGenes = b0)
a1_go <- getGo(listedGenes = a1)
b1_go <- getGo(listedGenes = b1)
##################
#rename GO terms
# a0_go <- descriptions(s4 = a0_go)
# b0_go <- descriptions(s4 = b0_go)
# a1_go <- descriptions(s4 = a1_go)
# b1_go <- descriptions(s4 = b1_go)
##################
#plot go enrichment
#setwd("../GO/")
#gogo(go = a0_go, titlestring = "hpmiR-199a/hpCon_D0 GO terms", namestring = "a0.jpg")
#gogo(go = a1_go, titlestring = "hpmiR-199a/hpCon_D1 GO terms", namestring = "a1.jpg")
#gogo(go = b0_go, titlestring = "hpmiR-199b/hpCon_D0 GO terms", namestring = "b0.jpg")
#gogo(go = b1_go, titlestring = "hpmiR-199b/hpCon_D1 GO terms", namestring = "b1.jpg")

################################################################################
#PCA
MSC <- (counts(dds))
MSC <- MSC[,order(colnames(MSC))]
#HAC <- (counts(dds2))
#HAC <- HAC[,order(colnames(HAC))]

# MSC_cols <- rep(c("hpmiR-199a_day0", "hpmiR-199a_day1", 
#                   "hpmiR-199b_day0", "hpmiR-199b_day1", 
#                   "hpCon_day0", "hpCon_day1"), each=3)
#HAC_cols <- rep(c("199a", "199b", "hpc"), each=4)
setwd("../pca/")
# PCAplease(data = dds, fileString = "MSCPCA.tiff", 
#           titleString = "PCA plot of Normalized Counts", 
#           colourCode = MSC_cols)
#PCAplease(data = dds2, fileString = "HACPCA.tiff", 
#          titleString = "PCA plot of Normalized HAC Counts", 
#          colourCode = HAC_cols)
################################################################################
#volcano
setwd("../volcano/")
# Volc2(de = hpc0_199a0, fileString = "199a0_volc", 
#      titleString = "hpmiR-199a/hpCon_D0 Volcano", t = MSC_199a_targets)
# Volc(de = hpc1_199a1, fileString = "199a1_volc", 
#      titleString = "hpmiR-199a/hpCon_D1 Volcano", t = MSC_199a_targets)
# Volc2(de = hpc0_199b0, fileString = "199b0_volc", 
#      titleString = "hpmiR-199b/hpCon_D0 Volcano", t = MSC_199b_targets)
# Volc(de = hpc1_199b1, fileString = "199b1_volc", 
#      titleString = "hpmiR-199b/hpCon_D1 Volcano", t = MSC_199b_targets)
#     
################################################################################
#scatter
setwd()