libs<-c("DESeq2", "tximport", "readr", "biomaRt", "miRNAtap", "miRNAtap.db",
        "org.Hs.eg.db", "topGO", "ggplot2", "biomaRt", "RColorBrewer", 
        "clusterProfiler", "dplyr", "ggrepel", "ggbiplot")
lapply(libs, require, character.only = TRUE)
#setwd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load functions.R
source("../codes/Functions.R")
load('../TPM/HAC_mir199_txiScaledTPM.RData')
txi2 <- txiScaledTPM
# load('HAC_mir199_txi.RData')
SampleTable2 <- read.table('../TPM/HAC_mir199_sampleTable.txt', row.names = 1, 
                           header = TRUE)
ddsTxi2 <- DESeqDataSetFromTximport(txi2,
                                    colData = SampleTable2,
                                    design = ~ Condition)
keep2 <- rowSums(counts(ddsTxi2)) >= 10
dds2 <- ddsTxi2[keep2,]

dds2$Condition <- factor(dds2$Condition, levels = c("hpC","hp199a","hp199b"))
dds2$Condition
#Deseq
dds2 <- DESeq(dds2)
################################################################################
# DE and get gene names
hpc_199a <- getDE(dds = dds2, Cont1 = "hp199a", Cont2 = "hpC")
hpc_199b <- getDE(dds = dds2, Cont1 = "hp199b", Cont2 = "hpC") 
################################################################################
# Which Target genes do we have confidence in, when both 199a and 199b are KO?
HAC_targets <- confidentTargets(DEres1 = hpc_199a, DEres2 = hpc_199b,
                                targets1 = miR199a5p_Targets,
                                targets2 = miR199b5p_Targets)