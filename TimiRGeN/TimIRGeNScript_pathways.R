# Load libraries. Using TimiRGeN v1.2.0 (BioManager 3.12).
setwd("~/Documents/Writing_Papers/ROCK1-SOX9/DataAnalysis/TimiRGeN/newanalysis/")

log <- readRDS("listofgenes.rds")

logun <- as.data.frame(unlist(log))

names(log)
Paths <- c(rep("VEGFA", 8), rep("GASTRIN", 3), rep("END OSS SKEL DYSP", 3),
           rep("EDO OSS", 3), rep("TGFB", 4), rep("EGFR", 4), 
           rep("FOCAL ADH", 8))
logun <- cbind(logun, Paths)
ogun <- logun[order(logun$`unlist(log)`),]
rownames(ogun) <- NULL

U <- unique(ogun$`unlist(log)`)

VEGFA <- c(1, 0, 1, 1, 0, 0, 1, 0, 1, 0,
           0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
           0, 0)

GASTRIN <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
             0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
             0, 0)

SKELDYSP <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 1)

ENDOSS <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             1, 1)

TGFB <- c(0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
          0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
          0, 0)

EGFR <- c(0, 0, 1, 1, 0, 0, 0, 0, 1, 0,
          0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
          0, 0)

FOCALADH <- c(0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
              1, 1, 0, 1, 1, 0, 1, 0, 0, 0,
              0, 0)

X <- cbind(VEGFA, GASTRIN, SKELDYSP, ENDOSS, TGFB, EGFR, FOCALADH)
rownames(X) <- U
X
######################################
#tileplot
library(ggplot2)
library(reshape2)
library(dplyr)
Y <- melt(X)
names(Y) <- c("Genes", "Pathways", "Value")


Y$Value <- sub(Y$Value, pattern = 0, replacement = "No")
Y$Value <- sub(Y$Value, pattern = 1, replacement = "Yes")


# Heatmap 
pdf(file = "tileplot.pdf")
ggplot(Y, aes(Pathways, Genes, fill= Value)) + 
    geom_tile()+
    scale_y_discrete(limits = rev)+
    scale_x_discrete(limits = rev)+
    theme_classic() +
    labs(title = "miR-199b-5p Targets", x = "Pathways", y = "Genes") 
dev.off()

