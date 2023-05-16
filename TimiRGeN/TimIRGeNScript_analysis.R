# Load libraries. Using TimiRGeN v1.0.4 (BioManager 3.12).
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load functions.R
source("function.R")
#Libs
library(TimiRGeN)
library(org.Hs.eg.db)
# Load data
miR <- read.table("miRNA.txt", header = TRUE, sep = ',', row.names = 1)
mRNA <- read.table("mRNA.txt", header = TRUE, sep = ',', row.names = 1)
# Create MAE, filter for significance per time point and separate data 
# into nested dataframes by time point. Using combined analysis.
MAE <- startObject(miR = miR, mRNA = mRNA)
MAE <- getIdsMir(MAE = MAE, miR = assay(MAE, 1), orgDB = org.Hs.eg.db, 
                 miRPrefix = "hsa")
MAE <- getIdsMrna(MAE = MAE, mRNA = assay(MAE, 2), mirror = 'www', 
                  species = "hsapiens")
MAE <- combineGenes(MAE = MAE, miR_data = assay(MAE, 1), 
                    mRNA_data = assay(MAE, 2))
MAE <- genesList(MAE = MAE, method = 'c',
                 genetic_data = assay(MAE, 9),
                 timeString = "D")
MAE <- significantVals(MAE = MAE, method = 'c',
                       geneList =metadata(MAE)[[1]],
                       maxVal = 0.05, stringVal = "adjPval")
MAE <- addIds(MAE = MAE, method = 'c',
              filtered_genelist = metadata(MAE)[[2]],
              miR_IDs = assay(MAE, 3),
              mRNA_IDs = assay(MAE, 7))
MAE <- eNames(MAE = MAE, method = 'c',
              gene_IDs = metadata(MAE)[[3]])
# start MAE2
MAE2 <- MultiAssayExperiment()
MAE2 <- dloadGmt(MAE = MAE2, species = "Homo sapiens")
probes <- read.csv("platforms/microarray_genes.csv")
probes$Entrez_Gene_ID <- as.character(probes$Entrez_Gene_ID)
# Overrepresentation analysis
MAE2 <- enrichWiki(MAE = MAE2, method = "c",
                   ID_list = metadata(MAE)[[4]],
                   orgDB = org.Hs.eg.db,
                   path_gene = assay(MAE2, 1),
                   path_name = assay(MAE2, 2),
                   ID = "ENTREZID",
                   universe = probes$Entrez_Gene_ID)
########################
plotEnrich(MAE = MAE2, select = 1, day = "D1")
plotEnrich(MAE = MAE2, select = 2, day = "D3")
plotEnrich(MAE = MAE2, select = 3, day = "D6")
plotEnrich(MAE = MAE2, select = 4, day = "D10")
plotEnrich(MAE = MAE2, select = 5, day = "D14")
########################

# create dynamic dataframes
MAE2 <- diffExpressRes(MAE = MAE2, df = assay(MAE, 1), dataType = "log2fc",
                       genes_ID = assay(MAE, 3), idColumn = "GENENAME",
                       name = "miR_log2fc")
MAE2 <- diffExpressRes(MAE = MAE2, df = assay(MAE, 2), dataType = "log2fc",
                       genes_ID = assay(MAE, 7), idColumn = "GENENAME",
                       name = "mRNA_log2fc")
options(timeout=180)
MAE2 <- dloadTargetscan(MAE = MAE2, species = "hsa")
MAE2 <- dloadMirdb(MAE = MAE2, species = "hsa", orgDB = org.Hs.eg.db)
MAE2 <- dloadMirtarbase(MAE = MAE2, species = "hsa")

paths <- list("VEGFA-VEGFR2 signaling pathway", 
              "Endochondral ossification with skeletal dysplasias",
              "Endochondral ossification",
              "EGF/EGFR signaling pathway",
              "Metabolic reprogramming in colon cancer",
              "TGF-beta signaling pathway",
              "Gastrin signaling pathway",
              "Adipogenesis",
              "Clear cell renal cell carcinoma pathways")

names <- list("VEGFA-VEGFR2", "EOSD", "EndOss", "EFG_EGFR", "MRCC", "TGFB",
              "GASTRIN", "Adipogenesis", "CCRCCP")


for (i in seq_along(paths)){
    M <- Pathways(selectedpath = paths[[i]])
    plotPaths(M = M, imagename =  names[[i]])
}


