#Functions
# function to perform DE and retrieve gene names
getDE <- function(dds, Cont1 = "", Cont2= ""){
    res <- results(dds, contrast= c("Condition", Cont1, Cont2))
    res_B <- suppressMessages(as.data.frame(lfcShrink(dds=dds, 
                                                      contrast=c("Condition",
                                                                 Cont1,
                                                                 Cont2), 
                                                      res=res,
                                                      type = 'ashr')))
    res_B$Ens <- rownames(res_B)
    human <- suppressMessages(biomaRt::useEnsembl("ensembl", 
                                                  dataset = "hsapiens_gene_ensembl",
                                                  host = "useast.ensembl.org",
                                                  version = "94"))
    Gs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                filters = "ensembl_gene_id", values = rownames(res), 
                mart = human)
    X <- merge(Gs, res_B, by.x = "ensembl_gene_id", by.y = "Ens")
    O <- X[order(X$baseMean, decreasing = TRUE),]
    D <- O[! duplicated(O$external_gene_name),]
    rownames(D) <- D$ensembl_gene_id
    D$ensembl_gene_id <- NULL
    rownames(D) <- D$external_gene_name
    D$external_gene_name <- NULL
    return(D)
}
#################################################################################
# function to find microRNA targets
getTargets <- function(microRNA = ""){
    predictions = as.data.frame(getPredictedTargets("miR-199b-5p", 
                                                    species = 'hsa',
                                                    method = 'geom',
                                                    min_src = 2))
    predictions$ENTREZID <- rownames(predictions)
    mRNA_symb <- suppressMessages(suppressWarnings(clusterProfiler::bitr(geneID = predictions$ENTREZID, 
                                                                         fromType = "ENTREZID", 
                                                                         toType = "SYMBOL", 
                                                                         OrgDb = org.Hs.eg.db)))
    X <- merge(predictions, mRNA_symb, by ="ENTREZID")
    return(X)
}
################################################################################
# Function to extract genes which we are confident in from 2 DEseq2 results, and
# which genes are known targets.
confidentTargets <- function(DEres1, DEres2, confidence = 0.05, targets1, targets2){
    targets1_reduced <- targets1[which(targets1$rank_final >50),]
    targets2_reduced <- targets2[which(targets2$rank_final >50),]
    DEres1_C <- DEres1[which(DEres1$padj < confidence),]
    DEres2_C <- DEres2[which(DEres2$padj < confidence),]
    aTargets <- DEres1_C[which(rownames(DEres1_C) %in% targets1_reduced$SYMBOL == TRUE),]
    bTargets <- DEres2_C[which(rownames(DEres2_C) %in% targets2_reduced$SYMBOL == TRUE),]
    bothTargets1 <- aTargets[which(rownames(aTargets) %in% rownames(bTargets) == TRUE),]
    return(rownames(bothTargets1))
}
################################################################################
# Which genes are found in D1, D0 and in both?
listOfGenes <- function(D0, D1){
    OnlyDay1Targets <- D1[which(D1 %in% D0 == FALSE)]
    OnlyDay0Targets <- D0[which(D0 %in% D1 == FALSE)]
    BothTargets  <- D0[which(D0 %in% D1 == TRUE)]
    GeneNames <- c(BothTargets, OnlyDay1Targets, OnlyDay0Targets)
    return(GeneNames)
}
################################################################################
# Function to find which genes are increasing between D0 and D1.
getLogfc <- function(D0, D1, Targets){
    D0_T <- D0[which(rownames(D0) %in% Targets == TRUE),]
    D0_R <- data.frame("Names" = rownames(D0_T),
                       "D0_log2FC" = D0_T$log2FoldChange,
                       "D0_padj" = D0_T$padj)
    D1_T <- D1[which(rownames(D1) %in% Targets == TRUE),]
    D1_R <- data.frame("Names" = rownames(D1_T),
                       "D1" = D1_T$log2FoldChange,
                       "D1_padj" = D1_T$padj)
    DF <- merge(D0_R, D1_R, by = 'Names')
    rownames(DF) <- DF$Names
    DF$Names <- NULL
    return(DF)
}
################################################################################
# Create heatmaps of potential mRNA targets.
createHeatMap <- function(targets, microRNA){
    DM <- data.matrix(frame = targets)
    par(oma =c(5,0,0,0))
    heatmap.2(DM, main=paste(microRNA , "Targets" ), Rowv = NULL, 
              column = NULL, dendrogram="none",
              scale="none", col= colorRampPalette(brewer.pal(6, "YlGnBu")), 
              density.info="none", trace="none", key.title = 'LogFC', 
              key.ylab = 'LogFC')
}
################################################################################
# Make simple plots of genes of interest from time course RNAseq data
plotGene <- function(Data, GENE){
    x <- Data[which(Data$external_gene_name == GENE),]
    rownames(x) <- x$external_gene_name
    x$external_gene_name <- NULL
    plot(x)
}
################################################################################
# extract significant genes
lessthan <- function(X, val){
    Y <- X[which(X$padj < val),]
    Z <- Y[order(Y$log2FoldChange, decreasing = TRUE),]
    A <- Z$log2FoldChange
    names(A) <- rownames(Z)
    return(A)
}
################################################################################
#kegg preparation
keggprep <- function(X){
    df_X <- as.data.frame(X)
    df_X$SYMBOL <- rownames(df_X)
    entrez_X <- as.data.frame(bitr(geneID = df_X$SYMBOL, fromType = "SYMBOL",
                                   toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    merge_X <- merge(df_X, entrez_X, "SYMBOL")
    ordered_X <- merge_X[order(merge_X[,2], decreasing = TRUE),]
    X <- as.vector(ordered_X[,2])
    names(X) <- as.character(ordered_X$ENTREZID)
    return(X)
}
################################################################################
# Time for scatter plots
Time <- function(x, Time, KD){
    x <- x[,-c(3, 4)]
    x$Time <- Time
    x$Gene <- rownames(x)
    x$KD <- KD
    colnames(x)[1:5] <- c("Mean.Count", "Log2FC", "Adj.P.Value", "Time", "Gene")
    return(x)
}
################################################################################
#prepare data for scatter plot creation
PrepForScatt <- function(){
    TimeList <- rep(c("D0", "D1"), 2)
    KDList <- rep(c("miR-199a-5p KD", "miR-199b-5p KD"), 2)
    DElist <- list("hpc0_199a0" = hpc0_199a0,
                   "hpc1_199a1" = hpc1_199a1,
                   "hpc0_199b0" = hpc0_199b0,
                   "hpc1_199b1" = hpc1_199b1)
    Names <- rownames(MSC_199a_targets)
    EmptyList1 <- list()
    for (i in seq_along(1:4)) {
        EmptyList1[[i]] <- DElist[[i]][which(rownames(DElist[[i]]) %in% Names == TRUE),]
    }
    names(EmptyList1) <- names(DElist)
    EmptyList2 <- list()
    for (i in seq_along(1:4)) {
        EmptyList2[[i]] <- Time(x = EmptyList1[[i]], TimeList[[i]], KDList[[i]])
    }
    adata <- rbind(EmptyList2[[1]], EmptyList2[[2]])
    bdata <- rbind(EmptyList2[[3]], EmptyList2[[4]])
    data <- list("adata" = adata, "bdata" = bdata)
    return(data)
}
################################################################################
#Scatter/dot plot
Scatt <- function(L, titleString, fileString){
    if(L == "a") {
        data <- data[[1]]
    } else if(L == "b") {
        data <- data[[2]]
    }
    data_num<- data %>% 
        mutate(Adj.P.Value = cut(Adj.P.Value, breaks = c(0, 0.0001, 0.001, 
                                                         0.05, 100), 
                                 labels = c("<0.0001", "<0.001", "<0.01",
                                            ">0.05")))
    p <- ggplot(data_num, aes(x = Gene, y = Log2FC, group=Time)) + 
        geom_point(aes(color = Adj.P.Value, size = Mean.Count, shape = Time), 
                   alpha = 20)+
        coord_flip() +
        scale_x_discrete(limits = rev)+
        theme_classic() +
        labs(title = titleString, x = "Genes", y = "Log2FC") +
        theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
              axis.text.x = element_text(size = 8), 
              axis.text.y = element_text(size = 8), 
              axis.title.x = element_text(size = 10), 
              axis.title.y = element_text(size = 10), 
              legend.text = element_text(size = 8))
    jpeg(fileString, units="in", width=4.5, height=4, res=300)
    show(p)
    dev.off()
}
################################################################################
#Go data
getGo <- function(listedGenes){
    gse <- gseGO(geneList=listedGenes,
                 ont ="ALL",
                 keyType = "SYMBOL",
                 minGSSize = 3,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 OrgDb = org.Hs.eg.db,
                 pAdjustMethod = "BH")
    return(gse)
}
################################################################################
#Volcano
Volc <- function(de, fileString, titleString, t){
    X <-de[which(rownames(de) %in% rownames(t) == TRUE),]
    Y <-de[which(rownames(de) %in% rownames(t) == FALSE),]
    
    X$STATUS <- "invariant miR-199a/b-5p target"
    X$STATUS[X$log2FoldChange > 0 & X$padj < 0.05] <- "up-regulated miR-199a/b-5p target"
    X$STATUS[X$log2FoldChange < 0 & X$padj < 0.05] <- "down-regulated miR-199a/b-5p target"
    X$delabel[X$STATUS != "not/ invariant miR-199a/b-5p target"] <- rownames(X)[X$STATUS != "not/ invariant miR-199a/b-5p target"] 
    
    Y$STATUS <- "not miR-199a/b-5p target"
    Y$STATUS[Y$log2FoldChange > 0 & Y$padj < 0.05] <- "up-regulated gene"
    Y$STATUS[Y$log2FoldChange < 0 & Y$padj < 0.05] <- "down-regulated gene"
    Y$delabel <- NA
    Z <- rbind(X, Y)
    
    tiff(filename = paste0(fileString, ".tiff"), width = 5, height = 4,
         units = "in", res = 300)
    p <- ggplot(data=Z, aes(x=log2FoldChange, y=-log10(padj), 
                            col=STATUS, label=delabel)) +
        geom_point(
            data=Y, 
            aes(x=log2FoldChange, y=-log10(padj)), 
            size=2, alpha=0.25)+
        geom_point(
            data=X,
            aes(x=log2FoldChange, y=-log10(padj)), 
            size=2) +
        geom_text_repel(size = 4, position = "identity", force_pull=0, max.overlaps=18) +
        scale_color_manual(values=c("pink", "red", "black", "lightgrey", "lightblue", "blue"))+
        theme_classic() +
        labs(title = titleString, x = "Log2FC", y = "-log10(Adjusted P Values)")+
        theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
              axis.text.x = element_text(size = 15), 
              axis.text.y = element_text(size = 15), 
              axis.title.x = element_text(size = 15), 
              axis.title.y = element_text(size = 15), 
              legend.text = element_text(size = 15))+
        theme(legend.position = "none")
    show(p)
    dev.off()
}
######
#VOLC2
Volc2 <- function(de, fileString, titleString, t){
    X <-de[which(rownames(de) %in% rownames(t) == TRUE),]
    Y <-de[which(rownames(de) %in% rownames(t) == FALSE),]
    
    X$STATUS <- "invariant miR-199a/b-5p target"
    X$STATUS[X$log2FoldChange > 0 & X$padj < 0.05] <- "up-regulated miR-199a/b-5p target"
    X$STATUS[X$log2FoldChange < 0 & X$padj < 0.05] <- "down-regulated miR-199a/b-5p target"
    X$delabel[X$STATUS != "not/ invariant miR-199a/b-5p target"] <- rownames(X)[X$STATUS != "not/ invariant miR-199a/b-5p target"] 
    
    Y$STATUS <- "not miR-199a/b-5p target"
    Y$STATUS[Y$log2FoldChange > 0 & Y$padj < 0.05] <- "up-regulated gene"
    Y$STATUS[Y$log2FoldChange < 0 & Y$padj < 0.05] <- "down-regulated gene"
    Y$delabel <- NA
    Z <- rbind(X, Y)
    
    tiff(filename = paste0(fileString, ".tiff"), width = 5, height = 4,
         units = "in", res = 300)
    p <- ggplot(data=Z, aes(x=log2FoldChange, y=-log10(padj), 
                            col=STATUS, label=delabel)) +
        geom_point(
            data=Y, 
            aes(x=log2FoldChange, y=-log10(padj)), 
            size=2, alpha=0.25)+
        geom_point(
            data=X,
            aes(x=log2FoldChange, y=-log10(padj)), 
            size=2) +
        geom_text_repel(size = 4, position = "identity", force_pull=0, max.overlaps=18) +
        scale_color_manual(values=c("pink", "black", "lightgrey", "lightblue", "blue"))+
        theme_classic() +
        labs(title = titleString, x = "Log2FC", y = "-log10(Adjusted P Values)")+
        theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
              axis.text.x = element_text(size = 15), 
              axis.text.y = element_text(size = 15), 
              axis.title.x = element_text(size = 15), 
              axis.title.y = element_text(size = 15), 
              legend.text = element_text(size = 15))+
        theme(legend.position = "none")
    show(p)
    dev.off()
}
################################################################################
#pca plot
PCAplease <- function(data, fileString, titleString, colourCode){
    count <- (counts(data))
    count <- count[,order(colnames(count))]
    count.pca <- prcomp(t(count), center = TRUE,scale. = TRUE)
    p <- ggbiplot(count.pca, varname.size=0, var.axes = F, groups = colourCode,
                  scale=0, obs.scale = 1) +
        theme_classic() +
        labs(title = titleString) +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              axis.text.x = element_text(size = 8), 
              axis.text.y = element_text(size = 8), 
              axis.title.x = element_text(size = 8), 
              axis.title.y = element_text(size = 8), 
              legend.text = element_text(size = 10))+
        scale_color_discrete(name = '') +  
        geom_point(aes(colour=colourCode), size = 5)+
        theme(plot.margin=unit(c(0,0,0,0), 'cm'))
    tiff(fileString, units="in", width=5, height=4, res=300)
    show(p)
    dev.off()
}
################################################################################
gogo <- function(go, titlestring, namestring){
    go@result$log10padj <- -log10(go[,8])
    
    x <- which(go@result$Description == "chromatin assembly or disassembly")
    go@result$Description[[x]] <- "chr assembly/dissassembly"
    
    x <- which(go@result$Description == "protein-DNA complex subunit organization")
    go@result$Description[[x]] <- "prot-DNA complex subunit org"
    
    p <- dotplot(go, showCategory=5, split=".sign", x="log10padj") + 
        facet_grid(~.sign) +
        theme(strip.text = element_text(size=14),
              plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
              axis.text.x = element_text(size = 10), 
              axis.text.y = element_text(size = 18), 
              axis.title.x = element_text(size = 12), 
              axis.title.y = element_text(size = 15), 
              legend.text = element_text(size = 15))+
        labs(title = titlestring, x = "-log10(Adjusted P Values)", y = "GO terms")+
        guides(shape = guide_legend(order = 1))
    jpeg(namestring, units="in", width=10.5, height=6, res=300)
    show(p)
    dev.off()
}
################################################################################
#rename GO terms 
descriptions <- function(s4){
    s4@result$Description <- sub(s4@result$Description, 
                                 pattern = "extracellular matrix",
                                 replacement = "ECM")
    return(s4)
}
