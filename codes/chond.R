#Chond genes
#Chondrogenic genes
D0A <- hpc0_199a0[rownames(hpc0_199a0) %in% c("SOX9", "COL2A1", "ACAN"),]
D0B <- hpc0_199b0[rownames(hpc0_199b0) %in% c("SOX9", "COL2A1", "ACAN"),]
D1A <- hpc1_199a1[rownames(hpc1_199a1) %in% c("SOX9", "COL2A1", "ACAN"),]
D1B <- hpc1_199b1[rownames(hpc1_199b1) %in% c("SOX9", "COL2A1", "ACAN"),]
D0A <- D0A[,c(2, 5)]
D0B <- D0B[,c(2, 5)]
D1A <- D1A[,c(2, 5)]
D1B <- D1B[,c(2, 5)]
colnames(D0A)[[1]] <- "log2fc_199a_D0"
colnames(D0B)[[1]] <- "log2fc_199b_D0"
colnames(D1A)[[1]] <- "log2fc_199a_D1"
colnames(D1B)[[1]] <- "log2fc_199b_D1"
colnames(D0A)[[2]] <- "padj_199a_D0"
colnames(D0B)[[2]] <- "padj_199b_D0"
colnames(D1A)[[2]] <- "padj_199a_D1"
colnames(D1B)[[2]] <- "padj_199b_D1"
A <- cbind(D0A[1], D1A[1])
B <- cbind(D0B[1], D1B[1])
DM <- cbind(A[2], B[2])
library(gplots)
par(oma=c(3,3,0,3))
x11()
heatmap.2(as.matrix(DM), main=paste("Bio-markers Log2FC at Day 1"), Rowv = NULL,
      column = NULL, dendrogram="none", Colv = FALSE,
      scale="none", col= colorRampPalette(brewer.pal(11, "Spectral")),
      density.info="none", trace="none", key = TRUE)

