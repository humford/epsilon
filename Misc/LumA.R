setwd("~/Documents/BRCA")

exprMatrix <- as.matrix(read.table("BRCA_Processed"))

lumA <- as.matrix(read.table("LumAExpression.txt", header = T, sep = "\t"))
k <- which(is.na(lumA), arr.ind=TRUE)
rmeans <- rowMeans(lumA, na.rm=TRUE)
lumA[k] <- rmeans[k[,1]]

rownames(lumA) <- rownames(exprMatrix)
colnames(lumA) <- sapply(colnames(lumA), function(x) paste(strsplit(x, "[.]")[[1]][1:3], collapse = ".", sep = ""))

skews <- NULL
for(i in 1:dim(lumA)[[1]])
{
  skews[i] <-moment((lumA[i,]), 3)
}

names(skews) <- rownames(lumA)


controlDiffTable[["LumA"]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]

library(mclust)
out <- Mclust(controlDiffTable[["LumA"]], G = 3, modelNames = "V")

survivalTable <- read.table("survival_luma.txt", sep = "", header= TRUE)

survivalTable$Samples <- sapply(survivalTable$Samples, function(x) gsub("-", ".", x))

genes <- names(controlDiffTable[["LumA"]])[which(out$classification == 3)] 

library(genefilter)

#Dead: 1 Alive: 0

survivalTable$UPGenes <- sapply(colnames(lumA), function(p) sum(lumA[genes,p] > rowMeans(lumA)[genes] + 3*rowSds(lumA)[genes]))
survivalTable$DownGenes <- sapply(colnames(lumA), function(p) sum(lumA[genes,p] < rowMeans(lumA)[genes] - 3*rowSds(lumA)[genes]))

pvalues <- NULL

means <- rowMeans(lumA)

for(gene in rownames(lumA))
{
  upper <- which(lumA[gene, survivalTable$Samples] > means[gene])
  lower <- which(lumA[gene, survivalTable$Samples] < means[gene])
  m <- matrix(c(sum(survivalTable$Status[upper] == 1), sum(survivalTable$Status[upper] == 0), sum(survivalTable$Status[lower] == 1), sum(survivalTable$Status[lower] == 0)), ncol = 2)
  pvalues <- c(pvalues, fisher.test(m, alternative = "less")$p.value)
}