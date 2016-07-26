moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

err <- function(Y, par)
{
  return(sum((Y-(par[1]*exp(-par[2] * (1:length(Y))) + par[3]))^2))
}


cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")

exprTable <- matrix(1:3*(length(cancer.names) + 1), ncol = 3,nrow = length(cancer.names) + 1)
colnames(exprTable) <- c("Lower", "Upper", "n")
rownames(exprTable) <- c(cancer.names, "Control")

controlDiffTable <- list()  


errors <- NULL
converge <- NULL
ControlSkews <- NULL

#Analyze Control

setwd("~/Documents/Control")
controlMatrix <- as.matrix(read.table("RNASeqControlProcessed.txt"))
setwd("~/Documents/Graphs")


for(i in 1:dim(controlMatrix)[[1]])
{
  ControlSkews[i] <-moment((controlMatrix[i,]), 3)
}

names(ControlSkews) <- rownames(controlMatrix)

exprTable["Control", "Lower"] <- length(which(ControlSkews <= 0))/length(ControlSkews)
exprTable["Control", "Upper"] <- length(which(ControlSkews > 0))/length(ControlSkews)
exprTable["Control", "n"] <- dim(controlMatrix)[[2]]

pdf("RNASeq_Control_Skews.pdf")
plot(density(ControlSkews), main = "Control", xlab= "Skewness")
dev.off()

lower <- NULL

for(a in 1:dim(controlMatrix)[[2]])
{
  skews <- NULL
  for(i in 1:dim(controlMatrix)[[1]])
  {
    skews[i] <- moment(sample(controlMatrix[i,], a), 3)
  }
  lower[a] <- length(which(skews < 0))
}

lower[1:3] <- lower[4]

lower <- lower/length(skews)

range <- 0.05 *abs(max(lower) - min(lower))
accept <- which(abs(lower  - lower[length(lower)])  <= range)
reject <- which(abs(lower  - lower[length(lower)]) > range)
  
pdf("RNASeq_Control_Sampling.pdf")

plot(reject, lower[reject],col = "red", main= "RNASeq Control Gene Spliting by Sample Size", xlab = "Samples", ylab = "Fraction of Genes in Lower Mode", 
                                            xlim = c(1, length(lower)), ylim = c(min(lower) - range, max(lower) + range))
points(accept, lower[accept], col = "black")

abline(a = (lower[length(lower)] - range), b = 0)
abline(a = (lower[length(lower)] + range), b = 0)
text(min(accept), lower[min(accept)] + 2*range, as.character(min(accept)))
dev.off()
  
fit <- optim(par = c(0,0.5,0), err, Y = lower[-(1:20)])
errors["Control"] <- (lower[length(lower)] - fit$par[3])/fit$par[3] * 100
converge["Control"] <- 1/exp(fit$par[2])
  
  

#Analyze Cancers
  
for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
  setwd("~/Documents/Graphs")
  
  skews <- NULL
  for(i in 1:dim(exprMatrix)[[1]])
  {
    skews[i] <-moment((exprMatrix[i,]), 3)
  }

  skews <- skews[!is.na(skews)]
  names(skews) <- rownames(exprMatrix)
  controlDiffTable[[cancer]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]

  exprTable[cancer, "Lower"] <- length(which(skews <= 0))/length(skews)
  exprTable[cancer, "Upper"] <- length(which(skews > 0))/length(skews)
  exprTable[cancer, "n"] <- dim(exprMatrix)[[2]]

  pdf(paste(cancer, "_Skews.pdf", sep = ""))
  plot(density(skews), main = cancer, xlab= "Skewness")
  dev.off()
  
  lower <- NULL

  for(a in 1:dim(exprMatrix)[[2]])
  {
    skews <- NULL
    for(i in 1:dim(exprMatrix)[[1]])
    {
      skews[i] <- moment(sample(exprMatrix[i,], a), 3)
    }
    lower[a] <- length(which(skews < 0))
  }

  lower[1:3] <- lower[4]

  lower <- lower/length(skews)

  range <- 0.05 *abs(max(lower) - min(lower))
  accept <- which(abs(lower  - lower[length(lower)])  <= range)
  reject <- which(abs(lower  - lower[length(lower)]) > range)
  
  pdf(paste(cancer, "_Sampling.pdf", sep = ""))

  plot(reject, lower[reject],col = "red", main= paste(cancer, " Gene Spliting by Sample Size", sep = ""), xlab = "Samples", ylab = "Fraction of Genes in Lower Mode", 
                                            xlim = c(1, length(lower)), ylim = c(min(lower) - range, max(lower) + range))
  points(accept, lower[accept], col = "black")

  abline(a = (lower[length(lower)] - range), b = 0)
  abline(a = (lower[length(lower)] + range), b = 0)
  text(min(accept), lower[min(accept)] + 2*range, as.character(min(accept)))
  dev.off()
  
  fit <- optim(par = c(0,0.5,0), err, Y = lower[-(1:20)])
  errors[cancer] <- (lower[length(lower)] - fit$par[3])/fit$par[3] * 100
  converge[cancer] <- 1/exp(fit$par[2])
}

convergenceTable <- data.frame(errors, converge)
names(convergenceTable) <- c("Sample Name", "Limit Error", "Rate of Convergence")
rownames(convergenceTable) <- c("Control", cancer.names)


setwd("~/Documents/Graphs")
pdf("RNASeq_Splitting_Comparison.eps")
par(mai = c(1,2,1,1))
barplot(rbind(exprTable[,"Lower"], exprTable[,"Upper"]), main="Percentage of Genes in Upper and Lower Modes",
        legend = colnames(exprTable)[1:2], beside = FALSE, names = paste(rownames(exprTable), " (n = ", exprTable[,"n"], ")", sep = ""), las = 1, horiz = FALSE)
dev.off()

setwd("~/Documents/Tables")
write.table(convergenceTable, "RNASeq_Convergence")


library(GOstats)
library(KEGG.db)
library(org.Hs.eg.db)
library(mclust)

GeneLists <- NULL

for(cancer in names(controlDiffTable))
{
  GeneLists[[cancer]] <- list() 
  
  out <- Mclust(controlDiffTable[[cancer]], G = 3, modelNames = "V")
    
  GeneLists[[cancer]][["UP"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 3)] 
  GeneLists[[cancer]][["MIDDLE"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 2)] 
  GeneLists[[cancer]][["DOWN"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 1)]  
  GeneLists[[cancer]][["ALL"]] <- names(controlDiffTable[[cancer]])
}

GeneLists[["Intersect"]] <- GeneLists[["OV"]]

for(cancer in names(controlDiffTable))
{
  GeneLists[["Intersect"]][["UP"]] <- intersect(GeneLists[["Intersect"]][["UP"]], GeneLists[[cancer]][["UP"]])
  GeneLists[["Intersect"]][["MIDDLE"]] <- intersect(GeneLists[["Intersect"]][["MIDDLE"]], GeneLists[[cancer]][["MIDDLE"]]) 
  GeneLists[["Intersect"]][["DOWN"]] <- intersect(GeneLists[["Intersect"]][["DOWN"]], GeneLists[[cancer]][["DOWN"]]) 
  GeneLists[["Intersect"]][["ALL"]] <- intersect(GeneLists[["Intersect"]][["ALL"]], GeneLists[[cancer]][["ALL"]])
}

hgCutoff <- 1
pValueCutoff <- 0.01
GeneStatsResults <- NULL
AdjustedResults <- NULL

for(cancer in names(GeneLists))
{
  GeneStatsResults[[cancer]] <- list()
  AdjustedResults[[cancer]] <- list()
  
  universelist <- GeneLists[[cancer]][["ALL"]]
  
  universe <- unlist(mget(intersect(universelist, ls(org.Hs.egSYMBOL2EG)), org.Hs.egSYMBOL2EG))

  universe <- unique(universe[!is.na(universe)])
  
  for(group in c("UP", "MIDDLE", "DOWN"))
  {
    
    GeneStatsResults[[cancer]][[group]] <- list()
    AdjustedResults[[cancer]][[group]] <- list()

    genelist <- GeneLists[[cancer]][[group]]
      
    genelist.eg <- unlist(mget(intersect(genelist, ls(org.Hs.egSYMBOL2EG)), org.Hs.egSYMBOL2EG))
  
    genelist.eg <- unique(genelist.eg[!is.na(genelist.eg)])
    
    for(marker in c("BP", "CC", "MF", "KEGG"))
    {
    
      if(marker == "KEGG")
      {
          params <- new("KEGGHyperGParams", geneIds=genelist.eg, universeGeneIds=universe, 
                annotation="org.Hs.eg.db", pvalueCutoff=hgCutoff, testDirection="over") 
      
      }
      else
      {
          params <- new("GOHyperGParams", geneIds=genelist.eg, universeGeneIds=universe, 
                annotation="org.Hs.eg.db", ontology=marker, pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")               
      }
      hgOver <- hyperGTest(params)
 
      GeneStatsResults[[cancer]][[group]][[marker]] <- summary(hgOver)
      adjustedPvalues <- p.adjust(GeneStatsResults[[cancer]][[group]][[marker]]$Pvalue, method = "BH")
      AdjustedResults[[cancer]][[group]][[marker]] <- GeneStatsResults[[cancer]][[group]][[marker]][which(adjustedPvalues < pValueCutoff),]
      AdjustedResults[[cancer]][[group]][[marker]]$Pvalues <- adjustedPvalues[which(adjustedPvalues < pValueCutoff)]
    }
  }
}

for(marker in c("BP", "CC", "MF", "KEGG"))
{
  for(group in c("UP", "MIDDLE", "DOWN"))
  {
    for(cancer in cancer.names) 
    {
      stats <- AdjustedResults[[cancer]][[group]][[marker]]
      if(dim(stats)[1] > 0)
      {
        write.table(t(c(marker, group, cancer)), "RNASeq_toControl_TMP.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(stats[1:min(length(stats[[1]]), 10), c(1,8,5,6,7)], "RNASeq_toControl_TMP.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        
        write.table(t(c(marker, group, cancer)), "RNASeq_toControl_TMP (NO PVAL).txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(stats[1:min(length(stats[[1]]), 10), c(1,6,7)], "RNASeq_toControl_TMP (NO PVAL).txt", append = TRUE, row.names = FALSE, col.names = FALSE) 
        
        if(marker != "KEGG")
        {
          setwd("~/Documents/Tables/GOSlim")
          write.table(stats[, c(1,8)], paste(group, cancer, "to_Control_TPM.txt", sep = "_"), append = TRUE, row.names = FALSE, col.names = FALSE)
        }
      }
    }
  }
}



      

      

