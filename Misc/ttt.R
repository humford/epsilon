moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")

exprTableRNA <- matrix(1:3*(length(cancer.names) + 1), ncol = 3,nrow = length(cancer.names) + 1)
colnames(exprTableRNA) <- c("Lower", "Upper", "n")
rownames(exprTableRNA) <- c(cancer.names, "Control")

controlDiffTable <- list()  

#Analyze Control

setwd("~/Documents/Control")
controlMatrix <- as.matrix(read.table("RNASeqControlProcessed.txt"))
setwd("~/Documents/Graphs")

ControlSkews <- NULL

for(i in 1:dim(controlMatrix)[[1]])
{
  ControlSkews[i] <-moment((controlMatrix[i,]), 3)
}

names(ControlSkews) <- rownames(controlMatrix)

exprTableRNA["Control", "Lower"] <- length(which(ControlSkews <= 0))/length(ControlSkews)
exprTableRNA["Control", "Upper"] <- length(which(ControlSkews > 0))/length(ControlSkews)
exprTableRNA["Control", "n"] <- dim(controlMatrix)[[2]]

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

  names(skews) <- rownames(exprMatrix)
  skews <- skews[!is.na(skews)]
  controlDiffTable[[cancer]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]

  exprTableRNA[cancer, "Lower"] <- length(which(skews <= 0))/length(skews)
  exprTableRNA[cancer, "Upper"] <- length(which(skews > 0))/length(skews)
  exprTableRNA[cancer, "n"] <- dim(exprMatrix)[[2]]

}





barplot(rbind(exprTableRNA[,"Lower"], exprTableRNA[,"Upper"]), main = "RNA-Seq",
         beside = FALSE, names = rownames(exprTableRNA), las = 1, horiz = TRUE)

      

      

