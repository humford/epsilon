#Load Files
setwd("~/Documents/LUAD")

BatchCorr <- as.matrix(read.table("LUAD_BC.txt", header = T, sep = "\t"))

exprMatrix <- as.matrix(read.table("LUAD_Processed"))



#Analysis of Moment Distributions
 
moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}


skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment((exprMatrix[i,]), 3)
}

corrskews <- NULL
for(i in 1:dim(BatchCorr)[[1]])
{
  corrskews[i] <-moment((BatchCorr[i,]), 3)
}

exprTable["LUAD", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["LUAD", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["BC LUAD", "Lower"] <- length(which(corrskews <= 0))/length(corrskews)
exprTable["BC LUAD", "Upper"] <- length(which(corrskews > 0))/length(corrskews)
exprTable["LUAD", "n"] <- dim(exprMatrix)[[2]]
exprTable["BC LUAD", "n"] <- dim(BatchCorr)[[2]]


old.par <- par(mfrow=c(1, 2))
plot(density(skews), main = "LUAD", xlab= "Skewness")
plot(density(corrskews), main = "LUAD Batch Corrected", xlab= "Skewness")
