#Load Files

setwd("~/Documents/OV")

BatchCorr <- read.table("ovarian_BC.txt", header = T, sep = ",")

rownames(BatchCorr) <- BatchCorr[,1]

BatchCorr <- BatchCorr[,-1]

BatchCorr <- as.matrix(BatchCorr)

exprMatrix <- as.matrix(read.table("OV_Processed"))








#Analysis of Moment Distributions

skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment((exprMatrix[i, ]), 3)
}
summary(skews)


corrskews <- NULL
for(i in 1:dim(BatchCorr)[[1]])
{
  corrskews[i] <-moment((BatchCorr[i,]), 3)
}

names(skews) <- rownames(exprMatrix)
names(corrskews) <- rownames(BatchCorr)

controlDiffTable[["OV"]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]
controlDiffTable[["BC OV"]] <- corrskews[intersect(names(corrskews), names(ControlSkews))] - ControlSkews[intersect(names(corrskews), names(ControlSkews))]

exprTable["Ovarian", "Lower"] <- length(which(corrskews <= 0))/length(corrskews)
exprTable["Ovarian", "Upper"] <- length(which(corrskews > 0))/length(corrskews)
exprTable["Ovarian", "n"] <- dim(BatchCorr)[[2]]

old.par <- par(mfrow=c(1, 2))
plot(density(skews), main = "OV", xlab= "Skewness")
plot(density(corrskews), main = "OV Batch Corrected", xlab= "Skewness")





