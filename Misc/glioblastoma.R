#Load Files

setwd("~/Documents/Glioblastoma")


BatchCorr <- read.table("glio_BC.txt", header = T, sep = ",")

rownames(BatchCorr) <- BatchCorr[,1]

BatchCorr <- BatchCorr[,-1]

BatchCorr <- as.matrix(BatchCorr)


exprMatrix <- as.matrix(read.table("Glioblastoma_Processed"))

#Analysis of Moment Distributions

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

names(skews) <- rownames(exprMatrix)
names(corrskews) <- rownames(BatchCorr)

controlDiffTable[["Glio"]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]
controlDiffTable[["BC Glio"]] <- corrskews[intersect(names(corrskews), names(ControlSkews))] - ControlSkews[intersect(names(corrskews), names(ControlSkews))]


exprTable["Glioblastoma", "Lower"] <- length(which(corrskews <= 0))/length(corrskews)
exprTable["Glioblastoma", "Upper"] <- length(which(corrskews > 0))/length(corrskews)
exprTable["Glioblastoma", "n"] <- dim(BatchCorr)[[2]]

old.par <- par(mfrow=c(1, 2))
plot(density(skews), main = "Glioblastoma", xlab= "Skewness")
plot(density(corrskews), main = "Glioblastoma Batch Corrected", xlab= "Skewness")





attach(mtcars)
par(mfrow=c(2,2))
plot(density(exprMatrix[which(skews < quantile(skews, 1:100/100)[1])[3], ]))
plot(density(exprMatrix[which(skews < quantile(skews, 1:100/100)[1])[4], ]))
plot(density(exprMatrix[which(skews < quantile(skews, 1:100/100)[1])[5], ]))
plot(density(exprMatrix[which(skews < quantile(skews, 1:100/100)[1])[6], ]))
