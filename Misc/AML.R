setwd("~/Documents/AML")

skewTable <- list()
#Analysis of Moment Distributions

exprMatrix <- as.matrix(read.table("AML_Processed"))


skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment((exprMatrix[i,]), 3)
}

 
names(skews) <- rownames(exprMatrix)
skewTable[["tgca"]] <- skews

plot(density(skews), main = "AML n = 196", xlab= "Skewness")




exprMatrix <- as.matrix(read.table("GSE17855_Processed"))

skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment((exprMatrix[i,]), 3)
}



names(skews) <- rownames(exprMatrix)
skewTable[["GSE17855"]] <- skews

plot(density(skews), main = "AML n = 237", xlab= "Skewness")



exprMatrix <- as.matrix(read.table("GSE6891_Processed"))

skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment((exprMatrix[i,]), 3)
}

names(skews) <- rownames(exprMatrix)
skews <- skews[!is.na(skews)]
controlDiffTable[["AML_6891"]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]

plot(density(skews), main = "AML n = 461", xlab= "Skewness")



exprTable["AML_1", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["AML_1", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["AML_1", "n"] <- dim(exprMatrix)[[2]]

names(skews) <- rownames(exprMatrix)
skewTable[["GSE6891"]] <- skews




exprMatrix <- as.matrix(read.table("GSE15434_Processed"))

skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment((exprMatrix[i,]), 3)
}

names(skews) <- rownames(exprMatrix)

controlDiffTable[["AML_15434"]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]

exprTable["AML_2", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["AML_2", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["AML_2", "n"] <- dim(exprMatrix)[[2]]


names(skews) <- rownames(exprMatrix)
skewTable[["GSE15434"]] <- skews


plot(density(skews), main = "AML n = 251", xlab= "Skewness")



