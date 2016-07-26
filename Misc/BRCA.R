setwd("~/Documents/BRCA")

#Import mRNA Expression Data

her2 <- as.matrix(read.table("Her2Expression.txt", header = T, sep = "\t"))
k <- which(is.na(her2), arr.ind=TRUE)
rmeans <- rowMeans(her2, na.rm=TRUE)
her2[k] <- rmeans[k[,1]]

basal <- as.matrix(read.table("BasalExpression.txt", header = T, sep = "\t"))
k <- which(is.na(basal), arr.ind=TRUE)
rmeans <- rowMeans(basal, na.rm=TRUE)
basal[k] <- rmeans[k[,1]]

lumA <- as.matrix(read.table("LumAExpression.txt", header = T, sep = "\t"))
k <- which(is.na(lumA), arr.ind=TRUE)
rmeans <- rowMeans(lumA, na.rm=TRUE)
lumA[k] <- rmeans[k[,1]]

lumB <- as.matrix(read.table("LumBExpression.txt", header = T, sep = "\t"))
k <- which(is.na(lumB), arr.ind=TRUE)
rmeans <- rowMeans(lumB, na.rm=TRUE)
lumB[k] <- rmeans[k[,1]]


subcategory_map <- read.table("classification.csv", header = T, sep = ",")

subcategory_map <- lapply(subcategory_map, function(x) gsub("-", ".", x, fixed = TRUE))

exprMatrix <- as.matrix(read.table("BRCA_Processed"))




#Analysis of Moment Distributions


skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment((exprMatrix[i,]), 3)
}



plot(density(skews), main = "BRCA", xlab= "Skewness")





par(mfrow=c(2,2),oma = c(0, 0, 2, 0))


skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment(exprMatrix[i, intersect(colnames(exprMatrix), subcategory_map$Her2)], 3)
}

exprTable["Her2", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["Her2", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["Her2", "n"] <- length(intersect(colnames(exprMatrix), subcategory_map$Her2))

plot(density(skews), main = "Her2", xlab= "Skewness")

skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment(exprMatrix[i, intersect(colnames(exprMatrix), subcategory_map$Basal)], 3)
}

exprTable["Basal", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["Basal", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["Basal", "n"] <- length(intersect(colnames(exprMatrix), subcategory_map$Basal))

plot(density(skews), main = "Basal", xlab= "Skewness")

skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment(exprMatrix[i, intersect(colnames(exprMatrix), subcategory_map$LumA)], 3)
}

exprTable["LumA", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["LumA", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["LumA", "n"] <- length(intersect(colnames(exprMatrix), subcategory_map$LumA))

plot(density(skews), main = "LumA", xlab= "Skewness")

skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment(exprMatrix[i, intersect(colnames(exprMatrix), subcategory_map$LumB)], 3)
}

exprTable["LumB", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["LumB", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["LumB", "n"] <- length(intersect(colnames(exprMatrix), subcategory_map$LumB))

plot(density(skews), main = "LumB", xlab= "Skewness")

mtext("Batch Uncorrected BRCA", outer = TRUE, cex = 1.5)







par(mfrow=c(2,2),oma = c(0, 0, 2, 0))


skews <- NULL
for(i in 1:dim(her2)[[1]])
{
  skews[i] <-moment((her2[i,]), 3)
}


exprTable["BC Her2", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["BC Her2", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["BC Her2", "n"] <- dim(her2)[[2]]

plot(density(skews), main = "Her2", xlab= "Skewness")

skews <- NULL
for(i in 1:dim(basal)[[1]])
{
  skews[i] <-moment((basal[i,]), 3)
}

exprTable["BC Basal", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["BC Basal", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["BC Basal", "n"] <- dim(basal)[[2]]

plot(density(skews), main = "Basal", xlab= "Skewness")

skews <- NULL
for(i in 1:dim(lumA)[[1]])
{
  skews[i] <-moment((lumA[i,]), 3)
}

exprTable["BC LumA", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["BC LumA", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["BC LumA", "n"] <- dim(lumA)[[2]]

plot(density(skews), main = "LumA", xlab= "Skewness")

skews <- NULL
for(i in 1:dim(lumB)[[1]])
{
  skews[i] <-moment((lumB[i,]), 3)
}

exprTable["BC LumB", "Lower"] <- length(which(skews <= 0))/length(skews)
exprTable["BC LumB", "Upper"] <- length(which(skews > 0))/length(skews)
exprTable["BC LumB", "n"] <- dim(lumB)[[2]]

plot(density(skews), main = "LumB", xlab= "Skewness")

mtext("Batch Corrected BRCA", outer = TRUE, cex = 1.5)
