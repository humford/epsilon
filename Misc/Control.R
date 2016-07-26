#Load Data

setwd("~/Documents/Control")

exprMatrix <- as.matrix(read.table("Control_Processed"))

phenoData <- read.table("Control_Clinical")

sex <- lapply(phenoData$characteristics_ch1, function(y) strsplit(as.character(y), " ")[[1]][[3]])
tempcategory <- lapply(phenoData$characteristics_ch1, function(y) strsplit(as.character(y), " ")[[1]][[2]])

category <- NULL

for(i in 1:length(tempcategory))
{

  category <- append(category, tempcategory[i][[1]])
}



exprMatrix.male <- exprMatrix[, which(sex == "Male")]
exprMatrix.female <- exprMatrix[, which(sex == "Female")]

#Analysis of Moment Distributions


ControlSkews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  ControlSkews[i] <-moment((exprMatrix[i,]), 3)
}


plot(density(ControlSkews), main = "Control", xlab= "Skewness")


names(ControlSkews) <- rownames(exprMatrix)





exprTable["Control", "Lower"] <- length(which(ControlSkews <= 0))/length(ControlSkews)
exprTable["Control", "Upper"] <- length(which(ControlSkews > 0))/length(ControlSkews)
exprTable["Control", "n"] <- dim(exprMatrix)[[2]]

table <- NULL
old.par <- par(mfrow=c(2, 2))
for(cat in unique(category))
{
  skews <- NULL
  for(i in 1:dim(exprMatrix[,which(category == cat)])[[1]])
  {
    skews[i] <-moment((exprMatrix[i, which(category == cat)]), 3)
  }
  plot(density(skews), main = cat, xlab= "Skewness")
  table <- rbind(table, c(length(which(skews <= 0))/length(skews), length(which(skews > 0))/ length(skews)))
}

barplot(t(table), main="Gene Counts in Upper and Lower Modes",
  xlab="Group",
 	legend = colnames(table), beside = FALSE)