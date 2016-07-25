moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

cancer.names <- c("Glioblastoma", "OV", "LumA", "AML_6891", "AML_15434")

exprTable <- matrix(1:3*(length(cancer.names) + 1), ncol = 3,nrow = length(cancer.names) + 1)
colnames(exprTable) <- c("Lower", "Upper", "n")
rownames(exprTable) <- c(cancer.names, "Control")

controlDiffTable <- list()  

#Analyze Control

setwd("~/Documents/Control")
controlMatrix <- as.matrix(read.table("Control_Processed"))
setwd("~/Documents/Graphs")

ControlSkews <- NULL

for(i in 1:dim(controlMatrix)[[1]])
{
  ControlSkews[i] <-moment((controlMatrix[i,]), 3)
}

names(ControlSkews) <- rownames(controlMatrix)

exprTable["Control", "Lower"] <- length(which(ControlSkews <= 0))/length(ControlSkews)
exprTable["Control", "Upper"] <- length(which(ControlSkews > 0))/length(ControlSkews)
exprTable["Control", "n"] <- dim(controlMatrix)[[2]]


#Analyze Cancers
  
for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancer, "ProcessedBC", sep = "_")))
  setwd("~/Documents/Graphs")
  
  skews <- NULL
  for(i in 1:dim(exprMatrix)[[1]])
  {
    skews[i] <-moment((exprMatrix[i,]), 3)
  }

  names(skews) <- rownames(exprMatrix)
  skews <- skews[!is.na(skews)]
  controlDiffTable[[cancer]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]

  exprTable[cancer, "Lower"] <- length(which(skews <= 0))/length(skews)
  exprTable[cancer, "Upper"] <- length(which(skews > 0))/length(skews)
  exprTable[cancer, "n"] <- dim(exprMatrix)[[2]]

}



old.par <- par(mfrow=c(1, 2), mai = c(1,1.15,1,1))
barplot(rbind(exprTable[,"Lower"], exprTable[,"Upper"]), main = "Microarray",
         beside = FALSE, names = rownames(exprTable), las = 1, horiz = TRUE)


      

      

