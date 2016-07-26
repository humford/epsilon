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



errors <- NULL
converge <- NULL


#Analyze Control

setwd("~/Documents/Control")
controlMatrix <- as.matrix(read.table("RNASeqControlProcessed.txt"))
setwd("~/Documents/")



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

  
fit <- optim(par = c(0,0.5,0), err, Y = lower[-(1:20)])
errors["Control"] <- (lower[length(lower)] - fit$par[3])/fit$par[3] * 100
converge["Control"] <- 1/exp(fit$par[2])
  
  

#Analyze Cancers
  
for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
  setwd("~/Documents/")
  
  
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


  
  fit <- optim(par = c(0,0.5,0), err, Y = lower[-(1:20)])
  errors[cancer] <- (lower[length(lower)] - fit$par[3])/fit$par[3] * 100
  converge[cancer] <- 1/exp(fit$par[2])
}

convergenceTable <- data.frame(errors, converge)
names(convergenceTable) <- c("Limit Error", "Rate of Convergence")
rownames(convergenceTable) <- c("Control", cancer.names)

setwd("~/Documents/Tables")
write.table(convergenceTable, "RNASeq_Convergence")

      

