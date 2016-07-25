setwd("~/Documents/LUSC")
exprMatrix <- as.matrix(read.table("LUSC_Processed"))

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
  

err <- function(Y, par)
{
  return(sum((Y-(par[1]*exp(-par[2] * (1:length(Y))) + par[3]))^2))
}

fit <- optim(par = c(0,0.5,0), err, Y = lower[-(1:20)])
plot(lower[-(1:20)])
lines(fit$par[1]*exp(-fit$par[2] * (1:length(lower[-(1:20)])) ) + fit$par[3])

fit

(lower[length(lower)] - fit$par[3])/fit$par[3] * 100