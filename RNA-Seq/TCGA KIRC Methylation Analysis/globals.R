#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#GLOBALS

library(mclust)

cutoff <- 0.01

moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

splitter <- function(exprs)
{
  mixmdl <- Mclust(exprs)
  mean <- mixmdl$parameters$mean
  sd <- sqrt(mixmdl$parameters$variance$sigmasq)
  if(moment(exprs, 3) > 0) 
  {
    return(1 - pnorm(exprs, mean, sd) < cutoff)
  }
  else if(moment(exprs, 3) < 0)
  {
    return(pnorm(exprs, mean, sd) < cutoff)
  }
}

TCGABarcode <- function(fileName)
{
  return(paste(as.list(strsplit(strsplit(fileName, "lvl-3.")[[1]][2], "-")[[1]][1:3]), sep = "", collapse = "-"))
}

Mvalue <- function(beta_value)
{
  return(log2(beta_value/(1-beta_value)))
}
