#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#GLOBALS

library(mclust)

cutoff <- 0.01
splitter_cutoff <- 0.05

moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

splitter <- function(exprs)
{
  mixmdl <- Mclust(exprs, G = 1)
  mean <- mixmdl$parameters$mean
  sd <- sqrt(mixmdl$parameters$variance$sigmasq)
  if(moment(exprs, 3) > 0) 
  {
    return(1 - pnorm(exprs, mean, sd) < splitter_cutoff)
  }
  else if(moment(exprs, 3) < 0)
  {
    return(pnorm(exprs, mean, sd) < splitter_cutoff)
  }
}

mkplot <- function(exprs)
{
  hist(exprs, col=rgb(1,0,0,0.5), main="Tail Splitting", xlab="log2(Expression)", breaks = seq(0, max(exprs) + 0.5, by = 0.25))
  hist(exprs[splitter(exprs)], col=rgb(0,0,1,0.5), add=T, breaks = seq(0, max(exprs) + 0.5, by = 0.25))
  box()
}

addRow <- function(DF, newRow)
{
  DF[dim(DF)[1] + 1, ] <- newRow
  return(DF)
}



TCGABarcode <- function(fileName)
{
  return(paste(as.list(strsplit(strsplit(fileName, "lvl-3.")[[1]][2], "-")[[1]][1:3]), sep = "", collapse = "-"))
}

Mvalue <- function(beta_value)
{
  return(log2(beta_value/(1-beta_value)))
}
