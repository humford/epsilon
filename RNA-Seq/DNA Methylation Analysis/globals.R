#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#GLOBALS

library(mclust)

cutoff <- 0.01
splitter_cutoff <- 0.1

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

mkplot <- function(exprs, symbol)
{
  hist(exprs, col=rgb(1,0,0,0.5), main=paste(symbol ,"Tail Splitting"), xlab="log2(Expression)", breaks = seq(0, max(exprs) + 0.5, by = 0.25))
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

function_mapper <- read.table("~/Documents/Methylation_Map", header = TRUE, nrow = 3656001, colClasses = c("character", "character"))

probe_function <- function(probe, symbol)
{
  s <- strsplit(as.character(function_mapper[probe, ][1]), ";")[[1]]
  group <- strsplit(as.character(function_mapper[probe, ][2]), ";")[[1]][which(s == symbol)]
  if(all(group == "TSS200" | group == "TSS1500")) return("promoter")
  else if(all(group == "3'UTR" | group == "5'UTR")) return("UTR")
  else if(all(group == "Body" | group == "1stExon")) return("body")
  else return(NULL)
}

quadrant <- function(x)
{
  if(x[1] > 0 && x[2] > 0) return(1)
  else if(x[1] <= 0 && x[2] > 0) return(2)
  else if(x[1] <= 0 && x[2] <= 0) return(3)
  else return(4)
}
