#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#GLOBALS

moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

arbitrary_split <- function(exprs)
{
  if (moment(exprs, 3) > 0) 
  {
  	return(exprs > mean(exprs) + sd(exprs))
  }
  else if (moment(expres, 3) < 0)
  {
  	return(exprs < mean(exprs) - sd(exprs))
  }
}

arbitraryplus_split <- function(exprs)
{
	return 0
}

gaussian_split <- function(exprs)
{
	return 0
}

probes <- function(gene, map = mapper)
{
  return(mapper[which(mapper[,"symbol"] == gene), "probe"])
}

TCGABarcode <- function(fileName)
{
	return(paste(as.list(strsplit(strsplit(fileName, "lvl-3.")[[1]][2], "-")[[1]][1:3]), sep = "", collapse = "-"))
}

Mvalue <- function(beta_value)
{
  return(log2(beta_value/(1-beta_value)))
}
