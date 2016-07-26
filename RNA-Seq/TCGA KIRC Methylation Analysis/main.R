#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#MAIN

source("globals.R")
source("import.R")

main <- function() 
{
	import()
	density(exprMatrix[, symbol])

	for (symbol in XXX) 
	{
		splitMatrix[, symbol] <- split(exprMatrix[, symbol])
	}

	for (gene in XXXX)
	{
		for (patient in gene) 
		{
			if (splitMatrix[patient, symbol]  == 1) 
			{
				TMatrix <- cbind(TMatrix, bvalues[probes(gene), patient])
			} 
			else if (splitMatrix[patient, symbol] == 0)
			{
				NTMatrix <- cbind(NTMatrix, bvalues[probes(gene), patient])
			}
		}

		for (probe in rowNames(TMatrix)) a
		{
			pvalues[probe] <- wilcox.test(TMatrix[probe], NTMatrix[probe], correct = FALSE)
			adjpvalues(pvalue[probe], method = "BH")
			sigprobes <- names(pvalues[which(adjpvalue < cutoff)])
			write.table(c(gene, length/length(pvalues), sigprobes), append = TRUE)
		}
	}

}