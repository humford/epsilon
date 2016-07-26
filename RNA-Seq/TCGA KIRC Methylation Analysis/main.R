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
			TMatrix <- NULL
			NTMatrix <- NULL
			mvalues <- MValues(bvalues[probes(gene), patient])
			
			if (splitMatrix[patient, symbol]  == 1) 
			{
				TMatrix <- cbind(TMatrix, mvalues)
			} 
			else if (splitMatrix[patient, symbol] == 0)
			{
				NTMatrix <- cbind(NTMatrix, mvalues)
			}
		}

		for (probe in rowNames(TMatrix))
		{
			pvalues[probe] <- wilcox.test(TMatrix[probe], NTMatrix[probe], correct = FALSE)
			adjpvalues <- p.adjust(pvalue[probe], method = "BH")
			sigprobes <- names(pvalues[which(adjpvalues < cutoff)])
			
			write.table(c(gene, length/length(pvalues), sigprobes), append = TRUE)
		}
	}

}