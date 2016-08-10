TCGABarcode <- function(fileName)
{
  return(paste(as.list(strsplit(strsplit(fileName, "lvl-3.")[[1]][2], "-")[[1]][1:3]), sep = "", collapse = "-"))
}

Mvalue <- function(beta_value)
{
  return(log2(beta_value/(1-beta_value)))
}


#cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
cancer.names <- c("SKCM")

for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  
  exprMatrix <- read.table(paste(cancer, "Processed", sep = "_"), check.names = FALSE)
  
  setwd("Methylation")
  
  MvalMatrix <- NULL
  
  library(data.table)

  for(d in dir())
  {
    setwd(d)
    file <- dir()[grepl("TCGA", dir())]
    if((TCGABarcode(file) %in% colnames(exprMatrix)) && !(TCGABarcode(file) %in% colnames(MvalMatrix)))
    {
      edat <- fread(file, header = FALSE, skip = 2, sep = "\t")
      MvalMatrix <- cbind(MvalMatrix, Mvalue(as.numeric(as.character(edat[[2]]))))
      colnames(MvalMatrix) <- c(colnames(MvalMatrix)[-length(colnames(MvalMatrix))], TCGABarcode(file))
    }
    setwd("..")
  }
  
  edat <- as.data.frame(edat)
  
  setwd(paste("~/Documents/", cancer, sep = ""))
  
  rownames(MvalMatrix) <- edat[,1]
  
  MvalMatrix <- MvalMatrix[!is.na(rowSums(MvalMatrix)),]
  
  MvalMatrix <- MvalMatrix[,!duplicated(colnames(MvalMatrix))]
  
  mapper <- edat[ edat[,1] %in% rownames(MvalMatrix),c(1,3)]
  colnames(mapper) <- c("probe", "symbol")
  
  write.table(mapper, "Methylation_Map", quote = FALSE)
  
  MvalMatrix <- MvalMatrix[, intersect(colnames(MvalMatrix), colnames(exprMatrix))]
  
  setwd("Gene_Methylation")
  
  m <- strsplit(mapper[, "symbol"], ";")
  lens <- lapply(m, length)
  toCheck <- which(lens > 1)
  v <- NULL
  for(i in toCheck)
  {
    for(s in m[[i]])
    {
      v <- rbind(v, c(mapper[i, "probe"], s))
    }
  }
  
  colnames(v) <- colnames(mapper)
  map <- rbind(mapper, v)
  
  genes <- intersect(edat[,3], rownames(exprMatrix))
  
  for(symbol in genes)
  {
    whichProbes <- map[which(map[,"symbol"] == symbol),"probe"]
    if(length(whichProbes) > 0)
    {
      if(length(whichProbes) > 1)write.table(MvalMatrix[whichProbes, ], symbol, quote = FALSE)
      else
      {
        m <- t(MvalMatrix[whichProbes, ])
        rownames(m) <- whichProbes
        write.table(m, symbol, quote = FALSE, )
      }
    }
  }
}
