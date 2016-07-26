TCGABarcode <- function(fileName)
{
  return(paste(as.list(strsplit(strsplit(fileName, "lvl-3.")[[1]][2], "-")[[1]][1:3]), sep = "", collapse = "-"))
}

Mvalue <- function(beta_value)
{
  return(log2(beta_value/(1-beta_value)))
}


#cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
cancer.names <- c("KIRC")

for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  
  setwd("Methylation")
  
  MvalMatrix <- NULL

  
  for(d in dir())
  {
    setwd(d)
    file <- dir()[grepl("TCGA", dir())]
    edat <- fread(file, header = FALSE, skip = 2, sep = "\t")
    MvalMatrix <- cbind(MvalMatrix, Mvalue(as.numeric(as.character(edat[[2]]))))
    colnames(MvalMatrix) <- c(colnames(MvalMatrix)[-length(colnames(MvalMatrix))], TCGABarcode(file))
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
  
  setwd("Gene_Methylation")
  
  for(symbol in unique(edat[which(edat[,3] != ""),3]))
  {
    write.table(MvalMatrix[which(mapper[,"symbol"] == symbol), ], symbol, quote = FALSE)
  }
}