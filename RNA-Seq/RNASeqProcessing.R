library("AnnotationDbi")
library(org.Hs.eg.db)


convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


#control

setwd("~/Documents/Control/")

controlMatrix <- read.table("RNASeqControlRPKM.txt", sep = "\t", header = TRUE)

rownames(controlMatrix) <- lapply(controlMatrix[,1], function(x) strsplit(as.character(x), "[.]")[[1]][[1]])
controlMatrix <- controlMatrix[,-c(1:4)]

controlMatrix <- as.matrix(controlMatrix)

symbols <- convertIDs(rownames(controlMatrix), "ENSEMBL", "SYMBOL", org.Hs.eg.db)

rownames(controlMatrix) <- symbols

controlMatrix <- controlMatrix[-which(is.na(rownames(controlMatrix))), ]


dups <- duplicated(rownames(controlMatrix), fromLast = FALSE) | duplicated(rownames(controlMatrix), fromLast = TRUE)

d <- t(sapply(unique(rownames(controlMatrix)[which(dups)]), function(x) colMeans(controlMatrix[grep(x, rownames(controlMatrix)), ])))

controlMatrix <- rbind(d, controlMatrix[which(!dups), ])

controlMatrix <- controlMatrix[-which(rowMeans(controlMatrix) < 1), ]

controlMatrix <- apply(controlMatrix, c(2), fpkmToTpm)

controlMatrix <- log2(controlMatrix + 1)

write.table(controlMatrix, "RNASeqControlProcessed.txt")





cancer.names <- c("LGG", "KIRC", "LUSC", "HNSC", "SKCM")


for(cancer in cancer.names)
{

    setwd(paste("~/Documents/", cancer, sep = ""))

    #Import Clinical Data

    setwd("Clinical/Clinical/Biotab")

    clindat <- read.table(paste("nationwidechildrens.org_clinical_patient_", tolower(cancer), ".txt", sep = ""), sep="\t", header=T)
    clindat.hd <- clindat[(1:3),]
    clindat <- clindat[-(1:3),]

    setwd(paste("~/Documents/", cancer, sep = ""))


    #Import mRNA Expression Data

    setwd("Expression")


    map <- read.table("FILE_SAMPLE_MAP.txt", header=T, sep="\t")

    setwd("RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3")

    files <- dir()

    files <- files[grepl("rsem.genes.results", files)]

    map$batch <- as.factor(sapply(map$barcode.s, function(x) strsplit(as.character(x), "-")[[1]][6]))

    exprMatrix <- NULL

    for (i in 1:length(files))
    {
      file <- files[i]
      edat <- read.table(file, header = FALSE, skip = 2, sep = "", fill = TRUE)
      v <- edat[,3] #select scaled estimate not raw count
      v <- as.numeric(as.character(v)) * 10**6 #scaled estimate to TPM 
      exprMatrix <- cbind(exprMatrix, v)
    }

    rownames(exprMatrix) <- lapply(as.character(edat[,1]), function(x) paste(strsplit(x, "\\|")[[1]][[1]])) 

    exprMatrix <- exprMatrix[-which(rownames(exprMatrix) == "?"), ]
    exprMatrix <- exprMatrix[-which(duplicated(rownames(exprMatrix), fromlast = TRUE)),]

    setwd(paste("~/Documents/", cancer, sep = ""))


    #Remove Duplicates and Patients not represented in both Datasets

    file.barcode <- NULL
    for(i in 1:length(files))
    {
      file <- files[i]
      barcodes <- map[map[,"filename"] %in% file, "barcode.s."] 
      barcodes <- as.character(barcodes)
      barcodes <- paste(strsplit(barcodes, "-")[[1]][1:3], collapse="-", sep="")
      file.barcode <- c(file.barcode, barcodes)
    }

    colnames(exprMatrix) <- file.barcode

    map <- map[map$filename %in% files,]

    map$barcode <- file.barcode

    exprMatrix.uni <- NULL
    for(i in 1:length(unique(file.barcode)))
    {
      barcodes <- unique(file.barcode)[i]
      e <- cbind(exprMatrix[,colnames(exprMatrix) %in% barcodes])
      if(sum(colnames(exprMatrix) %in% barcodes) > 1 )
      {
        e <- apply(e, 1, mean, na.rm=T)
      }
      exprMatrix.uni <- cbind(exprMatrix.uni, e)
    }

    colnames(exprMatrix.uni) <- unique(file.barcode)

    exprMatrix <- exprMatrix.uni[,intersect(colnames(exprMatrix.uni), clindat[,"bcr_patient_barcode"])]
    clindata <- clindat[intersect(colnames(exprMatrix.uni), clindat[,"bcr_patient_barcode"]),]


    #Remove NAs and mostly zero count genes

    k <- which(is.na(exprMatrix), arr.ind=TRUE)
    rmeans <- rowMeans(exprMatrix, na.rm=TRUE)
    exprMatrix[k] <- rmeans[k[,1]]

    exprMatrix <- exprMatrix[-which(rowMeans(exprMatrix) < 1), ]
    
    rlogMatrix <- log2(1 + exprMatrix)


    #Batch Correction 
    BATCH_CORRECTING = FALSE
    
    if(BATCH_CORRECTING)
    {
      library(reshape2)
      library(ggplot2)
      library(DESeq2)
      library(limma)


      batches <- map$batch[which(map$barcode %in% colnames(exprMatrix) & !duplicated(map$barcode))]



      pca <- prcomp(t(rlogMatrix), scale = TRUE)

      scores <- data.frame(batches, pca$x[,1:2])

      pc1 <- qplot(x = PC1, y = PC2, data=scores, colour = batches, main = paste(cancer, "Pre-Correction")) + theme(legend.position = "none") 


      logCMP <- removeBatchEffect(rlogMatrix, batch = batches)



      pca <- prcomp(t(logCMP), scale = TRUE)

      scores <- data.frame(batches, pca$x[,1:2])

      pc2 <- qplot(x = PC1, y = PC2, data=scores, colour = batches, main = paste(cancer, "Post-Correction")) + theme(legend.position = "none") 



      genes <- intersect(rownames(exprMatrix), rownames(controlMatrix))

      compMatrix <- cbind(logCMP[genes, ], controlMatrix[genes, ])

      batches <- c(1:dim(logCMP)[[2]]*0  + 1, 1:dim(controlMatrix)[[2]]*0)



      pca <- prcomp(t(compMatrix))
      scores <- data.frame(batches, pca$x[,1:2])

      pcComp <- qplot(x = PC1, y = PC2, data=scores, colour = factor(batches), main = paste("Batch Corrected", cancer, "to Control Comparison")) + theme(legend.position = "none") 

      pdf(paste(cancer, "(Unorm) Batch Correction"))
      plot(pc1)
      plot(pc2)
      plot(pcComp)
      dev.off()
      
      write.table(logCMP, paste(cancer, "_Processed", sep = ""))
    }

    else
    {
      write.table(rlogMatrix, paste(cancer, "_Processed", sep = ""))
    }
    
    write.table(clindat, paste(cancer, "_Clinical", sep = ""))
}


#Methylation

for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  setwd("Methylation")
  
  directories <- dir()
  
  for(d in directories)
  {
    setwd(d)
    file <- dir()[grepl("TCGA", dir())]
    edat <- read.table(file, header = FALSE, skip = 2, sep = "\t", fill = TRUE)
    colnames(edat)
  }
}






