moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}





#AML

setwd("~/Documents/AML")

#Import Clinical Data

setwd("Clinical/Clinical/Biotab")

clindat <- read.table("nationwidechildrens.org_clinical_patient_laml.txt", sep="\t", header=T)
clindat.hd <- clindat[(1:3),]
clindat <- clindat[-(1:3),]

setwd("~/Documents/AML")


#Import mRNA Expression Data

setwd("Expression")

map <- read.table("FILE_SAMPLE_MAP.txt", header=T, sep="\t")


setwd("Expression-Genes/WUSM__HG-U133_Plus_2/Level_3")

files <- dir()

exprMatrix <- NULL

for (i in 1:length(files))
{
  file <- files[i]
  edat <- read.table(file, header = FALSE, skip = 2, sep = "", fill = TRUE)
  v <- edat[,2]
  v <- as.numeric(as.character(v))
  exprMatrix <- cbind(exprMatrix, v)
}

rownames(exprMatrix) <- as.character(edat[,1]) 

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

setwd("~/Documents/AML")

#Remove Duplicates and Patients not represented in both Datasets (Also log2 Normalize the data)

exprMatrix.both <- log2(exprMatrix[,intersect(colnames(exprMatrix), clindat[,"bcr_patient_barcode"])])
clindata.both <- clindat[intersect(colnames(exprMatrix), clindat[,"bcr_patient_barcode"]),]

exprMatrix <- exprMatrix.both

write.table(exprMatrix, "AML_Processed")
write.table(clindat, "AML_Clinical")




#AML GSE17855


library(GEOquery) 
GSE17855 <- getGEO("GSE17855") 
GSE17855 <- GSE17855[[1]] 

AMLexprMatrixProbe <- exprs(GSE17855)
AMLphenoData <- pData(GSE17855)

library(hgu133plus2.db)

probes <- ls(hgu133plus2SYMBOL)

symbols <- unlist(mget(probes, hgu133plus2SYMBOL))

rownames(AMLexprMatrixProbe) <- symbols[rownames(AMLexprMatrixProbe)]

AMLexprMatrixProbe <- AMLexprMatrixProbe[which(rownames(AMLexprMatrixProbe) != "NA"), ]


dups <- duplicated(rownames(AMLexprMatrixProbe), fromLast = FALSE) | duplicated(rownames(AMLexprMatrixProbe), fromLast = TRUE)

d <- t(sapply(unique(rownames(AMLexprMatrixProbe)[which(dups)]), function(x) colMeans(AMLexprMatrixProbe[grep(x, rownames(AMLexprMatrixProbe)), ])))

exprMatrix <- rbind(d, AMLexprMatrixProbe[which(!dups), ])

write.table(exprMatrix, "GSE17855_Processed")
write.table(phenoData, "GSE17855_Clinical")



#AML GSE6891

library(GEOquery) 
GSE6891 <- getGEO("GSE6891") 
GSE6891 <- GSE6891[[1]] 

AMLexprMatrixProbe <- exprs(GSE6891)
AMLphenoData <- pData(GSE6891)

library(hgu133plus2.db)

probes <- ls(hgu133plus2SYMBOL)

symbols <- unlist(mget(probes, hgu133plus2SYMBOL))

rownames(AMLexprMatrixProbe) <- symbols[rownames(AMLexprMatrixProbe)]

AMLexprMatrixProbe <- AMLexprMatrixProbe[which(rownames(AMLexprMatrixProbe) != "NA"), ]


dups <- duplicated(rownames(AMLexprMatrixProbe), fromLast = FALSE) | duplicated(rownames(AMLexprMatrixProbe), fromLast = TRUE)

d <- t(sapply(unique(rownames(AMLexprMatrixProbe)[which(dups)]), function(x) colMeans(AMLexprMatrixProbe[grep(x, rownames(AMLexprMatrixProbe)), ])))

exprMatrix <- rbind(d, AMLexprMatrixProbe[which(!dups), ])


k <- which(is.na(exprMatrix), arr.ind=TRUE)
rmeans <- rowMeans(exprMatrix, na.rm=TRUE)
exprMatrix[k] <- rmeans[k[,1]]

means <- NULL
for(i in 1:dim(exprMatrix)[[2]])
{
  means[i] <- mean((exprMatrix[,i]))
}

exprMatrix <- exprMatrix[, -which(means < 5.8)]

write.table(exprMatrix, "GSE6891_Processed")
write.table(phenoData, "GSE6891_Clinical")



#AML GSE15434

library(GEOquery) 
GSE15434 <- getGEO("GSE15434") 
GSE15434 <- GSE15434[[1]] 

AMLexprMatrixProbe <- exprs(GSE15434)
AMLphenoData <- pData(GSE15434)

library(hgu133plus2.db)

probes <- ls(hgu133plus2SYMBOL)

symbols <- unlist(mget(probes, hgu133plus2SYMBOL))

rownames(AMLexprMatrixProbe) <- symbols[rownames(AMLexprMatrixProbe)]

AMLexprMatrixProbe <- AMLexprMatrixProbe[which(rownames(AMLexprMatrixProbe) != "NA"), ]


dups <- duplicated(rownames(AMLexprMatrixProbe), fromLast = FALSE) | duplicated(rownames(AMLexprMatrixProbe), fromLast = TRUE)

d <- t(sapply(unique(rownames(AMLexprMatrixProbe)[which(dups)]), function(x) colMeans(AMLexprMatrixProbe[grep(x, rownames(AMLexprMatrixProbe)), ])))

exprMatrix <- rbind(d, AMLexprMatrixProbe[which(!dups), ])

k <- which(is.na(exprMatrix), arr.ind=TRUE)
rmeans <- rowMeans(exprMatrix, na.rm=TRUE)
exprMatrix[k] <- rmeans[k[,1]]

write.table(exprMatrix, "GSE15434_Processed")
write.table(phenoData, "GSE15434_Clinical")







#OV

setwd("~/Documents/OV")

#Import Clinical Data

setwd("Clinical/Clinical/Biotab")

clindat <- read.table("nationwidechildrens.org_clinical_patient_ov.txt", sep="\t", header=T)
clindat.hd <- clindat[(1:3),]
clindat <- clindat[-(1:3),]

setwd("~/Documents/OV")


#Import mRNA Expression Data

setwd("Expression")


map <- read.table("FILE_SAMPLE_MAP.txt", header=T, sep="\t")

setwd("Expression-Genes/UNC__AgilentG4502A_07_3/Level_3")

files <- dir()

exprMatrix <- NULL

for (i in 1:length(files))
{
  file <- files[i]
  edat <- read.table(file, header = FALSE, skip = 2, sep = "", fill = TRUE)
  v <- edat[,2]
  v <- as.numeric(as.character(v))
  exprMatrix <- cbind(exprMatrix, v)
}

rownames(exprMatrix) <- as.character(edat[,1]) 

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

setwd("~/Documents/OV")

#Remove Duplicates and Patients not represented in both Datasets

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
clindata.both <- clindat[intersect(colnames(exprMatrix.uni), clindat[,"bcr_patient_barcode"]),]


k <- which(is.na(exprMatrix), arr.ind=TRUE)
rmeans <- rowMeans(exprMatrix, na.rm=TRUE)
exprMatrix[k] <- rmeans[k[,1]]

write.table(exprMatrix, "OV_Processed")
write.table(clindat, "OV_Clinical")




#BRCA

setwd("~/Documents/BRCA/Expression")

map <- read.table("FILE_SAMPLE_MAP.txt", header=T, sep="\t")

setwd("Expression-Genes/UNC__AgilentG4502A_07_3/Level_3")

files <- dir()

exprMatrix <- NULL

for (i in 1:length(files))
{
  file <- files[i]
  edat <- read.table(file, header = FALSE, skip = 2, sep = "", fill = TRUE)
  v <- edat[,2]
  v <- as.numeric(as.character(v))
  exprMatrix <- cbind(exprMatrix, v)
}

rownames(exprMatrix) <- as.character(edat[,1]) 


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

setwd("~/Documents/BRCA")

#Remove Duplicates and Patients not represented in both Datasets

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


k <- which(is.na(exprMatrix), arr.ind=TRUE)
rmeans <- rowMeans(exprMatrix, na.rm=TRUE)
exprMatrix[k] <- rmeans[k[,1]]

write.table(exprMatrix, "BRCA_Processed")



#Control


setwd("~/Documents/Control")

library(GEOquery) 
GSE6536 <- getGEO("GSE6536") 
GSE6536 <- GSE6536[[1]] 
exprMatrixProbe <- exprs(GSE6536)
phenoData <- pData(GSE6536)


map <- read.table("GPL2507_Human_WG-6.csv", header= TRUE, sep = ",")

exprMatrix <- exprMatrixProbe[map$Target,]

rownames(exprMatrix) <- map$Symbol

dups <- duplicated(rownames(exprMatrix), fromLast = FALSE) | duplicated(rownames(exprMatrix), fromLast = TRUE)

d <- t(sapply(unique(rownames(exprMatrix)[which(dups)]), function(x) colMeans(exprMatrix[grep(x, rownames(exprMatrix)), ])))

exprMatrix <- rbind(d, exprMatrix[which(!dups), ])

write.table(exprMatrix, "Control_Processed")
write.table(phenoData, "Control_Clinical")



#Glioblastoma

setwd("~/Documents/Glioblastoma")

#Import Clinical Data

setwd("Clinical/Clinical/Biotab")

clindat <- read.table("nationwidechildrens.org_clinical_patient_gbm.txt", sep="\t", header=T)
clindat.hd <- clindat[(1:3),]
clindat <- clindat[-(1:3),]

setwd("~/Documents/Glioblastoma")


#Import mRNA Expression Data

setwd("Expression")

map <- read.table("FILE_SAMPLE_MAP.txt", header=T, sep="\t")

setwd("Expression-Genes/BI__HT_HG-U133A/Level_3")

files <- dir()

exprMatrix <- NULL

for (i in 1:length(files))
{
  file <- files[i]
  edat <- read.table(file, header = FALSE, skip = 2, sep = "", fill = TRUE)
  v <- edat[,2]
  v <- as.numeric(as.character(v))
  exprMatrix <- cbind(exprMatrix, v)
}

rownames(exprMatrix) <- as.character(edat[,1]) 

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

setwd("~/Documents/Glioblastoma")

#Remove Duplicates and Patients not represented in both Datasets

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

write.table(exprMatrix, "Glioblastoma_Processed")
write.table(clindat, "Glioblastoma_Clinical")




#LUAD

setwd("~/Documents/LUAD")

#Import Clinical Data

setwd("Clinical/Clinical/Biotab")

clindat <- read.table("nationwidechildrens.org_clinical_patient_luad.txt", sep="\t", header=T)
clindat.hd <- clindat[(1:3),]
clindat <- clindat[-(1:3),]

setwd("~/Documents/LUAD")


#Import mRNA Expression Data

setwd("Expression")


map <- read.table("FILE_SAMPLE_MAP.txt", header=T, sep="\t")

setwd("Expression-Genes/UNC__AgilentG4502A_07_3/Level_3")

files <- dir()

exprMatrix <- NULL

for (i in 1:length(files))
{
  file <- files[i]
  edat <- read.table(file, header = FALSE, skip = 2, sep = "", fill = TRUE)
  v <- edat[,2]
  v <- as.numeric(as.character(v))
  exprMatrix <- cbind(exprMatrix, v)
}

rownames(exprMatrix) <- as.character(edat[,1]) 

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

setwd("~/Documents/LUAD")

#Remove Duplicates and Patients not represented in both Datasets

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



k <- which(is.na(exprMatrix), arr.ind=TRUE)
rmeans <- rowMeans(exprMatrix, na.rm=TRUE)
exprMatrix[k] <- rmeans[k[,1]]

write.table(exprMatrix, "LUAD_Processed")
write.table(clindat, "LUAD_Clinical")


