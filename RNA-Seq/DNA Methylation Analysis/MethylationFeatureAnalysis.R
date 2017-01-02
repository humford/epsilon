cancer.names <- c("KIRC", "LGG", "LUSC", "SKCM")

KIRC_Sig_Genes <- read.table("~/Documents/KIRC/KIRC_Significant_Methylated_Genes", header = TRUE)
KIRC_Sig_Probes <- read.table("~/Documents/KIRC/KIRC_Significant_Methylated_Probes", header = TRUE, fill = TRUE)

LGG_Sig_Genes <- read.table("~/Documents/LGG/LGG_Significant_Methylated_Genes", header = TRUE)
LGG_Sig_Probes <- read.table("~/Documents/LGG/LGG_Significant_Methylated_Probes", header = TRUE, fill = TRUE)

LUSC_Sig_Genes <- read.table("~/Documents/LUSC/LUSC_Significant_Methylated_Genes", header = TRUE)
LUSC_Sig_Probes <- read.table("~/Documents/LUSC/LUSC_Significant_Methylated_Probes", header = TRUE, fill = TRUE)

SKCM_Sig_Genes <- read.table("~/Documents/SKCM/SKCM_Significant_Methylated_Genes", header = TRUE)
SKCM_Sig_Probes <- read.table("~/Documents/SKCM/SKCM_Significant_Methylated_Probes", header = TRUE, fill = TRUE)

x <- intersect(KIRC_Sig_Genes[,1], LGG_Sig_Genes[,1])
y <- intersect(LUSC_Sig_Genes[,1], SKCM_Sig_Genes[,1])
Sig_Genes <- intersect(x,y)

KIRC_Only_Sig_Genes <- setdiff(KIRC_Sig_Genes[,1], Sig_Genes)
LGG_Only_Sig_Genes <- setdiff(LGG_Sig_Genes[,1], Sig_Genes)
LUSC_Only_Sig_Genes <- setdiff(LUSC_Sig_Genes[,1], Sig_Genes)
SKCM_Only_Sig_Genes <- setdiff(SKCM_Sig_Genes[,1], Sig_Genes)

write.table(KIRC_Only_Sig_Genes, file = "~/Documents/KIRC_Only_Signifigant_Methylated_Genes", quote = FALSE, row.names = FALSE)


write.table(Sig_Genes, file = "~/Documents/Signifigant_Methylated_Genes", quote = FALSE, row.names = FALSE)

background <- dir("~/Documents/Gene_Methylation/")
write.table(background, file = "~/Documents/Background_Methylated_Genes", quote = FALSE, row.names = FALSE)
