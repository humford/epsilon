immune <- read.table("~/Documents/Tables/Immune")
colnames(immune) <- c("GO ID", "Components", "Term")

immune.unique <- cbind(immune[!duplicated(immune[,"GO ID"]),], table(immune["GO ID"])[unique(immune[,"GO ID"])])
colnames(immune.unique) <- c("GO ID", "Components", "Term", "Frequency")
immune.unique <- immune.unique[order(immune.unique["Frequency"], decreasing = TRUE), c(1,4,2,3)]
write.table(immune.unique, "~/Documents/Tables/ImmuneTerms", col.names = TRUE, row.name = FALSE)



translation <- read.table("~/Documents/Tables/Translation")
colnames(translation) <- c("GO ID", "Components", "Term")

translation.unique <- cbind(translation[!duplicated(translation[,"GO ID"]),], table(translation["GO ID"])[unique(translation[,"GO ID"])])
colnames(translation.unique) <- c("GO ID", "Components", "Term", "Frequency")
translation.unique <- translation.unique[order(translation.unique["Frequency"], decreasing = TRUE), c(1,4,2,3)]
write.table(translation.unique, "~/Documents/Tables/TranslationTerms", col.names = TRUE, row.name = FALSE)




neural <- read.table("~/Documents/Tables/Neural")
colnames(neural) <- c("GO ID", "Components", "Term")

neural.unique <- cbind(neural[!duplicated(neural[,"GO ID"]),], table(neural["GO ID"])[unique(neural[,"GO ID"])])
colnames(neural.unique) <- c("GO ID", "Components", "Term", "Frequency")
neural.unique <- neural.unique[order(neural.unique["Frequency"], decreasing = TRUE), c(1,4,2,3)]
write.table(neural.unique, "~/Documents/Tables/NeuralTerms", col.names = TRUE, row.name = FALSE)



metabolism <- read.table("~/Documents/Tables/Metabolism")
colnames(metabolism) <- c("GO ID", "Components", "Term")

metabolism.unique <- cbind(metabolism[!duplicated(metabolism[,"GO ID"]),], table(metabolism["GO ID"])[unique(metabolism[,"GO ID"])])
colnames(metabolism.unique) <- c("GO ID", "Components", "Term", "Frequency")
metabolism.unique <- metabolism.unique[order(metabolism.unique["Frequency"], decreasing = TRUE), c(1,4,2,3)]
write.table(metabolism.unique, "~/Documents/Tables/MetabolismTerms", col.names = TRUE, row.name = FALSE)