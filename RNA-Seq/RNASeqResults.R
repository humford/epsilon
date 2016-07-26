par(mfrow = c(2,2), mai = c(1,1,1,1), oma = c(2, 1,1,1))


RNATranslation <- matrix(c(0,0,9,0,0,10,0,0,5,0,0,7,0,0,7), nrow = 3)
colnames(RNATranslation) <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
rownames(RNATranslation) <- c("DOWN", "MIDDLE", "UP")



plot(c(max(RNATranslation[3,]) - max(colSums(RNATranslation)), max(RNATranslation[3,])), c(0, 5),type = "n", axes = FALSE, ylab = "", xlab = "", main = "Translation: Cancer to Control")
barplot(height =  RNATranslation[-1,], add = TRUE,axes = TRUE, horiz = TRUE, las = 1, col = c("blue", "green"))
barplot(height =  -RNATranslation[1,],add = TRUE,axes = FALSE, horiz = TRUE, las = 1, col = "red")



MiddleTranslation <- matrix(c(0,7,0,0,4,0,0,6,0,0,7,0,0,0,0,0,4,2,0,6,0,0,5,0,0,2,0,0,7,0), nrow = 3)
colnames(MiddleTranslation) <- c("HNSC-SKCM", "LGG-SKCM", "LUSC-SKCM", "KIRC-SKCM", "LGG-HNSC", "LUSC-HNSC", "KIRC-HNSC", "LUSC-LGG", "KIRC-LGG", "KIRC-LUSC")
rownames(MiddleTranslation) <- c("DOWN", "MIDDLE", "UP")



plot(c(-max(colSums(MiddleTranslation)) + max(MiddleTranslation[3]), max(MiddleTranslation[3,])), c(0,10), type = "n", axes = FALSE, ylab = "", xlab = "", main = "Translation: Cancer to Cancer")
barplot(height =  MiddleTranslation[3,], add = TRUE,axes = TRUE, horiz = TRUE, las = 1, col = "green")
barplot(height =  -MiddleTranslation[-3,],add = TRUE,axes = FALSE, horiz = TRUE, las = 1, col = c("red", "blue"))



LGGNeuron <- matrix(c(7,0,0,5,0,0,9,0,0,6,0,0), nrow = 3)
colnames(LGGNeuron) <- c("SKCM", "HNSC", "LUSC", "KIRC")
rownames(LGGNeuron) <- c("DOWN", "MIDDLE", "UP")



plot(c(-max(LGGNeuron[1,]), max(colSums(LGGNeuron)) - max(LGGNeuron[1,])), c(0,4), type = "n", axes = FALSE, ylab = "", xlab = "", main = "LGG Neural")
barplot(height =  LGGNeuron[-1,], add = TRUE,axes = TRUE, horiz = TRUE, las = 1, col = c("blue", "green"))
barplot(height =  -LGGNeuron[1,],add = TRUE,axes = FALSE, horiz = TRUE, las = 1, col = "red")



LGGImmune <- matrix(c(0,0,3,0,0,8,0,0,8,0,0,9), nrow = 3)
colnames(LGGImmune) <- c("SKCM", "HNSC", "LUSC", "KIRC")
rownames(LGGImmune) <- c("DOWN", "MIDDLE", "UP")



plot(c(-max(LGGImmune[1,]), max(colSums(LGGImmune)) - max(LGGImmune[1,])), c(0,5), type = "n", axes = FALSE, ylab = "", xlab = "", main = "LGG Immune System")
barplot(height =  LGGImmune[-1,], add = TRUE,axes = TRUE, horiz = TRUE, las = 1, col = c("blue", "green"))
barplot(height =  -LGGImmune[1,],add = TRUE,axes = FALSE, horiz = TRUE, las = 1, col = "red")




par(fig = c(0, 1, 0, 1), oma = c(0.5, 0.5, 0.5, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend(x = "topleft", c(">", "~", "<"), lty=c(5,5), lwd=c(5,5),col=c("green", "blue", "red"))


legend(x = "left", c("LGG >", "LGG ~", "LGG <"), lty=c(5,5), lwd=c(5,5),col=c("green", "blue", "red"))
