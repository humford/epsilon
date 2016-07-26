par(mfrow = c(2,2), mai = c(1,1,1,1), oma = c(2, 1,1,1))

MicroImmune <- matrix(c(0,0,10,0,0,5,0,0,1,0,0,0,0,0,8), nrow = 3)
colnames(MicroImmune) <- c("GBM", "OV", "LumA", "AML", "AML NK")
rownames(MicroImmune) <- c("DOWN", "MIDDLE", "UP")




plot(c(max(MicroImmune[3,]) - max(colSums(MicroImmune)), max(MicroImmune[3,])), c(0, 5),type = "n", axes = FALSE, ylab = "", xlab = "", main = "Immune System")
barplot(height =  MicroImmune[-1,], add = TRUE,axes = TRUE, horiz = TRUE, las = 1, col = c("blue", "green"))
barplot(height =  -MicroImmune[1,],add = TRUE,axes = FALSE, horiz = TRUE, las = 1, col = "red")


MicroMetabolism <- matrix(c(5,0,0,6,0,0,0,0,0,7,1,0,6,0,0), nrow = 3)
colnames(MicroMetabolism) <- c("GBM", "OV", "LumA", "AML", "AML NK")
rownames(MicroMetabolism) <- c("DOWN", "MIDDLE", "UP")




plot(c(-max(MicroMetabolism[1,]), max(colSums(MicroMetabolism)) - max(MicroMetabolism[1,])), c(0,5), type = "n", axes = FALSE, ylab = "", xlab = "", main = "Metabolism")
barplot(height =  MicroMetabolism[-1,], add = TRUE,axes = TRUE, horiz = TRUE, las = 1, col = c("blue", "green"))
barplot(height =  -MicroMetabolism[1,],add = TRUE,axes = FALSE, horiz = TRUE, las = 1, col = "red")



AML_6891Metabolism <- matrix(c(4,5,0,6,3,0,4,5,0,6,4,0), nrow = 3)
colnames(AML_6891Metabolism) <- c("GBM", "OV", "LumA", "AML NK")
rownames(AML_6891Metabolism) <- c("DOWN", "MIDDLE", "UP")



AML_15434MicroMetabolism <- matrix(c(3,0,0,5,0,0,3,0,0,0,4,6), nrow = 3)
colnames(AML_15434MicroMetabolism) <- c("GBM", "OV", "LumA", "AML")
rownames(AML_15434MicroMetabolism) <- c("DOWN", "MIDDLE", "UP")




plot(c(-max(AML_6891Metabolism[1,]), max(colSums(AML_6891Metabolism)) - max(AML_6891Metabolism[1,])), c(0,4), type = "n", axes = FALSE, ylab = "", xlab = "", main = "AML Metabolism")
barplot(height =  AML_6891Metabolism[-1,], add = TRUE,axes = TRUE, horiz = TRUE, las = 1, col = c("blue", "green"))
barplot(height =  -AML_6891Metabolism[1,],add = TRUE,axes = FALSE, horiz = TRUE, las = 1, col = "red")



plot(c(-max(colSums(AML_15434MicroMetabolism)) + max(AML_15434MicroMetabolism[3]), max(AML_15434MicroMetabolism[3,])), c(0,4), type = "n", axes = FALSE, ylab = "", xlab = "", main = "AML NK Metabolism")
barplot(height =  AML_15434MicroMetabolism[-1,], add = TRUE,axes = TRUE, horiz = TRUE, las = 1, col = c("blue", "green"))
barplot(height =  -AML_15434MicroMetabolism[1,],add = TRUE,axes = FALSE, horiz = TRUE, las = 1, col = "red")


par(fig = c(0, 1, 0, 1), oma = c(0.5, 0.5, 0.5, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend(x = "topleft", c("> Control", "~ Control", "< Control"), lty=c(5,5), lwd=c(5,5),col=c("green", "blue", "red"))

legend(x = "left", c("AML >", "AML ~", "AML <"), lty=c(5,5), lwd=c(5,5),col=c("green", "blue", "red"))
