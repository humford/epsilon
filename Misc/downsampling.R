lower <- NULL

for(a in 1:dim(exprMatrix)[[2]])
{
  skews <- NULL
  for(i in 1:dim(exprMatrix)[[1]])
  {
    skews[i] <- moment(sample(exprMatrix[i,], a), 3)
  }
  lower[a] <- length(which(skews < 0))
}

lower[1:3] <- lower[4]

lower <- lower/length(skews)

range <- 0.05 *abs(max(lower) - min(lower))
accept <- which(abs(lower  - lower[length(lower)])  <= range)
reject <- which(abs(lower  - lower[length(lower)]) > range)

plot(reject, lower[reject],col = "red", main= "LumA Gene Spliting by Sample Size", xlab = "Samples", ylab = "Fraction of Genes in Lower Mode", 
                                          xlim = c(1, length(lower)), ylim = c(min(lower) - range, max(lower) + range))
points(accept, lower[accept], col = "black")

abline(a = (lower[length(lower)] - range), b = 0)
abline(a = (lower[length(lower)] + range), b = 0)
text(min(accept), lower[min(accept)] + 2*range, as.character(min(accept)))

