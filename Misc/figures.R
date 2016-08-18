skew <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s))
}

x <- rbeta(1000000,2,14.5)
y <- rbeta(1000000,4,5)

windows(width=10, height=6)
par(mfrow = c(1,2))
plot(density(x), xlim = c(0,1), main = "Beta Distribution: a = 2,  b = 14.5", xlab = paste(paste("skew =", round(skew(x,3), 4)), paste("     moment =", round(moment(x,3), 4))))
plot(density(y), xlim = c(0,1), main = "Beta Distribution: a = 4,  b = 5", xlab = paste(paste("skew =", round(skew(y,3), 4)), paste("     moment =", round(moment(y,3), 4))))



skew(x,3)
skew(y,3)

moment(x,3)
moment(y,3)




setwd("~/Documents/Control")
controlMatrix <- as.matrix(read.table("RNASeqControlProcessed.txt"))
setwd("~/Documents/Graphs")

ControlSkews <- NULL

for(i in 1:dim(controlMatrix)[[1]])
{ 
  ControlSkews[i] <-skew((controlMatrix[i,]), 3)
}

names(ControlSkews) <- rownames(controlMatrix)

dens <- density(ControlSkews)
plot(dens, main = "", xlab = "Skewness")
x1 <- min(which(dens$x >= 0))  
x2 <- max(which(dens$x <  100))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="lightgreen"))
x1 <- max(which(dens$x <= 0))  
x2 <- min(which(dens$x >  -100))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="gold"))





err <- function(Y, par)
{
  return(sum((Y-(par[1]*exp(-par[2] * (1:length(Y))) + par[3]))^2))
}



lower <- NULL

for(a in 1:dim(controlMatrix)[[2]])
{
  skews <- NULL
  for(i in 1:dim(controlMatrix)[[1]])
  {
    skews[i] <- moment(sample(controlMatrix[i,], a), 3)
  }
  lower[a] <- length(which(skews < 0))
}

lower[1:3] <- lower[4]

lower <- lower/length(skews)

range <- 0.05 *abs(max(lower) - min(lower))

plot(lower[-(1:20)], col = "black", xlim = c(1, length(lower)), ylim = c(min(lower) - range, max(lower) + range), lwd = 6, ann = FALSE, xaxt = 'n', yaxt = 'n', type = "l")

  
fit <- optim(par = c(0,0.5,0), err, Y = lower[-(1:20)])
errors <- (lower[length(lower)] - fit$par[3])/fit$par[3] * 100
converge <- 1/exp(fit$par[2])

lines(1:(length(lower)-20), fit$par[1]*exp(-fit$par[2] * (1:(length(lower)-20))) + fit$par[3], col = "red", lwd = 5) 

cancer <- "GBM"

setwd(paste("~/Documents/", cancer, sep = ""))
exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
setwd("~/Documents/Graphs")
skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <- skew((exprMatrix[i,]), 3)
}

skews <- skews[!is.na(skews)]
names(skews) <- rownames(exprMatrix)

controlDiff <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]
  
library(mclust)
out <- densityMclust(controlDiff, G = 3, modelNames = "V")

dens <- density(controlDiff)
plot(dens, main = "", xlab = "Skewness")

x1 <- min(which(dens$x >= min(controlDiff[which(out$classification == 1)])))
x2 <- max(which(dens$x < max(controlDiff[which(out$classification == 1)])))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="red"))
x1 <- min(which(dens$x >= min(controlDiff[which(out$classification == 2)])))
x2 <- max(which(dens$x < max(controlDiff[which(out$classification == 2)])))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="blue"))
x1 <- min(which(dens$x >= min(controlDiff[which(out$classification == 3)])))
x2 <- max(which(dens$x < max(controlDiff[which(out$classification == 3)])))  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="green"))

y <- dbeta(seq(0,1, by=0.01), 3,7)
x <- seq(0,1, by=0.01)
dens <- plot(x,y, main = "", xlab = "Skewness", type = 'l', xaxt='n', yaxt='n', ann=FALSE)
  
x1 <- min(which(x >= 0.0))
x2 <- max(which(x < 0.6))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="red"))
x1 <- min(which(x >= 0.59))
x2 <- max(which(x < 1.0))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="blue"))