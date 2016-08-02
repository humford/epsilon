addRow <- function(DF, newRow)
{
  DF[dim(DF)[1] + 1, ] <- newRow
  return(DF)
}

mkplot <- function(exprs, symbol, splitter, param, area, z)
{
  hist(exprs, col=rgb(1,0,0,0.5), main=paste(symbol ,"Tail Splitting"), xlab=paste("Value", "(", "Z_SD=", format(z[1], digits = 2), "Z_MEAN=", format(z[2], digits = 2), ")"), breaks = seq(0, max(exprs) + 0.5, by = 0.05))
  hist(exprs[splitter(exprs, param)], col=rgb(0,0,1,0.5), add=T, breaks = seq(0, max(exprs) + 0.5, by = 0.05))
  box()
}

mkplot2 <- function(exprs, symbol, splitter, param)
{
  hist(exprs, col=rgb(1,0,0,0.5), main=paste(symbol), xlab=paste("Value"), breaks = seq(0, max(exprs) + 0.5, by = 0.25))
  hist(exprs[splitter(exprs, param)], col=rgb(0,0,1,0.5), add=T, breaks = seq(0, max(exprs) + 0.5, by = 0.25))
  box()
}

test_runs = 1

test <- function(splitter, param) 
{
  stats <- data.frame(tail.area = numeric(0), tail.start = numeric(0), tail.z = numeric(0))
  for (x in 1:test_runs) 
  {
    for (shape1 in 1:20) {
      for (shape2 in 1:20) {
        dist <- rbeta(1000, shape1, shape2)
        tail <- dist[splitter(dist, param)]
        if(length(tail) > 0) 
        {
          stats <- addRow(stats, c(length(tail), min(tail), abs(min(tail) - mean(dist))/sd(dist)))
        } else {
          stats <- addRow(stats, c(0, 0, 0))
        }
      }
    }
  }
  stats
}

quantile_tests <- list()
gaussian_tests <- list()

area <- function(i, tests) 
{
  c( sd(tests[[i]]$tail.area), mean(tests[[i]]$tail.area), range(tests[[i]]$tail.area)[2] )
}

z <- function(i, tests) 
{
  c( sd(tests[[i]]$tail.z), mean(tests[[i]]$tail.z), range(tests[[i]]$tail.z)[2] )
}

for (i in c(0.70, 0.75, 0.80, 0.85, 0.90, 0.95)) 
{
  quantile_tests[[as.character(i)]] <- test(quantile_splitter, i)
}

for (i in c(0.005, 0.01, 0.025, 0.05, 0.075, 0.10))
{
  gaussian_tests[[as.character(i)]] <- test(gaussian_splitter, i)
}

for (i in names(quantile_tests)) 
{
  aval <- area(i, quantile_tests)
  zval <- z(i, quantile_tests)
  print(paste("AREA ", i, "SD", aval[1], "MEAN", aval[2], "RANGE", aval[3]))
  print(paste("Z ", i, "SD", zval[1], "MEAN", zval[2], "RANGE", zval[3]))
}

for (i in names(gaussian_tests)) 
{
  aval <- area(i, gaussian_tests)
  zval <- z(i, gaussian_tests)
  print(paste("AREA ", i, "SD", aval[1], "MEAN", aval[2], "RANGE", aval[3]))
  print(paste("Z ", i, "SD", zval[1], "MEAN", zval[2], "RANGE", zval[3]))
}

dist <- rbeta(1000, 10, 5)

par(mfrow = c(2, 3))
for (i in names(quantile_tests)) 
{ 
  mkplot(dist, paste("QUANTILE", i), quantile_splitter, as.numeric(i), area(i, quantile_tests), z(i, quantile_tests)) 
}

par(mfrow = c(2, 3))
for (i in names(gaussian_tests)) 
{
  mkplot(dist, paste("GAUSSIAN", i), gaussian_splitter, as.numeric(i), area(i, gaussian_tests), z(i, gaussian_tests)) 
}

test_signifigant <- function(i, splitter, param, name) 
{
  par(mfrow = c(4, 6))
  for (x in 2:i) 
  {
     data <- exprMatrix[as.character(sig[x,1]),]
     mkplot2(data, paste(as.character(sig[x,1]), name, param), splitter, param)
  }
}

G0.05GENES <- read.table("/Users/henrywilliams/Documents/KIRC/Results/G0.05/KIRC_Signifcant_Methylated_Genes")
G0.05PROBES <- read.table("/Users/henrywilliams/Documents/KIRC/Results/G0.05/KIRC_Signifcant_Methylated_Probes", fill = TRUE)

G0.1GENES <- read.table("/Users/henrywilliams/Documents/KIRC/Results/G0.1/KIRC_Signifcant_Methylated_Genes")
G0.1PROBES <- read.table("/Users/henrywilliams/Documents/KIRC/Results/G0.1/KIRC_Signifcant_Methylated_Probes", fill = TRUE)

Q90GENES <- read.table("/Users/henrywilliams/Documents/KIRC/Results/Q90/KIRC_Signifcant_Methylated_Genes")
Q90PROBES <- read.table("/Users/henrywilliams/Documents/KIRC/Results/Q90/KIRC_Signifcant_Methylated_Probes", fill = TRUE)

Q95GENES <- read.table("/Users/henrywilliams/Documents/KIRC/Results/Q95/KIRC_Signifcant_Methylated_Genes")
Q95PROBES <- read.table("/Users/henrywilliams/Documents/KIRC/Results/Q95/KIRC_Signifcant_Methylated_Probes", fill = TRUE)