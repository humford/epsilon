addRow <- function(DF, newRow)
{
  DF[dim(DF)[1] + 1, ] <- newRow
  return(DF)
}

mkplot <- function(exprs, symbol, splitter, param)
{
  hist(exprs, col=rgb(1,0,0,0.5), main=paste(symbol ,"Tail Splitting"), xlab="Value", breaks = seq(0, max(exprs) + 0.5, by = 0.05))
  hist(exprs[splitter(exprs, param)], col=rgb(0,0,1,0.5), add=T, breaks = seq(0, max(exprs) + 0.5, by = 0.05))
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

for (i in c(0.70, 0.75, 0.80, 0.85, 0.90, 0.95)) 
{
  quantile_tests[[as.character(i)]] <- test(quantile_splitter, i)
}

for (i in c(0.005, 0.01, 0.05, 0.025, 0.075, 0.10))
{
  gaussian_tests[[as.character(i)]] <- test(gaussian_splitter, i)
}

for (i in names(quantile_tests)) 
{
  print(paste("AREA ", i, "SD", sd(quantile_tests[[i]]$tail.area), "MEAN", mean(quantile_tests[[i]]$tail.area), "RANGE", range(quantile_tests[[i]]$tail.area)[2]))
  print(paste("Z ", i, "SD", sd(quantile_tests[[i]]$tail.z), "MEAN", mean(quantile_tests[[i]]$tail.z), "RANGE", range(quantile_tests[[i]]$tail.z)[2]))
}

for (i in names(gaussian_tests)) 
{
  print(paste("AREA ", i, "SD", sd(gaussian_tests[[i]]$tail.area), "MEAN", mean(gaussian_tests[[i]]$tail.area), "RANGE", range(gaussian_tests[[i]]$tail.area)[2]))
  print(paste("Z ", i, "SD", sd(gaussian_tests[[i]]$tail.z), "MEAN", mean(gaussian_tests[[i]]$tail.z), "RANGE", range(gaussian_tests[[i]]$tail.z)[2]))
}

dist <- rbeta(1000, 10, 5)

par(mfrow = c(2, 3))
for (i in names(quantile_tests)) { mkplot(dist, paste("QUANTILE", i), quantile_splitter, as.numeric(i)) }
par(mfrow = c(2, 3))
for (i in names(gaussian_tests)) { mkplot(dist, paste("GAUSSIAN", i), gaussian_splitter, as.numeric(i)) }