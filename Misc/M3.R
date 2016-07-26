trips <- read.csv("Documents/Ascii/Ascii/DAYV2PUB.CSV", header = TRUE)
important <- cbind(paste(c(trips["HOUSEID"])[[1]], c(trips["PERSONID"])[[1]], sep="-"), trips["TRVLCMIN"], trips["TRPMILES"])
important <- important[which(important["TRVLCMIN"] > 0 & important["TRPMILES"] > 0), ]
dailydriving <- aggregate(important[,-1], list(important[,1]), sum)
