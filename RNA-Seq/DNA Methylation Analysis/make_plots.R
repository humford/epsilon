
 setwd(paste("~/Documents/", cancer, "/Results", sep = ""))
 pdf(paste("~/Documents/Graphs/", cancer, "_Methylation_Results", sep = ""))
 
 clusterTable <- NULL
 
 par(mfrow = c(2,2))
 
#  m <- Mclust(geneStats[1:500, c(4,5)], G = 1:4)
#  reclass <- NULL
#  for(i in 1:4)
#  {
#     reclass <- c(reclass, quadrant(m$parameters$mean[,i]))
#  }
#  m$classification <- unlist(lapply(m$classification, function(x) reclass[x]))
#  plot(m, what = "classification", xlab = "Change in M-value", main = FALSE)
#  title(main = paste("all probes"))
   classification <- unlist(lapply(1:500, function(x) quadrant(c(geneStats$MvalDiff[x], geneStats$Skewness[x]))))
   
   plot(geneStats$MvalDiff[1:500], geneStats$Skewness[1:500], 
        xlab = "Change in M-value",
        ylab = "Skewness",
        main = paste(f, "probes"),
        col = c("blue", "red3", "purple", "green3")[classification], 
        pch = c(3, 0, 4, 2)[classification])
   abline(h = 0, lty = 2)
   abline(v = 0, lty = 2)
 
 for(f in functions)
 {
     corStats <- cor.test(functionStats[[f]]$Skewness[1:500], functionStats[[f]]$MvalDiff[1:500])  # Take 500 most significant genes
     write.table(t(c(f, corStats$conf.int[1], corStats$estimate, corStats$conf.int[2])), "corrTable", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
     
   #  m <- Mclust(functionStats[[f]][1:500, c(4,5)], G = 4)
   #  reclass <- NULL
   #  for(i in 1:4)
   #  {
   #     reclass <- c(reclass, quadrant(m$parameters$mean[,i]))
   #     
   #  }
   #  m$classification <- unlist(lapply(m$classification, function(x) reclass[x]))
   #  plot(m, what = "classification", xlab = "Change in M-value", main = FALSE)
   #  title(main = paste(f, "probes"))
   
     classification <- unlist(lapply(1:500, function(x) quadrant(c(functionStats[[f]]$MvalDiff[x], functionStats[[f]]$Skewness[x]))))
     
     plot(functionStats[[f]]$MvalDiff[1:500], functionStats[[f]]$Skewness[1:500], 
        xlab = "Change in M-value",
        ylab = "Skewness",
        main = "all probes",
        col = c("blue", "red3", "purple", "green3")[classification], 
        pch = c(3, 0, 4, 2)[classification])
     abline(h = 0, lty = 2)
     abline(v = 0, lty = 2)
     
     clusterTable <- rbind(clusterTable, c(length(which(classification == 1)), length(which(classification == 2)), length(which(classification == 3)), length(which(classification == 4))))
 }
 
 fTest <- fisher.test(clusterTable[c(1,2),])
 write.table(clusterTable, "clusterTable", col.names = FALSE, row.names = FALSE, quote = FALSE)
 write.table(format(fTest$p.value, digits = 5), "clusterTable", col.names = FALSE, append = TRUE, row.names = FALSE, quote = FALSE)

 corSampling <- NULL
 
 for(i in 100:min(1000, length(which(functionStats[[1]]$p.value < cutoff)), length(which(functionStats[[2]]$p.value < cutoff)), length(which(functionStats[[3]]$p.value < cutoff))))
 {
    row = i
    for(f in functions)
    {
      corStats <- cor.test(functionStats[[f]]$Skewness[1:i], functionStats[[f]]$MvalDiff[1:i])
      row <- c(row, corStats$estimate, corStats$conf.int[1], corStats$conf.int[2])
    }
    corSampling <- rbind(corSampling, row)
 }
 colnames(corSampling) <- c("num", "promoter_estimate", "promoter_lower", "promoter_upper", "UTR_estimate", "UTR_lower", "UTR_upper","body_estimate", "body_lower", "body_upper" )
 
 for(i in c(3,4,6,7,9,10))
 {
  corSampling[, i] <- smooth.spline(corSampling[, i], spar = 0.75)$y
 }
 
 library(ggplot2)
 
 chart<-
      ggplot(as.data.frame(corSampling)) +
      scale_colour_manual(name= "Probe Function", values=c("promoter" = "red", "UTR" = "blue", "body" = "green")) +
      geom_smooth(aes(x=num, y=promoter_estimate, color = "promoter"), size=1.0, span =  0.7) + 
      scale_x_continuous('Number of Genes',limits=c(100,1000)) +   
      scale_y_continuous('Correlation') +
      geom_ribbon(aes(x=num, ymin = promoter_lower, 
      ymax=promoter_upper), colour="red", fill="red",alpha=0.1) + 
      
      geom_smooth(aes(x=num, y=UTR_estimate, colour = "UTR"), size=1.0, span =  0.7) + 
      geom_ribbon(aes(x=num, ymin = UTR_lower, 
      ymax=UTR_upper), colour="blue", fill="blue",alpha=0.1)+
      
      geom_smooth(aes(x=num, y=body_estimate, colour = "body"), size=1.0, span =  0.7) + 
      geom_ribbon(aes(x=num, ymin = body_lower, 
      ymax=body_upper), colour="green", fill="green",alpha=0.1) +
      labs(title = paste("Skewness, M-value Correlation vs. Number of Included Genes in", cancer))
 print(chart) 
 
 dev.off()