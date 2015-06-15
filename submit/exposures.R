#! /usr/bin/Rscript 

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  
  library(plyr)
  library(ggplot2)
  library(grid)
  
  # LOAD NORM SAMPLES
  samps.files = list.files(pattern="+[[:digit:]]\\.Rdata")
  load.samps <- function(file)
  {
    load(file)
    cbind(data.frame(RunChannel = colnames(samps), ProteinID = strsplit(file, ".", fixed=T)[[1]][1]), t(samps))
  }
  samps <- as.data.frame(rbind.fill(lapply(samps.files, load.samps)))

  # CALCULATE EXPOSURES (SIMPLE MEDIANS LIKE PROTEIN PILOT)
  exposures <- melt(ddply(samps, .(RunChannel), function(x) sapply(x[,3:ncol(x)], median)))

  # EXTRACT SUMMARY STATISTICS AND OUTPUT TO DISK
  summaries <- ddply(exposures, .(RunChannel), function(x) data.frame(mean=mean(x$value), sd=sd(x$value)))
  save(summaries,file=paste0("exposures.Rdata"))
  
  # PLOTTING
  densities <- ddply(exposures, .(RunChannel), function(x)
  {
    dens <- density(x$value, n=4096)
    data.frame(x=dens$x, y=dens$y)     
  })  
  gaussians <- ddply(summaries, .(RunChannel), function(d)
  {
    x <- d$mean + seq(-4,4,length=100)*d$sd
    y <- dnorm(x, d$mean, d$sd)
    data.frame(x=x, y=y)  
  })  
  fcs <- ddply(summaries, .(RunChannel), function(x)
  {
    data.frame(mean=x$mean, fc=paste0("  ",ifelse(x$mean<0, format(2^x$mean,digits=3), format(2^x$mean,digits=4)),"fc"))
  })  
  g <- ggplot(summaries, aes(x=mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=2),
                 panel.grid.major=element_line(size=0.5),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank())
  g <- g + facet_wrap(~ RunChannel, ncol=1)
  g <- g + coord_cartesian(xlim=c(min(densities$x), max(densities$x)),ylim=c(-0.0,max(densities$y)*1.3))
  g <- g + xlab(expression('Log'[2]*' Ratio'))
  g <- g + ylab("Probability Density")
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=1,colour="red") 
  g <- g + geom_vline(aes(xintercept=mean),size=1,colour="blue") 
  g <- g + geom_text(data=fcs,aes(x=mean,label=fc),y=max(densities$y)*1.22,hjust=0,vjust=1,size=3,colour="blue")
  g
  ggsave("medians_posterior.png", g, height=1*length(levels(exposures$RunChannel)), width=6)
  g <- g + geom_line(data=gaussians,aes(x=x,y=y),size=1,colour="blue") 
  g
  ggsave("medians_gaussian.png", g, height=1*length(levels(exposures$RunChannel)), width=6)
}
