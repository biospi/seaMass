#! /usr/bin/Rscript 

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  
  library(plyr)
  library(ggplot2)
  library(grid)
  library(reshape2)
  
  load('parameters.Rdata')  
  load('design.Rdata')
  nsamps <- as.integer(ifelse("n_samps" %in% parameters$Key,parameters$Value[parameters$Key=="n_samps"],10000))
  files = list.files(pattern="+[[:digit:]]\\.Rdata")

  # read norm samples and calculate exposures (simple medians like protein pilot)
  samps <- mdply(levels(design$Run), function(r) {
    samps <- mdply(levels(design$Channel)[2:length(levels(design$Channel))], function(c) {
      runchannel <- paste0(r,c)
      message(paste0('Reading ', runchannel, '...'))
      samps <- matrix(NA, length(files), nsamps) 
      for (i in 1:1000) {
        load(files[i])
        if (runchannel %in% colnames(samps.runchannels)) samps[i,] <- samps.runchannels[,runchannel]
      }
      samps <- data.frame(value=aaply(samps, 2, function(x) median(x, na.rm=T)))
      samps$RunChannel <- runchannel
      samps
    })
    samps
  })
  samps$RunChannel <- factor(samps$RunChannel)

  # EXTRACT SUMMARY STATISTICS AND OUTPUT TO DISK
  exposures <- ddply(samps, .(RunChannel), function(x) data.frame(mean=mean(x$value), sd=sd(x$value)))
  save(exposures,file=paste0("exposures.Rdata"))
  
  # PLOTTING
  densities <- ddply(samps, .(RunChannel), function(x)
  {
    dens <- density(x$value, n=4096)
    data.frame(x=dens$x, y=dens$y)     
  })  
  y_range <- max(densities$y)*1.35
  x_range <- max(-min(densities$x[densities$y>y_range/100]),max(densities$x[densities$y>y_range/100]))
  
  gaussians <- ddply(exposures, .(RunChannel), function(d)
  {
    x <- d$mean + seq(-4,4,length=100)*d$sd
    y <- dnorm(x, d$mean, d$sd)
    data.frame(x=x, y=y)  
  })  
  fcs <- ddply(exposures, .(RunChannel), function(x)
  {
    data.frame(mean=x$mean, fc=paste0("  ",ifelse(x$mean<0, format(2^x$mean,digits=3), format(2^x$mean,digits=4))))
  })  
  g <- ggplot(exposures, aes(x=mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=2),
                 panel.grid.major=element_line(size=0.5),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank())
  g <- g + facet_wrap(~ RunChannel, ncol=1)
  g <- g + coord_cartesian(xlim=c(-x_range, x_range),ylim=c(-0.0,y_range))
  g <- g + xlab(expression('Log'[2]*' Ratio'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=1,colour="red") 
  g <- g + geom_vline(aes(xintercept=mean),size=1,colour="blue") 
  g <- g + geom_text(data=fcs,aes(x=mean,label=fc),y=max(densities$y)*1.22,hjust=0,vjust=1,size=3,colour="blue")
  g
  ggsave("medians_posterior.png", g, height=1*length(levels(samps$RunChannel)), width=6)
  g <- g + geom_line(data=gaussians,aes(x=x,y=y),size=1,colour="blue") 
  g
  ggsave("medians_gaussian.png", g, height=1*length(levels(samps$RunChannel)), width=6)
}
