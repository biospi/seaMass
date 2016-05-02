# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  print(paste(Sys.time(),"[Starting]"))
  
  load('index.Rdata')  
  load('parameters.Rdata')  
  load('design.Rdata')
  nsamp <- as.integer(ifelse("norm_nsamp" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nsamp"],1000))
  nchain <- as.integer(ifelse("norm_nchain" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nchain"],10))
  files = list.files(path="results",pattern="^[0-9]+\\.[0-9]+\\.Rdata")

  library(MCMCglmm)
  library(plyr)
  library(reshape2)
  library(ggplot2)
  
  # read norm samples and calculate exposures (simple medians like protein pilot)
  samps <- mdply(levels(design$Run), function(r) {
    samps <- mdply(levels(design$Channel)[2:length(levels(design$Channel))], function(c) {
      runchannel <- paste0(r,c)
      print(paste(Sys.time(),paste0("[Processing RunChannel ",runchannel,"]")))
      samps <- matrix(NA, nrow(data.index), nsamp) 
      for (f in files) {
        load(paste0("results/",f))
        
        if (runchannel %in% colnames(samps.runchannels)) {
          protein_id <- as.integer(gsub("\\.[0-9]+\\.Rdata$","",f))
          chain <- as.integer(gsub("\\.Rdata$","",gsub("^[0-9]+\\.","",f)))
          begin <- floor((chain-1)/nchain * nsamp) + 1
          rc.samps <- rev(samps.runchannels[max(1,nrow(samps.runchannels)-ceiling(nsamp/nchain)+1):nrow(samps.runchannels),runchannel])
          end <- min(begin+length(rc.samps)-1, nsamp)
          samps[protein_id+1,begin:end] <- rc.samps[1:(end-begin+1)]
        }
      }
      samps <- data.frame(value=aaply(samps, 2, function(x) median(x, na.rm=T)))
      samps$RunChannel <- runchannel
      samps$Channel <- c
      samps
    })
    samps$Run <- r
    samps
  })
  samps$RunChannel <- factor(samps$RunChannel)
  samps$Channel <- factor(samps$Channel)
  samps$Run <- factor(samps$Run)
  
  # EXTRACT SUMMARY STATISTICS AND OUTPUT TO DISK
  exposures <- ddply(samps, .(Run, Channel), function(x) data.frame(mean=mean(x$value), sd=sd(x$value)))
  save(exposures,file=paste0("exposures.Rdata"))
  
  # PLOTTING
  # densities <- ddply(samps, .(Run, Channel), function(x)
  # {
  #   dens <- density(x$value, n=4096)
  #   data.frame(x=dens$x, y=dens$y)     
  # })  
  # y_range <- max(densities$y)*1.35
  # x_range <- max(-min(densities$x[densities$y>y_range/100]),max(densities$x[densities$y>y_range/100]))*1.2
  # 
  # gaussians <- ddply(exposures, .(Run, Channel), function(d)
  # {
  #   x <- d$mean + seq(-4,4,length=100)*d$sd
  #   y <- dnorm(x, d$mean, d$sd)
  #   data.frame(x=x, y=y)  
  # })  
  # fcs <- ddply(exposures, .(Run, Channel), function(x)
  # {
  #   data.frame(mean=x$mean, fc=paste0("  ",ifelse(x$mean<0, format(2^x$mean,digits=3), format(2^x$mean,digits=4)),"fc"))
  # })  
  # g <- ggplot(exposures, aes(x=mean))
  # g <- g + theme_bw()
  # g <- g + theme(panel.border=element_rect(colour="black",size=2),
  #                panel.grid.major=element_line(size=0.5),
  #                axis.ticks=element_blank(),
  #                axis.text.y=element_blank(),
  #                plot.title=element_text(size=10),
  #                strip.background=element_blank())
  # g <- g + scale_x_continuous(expand = c(0,0))
  # g <- g + scale_y_continuous(expand = c(0,0))
  # g <- g + facet_wrap(Channel ~ Run, ncol=length(levels(samps$Run)))
  # g <- g + coord_cartesian(xlim=c(-x_range, x_range),ylim=c(-0.0,y_range))
  # g <- g + xlab(expression('Log'[2]*' Ratio'))
  # g <- g + ylab("Probability Density")
  # g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  # g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,size=1/2,alpha=0.3,fill='red')
  # g <- g + geom_line(data=densities,aes(x=x,y=y),size=1/2,colour="red") 
  # g <- g + geom_vline(aes(xintercept=mean),size=2/3,colour="red") 
  # g <- g + geom_text(data=fcs,aes(x=mean,label=fc),y=max(densities$y)*1.22,hjust=0,vjust=1,size=3,colour="red")
  # ggsave("medians_posterior.pdf", g, height=2*length(levels(samps$Channel)), width=6)
  # g <- g + geom_line(data=gaussians,aes(x=x,y=y),size=1/2,colour="black") 
  # ggsave("medians_gaussian.pdf", g, height=2*length(levels(samps$Channel)), width=6)
  
  print(paste(Sys.time(),"[Finished]"))
}
