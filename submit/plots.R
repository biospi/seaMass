library(plyr)
library(ggplot2)
library(grid)

wrapper <- function(x, ...) { paste(strwrap(x, ...), collapse = "\n") }

title <- function(meta) {
  paste0('[',
         'ProteinID ', meta$ProteinID, ', ',
         'N ', meta$N, ', ', 
         meta$Peptides, ' peptide', ifelse(meta$Peptides==1,'','s'), ', ',
         meta$Spectra, ' spectr', ifelse(meta$Spectra==1,'um','a'), ']')  
}

plot.runchannels <- function(samps, meta, dd, filename) {
  samps <- mdply(colnames(samps), function(i) {
    data.frame(Run = sub(paste0(paste(levels(dd$Channel),collapse='|'),'$'),'', i),
               Channel = sub(paste0('^',paste(levels(dd$Run),collapse='|')),'', i),
               Condition = dd$Condition[dd$RunChannel==i][1],
               value = samps[,i])
  })  
  stats <- ddply(samps, .(Run, Channel, Condition), function(x) {
    s <- data.frame(mean = mean(x$value), facet = ' ')
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps.trunc <- ddply(samps, .(Run, Channel, Condition), function(x) {
    lower = stats$lower[stats$Run == x$Run[1] & stats$Channel == x$Channel[1]]
    upper = stats$upper[stats$Run == x$Run[1] & stats$Channel == x$Channel[1]]
    x[x$value >= lower & x$value <= upper,]
  })
  
  stats$Channel = factor(stats$Channel, levels=sort(levels(stats$Channel)))
  samps.trunc$Channel = factor(samps.trunc$Channel, levels=sort(levels(samps.trunc$Channel)))
  
  y_range <- max(1,max(max(samps.trunc$value),-min(samps.trunc$value)))
  
  g <- ggplot(stats, aes(Channel, mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10),
                 axis.text.x=element_text(size=8),
                 legend.position="none")
  g <- g + facet_wrap(~ Run, ncol=1, scales="free_y")
  g <- g + ggtitle(title(meta))
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + coord_cartesian(ylim=c(-y_range,y_range))
  g <- g + geom_hline(yintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_violin(data = samps.trunc, aes(y = value, colour = Condition), alpha = 0.3, size = 2/3)
  g <- g + geom_segment(data = stats, aes(x = as.integer(Channel)-0.5, xend = as.integer(Channel) + 0.5, y = mean, yend = mean, colour = Condition),size = 2/3)
  g 
  ggsave(filename, g, height = 1 + 1*length(levels(stats$Run)), width=6, limitsize=F)
}

plot.conditions <- function(samps, stats, fc, meta, filename) {
  densities <- ddply(melt(samps), .(variable), function(x)
  {
    dens <- density(x$value, n=65536)
    data.frame(x=dens$x, y=dens$y)     
  })   
  densities$lower <- ifelse(densities$x<=log2(1/fc),densities$y,0) 
  densities$upper <- ifelse(densities$x>=log2(fc),densities$y,0) 
  y_range <- max(densities$y[densities$variable!=levels(densities$variable)[1]])*1.9
  x_range <- max(1.0,max(-min(densities$x[densities$y>y_range/100]),max(densities$x[densities$y>y_range/100])))
  
  stats$Up.text <- paste0("localFDR(up) = ",sapply(stats$Up, function(x) format(x,digits=2,scientific=F)))
  stats$Down.text <- paste0("localFDR(down) = ",sapply(stats$Down, function(x) format(x,digits=2,scientific=F)))
  stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc ")
  stats$mean.hjust <- ifelse(stats$mean<0,0,1)
  
  levels(stats$variable) = sub('^[0-9]+', '', levels(stats$variable))
  levels(densities$variable) = sub('^[0-9]+', '', levels(densities$variable))
  
  g <- ggplot(stats, aes(x=mean,colour=variable,fill=variable))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10),
                 legend.position="none")
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(-x_range, x_range),ylim=c(-0.0,y_range))
  g <- g + ggtitle(title(meta))
  g <- g + xlab(expression('Log'[2]*' Ratio'))
  g <- g + ylab("Probability Density")
  g <- g + geom_rect(xmin=log2(1/fc),xmax=log2(fc),ymin=-2^32,ymax=2^32,alpha=0.15,colour="lightgrey",fill="lightgrey") 
  g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=lower),ymin=0,alpha=0.3)    
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=upper),ymin=0,alpha=0.3)    
  g <- g + geom_vline(aes(xintercept=mean,colour=variable),size=2/3,) 
  g <- g + geom_vline(aes(xintercept=lower,colour=variable),size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper,colour=variable),size=1/2,lty=2)   
  g <- g + geom_text(aes(label=Up.text),x=x_range*0.98,y=y_range*0.94,hjust=1,vjust=1,,size=3,colour='black')
  g <- g + geom_text(aes(label=Down.text),x=-x_range*0.98,y=y_range*0.94,hjust=0,vjust=1,,size=3,colour='black')    
  g <- g + geom_text(aes(x=mean,label=mean.text,hjust=mean.hjust),y=y_range*0.7,vjust=1,size=3)
  g
  ggsave(filename, g, height = 1+ 1*length(levels(stats$variable)), width = 6, limitsize = F)
}

plot.conditions_sd <- function(samps, stats, meta, filename) {   
  densities <- ddply(melt(samps), .(variable), function(x)
  {
    dens <- density(x$value, n=4096)
    data.frame(x=dens$x, y=dens$y)     
  })   
  y_range <- max(densities$y)*1.4
  x_range <- max(1,max(densities$x[densities$y>y_range/100]))
  
  stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc ")
  stats$mean.hjust <- ifelse(stats$mean<0,0,1)
  
  levels(stats$variable) = sub('^[0-9]+', '', levels(stats$variable))
  levels(densities$variable) = sub('^[0-9]+', '', levels(densities$variable))
    
  g <- ggplot(stats, aes(x=mean, colour=variable, fill=variable))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10),
                 legend.position="none")
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
  g <- g + ggtitle(title(meta))
  g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
  g <- g + geom_vline(aes(xintercept=mean, colour=variable),size=2/3) 
  g <- g + geom_vline(aes(xintercept=lower, colour=variable),size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper, colour=variable),size=1/2,lty=2)   
  g <- g + geom_text(aes(x=mean,label=mean.text),y=y_range*0.9,hjust=0,vjust=1,size=3)
  g
  ggsave(filename, g, height = 1 + 1*length(levels(stats$variable)), width=6, limitsize=F)
}

plot.samples <- function(samps1, samps2, meta, dd, filename) {  
  samps1 <- mdply(colnames(samps1), function(i) {
    data.frame(Sample = i, Condition = dd$Condition[dd$Sample==i][1], value = samps1[,i])
  }) 
  samps2 <- mdply(colnames(samps2), function(i) {
    data.frame(Sample = i, Condition = dd$Condition[dd$Sample==i][1], value = samps2[,i])
  })   
  
  stats1 <- ddply(samps1, .(Sample, Condition), function(x) {
    s <- data.frame(mean = mean(x$value))
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps1.trunc <- ddply(samps1, .(Sample, Condition), function(x) {
    lower = stats1$lower[stats1$Sample == x$Sample[1]]
    upper = stats1$upper[stats1$Sample == x$Sample[1]]
    x[x$value >= lower & x$value <= upper,]
  })
  
  stats2 <- ddply(samps2, .(Sample, Condition), function(x) {
    s <- data.frame(mean = mean(x$value), facet = ' ')
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps2.trunc <- ddply(samps2, .(Sample, Condition), function(x) {
    lower = stats2$lower[stats2$Sample == x$Sample[1]]
    upper = stats2$upper[stats2$Sample == x$Sample[1]]
    x[x$value >= lower & x$value <= upper,]
  }) 
  
  stats1$Sample <- reorder(stats1$Sample,as.numeric(stats1$Condition))
  samps1.trunc$Sample <- reorder(samps1.trunc$Sample,as.numeric(samps1.trunc$Condition))
  stats2$Sample <- reorder(stats2$Sample,as.numeric(stats2$Condition))
  samps2.trunc$Sample <- reorder(samps2.trunc$Sample,as.numeric(samps2.trunc$Condition))
  
  g <- ggplot(stats2, aes(Sample, mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ facet, ncol=1)
  g <- g + ggtitle(title(meta))
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + geom_hline(yintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_violin(data = samps2.trunc, aes(y = value), alpha = 0.3, size = 2/3)
  g <- g + geom_segment(aes(x = as.integer(Sample)-0.5, xend = as.integer(Sample) + 0.5, mean, yend = mean),size = 2/3)
  g <- g + geom_violin(data = samps1.trunc, aes(y = value, colour = Condition, fill = Condition), alpha = 0.3, size = 2/3)
  g <- g + geom_segment(data = stats1, aes(x = as.integer(Sample)-0.5, xend = as.integer(Sample) + 0.5, y = mean, yend = mean, colour = Condition),size = 2/3)
  g 
  ggsave(filename, g, height=2, width=6, limitsize=F)
}

plot.peptides_sd <- function(samps, meta, dd, filename) { 
  stats <- data.frame(variable = colnames(samps),
                      N = daply(dd, .(Peptide), function(x) length(unique(x$Spectrum))),
                      mean = colMeans(samps))
  stats <- cbind(stats, HPDinterval(mcmc(samps)))    
  
  densities <- ddply(melt(samps), .(variable), function(x)
  {
    dens <- density(x$value, n=4096)
    data.frame(x=dens$x, y=dens$y)     
  })   
  y_range <- max(densities$y)*1.4
  x_range <- max(1,max(densities$x[densities$y>y_range/100]))
  
  stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc ")
  stats$mean.hjust <- ifelse(stats$mean<0,0,1)
  
  levels(stats$variable) <- levels(densities$variable) <- paste0(levels(stats$variable), ' [', stats$N, ' spectr', ifelse(stats$N==1,'um','a'), ']')
  
  g <- ggplot(stats, aes(x=mean, colour=variable, fill=variable))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
  g <- g + ggtitle(title(meta))
  g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
  g <- g + geom_vline(aes(xintercept=mean, colour=variable),size=2/3) 
  g <- g + geom_vline(aes(xintercept=lower, colour=variable),size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper, colour=variable),size=1/2,lty=2)   
  g <- g + geom_text(aes(x=mean,label=mean.text),y=y_range*0.9,hjust=0,vjust=1,size=3)
  g
  ggsave(filename, g, height=1+1*length(levels(stats$variable)), width=6, limitsize=F)
}

plot.peptides <- function(samps2, meta, dd, filename) {  
  N <- daply(samps2, .(Peptide), function(x) x$N[1])
  levels(samps2$Peptide) <- paste0(levels(samps2$Peptide), ' [', N, ' spectr', ifelse(N==1,'um','a'), ']')  
  
  stats2 <- ddply(samps2, .(Peptide, Sample, Condition), function(x) {
    s <- data.frame(N = x$N[1], mean = mean(x$value))
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps2.trunc <- ddply(samps2, .(Peptide, Sample, Condition), function(x) {
    lower = stats2$lower[stats2$Sample == x$Sample[1] & stats2$Peptide == x$Peptide[1]]
    upper = stats2$upper[stats2$Sample == x$Sample[1] & stats2$Peptide == x$Peptide[1]]
    x[x$value >= lower & x$value <= upper,]
  }) 
  
  stats2$Sample <- reorder(stats2$Sample,as.numeric(stats2$Condition))
  samps2.trunc$Sample <- reorder(samps2.trunc$Sample,as.numeric(samps2.trunc$Condition))  
  
  g <- ggplot(stats2, aes(Sample, mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ Peptide, ncol=1)
  g <- g + ggtitle(title(meta))
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + geom_hline(yintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_violin(data = samps2.trunc, aes(y = value, colour = Condition), alpha = 0.3, size = 2/3)
  g <- g + geom_segment(aes(x = as.integer(Sample)-0.5, xend = as.integer(Sample) + 0.5, mean, yend = mean, colour = Condition, fill = Condition),size = 2/3)
  g 
  ggsave(filename, g, height=1+1*length(levels(stats2$Peptide)), width=6, limitsize=F)
}

plot.spectra_sd <- function(samps, meta, dd, filename) { 
  stats <- mdply(colnames(samps), function(i) {
    s <- data.frame(Spectrum = i,
                    Peptide = dd$Peptide[dd$Spectrum==i][1],
                    mean = mean(samps[,i]),
                    Run=dd$Run[dd$Spectrum==i][1])
    s$variable <- paste0(s$Peptide, ' ', s$Run, '[', s$Spectrum, ']')
    cbind(s, HPDinterval(mcmc(samps[,i])))
  }) 
  
  densities <- mdply(colnames(samps), function(i) {
    dens <- density(samps[,i], n=4096)
    s <- data.frame(x=dens$x, y=dens$y, Peptide=dd$Peptide[dd$Spectrum==i][1], Spectrum=i, Run=dd$Run[dd$Spectrum==i][1])    
    s$variable <- paste0(s$Peptide, ' ', s$Run, '[', s$Spectrum, ']')
    s
  })

  y_range <- max(densities$y)*1.4
  x_range <- max(1,max(densities$x[densities$y>y_range/100]))
  
  stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc ")
  stats$mean.hjust <- ifelse(stats$mean<0,0,1)
  
  g <- ggplot(stats, aes(x=mean, colour=Peptide, fill=Peptide))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
  g <- g + ggtitle(title(meta))
  g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
  g <- g + geom_vline(aes(xintercept=mean,colour=Peptide),size=2/3) 
  g <- g + geom_vline(aes(xintercept=lower,colour=Peptide),size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper,colour=Peptide),size=1/2,lty=2)   
  g <- g + geom_text(aes(x=mean,label=mean.text),y=y_range*0.9,hjust=0,vjust=1,size=3)
  g
  ggsave(filename, g, height=1+0.7*length(levels(stats$Spectrum)), width=6, limitsize=F)
}

plot.spectra <- function(preds, meta, dd, filename) { 
  dd.plot <- dd
  dd.plot$Count.min <- dd$Count+0.5-sqrt(dd$Count+0.25) 
  dd.plot$Count.max <- dd$Count+0.5+sqrt(dd$Count+0.25)
  dd.plot$lwr <- preds$lwr
  dd.plot$fit <- ifelse(preds$fit > preds$lwr & preds$fit < preds$upr, preds$fit, preds$lwr)
  dd.plot$upr <- preds$upr
  dd.plot$PeptideRunSpectrum <- factor(paste0(dd$Peptide, ' ', dd$Run, '[', dd$Spectrum, ']')) 
  dd.plot$Channel <- factor(dd.plot$Channel, levels=rev(levels(dd.plot$Channel)))
  
  g <- ggplot(dd.plot, aes(x=Channel, y=Count))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_text(size=5),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ PeptideRunSpectrum, ncol=1)
  g <- g + scale_x_discrete(expand=c(0,1))
  g <- g + scale_y_continuous(trans="log2")
  g <- g + coord_flip()
  g <- g + ggtitle(title(meta))
  g <- g + ylab(expression('Ion Count'))
  g <- g + geom_crossbar(aes(ymin=lwr, y=fit, ymax=upr, colour=Peptide), width=0.7)
  g <- g + geom_errorbar(aes(ymin=Count.min, ymax=Count.max), width=0.4) 
  g <- g + geom_errorbar(aes(ymin=Count, ymax=Count), width=0.4) 
  g
  ggsave(filename, g, height=1+0.7*length(levels(dd.plot$Spectrum)), width=6, limitsize=F)
}
