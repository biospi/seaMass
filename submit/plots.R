
wrapper <- function(x, ...) { paste(strwrap(x, ...), collapse = "\n") }

plot.condition <- function(samps, stats, fc, meta, filename) {
  stats$variable <- factor(rownames(stats))
  stats$mean = colMeans(samps)
  stats <- cbind(stats, HPDinterval(mcmc(samps)))
  
  densities <- ddply(melt(samps), .(variable), function(x)
  {
    dens <- density(x$value, n=4096)
    data.frame(x=dens$x, y=dens$y)     
  })   
  densities$lower <- ifelse(densities$x<=log2(1/fc),densities$y,0) 
  densities$upper <- ifelse(densities$x>=log2(fc),densities$y,0) 
  y_range <- max(densities$y)*1.9
  x_range <- max(1.0,max(-min(densities$x[densities$y>y_range/100]),max(densities$x[densities$y>y_range/100])))
  
  stats$Up.text <- paste0("localFDR(up) = ",sapply(stats$Up, function(x) format(x,digits=2,scientific=F)))
  stats$Down.text <- paste0("localFDR(down) = ",sapply(stats$Down, function(x) format(x,digits=2,scientific=F)))
  stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc ")
  stats$mean.hjust <- ifelse(stats$mean<0,0,1)
  
  levels(stats$variable) = sub('^[0-9]+', '', levels(stats$variable))
  levels(densities$variable) = sub('^[0-9]+', '', levels(densities$variable))
  
  g <- ggplot(stats, aes(x=mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10))
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(-x_range, x_range),ylim=c(-0.0,y_range))
  g <- g + ggtitle(paste0('[', meta$Peptides, ' peptides, ', meta$Spectra, ' spectra]'))
  g <- g + xlab(expression('Log'[2]*' Ratio'))
  g <- g + ylab("Probability Density")
  g <- g + geom_rect(xmin=log2(1/fc),xmax=log2(fc),ymin=-2^32,ymax=2^32,alpha=0.15,colour="lightgrey") 
  g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3,colour="red") 
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=lower),ymin=0,alpha=0.3,fill="red")    
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=upper),ymin=0,alpha=0.3,fill="red")    
  g <- g + geom_vline(aes(xintercept=mean),size=2/3,colour="red") 
  g <- g + geom_vline(aes(xintercept=lower),colour="red",size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper),colour="red",size=1/2,lty=2)   
  g <- g + geom_text(aes(label=Up.text),x=x_range*0.98,y=y_range*0.9,hjust=1,vjust=1,,size=3)
  g <- g + geom_text(aes(label=Down.text),x=-x_range*0.98,y=y_range*0.9,hjust=0,vjust=1,,size=3)    
  g <- g + geom_text(aes(x=mean,label=mean.text,hjust=mean.hjust),y=y_range*0.7,vjust=1,size=3,colour="red")
  g
  ggsave(filename, g, height = 1+ 1*length(levels(stats$variable)), width = 6, limitsize = F)
}

plot.samples_sd <- function(samps, meta, filename) { 
  stats <- data.frame(variable = "Sample", mean = colMeans(samps))
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
  
  #levels(stats$variable) <- levels(densities$variable) <- paste0(levels(stats$variable), ' [', stats$N, ']')
  
  g <- ggplot(stats, aes(x=mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10))
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
  g <- g + ggtitle(paste0('[', meta$Peptides, ' peptides, ', meta$Spectra, ' spectra]'))
  g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3,colour="blue") 
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3,fill="blue")    
  g <- g + geom_vline(aes(xintercept=mean),size=2/3,colour="blue") 
  g <- g + geom_vline(aes(xintercept=lower),colour="blue",size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper),colour="blue",size=1/2,lty=2)   
  g <- g + geom_text(aes(x=mean,label=mean.text),y=y_range*0.9,hjust=0,vjust=1,size=3,colour="blue")
  g
  ggsave(filename, g, height = 1 + 1*length(levels(stats$variable)), width=6, limitsize=F)
}

plot.samples <- function(samps, meta, dd, filename) {
  samps <- mdply(colnames(samps), function(i) {
    data.frame(Sample = i, Condition = dd$Condition[dd$Sample==i][1], value = samps[,i])
  })  
  stats <- ddply(samps, .(Sample, Condition), function(x) {
    s <- data.frame(mean = mean(x$value), facet = ' ')
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps.trunc <- ddply(samps, .(Sample, Condition), function(x) {
    lower = stats$lower[stats$Sample == x$Sample[1]]
    upper = stats$upper[stats$Sample == x$Sample[1]]
    x[x$value >= lower & x$value <= upper,]
  })
  
  g <- ggplot(stats, aes(Sample, mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10),
                 legend.position="none")
  g <- g + facet_wrap(~ facet, ncol=1)
  g <- g + ggtitle(paste0('[', meta$Peptides, ' peptides, ', meta$Spectra, ' spectra]'))
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + geom_hline(yintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_violin(data = samps.trunc, aes(y = value, colour = Condition), alpha = 0.3, size = 2/3)
  g <- g + geom_segment(data = stats, aes(x = as.integer(Sample)-0.5, xend = as.integer(Sample) + 0.5, y = mean, yend = mean, colour = Condition),size = 2/3)
  g 
  ggsave(filename, g, height = 1 + 1*length(levels(stats$facet)), width=6, limitsize=F)
}

plot.peptides_sd <- function(samps, meta, dd, filename) { 
  stats <- data.frame(variable = colnames(samps),
                      N = daply(dd, .(Peptide), function(x) length(unique(x$Spectrum))),
                      mean = colMeans(samps.peptides_sd))
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
  
  levels(stats$variable) <- levels(densities$variable) <- paste0(levels(stats$variable), ' [', stats$N, ' spectra]')
  
  g <- ggplot(stats, aes(x=mean))
  g <- g + theme_bw()
  g <- g + theme(panel.margin=unit(0,"inches"),
                 panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6))
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
  g <- g + ggtitle(paste0('[', meta$Peptides, ' peptides, ', meta$Spectra, ' spectra]'))
  g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3,colour="blue") 
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3,fill="blue")    
  g <- g + geom_vline(aes(xintercept=mean),size=2/3,colour="blue") 
  g <- g + geom_vline(aes(xintercept=lower),colour="blue",size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper),colour="blue",size=1/2,lty=2)   
  g <- g + geom_text(aes(x=mean,label=mean.text),y=y_range*0.9,hjust=0,vjust=1,size=3,colour="blue")
  g
  ggsave(filename, g, height=1+1*length(levels(stats$variable)), width=6, limitsize=F)
}

plot.peptides <- function(samps1, samps2, meta, dd, filename) {  
  samps1 <- mdply(colnames(samps1), function(i) {
    data.frame(Sample = i, Condition = dd$Condition[dd$Sample==i][1], value = samps1[,i])
  }) 
  
  stats <- daply(samps2, .(Peptide), function(x) x$N[1])
  levels(samps2$Peptide) <- paste0(levels(samps2$Peptide), ' [', daply(samps2, .(Peptide), function(x) x$N[1]), ' spectra]')  
  
  stats1 <- ddply(samps1, .(Sample, Condition), function(x) {
    s <- data.frame(mean = mean(x$value))
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps1.trunc <- ddply(samps1, .(Sample, Condition), function(x) {
    lower = stats1$lower[stats1$Sample == x$Sample[1]]
    upper = stats1$upper[stats1$Sample == x$Sample[1]]
    x[x$value >= lower & x$value <= upper,]
  })
  
  stats2 <- ddply(samps2, .(Peptide, Sample), function(x) {
    s <- data.frame(N = x$N[1], mean = mean(x$value))
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps2.trunc <- ddply(samps2, .(Peptide, Sample), function(x) {
    lower = stats2$lower[stats2$Sample == x$Sample[1] & stats2$Peptide == x$Peptide[1]]
    upper = stats2$upper[stats2$Sample == x$Sample[1] & stats2$Peptide == x$Peptide[1]]
    x[x$value >= lower & x$value <= upper,]
  }) 
  
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
  g <- g + ggtitle(paste0('[', meta$Peptides, ' peptides, ', meta$Spectra, ' spectra]'))
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + geom_hline(yintercept=0,size=2/3,colour="darkgrey")          
  g <- g + geom_violin(data = samps1.trunc, aes(y = value, colour = Condition), alpha = 0.3, size = 2/3)
  g <- g + geom_segment(data = stats1, aes(x = as.integer(Sample)-0.5, xend = as.integer(Sample) + 0.5, y = mean, yend = mean, colour = Condition),size = 2/3)
  g <- g + geom_violin(data = samps2.trunc, aes(y = value), alpha = 0.3, size = 2/3)
  g <- g + geom_segment(aes(x = as.integer(Sample)-0.5, xend = as.integer(Sample) + 0.5, mean, yend = mean),size = 2/3)
  g 
  ggsave(filename, g, height=1+1*length(levels(stats2$Peptide)), width=6, limitsize=F)
}

plot.spectra_sd <- function(samps, meta, dd, filename) { 
  stats <- mdply(colnames(samps), function(i) {
    s <- data.frame(Spectrum = i,
                    Peptide = dd$Peptide[dd$Spectrum==i][1],
                    mean = mean(samps[,i]))
    s$variable <- paste0(s$Peptide, ' [', s$Spectrum, ']')
    cbind(s, HPDinterval(mcmc(samps[,i])))
  }) 
  
  densities <- mdply(colnames(samps), function(i) {
    dens <- density(samps[,i], n=4096)
    s <- data.frame(x=dens$x, y=dens$y, Peptide=dd$Peptide[dd$Spectrum==i][1], Spectrum=i)    
    s$variable <- paste0(s$Peptide, ' [', s$Spectrum, ']')
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
  g <- g + ggtitle(paste0('[', meta$Peptides, ' peptides, ', meta$Spectra, ' spectra]'))
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
