library(plyr)
library(ggplot2)
library(grid)
library(methods)
library(coda)

wrapper <- function(x, ...) { paste(strwrap(x, ...), collapse = "\n") }

title <- function(meta) {
  paste0('[',
         'ProteinID ', meta$ProteinID, ', ',
         'N ', meta$N, ', ', 
         meta$Peptides, ' peptide', ifelse(meta$Peptides==1,'','s'), ', ',
         meta$Spectra, ' spectr', ifelse(meta$Spectra==1,'um','a'), ']')  
}


plot.norm.runchannels <- function(samps.Sol, meta, dd, filename) {
  if (length(levels(dd$Run))==1) {
    runchannels <- colnames(samps.Sol) %in% c(maply(levels(dd$Run), function(x) paste0('Channel', levels(dd$Channel))))
    samps <- samps.Sol[,runchannels,drop=F]
    colnames(samps) <- paste0(dd$Run[1], sub('Channel', '', colnames(samps)))
  } else {
    runchannels <- colnames(samps.Sol) %in% c(maply(levels(dd$Run), function(x) paste0('Run', x, ':Channel', levels(dd$Channel))))
    samps <- samps.Sol[,runchannels,drop=F]
    colnames(samps) <- sub('^Run', '', colnames(samps))  
    colnames(samps) <- sub(':Channel', '', colnames(samps))  
  }  
  for (r in levels(dd$Run)) samps[paste0(r,"113")] <- 0
  samps <- mdply(colnames(samps), function(i) {
    data.frame(Run = sub(paste0(paste(levels(dd$Channel),collapse='|'),'$'),'', i),
               Channel = sub(paste0('^',paste(levels(dd$Run),collapse='|')),'', i),
               Condition = dd$Condition[dd$RunChannel==i][1],
               value = samps[,i])
  })  
  samps$Channel <- factor(samps$Channel, levels=sort(levels(samps$Channel)))
  stats <- ddply(samps, .(Run, Channel, Condition), function(x) {
    s <- data.frame(mean = mean(x$value))
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  
  samps.trunc <- ddply(samps, .(Run, Channel, Condition), function(x) {
    lower = stats$lower[stats$Run == x$Run[1] & as.character(stats$Channel) == as.character(x$Channel[1])]
    upper = stats$upper[stats$Run == x$Run[1] & as.character(stats$Channel) == as.character(x$Channel[1])]
    x[x$value >= lower & x$value <= upper,]
  })
  samps.violin <- samps.trunc[samps.trunc$Channel!="113",]
  
  y_range <- max(2,max(max(samps.trunc$value),-min(samps.trunc$value))*1.2)
  
  g <- ggplot(stats, aes(Channel, mean))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10),
                 axis.text.x=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ Run)
  g <- g + ggtitle(title(meta))
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + coord_cartesian(ylim=c(-y_range,y_range))
  g <- g + geom_hline(yintercept=0,size=1/2,colour="darkgrey")     
  g <- g + geom_boxplot(data = samps.trunc, aes(y = value), alpha = 0.3, size = 0, outlier.size = 0)
  g <- g + geom_violin(data = samps.violin, aes(y = value), alpha = 0.3, size = 1/2, trim=T)
  g <- g + geom_segment(data = stats, aes(x = as.integer(Channel)-0.5, xend = as.integer(Channel) + 0.5, y = mean, yend = mean),size = 1/2)
  ggsave(filename, g, height = 2, width=6, limitsize=F)
  
  y_range
}


plot.peptides_sd <- function(samps.sqrtVCV, meta, dd, filename) { 
  if (length(levels(dd$Peptide)) == 1) {
    samps <- samps.sqrtVCV[,"Digest",drop=F]
    colnames(samps) <- dd$Peptide[1]      
  } else {
    samps <- samps.sqrtVCV[,colnames(samps.sqrtVCV) %in% paste0(levels(dd$Peptide), ".Digest"),drop=F]
    colnames(samps) <- sub('\\.Digest$', '', colnames(samps))      
  }
  
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
  
  stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc ")
  stats$mean.hjust <- ifelse(stats$mean<0,0,1)
 
  stats$variable <- sub(":.*:","",stats$variable)
  densities$variable <- sub(":.*:","",densities$variable)
  stats$variable <- factor(stats$variable)
  densities$variable <- factor(densities$variable)
  levels(stats$variable) <- levels(densities$variable) <- paste0(levels(stats$variable), ' [', stats$N, ' spectr', ifelse(stats$N==1,'um','a'), ']')
  
  g <- ggplot(stats, aes(x=mean, colour=variable, fill=variable))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(0,1),ylim=c(-0.0,y_range))
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
  ggsave(filename, g, height=1+1*length(levels(stats$variable)), width=3.5, limitsize=F)
}


plot.norm.peptides <- function(samps.Sol, meta, dd, filename) {    
  samps2 <- mdply(levels(dd$Peptide), function(i) {
    s <- samps.Sol[,colnames(samps.Sol) %in% paste0('Digest.', i, '.Digest.', levels(dd$Digest))]
    colnames(s) <- sub(paste0('Digest.', i, '.Digest.'), '', colnames(s), fixed = T)
    # add Digest latent effects and melt
    mdply(colnames(s), function(j) {
      data.frame(Peptide = i,
                 Channel = dd$Channel[dd$Digest==j][1],
                 Run = dd$Run[dd$Digest==j][1],
                 value = s[,j])
    })
  })
  stats2 <- ddply(samps2, .(Peptide, Run, Channel), function(x) {
    s <- data.frame(mean = mean(x$value))
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps2.trunc <- ddply(samps2, .(Peptide, Run, Channel), function(x) {
    lower = stats2$lower[stats2$Run == x$Run[1] & stats2$Channel == x$Channel[1] & stats2$Peptide == x$Peptide[1]]
    upper = stats2$upper[stats2$Run == x$Run[1] & stats2$Channel == x$Channel[1] & stats2$Peptide == x$Peptide[1]]
    x[x$value >= lower & x$value <= upper,]
  })
  
  stats2$Channel <- factor(as.character(stats2$Channel))
  samps2.trunc$Channel <- factor(as.character(samps2.trunc$Channel))
  stats2$Peptide <- sub(":.*:","",stats2$Peptide)
  samps2.trunc$Peptide <- sub(":.*:","",samps2.trunc$Peptide)
  stats2$Peptide <- factor(stats2$Peptide)
  samps2.trunc$Peptide <- factor(samps2.trunc$Peptide)
  
  g <- ggplot(stats2, aes(Channel, mean))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 axis.text.x=element_text(size=6),
                 legend.position="none")
  g <- g + coord_cartesian(ylim=c(-2,2))
  g <- g + facet_wrap(Peptide~Run,drop=F,ncol=length(levels(dd$Run)))
  g <- g + ggtitle(title(meta))
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + geom_hline(yintercept=0,size=1/2,colour="darkgrey")          
  g <- g + geom_violin(data = samps2.trunc, aes(y = value, colour = Peptide), alpha = 0.3, size = 1/2, trim=T)
  g <- g + geom_segment(aes(x = as.integer(Channel)-0.5, xend = as.integer(Channel) + 0.5, mean, yend = mean, colour = Peptide),size = 1/2) 
  ggsave(filename, g, height=1+1*length(levels(stats2$Peptide)), width=6, limitsize=F)
}


plot.spectra_sd <- function(samps.sqrtVCV, meta, dd, filename) { 
  if (length(levels(dd$Spectrum)) == 1) {
    samps <- samps.sqrtVCV[,"units",drop=F]
    colnames(samps) <- dd$Spectrum[1] 
  } else {
    samps <- samps.sqrtVCV[,colnames(samps.sqrtVCV) %in% paste0(levels(dd$Spectrum), ".units"),drop=F]
    colnames(samps) <- sub('\\.units$', '', colnames(samps))    
  }
  
  stats <- mdply(colnames(samps), function(i) {
    s <- data.frame(Spectrum = i,
                    Peptide = dd$Peptide[dd$Spectrum==i][1],
                    mean = mean(samps[,i]),
                    Run=dd$Run[dd$Spectrum==i][1])
    s$variable <- paste0(sub(":.*:","",s$Peptide), ' ', s$Run, '[', s$Spectrum, ']')
    cbind(s, HPDinterval(mcmc(samps[,i])))
  }) 
  
  densities <- mdply(colnames(samps), function(i) {
    dens <- density(samps[,i], n=4096)
    s <- data.frame(x=dens$x, y=dens$y, Peptide=dd$Peptide[dd$Spectrum==i][1], Spectrum=i, Run=dd$Run[dd$Spectrum==i][1])    
    s$variable <- paste0(sub(":.*:","",s$Peptide), ' ', s$Run, '[', s$Spectrum, ']')
    s
  }) 
  y_range <- max(densities$y)*1.4
  
  stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc ")
  stats$mean.hjust <- ifelse(stats$mean<0,0,1)
  
  g <- ggplot(stats, aes(x=mean, colour=Peptide, fill=Peptide))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ variable, ncol=1)
  g <- g + coord_cartesian(xlim=c(0,1),ylim=c(-0.0,y_range))
  g <- g + ggtitle(title(meta))
  g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=1/2,colour="darkgrey")          
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
  g <- g + geom_vline(aes(xintercept=mean,colour=Peptide),size=2/3) 
  g <- g + geom_vline(aes(xintercept=lower,colour=Peptide),size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper,colour=Peptide),size=1/2,lty=2)   
  g <- g + geom_text(aes(x=mean,label=mean.text),y=y_range*0.9,hjust=0,vjust=1,size=3)
  g
  ggsave(filename, g, height=1+1*length(levels(stats$Spectrum)), width=3.5, limitsize=F)
}

  
plot.spectra <- function(results, marginal, meta, dd, filename) { 
  preds <- data.frame(predict(results[[1]],interval="confidence",marginal=marginal))
  preds50 <- data.frame(predict(results[[1]],interval="confidence",marginal=marginal,level=0))
  
  dd.plot <- dd
  dd.plot$Count.min <- dd.plot$Count+0.5-sqrt(dd.plot$Count+0.25) 
  dd.plot$Count.max <- dd.plot$Count+0.5+sqrt(dd.plot$Count+0.25)
  x_range=c(max(min(dd.plot$Count.min)/1.1,2^-0.5), max(dd.plot$Count.max)*1.1)
  dd.plot$Count.min <- pmax(x_range[1], dd.plot$Count.min)
  dd.plot$Count.max <- pmin(x_range[2], dd.plot$Count.max)
  
  dd.plot$lwr <- pmax(x_range[1], preds$lwr)
  dd.plot$upr <- pmin(x_range[2], preds$upr)
  dd.plot$med <- pmax(x_range[1], pmin(x_range[2], preds50$lwr))
  dd.plot$Spectrum <- factor(dd.plot$Spectrum)
  dd.plot$Peptide <- factor(dd.plot$Peptide)
  
  dd.plot$PeptideRunSpectrum <- factor(paste0(sub(":.*:","",dd.plot$Peptide), ' ', dd.plot$Run, '[', dd.plot$Spectrum, ']')) 
  dd.plot$Channel <- factor(dd.plot$Channel)
  dd.plot$lwr <- pmax(x_range[1],dd.plot$lwr)
  dd.plot$upr <- pmin(x_range[2],dd.plot$upr)
  
  g <- ggplot(dd.plot, aes(x=Channel, y=Count, colour=Peptide))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_text(size=5),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ PeptideRunSpectrum, ncol=1)
  g <- g + scale_x_discrete(expand=c(0,1))
  g <- g + scale_y_continuous(trans="log2", limits=x_range, expand=c(0,0))
  g <- g + coord_flip()
  g <- g + ggtitle(title(meta))
  g <- g + ylab(expression('Ion Count'))
  g <- g + geom_crossbar(aes(ymin=lwr, y=med, ymax=upr), width=0.8, colour="darkgrey")
  g <- g + geom_errorbar(aes(ymin=Count, ymax=Count), width=0.8, size=1) 
  g <- g + geom_errorbar(aes(ymin=Count.min, ymax=Count.max), width=0.8, size=1/2) 
  ggsave(filename, g, height=1+1*length(levels(dd.plot$Spectrum)), width=6, limitsize=F)
}


plot.conditions <- function(samps.Sol, fc, meta, dd, filename) {
  samps <- samps.Sol[,colnames(samps.Sol) %in% paste0('Condition', levels(dd$Condition)),drop=F]
  samps.baseline <- data.frame(x = rnorm(nrow(samps),1e-6,1e-6))
  colnames(samps.baseline)[1] <- paste0('Condition',levels(dd$Condition)[1])
  samps <- cbind(samps.baseline, samps)
  colnames(samps) <- sub('Condition', '', colnames(samps), fixed=T)  
  # stats for plot and csv output
  stats <- data.frame(variable = factor(colnames(samps), levels=levels(dd$Condition)),
                                 Up0 = 1 - colSums(samps > 0) / dim(samps.Sol)[1],
                                 Down0 = 1 - colSums(samps < 0) / dim(samps.Sol)[1],
                                 Up = 1 - colSums(samps > log2(fc)) / dim(samps.Sol)[1],
                                 Down = 1 - colSums(samps < -log2(fc)) / dim(samps.Sol)[1],
                                 Same = 1 - colSums(samps >= -log2(fc) & samps <= log2(fc)) / dim(samps.Sol)[1])
  stats$mean = colMeans(samps)
  stats <- cbind(stats, HPDinterval(mcmc(samps)))
  
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
  
  g <- ggplot(stats, aes(x=mean,fill=variable))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
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
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=lower),ymin=0)    
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=upper),ymin=0)    
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  g <- g + geom_vline(aes(xintercept=mean),size=2/3,) 
  g <- g + geom_vline(aes(xintercept=lower),size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper),size=1/2,lty=2)   
  g <- g + geom_text(aes(label=Up.text),x=x_range*0.98,y=y_range*0.94,hjust=1,vjust=1,,size=2.5,colour='black')
  g <- g + geom_text(aes(label=Down.text),x=-x_range*0.98,y=y_range*0.94,hjust=0,vjust=1,,size=2.5,colour='black')    
  g <- g + geom_text(aes(x=mean,label=mean.text,hjust=mean.hjust,colour=variable),y=y_range*0.7,vjust=1,size=2.5)
  ggsave(filename, g, height = 1+ 1*length(levels(stats$variable)), width = 6, limitsize = F)
}


plot.conditions_sd <- function(samps.sqrtVCV, meta, dd, filename) {   
  #if (var_equal) {
  #  samps <- samps.sqrtVCV[,"Sample",drop=F]
  #  stats <- data.frame(variable = "Sample", mean = colMeans(samps))
  #} else {
    samps <- samps.sqrtVCV[,colnames(samps.sqrtVCV) %in% paste0(levels(dd$Population), ".Sample"),drop=F]
    colnames(samps) <- sub('\\.Sample$', '', colnames(samps))    
    stats <- data.frame(variable = factor(colnames(samps), levels=levels(dd$Population)), mean = colMeans(samps))
  #}
  stats <- cbind(stats, HPDinterval(mcmc(samps)))  
  
  densities <- ddply(melt(samps), .(variable), function(x)
  {
    dens <- density(x$value, n=4096)
    data.frame(x=dens$x, y=dens$y)     
  })   
  y_range <- max(densities$y[densities$variable %in% levels(densities$variable)[ifelse(count(dd[!duplicated(dd$Sample),], "Population")$freq > 1, T, F)]])*1.4
  x_range <- max(1,max(densities$x[densities$y>y_range/100]))
  
  stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc")
  stats$mean.hjust <- ifelse(stats$mean<0,0,1)
  
  levels(stats$variable) = sub('^[0-9]+', '', levels(stats$variable))
  levels(densities$variable) = sub('^[0-9]+', '', levels(densities$variable))
    
  g <- ggplot(stats, aes(x=mean, fill=variable))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
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
  g <- g + geom_vline(xintercept=0,size=2/3)          
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0)    
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  g <- g + geom_vline(aes(xintercept=mean),size=2/3) 
  g <- g + geom_vline(aes(xintercept=lower),size=1/2,lty=2)      
  g <- g + geom_vline(aes(xintercept=upper),size=1/2,lty=2)   
  g <- g + geom_text(aes(x=mean,label=mean.text,colour=variable),y=y_range*0.9,hjust=0,vjust=1,size=3)
  ggsave(filename, g, height = 1 + 1*length(levels(stats$variable)), width=3.5, limitsize=F)
}


plot.samples <- function(samps.Sol, meta, dd, filename) {  
  samps.conditions <- samps.Sol[,colnames(samps.Sol) %in% paste0('Condition', levels(dd$Condition)),drop=F]
  samps.baseline <- data.frame(x = rnorm(nrow(samps.conditions),1e-6,1e-6))
  colnames(samps.baseline)[1] <- paste0('Condition',levels(dd$Condition)[1])
  samps.conditions <- cbind(samps.baseline, samps.conditions)
  colnames(samps.conditions) <- sub('Condition', '', colnames(samps.conditions), fixed=T)    
  
  #if (var_equal) {
  #  samples <- data.frame(Name=paste0('Sample.', levels(dd$Sample)),Sample=levels(dd$Sample))
  #} else {
    samples <- mdply(levels(dd$Sample), function(i) data.frame(Name=paste0('Sample.', dd$Condition[dd$Sample==i][1], '.Sample.', i),Sample=i))
  #}
  samps.samples <- samps.Sol[seq(1,nrow(samps.Sol),1),colnames(samps.Sol) %in% samples$Name]
  samps.conditions <- samps.conditions[seq(1,nrow(samps.Sol),1),]
  colnames(samps.samples) <- samples$Sample[match(colnames(samps.samples), samples$Name)]
  # add Condition fixed effects
  samps.samples_plus_conditions <- samps.samples
  for (i in colnames(samps.samples)) {
    if(as.character(dd$Condition[dd$Sample==i][1]) %in% colnames(samps.conditions))
    {
      samps.samples_plus_conditions[,i] <- samps.samples_plus_conditions[,i] + samps.conditions[,as.character(dd$Condition[dd$Sample==i][1])]
    }
  }  
  
  samps.samples_plus_conditions.melted <- mdply(colnames(samps.samples_plus_conditions), function(i) {
    data.frame(Sample = i, Condition = dd$Condition[dd$Sample==i][1], value = samps.samples_plus_conditions[,i])
  })  
  stats.samples_plus_conditions <- ddply(samps.samples_plus_conditions.melted, .(Sample, Condition), function(x) {
    s <- data.frame(mean = mean(x$value), facet = ' ')
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps.samples_plus_conditions.melted.trunc <- ddply(samps.samples_plus_conditions.melted, .(Sample, Condition), function(x) {
    lower = stats.samples_plus_conditions$lower[stats.samples_plus_conditions$Sample == x$Sample[1]]
    upper = stats.samples_plus_conditions$upper[stats.samples_plus_conditions$Sample == x$Sample[1]]
    x[x$value >= lower & x$value <= upper,]
  })   
  
  samps.conditions.melted <- mdply(colnames(samps.conditions), function(i) {
      data.frame(Condition = i, value = samps.conditions[,i])
  }) 
  samps.conditions.melted <- merge(samps.conditions.melted, data.frame(Sample=dd$Sample, Condition=dd$Condition))
  stats.conditions <- ddply(samps.conditions.melted, .(Sample, Condition), function(x) {
    s <- data.frame(mean = mean(x$value), facet = ' ')
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps.conditions.melted.trunc <- ddply(samps.conditions.melted, .(Sample, Condition), function(x) {
    lower = stats.conditions$lower[stats.conditions$Sample == x$Sample[1]]
    upper = stats.conditions$upper[stats.conditions$Sample == x$Sample[1]]
    x[x$value >= lower & x$value <= upper,]
  })
  
  stats.conditions$Sample <- reorder(stats.conditions$Sample,as.numeric(stats.conditions$Condition))
  samps.conditions.melted.trunc$Sample <- reorder(samps.conditions.melted.trunc$Sample,as.numeric(samps.conditions.melted.trunc$Condition))
  stats.samples_plus_conditions$Sample <- reorder(stats.samples_plus_conditions$Sample,as.numeric(stats.samples_plus_conditions$Condition))
  samps.samples_plus_conditions.melted.trunc$Sample <- reorder(samps.samples_plus_conditions.melted.trunc$Sample,as.numeric(samps.samples_plus_conditions.melted.trunc$Condition))
  
  samps.conditions.melted.trunc$X1 <- 1
  samps.samples_plus_conditions.melted.trunc$X1 <- 0
  samps <- rbind(samps.conditions.melted.trunc, samps.samples_plus_conditions.melted.trunc)
  samps$X1 <- factor(samps$X1)
  
  samps$value[samps$Condition==levels(samps$Condition)[1] & samps$X1=="1"] <- NA
  
  ylim <- c(min(samps$value,na.rm = T)*1.5, max(samps$value, na.rm = T)*1.5)
  
  g <- ggplot(stats.samples_plus_conditions, aes(Sample, mean))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 axis.text.x=element_text(size=6),
                 legend.position="none")
  g <- g + facet_wrap(~ facet, ncol=1)
  g <- g + ggtitle(title(meta))
  g <- g + coord_cartesian(ylim=ylim)
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + geom_hline(yintercept=0,size=1/2,colour="darkgrey")          
  g <- g + geom_boxplot(data = samps, aes(y = value), alpha = 0.3, size = 0, outlier.size = 0)
  g <- g + geom_segment(data = stats.samples_plus_conditions, aes(x = as.integer(Sample)-0.45, xend = as.integer(Sample) + 0.45, y = mean, yend = mean),size = 1/2,colour="darkgrey")
  g <- g + geom_violin(data = samps, aes(y = value, colour = X1, alpha = X1, fill = Condition), position="identity", trim=T, size = 1/2)
  g <- g + geom_segment(data = stats.conditions, aes(x = as.integer(Sample)-0.5, xend = as.integer(Sample) + 0.5, y = mean, yend = mean),size = 1/2)
  g <- g + scale_alpha_manual(values=c(0.0,1.0))
  g <- g + scale_colour_manual(values=c("darkgrey","black"))
  ggsave(filename, g, height=2, width=6, limitsize=F)
  
  ylim
}


plot.model.peptides <- function(samps.Sol, meta, dd, filename) {    
  samps2 <- mdply(levels(dd$Peptide), function(i) {
    if (length(levels(dd$Peptide))==1) {
      s <- samps.Sol[,colnames(samps.Sol) %in% paste0('Digest.', levels(dd$Digest))]
      colnames(s) <- sub('Digest.', '', colnames(s), fixed = T)
    } else {
      s <- samps.Sol[,colnames(samps.Sol) %in% paste0('Digest.', i, '.Digest.', levels(dd$Digest))]
      colnames(s) <- sub(paste0('Digest.', i, '.Digest.'), '', colnames(s), fixed = T)
    }    
    # add Sample latent effects and melt
    mdply(colnames(s), function(j) {
      data.frame(Peptide = i,
                 Digest = j,
                 Condition = dd$Condition[dd$Digest==j][1],
                 value = s[,j])
    })
  })
  stats2 <- ddply(samps2, .(Peptide, Digest, Condition), function(x) {
    s <- data.frame(mean = mean(x$value))
    cbind(s, HPDinterval(mcmc(x$value)))
  })
  samps2.trunc <- ddply(samps2, .(Peptide, Digest, Condition), function(x) {
    lower = stats2$lower[stats2$Digest == x$Digest[1] & stats2$Peptide == x$Peptide[1]]
    upper = stats2$upper[stats2$Digest == x$Digest[1] & stats2$Peptide == x$Peptide[1]]
    x[x$value >= lower & x$value <= upper,]
  })
  
  stats2$Digest <- factor(as.character(stats2$Digest))
  samps2.trunc$Digest <- factor(as.character(samps2.trunc$Digest))
  N = maply(levels(stats2$Peptide), function (p) length(unique(dd$Spectrum[dd$Peptide==p])))
  stats2$Peptide <- sub(":.*:","",stats2$Peptide)
  samps2.trunc$Peptide <- sub(":.*:","",samps2.trunc$Peptide)
  stats2$Peptide <- factor(stats2$Peptide)
  samps2.trunc$Peptide <- factor(samps2.trunc$Peptide)
  
  levels(stats2$Peptide) <- levels(samps2.trunc$Peptide) <- paste0(levels(stats2$Peptide), ' [', N, ' spectr', ifelse(N==1,'um','a'), ']')
  
  stats2$Digest <- reorder(stats2$Digest,as.numeric(stats2$Condition))
  samps2.trunc$Digest <- reorder(samps2.trunc$Digest,as.numeric(samps2.trunc$Condition))  
  
  g <- ggplot(stats2, aes(Digest, mean))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank(),
                 strip.text=element_text(size=6),
                 axis.text.x=element_text(size=6),
                 legend.position="none")
  g <- g + coord_cartesian(ylim=c(-2,2))
  g <- g + facet_wrap(~Peptide,drop=F,ncol=1)
  g <- g + ggtitle(title(meta))
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + geom_hline(yintercept=0,size=1/2,colour="darkgrey")          
  g <- g + geom_violin(data = samps2.trunc, aes(y = value, fill = Condition), size = 1/2, trim=T)
  g <- g + geom_segment(aes(x = as.integer(Digest)-0.5, xend = as.integer(Digest) + 0.5, mean, yend = mean),size = 1/2) 
  ggsave(filename, g, height=1+1*length(levels(stats2$Peptide)), width=6, limitsize=F)
}
