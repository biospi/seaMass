Sys.setlocale("LC_COLLATE","C")

# FOR EXECUTING UNDER HPC
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HPC")
{
  print(paste(Sys.time(),"[Starting]"))
  
  library(plyr)
  library(reshape2)
  library(ggplot2)
  library(coda)
  
  load('index.Rdata')  
  load('design.Rdata')
  load('parameters.Rdata')  

  files = list.files(path="stats",pattern="^[0-9]+\\.Rdata")
  
  # conditions
  
  test_conditions <- levels(design$Condition)[levels(design$Condition) != tolower(levels(design$Condition))]
  test_conditions <- test_conditions[2:length(test_conditions)]
  
  samps <- mdply(files, .id=NULL, function(f) {
    load(paste0("stats/",f))
    samps <- data.frame(t(colMeans(s.Sol[,colnames(s.Sol) %in% paste0('Condition', levels(design$Condition)),drop=F])))
    colnames(samps) <- sub('Condition', '', colnames(samps))    
    samps$ProteinID <- factor(as.integer(gsub("\\.Rdata","",f)))
    #samps$itt <- seq(1,nrow(samps))
    samps
  }) 
  
  densities <- ddply(melt(data.frame(samps), variable.name="Condition"), .(Condition), function(x)
  {
    dens <- density(x$value, n=65536, na.rm=T)
    data.frame(x=dens$x, y=dens$y)     
  })   
  y_range <- max(densities$y[densities$Condition %in% test_conditions])*1.4
  x_range <- max(1,max(densities$x[densities$y>y_range/100]))
  
  g <- ggplot(densities, aes(fill=Condition))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10),
                 legend.position="none")
  g <- g + scale_x_continuous(expand = c(0,0))
  g <- g + scale_y_continuous(expand = c(0,0))
  g <- g + facet_wrap(~ Condition, ncol=1)
  g <- g + coord_cartesian(xlim=c(-x_range,x_range),ylim=c(-0.0,y_range))
  g <- g + xlab(expression('Log'[2]*' Ratio'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3)          
  g <- g + geom_ribbon(aes(x=x,ymax=y),ymin=0,alpha=0.3)    
  g <- g + geom_line(aes(x=x,y=y),size=2/3) 
  ggsave("study_conditions.png", g, height = 1 + 1*length(levels(densities$Condition)), width=6, limitsize=F, device="png")
 
  
  
 
  # conditions_sd
  
  test_populations <- levels(design$Population)[levels(design$Population) != tolower(levels(design$Population))]
  nsample <- count(design[!duplicated(design$Sample),]$Population)
  test_populations <- test_populations[test_populations %in% nsample$x[nsample$freq>1]]
  
  samps <- mdply(files, .id=NULL, function(f) {
    load(paste0("stats/",f))
    if (length(levels(design$Population))==1) {
      samps <- data.frame(t(colMeans(sqrt(s.VCV[,"Sample",drop=F]))))
      colnames(samps) <- levels(design$Population)
    } else {
      samps <- data.frame(t(colMeans(sqrt(s.VCV[,colnames(s.VCV) %in% paste0("Population", levels(design$Population), ".Sample"),drop=F]))))
      colnames(samps) <- sub('\\.Sample$', '', colnames(samps))    
      colnames(samps) <- sub('Population', '', colnames(samps))    
    }
    samps
  }) 

  samps <- samps[,test_populations,drop=F]
 
  if (length(levels(design$Population))==1) {
    stats <- data.frame(Population = levels(design$Population), mean = colMeans(samps))
  } else {
    stats <- data.frame(Population = factor(colnames(samps), levels=levels(design$Population)), mean = colMeans(samps))
  }
  stats <- cbind(stats, HPDinterval(mcmc(samps)))  
  
  densities <- ddply(melt(data.frame(samps), variable.name="Population"), .(Population), function(x)
  {
    dens <- density(x$value, n=65536)
    data.frame(x=dens$x, y=dens$y)     
  })   
  y_range <- max(densities$y[densities$Population %in% test_populations])*0.5
  x_range <- max(1,max(densities$x[densities$y>y_range/100]))
  
  levels(densities$Population) = sub('^[0-9]+', '', levels(densities$Population))
  
  g <- ggplot(stats, aes(x=mean, fill=Population))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10),
                 legend.position="none")
  g <- g + scale_x_continuous(expand = c(0,0))
  g <- g + scale_y_continuous(expand = c(0,0))
  g <- g + facet_wrap(~ Population, ncol=1)
  g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
  g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3)          
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  #g <- g + geom_vline(aes(xintercept=mean),size=2/3) 
  #g <- g + geom_vline(aes(xintercept=lower),size=1/2,lty=2)      
  #g <- g + geom_vline(aes(xintercept=upper),size=1/2,lty=2)   
  ggsave("study_populations.png", g, height = 1 + 1*length(levels(densities$Population)), width=6, limitsize=F, device="png")
  
  
  # peptide sd
 
  samps <- mdply(files, .id=NULL, function(f) {
    load(paste0("stats/",f))
    samps <- data.frame(colMeans(sqrt(s.VCV[,grepl("^Peptide.*\\.Digest$", colnames(s.VCV)),drop=F])))
    #colnames(samps) <- "Peptides"
    #if (nrow(samps > 0)) samps$ProteinID <- factor(as.integer(gsub("\\.Rdata","",f)))
    samps
  })  
  
  densities <- ddply(melt(data.frame(samps), variable.name="Population"), .(Population), function(x)
  {
    dens <- density(x$value, n=65536)
    data.frame(x=dens$x, y=dens$y)     
  })   
  
  y_range <- max(densities$y)*1.4
  g <- ggplot(densities, aes(x=mean))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
                 panel.grid.major=element_line(size=0.2),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
                 strip.background=element_blank(),
                 strip.text=element_text(size=10),
                 legend.position="none")
  g <- g + scale_x_continuous(expand = c(0,0))
  g <- g + scale_y_continuous(expand = c(0,0))
  g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
  g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=2/3)          
  g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
  g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
  ggsave("study_peptides.png", g, height = 2, width=6, limitsize=F,device="png")
  
  
  # output csv
  
  stats <- mdply(files, .id=NULL, function(f) {
    load(paste0("stats/",f))
    stats$ProteinID <- as.integer(gsub("\\.Rdata","",f))
    stats
  })
  stats$Condition <- factor(stats$Condition)
  
  results <- merge(data.index, stats)
  for (con in levels(stats$Condition))
  {
    out <- results[results$Condition == con & (results$Test == "Up" | results$Test == "Down"),]
    out <- out[order(out$localFDR),]
    out <- out[!duplicated(out$ProteinID),]
    out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
    write.csv(out, paste0(parameters$Value[parameters$Key=="id"],"_",con,"_UpDown.csv"), row.names=F)
    
    out <- results[results$Condition == con & results$Test == "Same",]
    out <- out[order(out$localFDR),]
    out <- out[!duplicated(out$ProteinID),]
    out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
    write.csv(out, paste0(parameters$Value[parameters$Key=="id"],"_",con,"_Same.csv"), row.names=F)
  }

  print(paste(Sys.time(),"[Finished]"))
}
