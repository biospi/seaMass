invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[",Sys.time(), " Starting]"))

library(data.table)
library(ggplot2)
library(coda)
library(mcgibbsit)

# plot.conditions <- function(s.Sol, design, fc, filename) {
#   
#   test_conditions <- levels(design$Condition)[levels(design$Condition) != tolower(levels(design$Condition))]
#   test_conditions <- test_conditions[2:length(test_conditions)]
#   
#   samps <- s.Sol[,colnames(s.Sol) %in% paste0('Condition', levels(design$Condition)),drop=F]
#   samps.baseline <- data.frame(x = rnorm(nrow(samps),1e-6,1e-6))
#   colnames(samps.baseline)[1] <- paste0('Condition',levels(design$Condition)[1])
#   samps <- cbind(samps.baseline, samps)
#   colnames(samps) <- sub('Condition', '', colnames(samps), fixed=T)    
#   for (i in colnames(samps)) {
#     if (!(i %in% test_conditions))
#     {
#       samps[,i] <- rnorm(nrow(samps),1e-6,1e-6)
#     }
#   }
# 
#   samps <- samps[,colnames(samps) %in% test_conditions,drop=FALSE] 
# 
#   densities <- ddply(melt(samps, variable.name="Condition"), .(Condition), function(x)
#   {
#     dens <- density(x$value, n=65536)
#     data.frame(x=dens$x, y=dens$y)     
#   })   
#   densities$lower <- ifelse(densities$x<=log2(1/fc),densities$y,0) 
#   densities$upper <- ifelse(densities$x>=log2(fc),densities$y,0) 
#   y_range <- max(densities$y[densities$Condition %in% test_conditions])*1.9
#   x_range <- max(1.0,max(-min(densities$x[densities$y>y_range/100]),max(densities$x[densities$y>y_range/100])))
# 
#   stats <- data.frame(Condition = factor(colnames(samps), levels=levels(densities$Condition)),
#                       Up0 = 1 - colSums(samps > 0) / nrow(samps),
#                       Down0 = 1 - colSums(samps < 0) / nrow(samps),
#                       Up = 1 - colSums(samps > log2(fc)) / nrow(samps),
#                       Down = 1 - colSums(samps < -log2(fc)) / nrow(samps),
#                       Same = 1 - colSums(samps >= -log2(fc) & samps <= log2(fc)) / nrow(samps))
#   stats$mean = colMeans(samps)
#   stats <- cbind(stats, HPDinterval(mcmc(samps)))
#   stats$Up.text <- ifelse(stats$Condition %in% test_conditions, paste0("localFDR(up) = ",sapply(stats$Up, function(x) format(x,digits=2,scientific=F))), "") 
#   stats$Down.text <- ifelse(stats$Condition %in% test_conditions, paste0("localFDR(down) = ",sapply(stats$Down, function(x) format(x,digits=2,scientific=F))), "") 
#   stats$mean.text <- ifelse(stats$Condition %in% test_conditions, paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc "), "") 
#   stats$mean.hjust <- ifelse(stats$mean<0,0,1)
#   
#   g <- ggplot(stats,aes(mean,fill=Condition))
#   g <- g + theme_bw()
#   g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                  panel.grid.major=element_line(size=0.2),
#                  axis.ticks=element_blank(),
#                  axis.text.y=element_blank(),
#                  plot.title=element_text(size=10),
#                  plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                  strip.background=element_blank(),
#                  strip.text=element_text(size=10),
#                  legend.position="none")
#   g <- g + scale_x_continuous(expand = c(0,0))
#   g <- g + scale_y_continuous(expand = c(0,0))
#   g <- g + coord_cartesian(xlim=c(-x_range, x_range),ylim=c(-0.0,y_range))
#   g <- g + facet_wrap(~ Condition, ncol=1)
#   g <- g + xlab(expression('Log'[2]*' Ratio'))
#   g <- g + ylab("Probability Density")
#   g <- g + geom_rect(xmin=log2(1/fc),xmax=log2(fc),ymin=-2^32,ymax=2^32,alpha=0.15,colour="lightgrey",fill="lightgrey") 
#   g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
#   g <- g + geom_ribbon(data=densities,aes(x=x,ymax=lower),ymin=0)    
#   g <- g + geom_ribbon(data=densities,aes(x=x,ymax=upper),ymin=0)    
#   g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
#   g <- g + geom_vline(aes(xintercept=mean),size=2/3) 
#   g <- g + geom_linerange(aes(x=lower),ymin=0,ymax=y_range*0.7,size=2/3,lty=2)      
#   g <- g + geom_linerange(aes(x=upper),ymin=0,ymax=y_range*0.7,size=2/3,lty=2)   
#   g <- g + geom_text(aes(label=Up.text),x=x_range*0.98,y=y_range*0.94,hjust=1,vjust=1,size=3.5,colour='black')
#   g <- g + geom_text(aes(label=Down.text),x=-x_range*0.98,y=y_range*0.94,hjust=0,vjust=1,size=3.5,colour='black')    
#   g <- g + geom_text(aes(x=mean,label=mean.text,hjust=mean.hjust,colour=Condition),y=y_range*0.7,vjust=1,size=3.5)
#   ggsave(filename, g, height = 1+ 1*length(levels(stats$Condition)), width = 6, limitsize = F,device="png")
# }
# 
# 
# plot.conditions_sd <- function(s.VCV, design, filename) {   
# 
#   test_populations <- levels(design$Population)[levels(design$Population) != tolower(levels(design$Population))]
#   nsample <- count(design[!duplicated(design$Sample),]$Population)
#   test_populations <- test_populations[test_populations %in% nsample$x[nsample$freq>1]]
#   
#   if (length(levels(design$Population))==1) {
#     samps <- sqrt(s.VCV[,"Sample",drop=F])
#     colnames(samps) <- levels(design$Population)
#     stats <- data.frame(Population = levels(design$Population), mean = colMeans(samps))
#   } else {
#     samps <- sqrt(s.VCV[,colnames(s.VCV) %in% paste0("Population", levels(design$Population), ".Sample"),drop=F])
#     colnames(samps) <- sub('\\.Sample$', '', colnames(samps))    
#     colnames(samps) <- sub('Population', '', colnames(samps))    
#     stats <- data.frame(Population = factor(colnames(samps), levels=levels(design$Population)), mean = colMeans(samps))
#   }
#   stats <- cbind(stats, HPDinterval(mcmc(samps)))
# 
#   stats <- stats[stats$Population %in% test_populations,]
# 
#   densities <- ddply(melt(data.frame(samps), variable.name="Population"), .(Population), function(x)
#   {
#     dens <- density(x$value, n=65536)
#     data.frame(x=dens$x, y=dens$y)     
#   })   
#   y_range <- max(densities$y[densities$Population %in% test_populations])*1.4
#   x_range <- max(1,max(densities$x[densities$y>y_range/100]))
#   
#   stats$mean.text <- ifelse(stats$Population %in% test_populations, paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc"), "") 
#   stats$mean.hjust <- ifelse(stats$mean<0,0,1)
#   
#   levels(stats$Population) = sub('^[0-9]+', '', levels(stats$Population))
#   levels(densities$Population) = sub('^[0-9]+', '', levels(densities$Population))
#   
#   g <- ggplot(stats, aes(x=mean, fill=Population))
#   g <- g + theme_bw()
#   g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                  panel.grid.major=element_line(size=0.2),
#                  axis.ticks=element_blank(),
#                  axis.text.y=element_blank(),
#                  plot.title=element_text(size=10),
#                  plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                  strip.background=element_blank(),
#                  strip.text=element_text(size=10),
#                  legend.position="none")
#   g <- g + scale_x_continuous(expand = c(0,0))
#   g <- g + scale_y_continuous(expand = c(0,0))
#   g <- g + facet_wrap(~ Population, ncol=1)
#   g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
#   g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
#   g <- g + ylab("Probability Density")
#   g <- g + geom_vline(xintercept=0,size=2/3)          
#   g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0)    
#   g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
#   g <- g + geom_vline(aes(xintercept=mean),size=2/3) 
#   g <- g + geom_vline(aes(xintercept=lower),size=1/2,lty=2)      
#   g <- g + geom_vline(aes(xintercept=upper),size=1/2,lty=2)   
#   g <- g + geom_text(aes(x=mean,label=mean.text,colour=Population),y=y_range*0.9,hjust=0,vjust=1,size=3)
#   ggsave(filename, g, height = 1 + 1*length(levels(stats$Population)), width=3.5, limitsize=F,device="png")
# }
# 
# 
# plot.samples <- function(s.Sol, design, filename) {  
# 
#   samps.conditions <- s.Sol[,colnames(s.Sol) %in% paste0('Condition', levels(design$Condition)),drop=F]
#   samps.baseline <- data.frame(x = rnorm(nrow(samps.conditions),1e-6,1e-6))
#   colnames(samps.baseline)[1] <- paste0('Condition',levels(design$Condition)[1])
#   samps.conditions <- cbind(samps.baseline, samps.conditions)
#   colnames(samps.conditions) <- sub('Condition', '', colnames(samps.conditions), fixed=T)    
#   
#   if (length(levels(design$Population))==1) {
#     samples <- data.frame(Name=paste0('Sample.', levels(design$Sample)),Sample=levels(design$Sample))
#   } else {
#     samples <- mdply(levels(design$Sample), function(i) data.frame(Name=paste0('Population', design$Population[design$Sample==i][1], '.Sample.', i),Sample=i))
#   }
#   samps.samples <- s.Sol[seq(1,nrow(s.Sol),1),colnames(s.Sol) %in% samples$Name]
#   colnames(samps.samples) <- samples$Sample[match(colnames(samps.samples), samples$Name)]
#   samps.conditions <- samps.conditions[seq(1,nrow(s.Sol),1),]
#   
#   # add Condition fixed effects to samples
#   samps.samples_plus_conditions <- samps.samples
#   for (i in colnames(samps.samples)) {
#     if (as.character(design$Condition[design$Sample==i][1]) %in% colnames(samps.conditions))
#     {
#       samps.samples_plus_conditions[,i] <- samps.samples_plus_conditions[,i] + samps.conditions[,as.character(design$Condition[design$Sample==i][1])]
#     }
#   }  
#   
#   samps.samples_plus_conditions.melted <- mdply(colnames(samps.samples_plus_conditions), function(i) {
#     data.frame(Sample = i, Condition = design$Condition[design$Sample==i][1], value = samps.samples_plus_conditions[,i])
#   })  
#   stats.samples_plus_conditions <- ddply(samps.samples_plus_conditions.melted, .(Sample, Condition), function(x) {
#     s <- data.frame(mean = mean(x$value), facet = ' ')
#     cbind(s, HPDinterval(mcmc(x$value)))
#   })
#   samps.samples_plus_conditions.melted.trunc <- ddply(samps.samples_plus_conditions.melted, .(Sample, Condition), function(x) {
#     lower = stats.samples_plus_conditions$lower[stats.samples_plus_conditions$Sample == x$Sample[1]]
#     upper = stats.samples_plus_conditions$upper[stats.samples_plus_conditions$Sample == x$Sample[1]]
#     x[x$value >= lower & x$value <= upper,]
#   })   
#   
#   stats.samples_plus_conditions$Sample <- reorder(stats.samples_plus_conditions$Sample,as.numeric(stats.samples_plus_conditions$Condition))
#   samps.samples_plus_conditions.melted.trunc$Sample <- reorder(samps.samples_plus_conditions.melted.trunc$Sample,as.numeric(samps.samples_plus_conditions.melted.trunc$Condition))
#   
#   test_conditions <- levels(design$Condition)[levels(design$Condition) != tolower(levels(design$Condition))]
#   test_conditions <- test_conditions[2:length(test_conditions)]
# 
#   ymin <- min(samps.samples_plus_conditions.melted.trunc$value, na.rm = T)
#   ymax <- max(samps.samples_plus_conditions.melted.trunc$value, na.rm = T)
#   mid <- ymax - (ymax-ymin)/2
#   d <- (ymax-mid)*1.5
#   ylim <- c(mid+d,mid-d)
# 
#   g <- ggplot(stats.samples_plus_conditions, aes(Sample, mean))
#   g <- g + theme_bw()
#   g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                  panel.grid.major=element_line(size=0.2),
#                  axis.ticks=element_blank(),
#                  plot.title=element_text(size=10),
#                  plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                  strip.background=element_blank(),
#                  strip.text=element_text(size=6),
#                  axis.text.x=element_text(size=6,angle = 90, hjust = 1),
#                  legend.position="right")
#   g <- g + facet_wrap(~ facet, ncol=1)
#   g <- g + coord_cartesian(ylim=ylim)
#   g <- g + ylab(expression('Log'[2]*' Ratio'))
#   g <- g + geom_boxplot(data = samps.samples_plus_conditions.melted.trunc, aes(y = value), alpha = 0.0, weight = 0, colour = "white", size = 0, outlier.size = 0)
#   g <- g + geom_hline(yintercept=0,size=1/2,colour="darkgrey")          
#   g <- g + geom_violin(data = samps.samples_plus_conditions.melted.trunc, aes(y = value, fill = Condition), position="identity", trim=T, size = 1/2)
#   g <- g + geom_segment(data = stats.samples_plus_conditions, aes(x = as.integer(Sample)-0.45, xend = as.integer(Sample) + 0.45, y = mean, yend = mean),size = 1/2)
#   ggsave(filename, g, height=2, width=8, limitsize=F,device="png")
#   
#   f = gsub("samples\\/","",gsub("\\.png","",filename))
#   save(stats.samples_plus_conditions, file=paste0("samplestats/",f,".Rdata"))  
# 
#   ylim
# }
# 
# 
# plot.peptides_sd <- function(s.VCV, design, filename) {   
#   if (sum(grepl("^Peptide.*\\.Digest$", colnames(s.VCV))) > 0) {
# 
#     samps <- sqrt(s.VCV[,grepl("^Peptide.*\\.Digest$", colnames(s.VCV)),drop=F])
#     colnames(samps) <- sub("^Peptide", '', colnames(samps))      
#     colnames(samps) <- sub("\\.Digest$", '', colnames(samps))      
#     
#     stats <- data.frame(Peptide = colnames(samps), mean = colMeans(samps))
#     stats <- cbind(stats, HPDinterval(mcmc(samps)))    
#     
#     densities <- ddply(melt(data.frame(samps), variable.name="Peptide"), .(Peptide), function(x)
#     {
#       dens <- density(x$value, n=4096)
#       data.frame(x=dens$x, y=dens$y)     
#     })   
#     levels(densities$Peptide) <- colnames(samps)
#     y_range <- max(densities$y)*1.4
#     
#     stats$mean.text <- paste0(" ",sapply(stats$mean, function(x) format(ifelse(x<0,-1/2^x,2^x),digits=3,scientific=F)),"fc ")
#     stats$mean.hjust <- ifelse(stats$mean<0,0,1)
#     
#     #stats$Peptide <- sub(":.*:","",stats$Peptide)
#     #densities$Peptide <- sub(":.*:","",densities$Peptide)
#     #stats$Peptide <- factor(stats$Peptide)
#     #densities$Peptide <- factor(densities$Peptide)
#     #levels(stats$Peptide) <- paste0(levels(stats$Peptide), ' []')
#     #levels(densities$Peptide) <- paste0(levels(densities$Peptide), ' []')
#     
#     g <- ggplot(stats, aes(x=mean, colour=Peptide, fill=Peptide))
#     g <- g + theme_bw()
#     g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                    panel.grid.major=element_line(size=0.2),
#                    axis.ticks=element_blank(),
#                    axis.text.y=element_blank(),
#                    plot.title=element_text(size=10),
#                    plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                    strip.background=element_blank(),
#                    strip.text=element_text(size=4),
#                    legend.position="none")
#     g <- g + scale_x_continuous(expand = c(0,0))
#     g <- g + scale_y_continuous(expand = c(0,0))
#     g <- g + facet_wrap(~ Peptide, ncol=1)
#     g <- g + coord_cartesian(xlim=c(0,1),ylim=c(-0.0,y_range))
#     g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
#     g <- g + ylab("Probability Density")
#     g <- g + geom_vline(xintercept=0,size=2/3,colour="darkgrey")          
#     g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
#     g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
#     g <- g + geom_vline(aes(xintercept=mean),size=2/3) 
#     g <- g + geom_vline(aes(xintercept=lower),size=1/2,lty=2)      
#     g <- g + geom_vline(aes(xintercept=upper),size=1/2,lty=2)   
#     g <- g + geom_text(aes(x=mean,label=mean.text),y=y_range*0.9,hjust=0,vjust=1,size=3)
#     ggsave(filename, g, height=1+1*length(levels(stats$Peptide)), width=3.5, limitsize=F,device="png")  
#   }
# }
# 
# 
# plot.peptides <- function(s.Sol, design, filename, ylim) { 
#   if(sum(grepl("^Peptide.*\\.Digest.*$", colnames(s.Sol))) > 0) {
# 
#     samps.Sol <- as.data.frame(s.Sol[,grepl("^Peptide.*\\.Digest.*$", colnames(s.Sol)),drop=F])
#     samps <- mdply(colnames(samps.Sol), function(s) {
#       s2 <- sub("^Peptide", '', s) 
#       peptide <- sub("\\.Digest.*$", '', s2)   
#       digest <- sub("^.*\\.Digest.", '', s2)   
#       data.frame(
#         Peptide = peptide,
#         Digest = digest,
#         Condition = design$Condition[design$Digest==digest][1],
#         value = samps.Sol[,s]
#       )
#     })
#   
#     stats <- ddply(samps, .(Peptide, Digest, Condition), function(x) {
#       s <- data.frame(mean = mean(x$value))
#       cbind(s, HPDinterval(mcmc(x$value)))
#     })
#     
#     samps.trunc <- ddply(samps, .(Peptide, Digest, Condition), function(s) {
#       lower = stats$lower[stats$Digest == s$Digest[1] & stats$Peptide == s$Peptide[1]]
#       upper = stats$upper[stats$Digest == s$Digest[1] & stats$Peptide == s$Peptide[1]]
#       s[s$value >= lower & s$value <= upper,]
#     })  
#     
#     stats$Digest <- reorder(stats$Digest,as.numeric(stats$Condition))
#     samps.trunc$Digest <- reorder(samps.trunc$Digest,as.numeric(samps.trunc$Condition))
#     stats$Peptide <- factor(stats$Peptide, levels=levels(samps$Peptide)[order(levels(samps$Peptide))])
#     
#     g <- ggplot(stats, aes(Digest, mean))
#     g <- g + theme_bw()
#     g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                    panel.grid.major=element_line(size=0.2),
#                    axis.ticks=element_blank(),
#                    plot.title=element_text(size=10),
#                    plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                    strip.background=element_blank(),
#                    strip.text=element_text(size=6),
#                    axis.text.x=element_text(size=6,angle = 90, hjust = 1),
#                    legend.position="right")
#     g <- g + coord_cartesian(ylim=ylim)
#     g <- g + facet_wrap(~Peptide,drop=F,ncol=1)
#     g <- g + ylab(expression('Log'[2]*' Ratio'))
#     g <- g + geom_boxplot(data = samps.trunc, aes(y = value), alpha = 0.0, weight = 0, colour = "white", size = 0, outlier.size = 0)
#     g <- g + geom_hline(yintercept=0,size=1/2,colour="darkgrey")          
#     g <- g + geom_violin(data = samps.trunc, aes(y = value, fill = Condition), position="identity", size = 1/2, trim=T)
#     g <- g + geom_segment(aes(x = as.integer(Digest)-0.45, xend = as.integer(Digest) + 0.45, mean, yend = mean),size = 1/2) 
#     ggsave(filename, g, height=1+1*length(levels(stats$Peptide)), width=8, limitsize=F,device="png")
#   }
# }
# 
# 
# plot.digests_sd <- function(s.VCV, design, filename) {   
#   if (sum(grepl("^Digest.*\\.Peptide$", colnames(s.VCV))) > 0) {
# 
#     samps <- sqrt(s.VCV[,grepl("^Digest.*\\.Peptide$", colnames(s.VCV)),drop=F])
#     colnames(samps) <- sub("^Digest", '', colnames(samps))      
#     colnames(samps) <- sub("\\.Peptide$", '', colnames(samps))      
#     
#     samps.melted <- mdply(colnames(samps), function(i) {
#       data.frame(Digest = i, Condition = design$Condition[design$Digest==i][1], value = samps[,i])
#     })  
#     stats <- ddply(samps.melted, .(Digest, Condition), function(x) {
#       s <- data.frame(mean = mean(x$value), facet = ' ')
#       cbind(s, HPDinterval(mcmc(x$value)))
#     })
#     samps.melted.trunc <- ddply(samps.melted, .(Digest, Condition), function(x) {
#       lower = stats$lower[stats$Digest == x$Digest[1]]
#       upper = stats$upper[stats$Digest == x$Digest[1]]
#       x[x$value >= lower & x$value <= upper,]
#     })  
#  
#     stats$Digest <- reorder(stats$Digest,as.numeric(stats$Condition))
#     samps.melted.trunc$Digest <- reorder(samps.melted.trunc$Digest,as.numeric(samps.melted.trunc$Condition))
#     
#     g <- ggplot(stats, aes(Digest, mean))
#     g <- g + theme_bw()
#     g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                    panel.grid.major=element_line(size=0.2),
#                    axis.ticks=element_blank(),
#                    plot.title=element_text(size=10),
#                    plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                    strip.background=element_blank(),
#                    strip.text=element_text(size=6),
#                    axis.text.x=element_text(size=6),
#                    legend.position="none")
#     g <- g + facet_wrap(~ facet, ncol=1)
#     g <- g + scale_y_continuous(expand = c(0,0))
#     g <- g + ylab(expression('Log'[2]*' Std Dev'))
#     g <- g + geom_boxplot(data = samps.melted.trunc, aes(y = value), alpha = 0.0, weight = 0, colour = "white", size = 0, outlier.size = 0)
#     g <- g + geom_hline(yintercept=0,size=1/2,colour="darkgrey")          
#     g <- g + geom_violin(data = samps.melted.trunc, aes(y = value, fill = Condition), position="identity", trim=T, size = 1/2)
#     g <- g + geom_segment(aes(x = as.integer(Digest)-0.45, xend = as.integer(Digest) + 0.45, y = mean, yend = mean),size = 1/2)
#     ggsave(filename, g, height=2, width=6, limitsize=F,device="png")
#   }
# }

calc.stats.conditions <- function(iterations.Sol, samples.Sol, design, nitt, nburnin, nchain, tol, modeltime, dic) { 
    
  # one-sided statistical tests, checking precision
  test_conditions <- levels(design$Condition)[levels(design$Condition) != tolower(levels(design$Condition))]
  test_conditions <- test_conditions[2:length(test_conditions)]
  
  stats <- data.table(
    Condition=test_conditions
  )  
  stats.func <- function(con) {
    if (paste0('Condition', con) %in% names(samples.Sol[1,])) {
      
      iterations <- iterations.Sol[,paste0('Condition', con)]
      samples <- samples.Sol[,paste0('Condition', con)]
      samples.mean <- mean(samples)
      samples.hpdi <- HPDinterval(mcmc(samples))

      qs <- data.table(
        Test=c("Up","Down"),
        q=pmin(pmax(c(sum(samples > 0.0), sum(samples < 0.0))  ,1),length(samples)-1) / length(samples)
      )
      
      tests.func <- function(q) {
        res <- mcgibbsit(iterations, q, tol)
        
        nitt_pred <- NA
        try(nitt_pred <- max(ceiling(res$resmatrix[,'Total'] / nchain),1), silent=T)
        nburnin_pred <- NA
        try(nburnin_pred <- max(ceiling(res$resmatrix[,'M'] / nchain),1), silent=T)
        
        data.table(
          ModelTime = modeltime,
          DIC = dic,
          itt_pred = nitt_pred,
          itt_ok = ifelse(nitt>=nitt_pred,"Yes","No"),
          burnin_pred = nburnin_pred,
          burnin_ok = ifelse(nburnin>=nburnin_pred,"Yes","No"),
          mean = samples.mean,
          lower = samples.hpdi[1],
          upper = samples.hpdi[2],
          localFDR = 1-q
        )
      }
      qs <- qs[,as.list(tests.func(q)), by=Test]
      
    } else {
      data.table()
    }
  }
  
  stats <- stats[,as.list(stats.func(Condition)), by=Condition]
}


calc.stats.samples <- function(s.Sol, design) {

  samps.conditions <- s.Sol[,colnames(s.Sol) %in% paste0('Condition', levels(design$Condition)),drop=F]
  samps.baseline <- data.table(x = rnorm(nrow(samps.conditions),1e-6,1e-6))
  colnames(samps.baseline)[1] <- paste0('Condition',levels(design$Condition)[1])
  samps.conditions <- cbind(samps.baseline, samps.conditions)
  colnames(samps.conditions) <- sub('Condition', '', colnames(samps.conditions), fixed=T)

  samples <- data.table(Sample=levels(design$Sample))
  if (length(levels(design$Population))==1) {
    samples$Name = paste0('Sample.', levels(design$Sample))
  } else {
    samples$Name <- sapply(samples$Sample, function(i) paste0('Population', design$Population[design$Sample==i][1], '.Sample.', i))
  }
  samps.samples <- s.Sol[seq(1,nrow(s.Sol),1),colnames(s.Sol) %in% samples$Name]
  colnames(samps.samples) <- samples$Sample[match(colnames(samps.samples), samples$Name)]
  samps.conditions <- samps.conditions[seq(1,nrow(s.Sol),1),]

  # add Condition fixed effects to samples
  samps.samples_plus_conditions <- samps.samples
  for (i in colnames(samps.samples)) {
    if (as.character(design$Condition[design$Sample==i][1]) %in% colnames(samps.conditions)) {
      samps.samples_plus_conditions[,i] <- samps.samples_plus_conditions[,i] + samps.conditions[,get(as.character(design$Condition[design$Sample==i][1]))]
    }
  }
  
  stats.mean <-colMeans(samps.samples_plus_conditions)
  names(stats.mean) <- paste0(names(stats.mean), ".mean")
  
  stats.sd <- sapply(1:ncol(samps.samples_plus_conditions), function(i) sd(samps.samples_plus_conditions[,i]))
  names(stats.sd) <- paste0(names(stats.mean), ".sd")

  as.data.table(t(c(stats.mean, stats.sd)))
}


args <- commandArgs(T)
if (length(args) == 0) args <- c("3") # test case

# some tuning parameters (should come from parameters.Rdata with defaults given here)
prefix <- ifelse(file.exists("parameters.Rdata"),".",file.path("..","..","input"))
load(file.path(prefix,"design.Rdata"))
load(file.path(prefix,"parameters.Rdata"))
nitt <- as.integer(ifelse("model_iterations" %in% parameters$Key,parameters$Value[parameters$Key=="model_iterations"],13000))
nburnin <- as.integer(ifelse("model_warmup" %in% parameters$Key,parameters$Value[parameters$Key=="model_warmup"],3000))
nchain <- as.integer(ifelse("model_chains" %in% parameters$Key,parameters$Value[parameters$Key=="model_chains"],100))
tol <- as.double(ifelse("model_tol" %in% parameters$Key,parameters$Value[parameters$Key=="model_tol"],0.0125))
do_plots <- as.integer(ifelse("do_plots" %in% parameters$Key,parameters$Value[parameters$Key=="do_plots"],1))

# if(!dir.exists("stats")) {
#   dir.create("conditions")
#   dir.create("conditions_sd")
#   dir.create("samples")
#   dir.create("peptides")
#   dir.create("peptides_sd")
#   dir.create("digests_sd")
# }

# get files for this batch
prefix <- "results"
files <- list.files(prefix,pattern=paste0("^",args[1],"\\.[0-9]+\\.Rdata"))
if (length(files) < nchain) {
  prefix = file.path("..","..","model","results")
  files <- list.files(prefix,pattern=paste0("^",args[1],"\\.[0-9]+\\.Rdata"))
  
  if (length(files) < nchain) stop(paste0("ERROR: Missing model input"))
}

# load files
for (f in files) {
  message("[",paste0(Sys.time()," Processing ",f,"...]"))
  load(file.path(prefix,f))
  chain <- as.integer(gsub("\\.Rdata$","",gsub("^[0-9]+\\.","",f)))
  
  if (f == files[1]) {
    proteins <- vector("list", length(samples.Sol))
    names(proteins) <- names(samples.Sol)
    for (p in names(proteins)) proteins[[p]] <- vector("list", length(nchain))
  }

  for (p in names(proteins)) {
    proteins[[p]][[chain+1]] <- list(DIC=dic[[p]], ModelTime=time[[p]], samples.Sol=samples.Sol[[p]], samples.VCV=samples.VCV[[p]])
  }
}

# combine across chains, compute stats and construct plots for each protein
stats.conditions <- vector("list", length(proteins))
names(stats.conditions) <- names(proteins)
stats.samples <- vector("list", length(proteins))
names(stats.samples) <- names(proteins)
for (p in names(proteins)) {
  message(paste0("[", Sys.time()," Processing protein ", p,"]"))    
  
  dic <- mean(sapply(1:nchain, function(x) proteins[[p]][[x]]$DIC))
  time <- sum(sapply(1:nchain, function(x) proteins[[p]][[x]]$ModelTime))
  samples.Sol <- mcmc.list(lapply(1:nchain, function(x) proteins[[p]][[x]]$samples.Sol))
  samples.VCV <- mcmc.list(lapply(1:nchain, function(x) proteins[[p]][[x]]$samples.VCV))
  
  proteins[[p]] <- list(DIC=dic, ModelTime=time, samples.Sol=samples.Sol, samples.VCV=samples.VCV)
  
  # mcgibbsit needs burnin and samples across all chains
  iterations.Sol <- proteins[[p]][["samples.Sol"]] 
  
  # merge samples across all chains for stats and plots 
  samples.Sol <- as.matrix(window(proteins[[p]][["samples.Sol"]],nburnin+1,end(proteins[[p]][["samples.Sol"]])))
  samples.Sol <- samples.Sol/log(2) #transforming to log2
  samples.VCV <- as.matrix(window(proteins[[p]][["samples.VCV"]],nburnin+1,end(proteins[[p]][["samples.VCV"]])))
  samples.VCV <- samples.VCV/log(2) #transforming to log2  
  
  # conditions
  stats.conditions[[p]] <- calc.stats.conditions(iterations.Sol, samples.Sol, design, nitt, nburnin, nchain, tol, proteins[[p]][["ModelTime"]], proteins[[p]][["DIC"]])
  # samples
  stats.samples[[p]] <- calc.stats.samples(samples.Sol, design)
  
  proteins[[p]] <- NULL
}
save(stats.conditions, stats.samples, file=paste0(args[1],".Rdata"))

message(paste0("[",Sys.time()," Finished]"))

