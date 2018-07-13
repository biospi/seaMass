invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[",Sys.time(), " Starting]"))

library(methods)
library(data.table)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(coda)

# some tuning parameters (should come from parameters.Rdata with defaults given here)
prefix <- ifelse(file.exists("parameters.Rdata"),".",file.path("..","..","input"))
load(file.path(prefix,"design.Rdata"))
load(file.path(prefix,"parameters.Rdata"))
nbatch <- as.integer(ifelse("batches" %in% parameters$Key,parameters$Value[parameters$Key=="batches"],10))

prefix <- ifelse(file.exists("exposures.Rdata"),".",file.path("..","..","exposures","results"))
load(file.path(prefix,"exposures.Rdata"))

# get stats files
prefix <- "."
files <- list.files(prefix,pattern="^[0-9]+\\.Rdata")
if (length(files) < nbatch) {
  prefix = file.path("..","..","plots","results")
  files <- list.files(prefix,pattern="^[0-9]+\\.Rdata")
  
  if (length(files) < nbatch) stop(paste0("ERROR: Missing model input"))
}

# load stats
stats.samples.all <- vector("list", nbatch)
stats.conditions.all <- vector("list", nbatch)
for (i in 1:nbatch) {
  load(file.path(prefix, files[i]))
  stats.samples.all[[i]] <- rbindlist(stats.samples, idcol="ProteinID", fill=T)
  stats.conditions.all[[i]] <- rbindlist(stats.conditions, idcol="ProteinID")
}
stats.samples.all <- rbindlist(stats.samples.all)
stats.samples.all$ProteinID <- factor(stats.samples.all$ProteinID)
stats.samples.all <- merge(dd.index, stats.samples.all)
stats.conditions.all <- rbindlist(stats.conditions.all)
stats.conditions.all$ProteinID <- factor(stats.conditions.all$ProteinID)
stats.conditions.all$Condition <- factor(stats.conditions.all$Condition)
stats.conditions.all <- merge(dd.index, stats.conditions.all)

# write stats.samples
stats.samples.all <- stats.samples.all[order(as.numeric(as.character(stats.samples.all$ProteinID))),]
fwrite(stats.samples.all, paste0(parameters$Value[parameters$Key=="id"],"_Samples.csv"))

# write stats.conditions
test_conditions <- levels(design$Condition)[levels(design$Condition) != tolower(levels(design$Condition))]
test_conditions <- test_conditions[2:length(test_conditions)]
for (con in levels(stats.conditions.all$Condition))
{
  out <- stats.conditions.all[stats.conditions.all$Condition == con,]
  out$Condition <- NULL
  out <- out[order(out$localFDR),]
  out <- out[!duplicated(out$ProteinID),]
  out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
  
  fwrite(out, paste0(parameters$Value[parameters$Key=="id"],"_Condition_", con, ".csv"))
}

# diagnostic plots

# generate wide quant data.table with samples as row names
dd.samples <- as.matrix(stats.samples.all[,grepl("\\.mean$", colnames(stats.samples.all)), with=F])
samples <- gsub("\\.mean$", "", colnames(dd.samples))
rownames(dd.samples) <- as.character(stats.samples.all$ProteinID)
dd.samples <- as.data.table(t(na.omit(dd.samples)))
dd.samples$Sample <- samples

# PCA
pca.samples <- prcomp(dd.samples[, !"Sample", with = F], center = F, scale = F)
dd.pca.samples <- fortify(pca.samples)
dd.pca.samples$Sample <- dd.samples$Sample
dd.pca.samples <- merge(dd.pca.samples, design[!duplicated(Sample), list(Sample, Condition)])

g <- autoplot(pca.samples, data = dd.pca.samples, scale = 0, colour = "Condition")
g <- g + geom_label_repel(aes(label = Sample, colour = Condition))
g <- g + theme(aspect.ratio=1) + coord_equal() 

ggsave(paste0(parameters$Value[parameters$Key=="id"],"_PCA.pdf"), g, width=8, height=8)


# samps <- mdply(files, .id=NULL, function(f) {
#   load(paste0("stats/",f))
#   tmp <- colMeans(s.Sol[,colnames(s.Sol) %in% paste0('Condition', levels(design$Condition)),drop=F])
#   samps <- data.frame(t(tmp))
#   colnames(samps) <- names(tmp)
#   colnames(samps) <- sub('Condition', '', colnames(samps))    
#   samps$ProteinID <- factor(as.integer(gsub("\\.Rdata","",f)))
#   #samps$itt <- seq(1,nrow(samps))
#   samps
# }) 

# densities <- ddply(melt(samps, variable.name="Condition"), .(Condition), function(x)
# {
#   dens <- density(x$value, n=65536, na.rm=T)
#   data.frame(x=dens$x, y=dens$y)     
# })   
# y_range <- max(densities$y[densities$Condition %in% test_conditions])*1.4
# x_range <- max(1,max(densities$x[densities$y>y_range/100]))
# 
# g <- ggplot(densities, aes(fill=Condition))
# g <- g + theme_bw()
# g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                panel.grid.major=element_line(size=0.2),
#                axis.ticks=element_blank(),
#                axis.text.y=element_blank(),
#                plot.title=element_text(size=10),
#                plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                strip.background=element_blank(),
#                strip.text=element_text(size=10),
#                legend.position="none")
# g <- g + scale_x_continuous(expand = c(0,0))
# g <- g + scale_y_continuous(expand = c(0,0))
# g <- g + facet_wrap(~ Condition, ncol=1)
# g <- g + coord_cartesian(xlim=c(-x_range,x_range),ylim=c(-0.0,y_range))
# g <- g + xlab(expression('Log'[2]*' Ratio'))
# g <- g + ylab("Probability Density")
# g <- g + geom_vline(xintercept=0,size=2/3)          
# g <- g + geom_ribbon(aes(x=x,ymax=y),ymin=0,alpha=0.3)    
# g <- g + geom_line(aes(x=x,y=y),size=2/3) 
# ggsave("study_conditions.png", g, height = 1 + 1*length(levels(densities$Condition)), width=6, limitsize=F, device="png")


# # conditions_sd
# test_populations <- levels(design$Population)[levels(design$Population) != tolower(levels(design$Population))]
# nsample <- count(design[!duplicated(design$Sample),]$Population)
# test_populations <- test_populations[test_populations %in% nsample$x[nsample$freq>1]]
# 
# samps <- mdply(files, .id=NULL, function(f) {
#   load(paste0("stats/",f))
#   if (length(levels(design$Population))==1) {
#     samps <- data.frame(t(colMeans(sqrt(s.VCV[,"Sample",drop=F]))))
#     colnames(samps) <- levels(design$Population)
#   } else {
#     samps <- data.frame(t(colMeans(sqrt(s.VCV[,colnames(s.VCV) %in% paste0("Population", levels(design$Population), ".Sample"),drop=F]))))
#     colnames(samps) <- sub('\\.Sample$', '', colnames(samps))    
#     colnames(samps) <- sub('Population', '', colnames(samps))    
#   }
#   samps
# }) 
# 
# samps <- samps[,test_populations,drop=F]
# 
# if (length(levels(design$Population))==1) {
#   stats <- data.frame(Population = levels(design$Population), mean = colMeans(samps))
# } else {
#   stats <- data.frame(Population = factor(colnames(samps), levels=levels(design$Population)), mean = colMeans(samps))
# }
# stats <- cbind(stats, HPDinterval(mcmc(samps)))  
# 
# densities <- ddply(melt(data.frame(samps), variable.name="Population"), .(Population), function(x)
# {
#   dens <- density(x$value, n=65536)
#   data.frame(x=dens$x, y=dens$y)     
# })   
# y_range <- max(densities$y[densities$Population %in% test_populations])*0.5
# x_range <- max(1,max(densities$x[densities$y>y_range/100]))
# 
# levels(densities$Population) = sub('^[0-9]+', '', levels(densities$Population))
# 
# g <- ggplot(stats, aes(x=mean, fill=Population))
# g <- g + theme_bw()
# g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                panel.grid.major=element_line(size=0.2),
#                axis.ticks=element_blank(),
#                axis.text.y=element_blank(),
#                plot.title=element_text(size=10),
#                plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                strip.background=element_blank(),
#                strip.text=element_text(size=10),
#                legend.position="none")
# g <- g + scale_x_continuous(expand = c(0,0))
# g <- g + scale_y_continuous(expand = c(0,0))
# g <- g + facet_wrap(~ Population, ncol=1)
# g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
# g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
# g <- g + ylab("Probability Density")
# g <- g + geom_vline(xintercept=0,size=2/3)          
# g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
# g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
# #g <- g + geom_vline(aes(xintercept=mean),size=2/3) 
# #g <- g + geom_vline(aes(xintercept=lower),size=1/2,lty=2)      
# #g <- g + geom_vline(aes(xintercept=upper),size=1/2,lty=2)   
# ggsave("study_populations.png", g, height = 1 + 1*length(levels(densities$Population)), width=6, limitsize=F, device="png")
# 
# 
# # peptide sd
# 
# samps <- mdply(files, .id=NULL, function(f) {
#   load(paste0("stats/",f))
#   samps <- data.frame(colMeans(sqrt(s.VCV[,grepl("^Peptide.*\\.Digest$", colnames(s.VCV)),drop=F])))
#   #colnames(samps) <- "Peptides"
#   #if (nrow(samps > 0)) samps$ProteinID <- factor(as.integer(gsub("\\.Rdata","",f)))
#   samps
# })  
# 
# densities <- ddply(melt(data.frame(samps), variable.name="Population"), .(Population), function(x)
# {
#   dens <- density(x$value, n=65536)
#   data.frame(x=dens$x, y=dens$y)     
# })   
# 
# y_range <- max(densities$y)*1.4
# g <- ggplot(densities, aes(x=mean))
# g <- g + theme_bw()
# g <- g + theme(panel.border=element_rect(colour="black",size=1.5),
#                panel.grid.major=element_line(size=0.2),
#                axis.ticks=element_blank(),
#                axis.text.y=element_blank(),
#                plot.title=element_text(size=10),
#                plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"),
#                strip.background=element_blank(),
#                strip.text=element_text(size=10),
#                legend.position="none")
# g <- g + scale_x_continuous(expand = c(0,0))
# g <- g + scale_y_continuous(expand = c(0,0))
# g <- g + coord_cartesian(xlim=c(0,x_range),ylim=c(-0.0,y_range))
# g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
# g <- g + ylab("Probability Density")
# g <- g + geom_vline(xintercept=0,size=2/3)          
# g <- g + geom_ribbon(data=densities,aes(x=x,ymax=y),ymin=0,alpha=0.3)    
# g <- g + geom_line(data=densities,aes(x=x,y=y),size=2/3) 
# ggsave("study_peptides.png", g, height = 2, width=6, limitsize=F,device="png")

print(paste(Sys.time(),"[Finished]"))

