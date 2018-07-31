invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Starting]"))
  
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(coda))

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

nbatch <- as.integer(dd.params[Key=="nbatch", Value])
nchain <- as.integer(dd.params[Key=="quant.nchain", Value])
nitt <- as.integer(dd.params[Key=="quant.nitt", Value])
burnin <- as.integer(dd.params[Key=="quant.burnin", Value])
thin <- as.integer(dd.params[Key=="quant.thin", Value])
nsamp <- (nitt - burnin) / thin

nP <- length(levels(dd.proteins$ProteinID))
#nR <- length(levels(dd.assays$RunID))
#nL <- length(levels(dd.assays$LabelID))

prefix <- ifelse(file.exists("quants"), ".", file.path("..", "..", "norm", "results", "quants"))

# calculate exposures
mcmc.exposures <- matrix(0.0, nchain * nsamp, nrow(dd.assays))
for (j in 1:nchain) {
  for(k in 1:nrow(dd.assays)) {
    message("[", paste0(Sys.time(), " Calculating exposures for chain ", j, "/", nchain, ", run ", dd.assays$RunID[k], ", label ", dd.assays$LabelID[k], "...]"))
    
    if(dd.assays$LabelID[k] != levels(dd.assays$LabelID)[1]) {
      
      # read MCMC samps
      mcmc.quants <- vector("list", nbatch)
      for (i in 1:nbatch)  mcmc.quants[[i]] <- readRDS(file.path(prefix, paste0(i, ".", j, ".", dd.assays$RunID[k], ".", dd.assays$LabelID[k], ".rds")))
      mcmc.quants <- rbindlist(mcmc.quants)
      
      # calculate exposure
      mcmc.exposures[((j-1)*nsamp+1):(j*nsamp), k] <- mcmc.quants[, median(log2quant), by = sampID]$V1
    }
  }
}
saveRDS(mcmc.exposures, "exposures.rds")
    
# calculate normalised quants within each run
stats.quants <- vector("list", nbatch)
for (i in 1:nbatch) {
  for(r in levels(dd.assays$RunID)) {
    message("[", paste0(Sys.time(), " Calculating  normalised quants for batch ", i, "/", nbatch, ", run ", r, "...]"))
    
    labels <- dd.assays[RunID == r,]$LabelID
    labels <- labels[labels != levels(dd.assays$LabelID)[1]]
    
    mcmc.quants <- NULL
    for(l in labels) {
      for(j in 1:nchain) {
        dd <- readRDS(file.path(prefix, paste0(i, ".", j, ".", r, ".", l, ".rds"))) 
        
        if (is.null(mcmc.quants)) mcmc.quants <- array(0.0, c(nchain * nsamp, length(labels), length(levels(dd$ProteinID))))
 
        mcmc.quants[((i-1)*nsamp+1):(i*nsamp), j] <- mcmc.quants[, median(log2quant), by = sampID]$V1
        
        mcmc.quants[[k*nchain + l]] <- 
        mcmc.quants[[k*nchain + l]][, chainID := l]
        mcmc.quants[[k*nchain + l]][, RunID := j]
        mcmc.quants[[k*nchain + l]][, LabelID := k]
      }  
    }
    mcmc.quants <- merge(rbindlist(mcmc.quants), mcmc.exposures)
    
    mcmc.quants[, log2quant := log2quant - log2exposure]
    mcmc.quants[, log2exposure := NULL]
    mcmc.quants <- rbind(mcmc.quants, data.table(
      ProteinID = rep(levels(mcmc.quants$ProteinID), each = nsamp * nchain),
      sampID = rep(1:nsamp, nchain * length(levels(mcmc.quants$ProteinID))),
      chainID = rep(rep(1:nchain, each = nsamp), length(levels(mcmc.quants$ProteinID))),
      log2quant = 0.0,
      AssayID = 1
    ))   
    
    mean.centre <- function(x) x - mean(x)
    mcmc.quants[, log2quant := mean.centre(log2quant), by = c("ProteinID", "sampID", "chainID")]
    stats.quants[[i]] <- mcmc.quants[, list(log2mean = mean(log2quant), log2stdev = sd(log2quant)), by = c("ProteinID", "AssayID")]
  }
}
stats.quants <- rbindlist(stats.quants)
stats.quants$AssayID <- factor(stats.quants$AssayID, levels = levels(dd.assays$AssayID))
stats.quants <- merge(dd.assays[, list(AssayID, Assay)], stats.quants, by = "AssayID")
stats.quants$Assay.stdev <- factor(paste0(stats.quants$Assay, ".log2stdev"))
stats.quants$Assay <- factor(paste0(stats.quants$Assay, ".log2mean"))
stats.quants <- merge(
  dcast(stats.quants, ProteinID ~ Assay, value.var = "log2mean"),
  dcast(stats.quants, ProteinID ~ Assay.stdev, value.var = "log2stdev"),
  by = "ProteinID"
)
stats.quants <- merge(dd.proteins, stats.quants, by = "ProteinID")

fwrite(stats.quants, paste0(dd.params[Key == "bayesprot.id", Value], "_quants.csv") )

# stats.peptidoforms.all <- vector("list", nbatch)
# stats.quants.all <- vector("list", nbatch)
# mcmc.quants.all <- array(NA, c(nP, nR, nL, nsamp))
# time.mcmc.all <- rep(0, nbatch)
# for (i in 1:nbatch) {
#   message("[", paste0(Sys.time(), " Processing batch ", i,"/", nbatch, "...]"))
# 
#   mcmc.quants.chains <- vector("list", nchain)
#   mcmc.peptidoforms.chains <- vector("list", nchain)
#   for (j in 1:nchain) {
#     load(file.path(prefix, paste0(i-1, ".", j-1, ".Rdata")))
#     mcmc.quants.chains[[j]] <- as.mcmc(mcmc.quants)
#     mcmc.peptidoforms.chains[[j]] <- as.mcmc(mcmc.peptidoforms)
#     time.mcmc.all[i] <- time.mcmc.all[i] + time.mcmc["elapsed"]
#   }
#   
#   # stats for peptidoform variances
#   mcmc.peptidoforms.chains <- mcmc.list(mcmc.peptidoforms.chains)
#   stats.peptidoforms.all[[i]] <- summary(mcmc.peptidoforms.chains)$statistics
#   stats.peptidoforms.all[[i]] <- data.table(
#     Effect = rownames(stats.peptidoforms.all[[i]]),
#     Mean = stats.peptidoforms.all[[i]][, "Mean"],
#     SD = stats.peptidoforms.all[[i]][, "SD"],
#     Rhat = gelman.diag(mcmc.peptidoforms.chains)$psrf[, "Point est."]
#   )
#   mcmc.peptidoforms.chains <- NULL
#   
#   # stats for quant ratios
#   mcmc.quants.chains <- mcmc.list(mcmc.quants.chains)
#   stats.quants.all[[i]] <- summary(mcmc.quants.chains)$statistics
#   rhat.cols <- stats.quants.all[[i]][, "SD"] != 0
#   stats.quants.all[[i]] <- data.table(
#     Effect = rownames(stats.quants.all[[i]]),
#     Rhat = NA
#   )
#   stats.quants.all[[i]]$Rhat[rhat.cols] <- gelman.diag(mcmc.quants.chains[, rhat.cols])$psrf[, "Point est."]
#   
#   # 4D array of Protein, Run, Label, MCMC samps
#   mcmc.quants.chains <- t(as.matrix(mcmc.quants.chains))
#   for (j in 1:nrow(mcmc.quants.chains)) {
#     p <- as.integer(gsub("^([0-9]+):[0-9]+:[0-9]+$", "\\1", rownames(mcmc.quants.chains)[j]))
#     r <- as.integer(gsub("^[0-9]+:([0-9])+:[0-9]+$", "\\1", rownames(mcmc.quants.chains)[j]))
#     l <- as.integer(gsub("^[0-9]+:[0-9]+:([0-9])+$", "\\1", rownames(mcmc.quants.chains)[j]))
#     mcmc.quants.all[p, r, l,] <- mcmc.quants.chains[j,]
#   }
#   mcmc.quants.chains <- NULL
# }
# stats.peptidoforms.all <- rbindlist(stats.peptidoforms.all)
# stats.quants.all <- rbindlist(stats.quants.all)
# 
# # calculate medians to use as exposures
# mcmc.exposures <- array(NA, c(nR, nL, nsamp))
# for (r in 1:nR) {
#   for (l in 1:nL) {
#     print(paste(r,l))
#     for (i in 1:nsamp) {
#       mcmc.exposures[r, l, i] <- median(log2(exp(mcmc.quants.all[, r, l, i])))
#     }
#   }
# }
# 
# # save times for load balancing
# save(time.mcmc.all, file = "time.mcmc.all.Rdata")
# 
# 
# 
# 
# plot.mcmc.exposures <- function(mcmc.exposures)
# {
#   # construct metadata
#   mcmc.exposures.meta.func <- function(x) {
#     m <- mean(x,na.rm=T)
#     if (is.nan(m)) m <- NA
#     
#     data.table(mean=m, fc=paste0("  ", ifelse(m<0, format(-2^-m,digits=3), format(2^m,digits=3)),"fc"))
#   }
#   mcmc.exposures.meta <- mcmc.exposures[,as.list(mcmc.exposures.meta.func(Exposure)), by=list(Run,Channel)]
#   
#   # construct densities
#   mcmc.exposures.density.func <- function(x) {
#     if (all(x == 0.0)) {
#       data.table()
#     }
#     else {
#       dens <- density(x, n=4096, na.rm=T)
#       data.table(x=dens$x, y=dens$y)
#     }  
#   }
#   mcmc.exposures.density <- mcmc.exposures[,as.list(mcmc.exposures.density.func(Exposure)), by=list(Run,Channel)]
#   
#   y_range <- max(mcmc.exposures.density$y)*1.35
#   x_range <- max(-min(mcmc.exposures.density$x[mcmc.exposures.density$y>y_range/100]),max(mcmc.exposures.density$x[mcmc.exposures.density$y>y_range/100]))*1.2
#   
#   g <- ggplot(mcmc.exposures, aes(x=mean))
#   g <- g + theme_bw()
#   g <- g + theme(panel.border=element_rect(colour="black",size=1),
#                  panel.grid.major=element_line(size=0.5),
#                  axis.ticks=element_blank(),
#                  axis.text.y=element_blank(),
#                  plot.title=element_text(size=10),
#                  strip.background=element_blank())
#   g <- g + scale_x_continuous(expand = c(0,0))
#   g <- g + scale_y_continuous(expand = c(0,0))
#   g <- g + facet_grid(Channel ~ Run)
#   g <- g + coord_cartesian(xlim=c(-x_range, x_range),ylim=c(-0.0,y_range))
#   g <- g + xlab(expression('Log'[2]*' Ratio'))
#   g <- g + ylab("Probability Density")
#   g <- g + geom_vline(xintercept=0,size=1/2,colour="darkgrey")
#   g <- g + geom_ribbon(data=mcmc.exposures.density,aes(x=x,ymax=y),ymin=0,size=1/2,alpha=0.3,fill='red')
#   g <- g + geom_vline(data=mcmc.exposures.meta,aes(xintercept=mean),size=1/2,colour="red")
#   g <- g + geom_line(data=mcmc.exposures.density,aes(x=x,y=y),size=1/2)
#   g <- g + geom_text(data=mcmc.exposures.meta,aes(x=mean,label=fc),y=max(mcmc.exposures.density$y)*1.22,hjust=0,vjust=1,size=3,colour="red")
#   g  
# }
# 
# message("[",paste0(Sys.time()," Writing exposures.pdf...]"))
# g <- plot.exposures(exposures.plot)
# ggsave("exposures.pdf", g, width=4*length(levels(design$Run)), height=1*length(levels(design$Channel)))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# exposures <- matrix(nrow=nsamp,ncol=length(paste0(design$Run,design$Channel)))
# rownames(exposures) <- 1:nsamp
# colnames(exposures) <- paste0(design$Run,design$Channel)
# 
# for (rc in colnames(samples.all))
# {
#   message("[",paste0(Sys.time()," Processing ",rc,"...]"))
#   exposures[,rc] <- apply(samples.all[,rc,], 2, function(x) median(x, na.rm=T))
#   
#   if (all(is.na(exposures[,rc]))) exposures[,rc] <- 0.0
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mcmc.quants.chains <- mcmc.list(mcmc.quants.chains)
# stats.quants <- summary(mcmc.quants.chains)$statistics
# stats.quants <- data.table(
#   Effect = rownames(stats.quants),
#   Mean = stats.quants[, "Mean"],
#   SD = stats.quants[, "SD"],
#   Rhat = gelman.diag(mcmc.quants.chains)$psrf[, "Point est."]
# )
# 
# 
# 
# 
# 
# 
# 
# prefix <- "results"
# if (length(files) < nbatch * nchain) {
#   prefix <- file.path("..","..","norm","results")
#   files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
#   
#   if (length(files) < nbatch * nchain) stop(paste0("ERROR: Missing norm input"))
# }
# 
# 
# 
# 
# 
# 
# 
#   
# nbatch <- as.integer(ifelse("batches" %in% parameters$Key,parameters$Value[parameters$Key=="batches"],10))
# nsamp <- as.integer(ifelse("norm_samples" %in% parameters$Key,parameters$Value[parameters$Key=="norm_samples"],1000))
# nchain <- as.integer(ifelse("norm_chains" %in% parameters$Key,parameters$Value[parameters$Key=="norm_chains"],10))
# 
# prefix <- ifelse(file.exists("index.Rdata"),".",file.path("..","..","input"))
# load(file.path(prefix,"index.Rdata"))
# 
# prefix <- "results"
# files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
# if (length(files) < nbatch * nchain) {
#   prefix <- file.path("..","..","norm","results")
#   files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
#   
#   if (length(files) < nbatch * nchain) stop(paste0("ERROR: Missing norm input"))
# }
# 
# # load samples from norm output into samples.all, arranging chains
# samples.all <- array(dim=c(nrow(dd.index),length(paste0(design$Run,design$Channel)),col=nsamp))
# rownames(samples.all) <- 0:(nrow(dd.index)-1)
# colnames(samples.all) <- paste0(design$Run,design$Channel)
# dd.index$NormTime <- 0.0
# 
# for (f in files) {
#   message("[",paste0(Sys.time()," Processing ",f,"...]"))
#   
#   load(file.path(prefix,f))
#   chain <- as.integer(gsub("\\.Rdata$","",gsub("^[0-9]+\\.","",f)))
#   begin <- floor(chain/nchain * nsamp) + 1
#   
#   for (p in names(samples))
#   {
#     dd.index$NormTime[dd.index$ProteinID==p] <- dd.index$NormTime[dd.index$ProteinID==p] + time[[p]]["elapsed"]
#     
#     for (rc in colnames(samples[[p]]))
#     {
#       samps.rc <- rev(samples[[p]][max(1,nrow(samples[[p]])-ceiling(nsamp/nchain)+1):nrow(samples[[p]]),rc])
#       samples.all[as.integer(p),rc,begin:(begin+length(samps.rc)-1)] <- samps.rc
#     }
#   }
# }
# 
# # calculate and save exposures as medians
# exposures <- matrix(nrow=nsamp,ncol=length(paste0(design$Run,design$Channel)))
# rownames(exposures) <- 1:nsamp
# colnames(exposures) <- paste0(design$Run,design$Channel)
# 
# for (rc in colnames(samples.all))
# {
#   message("[",paste0(Sys.time()," Processing ",rc,"...]"))
#   exposures[,rc] <- apply(samples.all[,rc,], 2, function(x) median(x, na.rm=T))
#   
#   if (all(is.na(exposures[,rc]))) exposures[,rc] <- 0.0
# }
# 
# save(exposures,dd.index,file="exposures.Rdata")
# 
# # plot
# exposures.plot <- data.table(t(exposures))
# exposures.plot$Run <- design$Run
# exposures.plot$Channel <- design$Channel
# exposures.plot <- melt(exposures.plot,variable.name="sample",value.name="Exposure",id.vars=c("Run","Channel"))
# exposures.plot$Exposure <- exposures.plot$Exposure / log(2)
# 
# # contruct mean-centred exposures
# exposures.centred.func <- function(x) {
#   x$Exposure <- x$Exposure - mean(x$Exposure)
#   x
# }
# exposures.centred <- exposures.plot[,as.list(exposures.centred.func(.SD)), by=list(Run,sample)]
# 
# 
# 
# message("[",paste0(Sys.time()," Writing exposures.pdf...]"))
# g <- plot.exposures(exposures.plot)
# ggsave("exposures.pdf", g, width=4*length(levels(design$Run)), height=1*length(levels(design$Channel)))
# 
# message("[",paste0(Sys.time()," Writing exposures_mean-centred.pdf...]"))
# g <- plot.exposures(exposures.centred)
# ggsave("exposures_mean-centred.pdf", g, width=4*length(levels(design$Run)), height=1*length(levels(design$Channel)))
# 
# message("[",paste0(Sys.time()," Finished]"))
