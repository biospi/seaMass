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
nsamp <- nchain * (nitt - burnin) / thin

nproteingroup <- as.integer(dd.params[Key=="internal.nproteingroup", Value])
nrun <- as.integer(dd.params[Key=="internal.nrun", Value])
nlabel <- as.integer(dd.params[Key=="internal.nlabel", Value])

prefix <- ifelse(file.exists("0.0.Rdata"), ".", file.path("..", "..", "norm", "results"))

# load batches
stats.peptidoforms.all <- vector("list", nbatch)
stats.quants.all <- vector("list", nbatch)
mcmc.quants.all <- array(NA, c(nproteingroup, nrun, nlabel, nsamp))
for (i in 1:nbatch) {
  message("[", paste0(Sys.time(), " Processing batch ", i,"/", nbatch, "...]"))

  mcmc.quants.chains <- vector("list", nchain)
  mcmc.peptidoforms.chains <- vector("list", nchain)
  for (j in 1:nchain) {
    load(file.path(prefix, paste0(i-1, ".", j-1, ".Rdata")))
    mcmc.quants.chains[[j]] <- mcmc.quants
    mcmc.peptidoforms.chains[[j]] <- mcmc.peptidoforms
  }
  
  # stats for peptidoform variances
  mcmc.peptidoforms.chains <- mcmc.list(mcmc.peptidoforms.chains)
  stats.peptidoforms.all[[i]] <- summary(mcmc.peptidoforms.chains)$statistics
  stats.peptidoforms.all[[i]] <- data.table(
    Effect = rownames(stats.peptidoforms.all[[i]]),
    Mean = stats.peptidoforms.all[[i]][, "Mean"],
    SD = stats.peptidoforms.all[[i]][, "SD"],
    Rhat = gelman.diag(mcmc.peptidoforms.chains)$psrf[, "Point est."]
  )
  mcmc.peptidoforms.chains <- NULL
  
  # stats for quant ratios
  mcmc.quants.chains <- mcmc.list(mcmc.quants.chains)
  stats.quants.all[[i]] <- summary(mcmc.quants.chains)$statistics
  stats.quants.all[[i]] <- data.table(
    Effect = rownames(stats.quants.all[[i]]),
    Rhat = gelman.diag(mcmc.quants.chains)$psrf[, "Point est."]
  )
  
  # 4D array of ProteinGroup, Run, Label, MCMC samps
  mcmc.quants.chains <- t(as.matrix(mcmc.quants.chains))
  for (j in 1:nrow(mcmc.quants.chains)) {
    p <- as.integer(sub(":Run[0-9]+:Label[0-9]+", "", sub("^ProteinGroup", "", rownames(mcmc.quants.chains)[j])))
    r <- as.integer(sub(":Label[0-9]+", "", sub("^ProteinGroup[0-9]+:Run", "", rownames(mcmc.quants.chains)[j])))
    l <- as.integer(sub("^ProteinGroup[0-9]+:Run[0-9]+:Label", "", rownames(mcmc.quants.chains)[j]))
    print(paste(p, r, l))
    mcmc.quants.all[p, r, l,] <- mcmc.quants.chains[j,]
  }
  mcmc.quants.chains <- NULL
}
stats.peptidoforms.all <- rbindlist(stats.peptidoforms.all)
stats.quants.all <- rbindlist(stats.quants.all)

# calculate medians to use as exposures
mcmc.exposures <- array(NA, c(nrun, nlabel, nsamp))
for (r in 1:nrun) {
  for (l in 1:nlabel) {
    print(paste(r,l))
    for (i in 1:nsamp) {
      mcmc.exposures[r, l, i] <- median(mcmc.quants.all[, r, l, i])
    }
  }
}


plot.mcmc.exposures <- function(mcmc.exposures)
{
  # construct metadata
  mcmc.exposures.meta.func <- function(x) {
    m <- mean(x,na.rm=T)
    if (is.nan(m)) m <- NA
    
    data.table(mean=m, fc=paste0("  ", ifelse(m<0, format(-2^-m,digits=3), format(2^m,digits=3)),"fc"))
  }
  mcmc.exposures.meta <- mcmc.exposures[,as.list(mcmc.exposures.meta.func(Exposure)), by=list(Run,Channel)]
  
  # construct densities
  mcmc.exposures.density.func <- function(x) {
    if (all(x == 0.0)) {
      data.table()
    }
    else {
      dens <- density(x, n=4096, na.rm=T)
      data.table(x=dens$x, y=dens$y)
    }  
  }
  mcmc.exposures.density <- mcmc.exposures[,as.list(mcmc.exposures.density.func(Exposure)), by=list(Run,Channel)]
  
  y_range <- max(mcmc.exposures.density$y)*1.35
  x_range <- max(-min(mcmc.exposures.density$x[mcmc.exposures.density$y>y_range/100]),max(mcmc.exposures.density$x[mcmc.exposures.density$y>y_range/100]))*1.2
  
  g <- ggplot(mcmc.exposures, aes(x=mean))
  g <- g + theme_bw()
  g <- g + theme(panel.border=element_rect(colour="black",size=1),
                 panel.grid.major=element_line(size=0.5),
                 axis.ticks=element_blank(),
                 axis.text.y=element_blank(),
                 plot.title=element_text(size=10),
                 strip.background=element_blank())
  g <- g + scale_x_continuous(expand = c(0,0))
  g <- g + scale_y_continuous(expand = c(0,0))
  g <- g + facet_grid(Channel ~ Run)
  g <- g + coord_cartesian(xlim=c(-x_range, x_range),ylim=c(-0.0,y_range))
  g <- g + xlab(expression('Log'[2]*' Ratio'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept=0,size=1/2,colour="darkgrey")
  g <- g + geom_ribbon(data=mcmc.exposures.density,aes(x=x,ymax=y),ymin=0,size=1/2,alpha=0.3,fill='red')
  g <- g + geom_vline(data=mcmc.exposures.meta,aes(xintercept=mean),size=1/2,colour="red")
  g <- g + geom_line(data=mcmc.exposures.density,aes(x=x,y=y),size=1/2)
  g <- g + geom_text(data=mcmc.exposures.meta,aes(x=mean,label=fc),y=max(mcmc.exposures.density$y)*1.22,hjust=0,vjust=1,size=3,colour="red")
  g  
}

message("[",paste0(Sys.time()," Writing exposures.pdf...]"))
g <- plot.exposures(exposures.plot)
ggsave("exposures.pdf", g, width=4*length(levels(design$Run)), height=1*length(levels(design$Channel)))












exposures <- matrix(nrow=nsamp,ncol=length(paste0(design$Run,design$Channel)))
rownames(exposures) <- 1:nsamp
colnames(exposures) <- paste0(design$Run,design$Channel)

for (rc in colnames(samples.all))
{
  message("[",paste0(Sys.time()," Processing ",rc,"...]"))
  exposures[,rc] <- apply(samples.all[,rc,], 2, function(x) median(x, na.rm=T))
  
  if (all(is.na(exposures[,rc]))) exposures[,rc] <- 0.0
}












mcmc.quants.chains <- mcmc.list(mcmc.quants.chains)
stats.quants <- summary(mcmc.quants.chains)$statistics
stats.quants <- data.table(
  Effect = rownames(stats.quants),
  Mean = stats.quants[, "Mean"],
  SD = stats.quants[, "SD"],
  Rhat = gelman.diag(mcmc.quants.chains)$psrf[, "Point est."]
)







prefix <- "results"
if (length(files) < nbatch * nchain) {
  prefix <- file.path("..","..","norm","results")
  files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
  
  if (length(files) < nbatch * nchain) stop(paste0("ERROR: Missing norm input"))
}







  
nbatch <- as.integer(ifelse("batches" %in% parameters$Key,parameters$Value[parameters$Key=="batches"],10))
nsamp <- as.integer(ifelse("norm_samples" %in% parameters$Key,parameters$Value[parameters$Key=="norm_samples"],1000))
nchain <- as.integer(ifelse("norm_chains" %in% parameters$Key,parameters$Value[parameters$Key=="norm_chains"],10))

prefix <- ifelse(file.exists("index.Rdata"),".",file.path("..","..","input"))
load(file.path(prefix,"index.Rdata"))

prefix <- "results"
files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
if (length(files) < nbatch * nchain) {
  prefix <- file.path("..","..","norm","results")
  files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
  
  if (length(files) < nbatch * nchain) stop(paste0("ERROR: Missing norm input"))
}

# load samples from norm output into samples.all, arranging chains
samples.all <- array(dim=c(nrow(dd.index),length(paste0(design$Run,design$Channel)),col=nsamp))
rownames(samples.all) <- 0:(nrow(dd.index)-1)
colnames(samples.all) <- paste0(design$Run,design$Channel)
dd.index$NormTime <- 0.0

for (f in files) {
  message("[",paste0(Sys.time()," Processing ",f,"...]"))
  
  load(file.path(prefix,f))
  chain <- as.integer(gsub("\\.Rdata$","",gsub("^[0-9]+\\.","",f)))
  begin <- floor(chain/nchain * nsamp) + 1
  
  for (p in names(samples))
  {
    dd.index$NormTime[dd.index$ProteinID==p] <- dd.index$NormTime[dd.index$ProteinID==p] + time[[p]]["elapsed"]
    
    for (rc in colnames(samples[[p]]))
    {
      samps.rc <- rev(samples[[p]][max(1,nrow(samples[[p]])-ceiling(nsamp/nchain)+1):nrow(samples[[p]]),rc])
      samples.all[as.integer(p),rc,begin:(begin+length(samps.rc)-1)] <- samps.rc
    }
  }
}

# calculate and save exposures as medians
exposures <- matrix(nrow=nsamp,ncol=length(paste0(design$Run,design$Channel)))
rownames(exposures) <- 1:nsamp
colnames(exposures) <- paste0(design$Run,design$Channel)

for (rc in colnames(samples.all))
{
  message("[",paste0(Sys.time()," Processing ",rc,"...]"))
  exposures[,rc] <- apply(samples.all[,rc,], 2, function(x) median(x, na.rm=T))
  
  if (all(is.na(exposures[,rc]))) exposures[,rc] <- 0.0
}

save(exposures,dd.index,file="exposures.Rdata")

# plot
exposures.plot <- data.table(t(exposures))
exposures.plot$Run <- design$Run
exposures.plot$Channel <- design$Channel
exposures.plot <- melt(exposures.plot,variable.name="sample",value.name="Exposure",id.vars=c("Run","Channel"))
exposures.plot$Exposure <- exposures.plot$Exposure / log(2)

# contruct mean-centred exposures
exposures.centred.func <- function(x) {
  x$Exposure <- x$Exposure - mean(x$Exposure)
  x
}
exposures.centred <- exposures.plot[,as.list(exposures.centred.func(.SD)), by=list(Run,sample)]



message("[",paste0(Sys.time()," Writing exposures.pdf...]"))
g <- plot.exposures(exposures.plot)
ggsave("exposures.pdf", g, width=4*length(levels(design$Run)), height=1*length(levels(design$Channel)))

message("[",paste0(Sys.time()," Writing exposures_mean-centred.pdf...]"))
g <- plot.exposures(exposures.centred)
ggsave("exposures_mean-centred.pdf", g, width=4*length(levels(design$Run)), height=1*length(levels(design$Channel)))

message("[",paste0(Sys.time()," Finished]"))
