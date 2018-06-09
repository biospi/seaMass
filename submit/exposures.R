invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Starting]"))
  
library(data.table)
library(ggplot2)
  
prefix <- ifelse(file.exists("parameters.Rdata"),".",file.path("..","..","input"))
load(file.path(prefix,"parameters.Rdata"))
load(file.path(prefix,"design.Rdata"))
  
nbatch <- as.integer(ifelse("batches" %in% parameters$Key,parameters$Value[parameters$Key=="batches"],100))
nsamp <- as.integer(ifelse("norm_samples" %in% parameters$Key,parameters$Value[parameters$Key=="norm_samples"],1000))
nchain <- as.integer(ifelse("norm_chains" %in% parameters$Key,parameters$Value[parameters$Key=="norm_chains"],10))

prefix <- ifelse(file.exists("index.Rdata"),".",file.path("..","..","input"))
load(file.path(prefix,"index.Rdata"))

prefix <- "results"
files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
if (length(files) < nbatch * nchain) {
  prefix <- file.path("..","..","norm","results")
  files <- list.files(path=file.path("..","..","norm","results"),pattern="^[0-9]+\\.[0-9]+\\.Rdata")
  
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

# contruct mean-shifted exposures
exposures.shifted.func <- function(x) {
  x$Exposure <- x$Exposure - mean(x$Exposure)
  x
}
exposures.shifted <- exposures.plot[,as.list(exposures.shifted.func(.SD)), by=list(Run,sample)]

plot.exposures <- function(exposures)
{
  # construct metadata
  exposures.meta.func <- function(x) {
    m <- mean(x,na.rm=T)
    if (is.nan(m)) m <- NA
    
    data.table(mean=m, fc=paste0("  ", ifelse(m<0, format(-2^-m,digits=2), format(2^m,digits=3)),"fc"))
  }
  exposures.meta <- exposures[,as.list(exposures.meta.func(Exposure)), by=list(Run,Channel)]
  
  # construct densities
  exposures.density.func <- function(x) {
    if (all(x == 0.0)) {
      data.table()
    }
    else {
      dens <- density(x, n=4096, na.rm=T)
      data.table(x=dens$x, y=dens$y)
    }  
  }
  exposures.density <- exposures[,as.list(exposures.density.func(Exposure)), by=list(Run,Channel)]
  
  y_range <- max(exposures.density$y)*1.35
  x_range <- max(-min(exposures.density$x[exposures.density$y>y_range/100]),max(exposures.density$x[exposures.density$y>y_range/100]))*1.2
  
  g <- ggplot(exposures, aes(x=mean))
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
  g <- g + geom_ribbon(data=exposures.density,aes(x=x,ymax=y),ymin=0,size=1/2,alpha=0.3,fill='red')
  g <- g + geom_vline(data=exposures.meta,aes(xintercept=mean),size=1/2,colour="red")
  g <- g + geom_line(data=exposures.density,aes(x=x,y=y),size=1/2)
  g <- g + geom_text(data=exposures.meta,aes(x=mean,label=fc),y=max(exposures.density$y)*1.22,hjust=0,vjust=1,size=3,colour="red")
  g  
}

message("[",paste0(Sys.time()," Writing exposures.pdf...]"))
g <- plot.exposures(exposures.plot)
ggsave("exposures.pdf", g, width=4*length(levels(design$Run)), height=1*length(levels(design$Channel)))

message("[",paste0(Sys.time()," Writing exposures_mean-shifted.pdf...]"))
g <- plot.exposures(exposures.shifted)
ggsave("exposures_mean-shifted.pdf", g, width=4*length(levels(design$Run)), height=1*length(levels(design$Channel)))

message("[",paste0(Sys.time()," Finished]"))
