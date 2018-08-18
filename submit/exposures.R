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
nL <- length(levels(dd.assays$LabelID))
nA <- length(levels(dd.assays$AssayID))

prefix <- ifelse(file.exists("quants"), ".", file.path("..", "..", "quant", "results", "quants"))

# calculate exposures
mcmc.exposures.label <- matrix(NA, nchain * nsamp, nA)
mcmc.exposures.run <- matrix(NA, nchain * nsamp, nA)
baseline.exposures.label <- rep(NA, nA)
baseline.exposures.run <- rep(NA, nA)
for (j in 1:nchain) {
  for (l in 1:nL) {
    baseline.run <- NA
    
    for (r in dd.assays[, .N, by = RunID][order(-N), RunID]) {
      a <- dd.assays[RunID == r & LabelID == l, AssayID]
      if (length(a) > 0) {
        files <- list.files(prefix, paste0("^[0-9]+\\.", j, "\\.", a, "\\.rds$"))
        
        if (length(files) > 0) {
          if (length(files) < nbatch) stop("ERROR: Some quant output is missing")
          
          message("[", paste0(Sys.time(), " Calculating exposures for chain ", j, "/", nchain, ", run ", r, ", label ", l, "...]"))
          
          # read MCMC samps
          mcmc.quants <- matrix(NA, nsamp, nP)
          baseline.quants <- array(NA, nP)
          for (f in files) {
            mcmc.quants.batch <- readRDS(file.path(prefix, f))
            for (k in 1:ncol(mcmc.quants.batch)) {
              proteinID <- as.integer(sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.quants.batch)[k]))
              baseline.quants[proteinID] <- as.integer(sub("^[0-9]+\\.([0-9]+)$", "\\1", colnames(mcmc.quants.batch)[k]))
              mcmc.quants[, proteinID] <- mcmc.quants.batch[, k]
            }
          }
          
          # calculate intra-run exposure
          if (is.na(baseline.exposures.label[a])) baseline.exposures.label[a] <- names(sort(table(baseline.quants), decreasing = T))[1]
          mcmc.exposures.label[((j-1)*nsamp+1):(j*nsamp), a] <- apply(mcmc.quants[, which(baseline.quants == baseline.exposures.label[a])], 1, function(x) median(x, na.rm = T))
          
          # # calculate inter-run exposure
          # mcmc.quants <- apply(mcmc.quants, 2, '-', mcmc.exposures.label[((j-1)*nsamp+1):(j*nsamp), a])
          # 
          # if (r == dd.assays[, .N, by = RunID][order(-N), RunID][1]) {
          #   if (is.na(baseline.run)) baseline.run <- r
          #   mcmc.quants.baseline <- mcmc.quants
          #   mcmc.exposures.run[((j-1)*nsamp+1):(j*nsamp), a] <- 0.0
          # } else {
          #   baseline.exposures.run[a] <- dd.assays[, .N, by = RunID][order(-N), RunID][1]
          #   mcmc.quants <- mcmc.quants - mcmc.quants.baseline
          #   mcmc.exposures.run[((j-1)*nsamp+1):(j*nsamp), a] <- apply(mcmc.quants[, which(baseline.quants == baseline.exposures.label[a])], 1, function(x) median(x, na.rm = T))
          # }
          # baseline.exposures.run[a] <- r
        }        
      }
    }
  }
}

# fill denominators with zeros
baselines <- as.integer(unique(baseline.exposures.label[!is.na(baseline.exposures.label)]))
mcmc.exposures.label[, baselines] <- 0.0
baseline.exposures.label[baselines] <- baselines
baseline.exposures.label <- factor(baseline.exposures.label, levels = unique(baseline.exposures.label))

# ploting function for exposures
plot.exposures <- function(mcmc.exposures)
{
  dd.exposures <- data.table(t(mcmc.exposures))
  dd.exposures$Assay <- dd.assays$Assay
  dd.exposures <- melt(dd.exposures, variable.name="mcmc", value.name="Exposure", id.vars = "Assay")
  dd.exposures <- dd.exposures[complete.cases(dd.exposures),]
  
  # construct metadata
  dd.exposures.meta.func <- function(x) {
    m <- mean(x, na.rm=T)
    if (is.nan(m)) m <- NA

    data.table(mean = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
  }
  dd.exposures.meta <- dd.exposures[, as.list(dd.exposures.meta.func(Exposure)), by = Assay]

  # construct densities
  dd.exposures.density.func <- function(x) {
    if (all(x == 0.0)) {
      data.table()
    }
    else {
      dens <- density(x, n = 4096, na.rm = T)
      data.table(x = dens$x, y = dens$y)
    }
  }
  dd.exposures.density <- dd.exposures[, as.list(dd.exposures.density.func(Exposure)), by = Assay]

  y_range <- max(dd.exposures.density$y) * 1.35
  x_range <- max(-min(dd.exposures.density$x[dd.exposures.density$y > y_range/100]), max(dd.exposures.density$x[dd.exposures.density$y > y_range/100])) * 1.2

  g <- ggplot(dd.exposures, aes(x = mean))
  g <- g + theme_bw()
  g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
                 panel.grid.major = element_line(size = 0.5),
                 axis.ticks = element_blank(),
                 axis.text.y = element_blank(),
                 plot.title = element_text(size = 10),
                 strip.background=element_blank())
  g <- g + scale_x_continuous(expand = c(0, 0))
  g <- g + scale_y_continuous(expand = c(0, 0))
  g <- g + facet_grid(Assay ~ .)
  g <- g + coord_cartesian(xlim = c(-x_range, x_range), ylim = c(-0.0, y_range))
  g <- g + xlab(expression('Log'[2]*' Ratio'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
  g <- g + geom_ribbon(data = dd.exposures.density,aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3, fill = "red")
  g <- g + geom_vline(data = dd.exposures.meta,aes(xintercept = mean), size = 1/2, colour = "red")
  g <- g + geom_line(data = dd.exposures.density, aes(x = x,y = y), size = 1/2)
  g <- g + geom_text(data = dd.exposures.meta, aes(x = mean, label = fc), y = max(dd.exposures.density$y) * 1.22, hjust = 0, vjust = 1, size = 3, colour="red")
  g
}

# plot label exposures
message("[", paste0(Sys.time(), " Writing exposures_intrarun-alr.pdf...]"))
g <- plot.exposures(mcmc.exposures.label)
ggsave("medians.pdf", g, width = 8, height = 0.5 * nA)
saveRDS(mcmc.exposures.label, "medians.rds")

# # plot run exposures
# message("[",paste0(Sys.time()," Writing exposures-run.pdf...]"))
# g <- plot.exposures(mcmc.exposures.run)
# ggsave("exposures-run.pdf", g, width = 8, height = 0.5 * nA)
# 
# test <- mcmc.exposures.run - mcmc.exposures.label
# plot.exposures(mcmc.exposures.label)

# # contruct label mean-centred exposures
# for (i in levels(baseline.exposures.label)) {
#  cols <- which(baseline.exposures.label == i)
#  mcmc.exposures.label[, cols] <- mcmc.exposures.label[, cols] - rowMeans(mcmc.exposures.label[, cols])
#}
 
# # plot label mean-centred exposures
# message("[",paste0(Sys.time()," Writing exposures_intrarun-clr.pdf...]"))
# g <- plot.exposures(mcmc.exposures.label)
# ggsave("exposures_intrarun-clr.pdf", g, width = 8, height = 0.5 * nA)
# saveRDS(mcmc.exposures.label, "exposures_intrarun-clr.rds")

# calculate normalised quants within each run
dd.quants <- vector("list", nbatch * nA)
for (i in 1:nbatch) {
  
  mcmc.quants <- array(NA, c(nchain * nsamp, nrow(dd.proteins[batchID == i,]), nA))
  colnames(mcmc.quants) <- dd.proteins[batchID == i, ProteinID]
  baseline.quants <- matrix(NA, nrow(dd.proteins[batchID == i,]), nA)
  rownames(baseline.quants) <- dd.proteins[batchID == i, ProteinID]
  for (a in 1:nA) {
    files <- list.files(prefix, paste0("^", i, "\\.[0-9]+\\.", a, "\\.rds"))
    if (length(files) > 0) {
      if (length(files) < nchain) stop("ERROR: Some quant output is missing")
      
      message("[", paste0(Sys.time(), " Determining quants for batch ", i, "/", nbatch, ", assay ", a, "...]"))
      
      # read MCMC samps and normalise
      for (j in 1:nchain) {
        mcmc.quants.chain <- readRDS(file.path(prefix, files[j]))
        for (k in 1:ncol(mcmc.quants.chain)) {
          proteinID <-  sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.quants.chain)[k])
          baseline.quants[proteinID, a] <- as.integer(sub("^[0-9]+\\.([0-9]+)$", "\\1", colnames(mcmc.quants.chain)[k]))
          if (is.na(baseline.quants[proteinID, baseline.quants[proteinID, a]])) {
            baseline.quants[proteinID, baseline.quants[proteinID, a]] <- baseline.quants[proteinID, a]
            mcmc.quants[, proteinID, baseline.quants[proteinID, a]] <- 0.0
          }
          mcmc.quants[, proteinID, a] <- mcmc.quants.chain[, k] - mcmc.exposures.label[, a]
        }
      }
    }
  }
  
  # clr & ilr transform
   for (l in levels(baseline.exposures.label)) {
     cols <- which(baseline.exposures.label == l)
     for (p in 1:ncol(mcmc.quants)) {
       mcmc.quants[, p, cols] <- mcmc.quants[, p, cols] - rowMeans(mcmc.quants[, p, cols])
       #mcmc.quants[, p, cols] <- cbind(clr2ilr(mcmc.quants[, p, cols], ilrBase(D = nL, method = "balanced")), 0)
     }
   }
  
  for (a in 1:nA) {
    dd.quants[[(i-1)*nA + a]] <- data.table(
      ProteinID = colnames(mcmc.quants),
      AssayID = a,
      log2mean = colMeans(mcmc.quants[,, a]),
      log2sd = apply(mcmc.quants[,, a], 2, sd)
    )  
  }
}
dd.quants <- rbindlist(dd.quants)
dd.quants$ProteinID <- factor(dd.quants$ProteinID, levels = levels(dd.proteins$ProteinID))
dd.quants$AssayID <- factor(dd.quants$AssayID, levels = levels(dd.assays$AssayID))
dd.quants <- merge(dd.assays[, list(AssayID, Assay)], dd.quants, by = "AssayID")
dd.quants <- dd.quants[order(ProteinID, AssayID),]

# write out quants
dd.quants.out <- dd.quants
dd.quants.out$Assay.sd <- factor(paste0(dd.quants.out$Assay, ".log2sd"))
dd.quants.out$Assay <- factor(paste0(dd.quants.out$Assay, ".log2mean"))
dd.quants.out <- merge(
  dcast(dd.quants.out, ProteinID ~ Assay, value.var = "log2mean"),
  dcast(dd.quants.out, ProteinID ~ Assay.sd, value.var = "log2sd"),
  by = "ProteinID"
)
#dd.quants.out <- merge(dd.proteins[, list(ProteinID, Protein, nPeptide, nFeature, nCount)], dd.quants.out, by = "ProteinID")
dd.quants.out <- merge(dd.proteins, dd.quants.out, by = "ProteinID")
fwrite(dd.quants.out, paste0(dd.params[Key == "bayesprot.id", Value], "_quants.csv") )

# plot quant distributions
dd.quants.meta <- dd.quants[, list(
  log2median = median(log2mean, na.rm = T),
  log2lower = quantile(log2mean, probs = 0.025, na.rm = T),
  log2upper = quantile(log2mean, probs = 0.975, na.rm = T)
), by = AssayID]

dd.quants.plot <- merge(dd.quants, dd.quants.meta)
dd.quants.plot <- dd.quants.plot[log2mean >= log2lower & log2mean <= log2upper,]

g <- ggplot(dd.quants.plot, aes(x = AssayID, y = log2mean))
g <- g + scale_y_continuous(expand = c(0, 0))
#g <- g + coord_cartesian(ylim = c(-0.03, 0.03))
g <- g + geom_violin()
g <- g + geom_segment(data = dd.quants.meta, aes(x = as.integer(AssayID) - 0.45, xend = as.integer(AssayID) + 0.45, y = log2median, yend = log2median),size = 1/2)
g





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
#   prefix <- file.path("..","..","quant","results")
#   files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
#   
#   if (length(files) < nbatch * nchain) stop(paste0("ERROR: Missing quant input"))
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
# nsamp <- as.integer(ifelse("quant_samples" %in% parameters$Key,parameters$Value[parameters$Key=="quant_samples"],1000))
# nchain <- as.integer(ifelse("quant_chains" %in% parameters$Key,parameters$Value[parameters$Key=="quant_chains"],10))
# 
# prefix <- ifelse(file.exists("index.Rdata"),".",file.path("..","..","input"))
# load(file.path(prefix,"index.Rdata"))
# 
# prefix <- "results"
# files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
# if (length(files) < nbatch * nchain) {
#   prefix <- file.path("..","..","quant","results")
#   files <- list.files(path=prefix,pattern="^[0-9]+\\.[0-9]+\\.Rdata")
#   
#   if (length(files) < nbatch * nchain) stop(paste0("ERROR: Missing quant input"))
# }
# 
# # load samples from quant output into samples.all, arranging chains
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

message("[",paste0(Sys.time()," Finished]"))
