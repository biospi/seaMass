plot.samples <- function(s.Sol, meta, design, filename) {  
  library(plyr)
  library(ggplot2)
  
  samps.conditions <- s.Sol[,colnames(s.Sol) %in% paste0('Condition', levels(design$Condition)),drop=F]
  samps.baseline <- data.frame(x = rnorm(nrow(samps.conditions),1e-6,1e-6))
  colnames(samps.baseline)[1] <- paste0('Condition',levels(design$Condition)[1])
  samps.conditions <- cbind(samps.baseline, samps.conditions)
  colnames(samps.conditions) <- sub('Condition', '', colnames(samps.conditions), fixed=T)    
  
  if (length(levels(design$Population))==1) {
    samples <- data.frame(Name=paste0('Population', levels(design$Sample)),Sample=levels(design$Sample))
  } else {
    samples <- mdply(levels(design$Sample), function(i) data.frame(Name=paste0('Population', design$Population[design$Sample==i][1], '.Sample.', i),Sample=i))
  }
  samps.samples <- s.Sol[seq(1,nrow(s.Sol),1),colnames(s.Sol) %in% samples$Name]
  colnames(samps.samples) <- samples$Sample[match(colnames(samps.samples), samples$Name)]
  samps.conditions <- samps.conditions[seq(1,nrow(s.Sol),1),]
  
  # add Condition fixed effects to samples
  samps.samples_plus_conditions <- samps.samples
  for (i in colnames(samps.samples)) {
    if (as.character(design$Condition[design$Sample==i][1]) %in% colnames(samps.conditions))
    {
      samps.samples_plus_conditions[,i] <- samps.samples_plus_conditions[,i] + samps.conditions[,as.character(design$Condition[design$Sample==i][1])]
    }
  }  
  
  samps.samples_plus_conditions.melted <- mdply(colnames(samps.samples_plus_conditions), function(i) {
    data.frame(Sample = i, Condition = design$Condition[design$Sample==i][1], value = samps.samples_plus_conditions[,i])
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
  samps.conditions.melted <- merge(samps.conditions.melted, data.frame(Sample=design$Sample, Condition=design$Condition))
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
  
  for (i in colnames(samps.conditions)) {
    if (i == tolower(i))
    {
      samps.conditions[,i] <- mean(samps.conditions[,i]) + rnorm(nrow(samps.conditions),1e-6,1e-6)
    }
  }   
  test_conditions <- levels(design$Condition)[levels(design$Condition) != tolower(levels(design$Condition))]
  test_conditions <- test_conditions[2:length(test_conditions)]
  samps$value[!(samps$Condition %in% test_conditions) & samps$X1=="1"] <- NA
  
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
  g <- g + ggtitle(meta)
  g <- g + coord_cartesian(ylim=ylim)
  g <- g + ylab(expression('Log'[2]*' Ratio'))
  g <- g + geom_boxplot(data = samps, aes(y = value), alpha = 0.0, weight = 0, colour = "white", size = 0, outlier.size = 0)
  g <- g + geom_hline(yintercept=0,size=1/2,colour="darkgrey")          
  g <- g + geom_segment(data = stats.samples_plus_conditions, aes(x = as.integer(Sample)-0.45, xend = as.integer(Sample) + 0.45, y = mean, yend = mean),size = 1/2,colour="darkgrey")
  g <- g + geom_violin(data = samps, aes(y = value, colour = X1, alpha = X1, fill = Condition), position="identity", trim=T, size = 1/2)
  g <- g + geom_segment(data = stats.conditions, aes(x = as.integer(Sample)-0.5, xend = as.integer(Sample) + 0.5, y = mean, yend = mean),size = 1/2)
  g <- g + scale_alpha_manual(values=c(0.0,1.0))
  g <- g + scale_colour_manual(values=c("darkgrey","black"))
  ggsave(filename, g, height=2, width=6, limitsize=F)
  
  ylim
}


plots <- function(protein_id,design,nburnin,nsamp,nchain,thin,fc,tol) { 
  print(paste(Sys.time(),"[Starting]"))
  
  library(coda)
  library(mcgibbsit)
  library(plyr)
  library(reshape2)

  test_conditions <- levels(design$Condition)[levels(design$Condition) != tolower(levels(design$Condition))]
  test_conditions <- test_conditions[2:length(test_conditions)]
  
  print(paste0(Sys.time()," [output() Processing protein ",protein_id,"]"))    
    
  files = list.files(path=paste0(protein_id),pattern=paste0("^[0-9]+\\.Rdata"))
  dics <- rep(NA,length(files))
  samps.Sol <- mcmc.list(mlply(files, function(f) {
    load(paste0(protein_id,"/",f))
    samps.Sol
  }))
  end <- summary(samps.Sol[[1]])$end
  s.Sol <- as.matrix(window(samps.Sol,end-nsamp/nchain,end))
  samps.VCV <- mcmc.list(mlply(files, function(f) {
    load(paste0(protein_id,"/",f))
    samps.VCV
  }))
  s.VCV <- as.matrix(window(samps.VCV,end-nsamp/nchain,end))
  dics <- mdply(files, function(f) {
    load(paste0(protein_id,"/",f))
    dic
  })
    
  # one-sided statistical tests, checking precision
  stats <- mdply(test_conditions, function(con) {
    samps <- samps.Sol[,paste0('Condition', con)]
    s <- s.Sol[,paste0('Condition', con)]
    s.mean <- mean(s)
    s.hpdi <- HPDinterval(mcmc(s))
    
    qs <- c(sum(s > log2(fc)), sum(s < -log2(fc)), sum(s >= -log2(fc) & s <= log2(fc)))  
    qs <- pmin(pmax(qs,1),length(s)-1) / length(s)
    
    b <- mdply(qs, function(q) {
      res <- mcgibbsit(samps,q,tol)
      
      nburnin_pred <- NA
      try(nburnin_pred <- ceiling(res$resmatrix[,'M']/thin(samps)), silent=T)
      nsamp_pred <- NA
      try(nsamp_pred <- ceiling(res$resmatrix[,'N']/thin(samps)), silent=T)
      
      data.frame(burnin_pred=nburnin_pred,
                 burnin_ok=ifelse(nburnin>=nburnin_pred,"Yes","No"),
                 samp_pred=nsamp_pred,
                 samp_ok=ifelse(nsamp>=nsamp_pred,"Yes","No"),
                 dic=mean(dics$V1),
                 mean=s.mean,
                 lower=s.hpdi[1],
                 upper=s.hpdi[2],
                 localFDR=1-q)
    })
    b$X1 <- c("Up","Down","Same")
    colnames(b)[1] <- "Test"
    b
  })  
  stats$X1 <- test_conditions[stats$X1]
  colnames(stats)[1] <- "Condition"
  
  # save stats
  dic <- mean(dics$V1)
  save(stats, dic, file=paste0("stats/",protein_id,".Rdata"))   
  
  # plots
  plot.samples(s.Sol, protein_id, design, paste0("samples/",protein_id,".pdf"))
}


# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  print(paste(Sys.time(),"[Starting]"))

  # some tuning parameters (should come from parameters.Rdata with defaults given here)
  load("parameters.Rdata")  
  nburnin <- as.integer(ifelse("model_nburnin" %in% parameters$Key,parameters$Value[parameters$Key=="model_nburnin"],1000000))
  nsamp <- as.integer(ifelse("model_nsamp" %in% parameters$Key,parameters$Value[parameters$Key=="model_nsamp"],1000000))
  nchain <- as.integer(ifelse("model_nchain" %in% parameters$Key,parameters$Value[parameters$Key=="model_nchain"],100))
  thin <- as.integer(ifelse("model_thin" %in% parameters$Key,parameters$Value[parameters$Key=="model_thin"],10))
  fc <- as.double(ifelse("model_fc" %in% parameters$Key,parameters$Value[parameters$Key=="model_fc"],1.05))
  tol <- as.double(ifelse("model_tol" %in% parameters$Key,parameters$Value[parameters$Key=="model_tol"],0.025))

  dir.create("stats")
  dir.create("samples")
  
  # run jobs
  protein_ids <- commandArgs(T)[3:length(commandArgs(T))]
  
  load("design.Rdata")  
  devnull <- sapply(protein_ids, function(protein_id) {
    print(paste(Sys.time(),paste0("[Processing job ",protein_id,"]")))
    plots(protein_id,design,nburnin,nsamp,nchain,thin,fc,tol)
  })
  
  print(paste(Sys.time(),"[Finished]"))  
}
