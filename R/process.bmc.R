#' process.bmc (BayesProt internal function)
#'
#' @return .
#' @import data.table
#' @export

process.bmc <- function(chain) {
  message(paste0("[", Sys.time(), "] BMC started, chain=", chain))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
  nsamp <- (params$quant.nitt - params$quant.burnin) / params$quant.thin

  nP <- length(levels(dd.proteins$ProteinID))
  nA <- length(levels(dd.assays$AssayID))

  prefix <- ifelse(file.exists("1.1.rds"), ".", file.path("..", "..", "quant", "results", "quants"))

  # process design
  cts <- combn(levels(dd.assays$Condition), 2)

  # read MCMC samps
  mcmc.quants.all <- array(NA, c(nsamp, nP, nA))
  for (a in 1:nA) {
    mcmc.quants <- readRDS(file.path(prefix, paste0(a, ".", chain, ".rds")))

    # todo: check baseline
    mcmc.quants.all[, as.integer(sub("\\.[0-9]+$", "", colnames(mcmc.quants))), a] <- mcmc.quants
  }

  # parallel processing - modelComparisonBatch
  message(paste0("[", Sys.time(), "]  modelComparisonBatch..."))
  doParallel::registerDoParallel(params$nthread)
  `%dopar%` <- foreach::`%dopar%`
  dd.output <- foreach::foreach(ct = rep(1:ncol(cts), each = nsamp), s = rep(1:nsamp, ncol(cts)), .combine = rbind, .options.multicore = list(preschedule = F, silent = T)) %dopar% {
    dd <- as.data.table(mcmc.quants.all[s,,])
    colnames(dd) <- dd.assays[, Assay]
    dd.bmc <- as.data.table(bayesmodelquant::modelComparisonBatch(dd, list(dd.assays[Condition == cts[1, ct], Assay], dd.assays[Condition == cts[2, ct], Assay])))
    dd.bmc <- cbind(dd.proteins[, .(ProteinID)], dd.bmc[, .(log2fc.lower = lower, log2fc.mean = mean, log2fc.upper = upper, PEP, Baseline = cts[1, ct], Condition = cts[2, ct], samp = s)])
    setorder(dd.bmc, PEP)
    dd.bmc[, Discoveries := 1:nrow(dd.bmc)]
    dd.bmc[, FDR := cumsum(PEP) / Discoveries]
    dd.bmc
  }
  doParallel::stopImplicitCluster()
  dd.output[, ProteinID := factor(ProteinID)]
  dd.output[, Baseline := factor(Baseline)]
  dd.output[, Condition := factor(Condition)]
  dd.output[, samp := factor(samp)]
  fst::write.fst(dd.output, paste0(chain, ".fst"))

  # parallel processing - populationLevel K = 3
  message(paste0("[", Sys.time(), "]  populationLevel K = 3 ..."))
  doParallel::registerDoParallel(params$nthread)
  `%dopar%` <- foreach::`%dopar%`
  dd.output <- foreach::foreach(ct = rep(1:ncol(cts), each = nsamp), s = rep(1:nsamp, ncol(cts)), .combine = rbind, .options.multicore = list(preschedule = F, silent = T)) %dopar% {
    dd <- as.data.table(mcmc.quants.all[s,,])
    colnames(dd) <- dd.assays[, Assay]
    bmc <- bayesmodelquant::populationLevel(dd, list(dd.assays[Condition == cts[1, ct], Assay], dd.assays[Condition == cts[2, ct], Assay]))
    dd.bmc <- cbind(dd.proteins[, .(ProteinID)], data.table(log2fc.mean = bmc$mean, PEP = bmc$PEP, Baseline = cts[1, ct], Condition = cts[2, ct], samp = s))
    setorder(dd.bmc, PEP)
    dd.bmc[, Discoveries := 1:nrow(dd.bmc)]
    dd.bmc[, FDR := cumsum(PEP) / 1:nrow(dd.bmc)]
    dd.bmc
  }
  doParallel::stopImplicitCluster()
  dd.output[, ProteinID := factor(ProteinID)]
  dd.output[, Baseline := factor(Baseline)]
  dd.output[, Condition := factor(Condition)]
  dd.output[, samp := factor(samp)]
  fst::write.fst(dd.output, paste0(chain, ".3.fst"))

  # parallel processing - populationLevel K = 11
  message(paste0("[", Sys.time(), "]  populationLevel K = 11 ..."))
  doParallel::registerDoParallel(params$nthread)
  `%dopar%` <- foreach::`%dopar%`
  dd.output <- foreach::foreach(ct = rep(1:ncol(cts), each = nsamp), s = rep(1:nsamp, ncol(cts)), .combine = rbind, .options.multicore = list(preschedule = F, silent = T)) %dopar% {
    dd <- as.data.table(mcmc.quants.all[s,,])
    colnames(dd) <- dd.assays[, Assay]
    bmc <- bayesmodelquant::populationLevel(dd, list(dd.assays[Condition == cts[1, ct], Assay], dd.assays[Condition == cts[2, ct], Assay]), K = 11)
    dd.bmc <- cbind(dd.proteins[, .(ProteinID)], data.table(log2fc.mean = bmc$mean, PEP = bmc$PEP, Baseline = cts[1, ct], Condition = cts[2, ct], samp = s))
    setorder(dd.bmc, PEP)
    dd.bmc[, Discoveries := 1:nrow(dd.bmc)]
    dd.bmc[, FDR := cumsum(PEP) / 1:nrow(dd.bmc)]
    dd.bmc
  }
  doParallel::stopImplicitCluster()
  dd.output[, ProteinID := factor(ProteinID)]
  dd.output[, Baseline := factor(Baseline)]
  dd.output[, Condition := factor(Condition)]
  dd.output[, samp := factor(samp)]
  fst::write.fst(dd.output, paste0(chain, ".11.fst"))

  if (params$qprot) {
    # set up parallel processing, seed and go
    message(paste0("[", Sys.time(), "]  qprot..."))
    doParallel::registerDoParallel(params$nthread)
    set.seed(params$qprot.seed * params$quant.nchain + chain - 1)
    `%dopar%` <- foreach::`%dopar%`
    `%dorng%` <- doRNG::`%dorng%`
    dd.output <- foreach::foreach(ct = rep(1:ncol(cts), each = nsamp), s = rep(1:nsamp, ncol(cts)), .combine = rbind, .options.multicore = list(preschedule = F, silent = T)) %dorng% {
      # process samp s
      dd.0 <- as.data.table(mcmc.quants.all[s,, dd.assays[Condition == cts[1, ct], AssayID]])
      colnames(dd.0) <- rep("0", ncol(dd.0))

      dd.1 <- as.data.table(mcmc.quants.all[s,, dd.assays[Condition == cts[2, ct], AssayID]])
      colnames(dd.1) <- rep("1", ncol(dd.1))

      dd.qprot <- cbind(dd.proteins$ProteinID, dd.0, dd.1)
      colnames(dd.qprot)[1] <- "Protein"

      # exponent as qprot needs intensities, not log ratios
      for (j in 2:ncol(dd.qprot)) dd.qprot[[j]] <- 2^dd.qprot[[j]]

      # because we can't pass seed to qprot, randomise the rows to get the desired effect
      dd.qprot <- dd.qprot[sample(1:nrow(dd.qprot), nrow(dd.qprot)),]

      # run qprot
      filename.qprot <- paste0("_", ct, ".", chain, ".", s, ".tsv")
      fwrite(dd.qprot, filename.qprot, sep = "\t")
      if (params$de.paired) {
        system2(paste0(params$qprot.path, "qprot-paired"),
                args = c(filename.qprot, format(params$qprot.burnin, scientific = F), format(params$qprot.nitt - params$qprot.burnin, scientific = F), "0"),
                stdout = NULL, stderr = NULL)
      } else {
        system2(paste0(params$qprot.path, "qprot-param"),
                args = c(filename.qprot, format(params$qprot.burnin, scientific = F), format(params$qprot.nitt - params$qprot.burnin, scientific = F), "0"),
                stdout = NULL, stderr = NULL)
      }
      system2(paste0(params$qprot.path, "getfdr"), arg = c(paste0(filename.qprot, "_qprot")), stdout = NULL, stderr = NULL)
      dd.qprot <- fread(paste0(filename.qprot, "_qprot_fdr"))[, .(ProteinID = Protein, log2fc.mean = LogFoldChange, Z = Zstatistic, PEP = fdr, Baseline = cts[1, ct], Condition = cts[2, ct], samp = s)]

      file.remove(filename.qprot)
      file.remove(paste0(filename.qprot, "_qprot"))
      file.remove(paste0(filename.qprot, "_qprot_density"))
      file.remove(paste0(filename.qprot, "_qprot_fdr"))

      setorder(dd.qprot, PEP)
      dd.qprot
    }
    doParallel::stopImplicitCluster()
    dd.output[, ProteinID := factor(ProteinID)]
    dd.output[, Baseline := factor(Baseline)]
    dd.output[, Condition := factor(Condition)]
    dd.output[, samp := factor(samp)]
    fst::write.fst(dd.output, paste0(chain, ".qprot.fst"))
  }

  message(paste0("[", Sys.time(), "] BMC finished]"))
}


