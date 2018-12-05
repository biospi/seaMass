#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @importFrom foreach %dopar%
#' @importFrom doRNG %dorng%
#' @export

process.qprot <- function(chain) {
  message(paste0("[", Sys.time(), "] QPROT started, chain=", chain))

  # load parameters
  prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
  load(file.path(prefix, "metadata.Rdata"))
  nsamp <- (params$quant.nitt - params$quant.burnin) / params$quant.thin

  nP <- length(levels(dd.proteins$ProteinID))
  nA <- length(levels(dd.assays$AssayID))

  prefix <- ifelse(file.exists("1.1.rds"), ".", file.path("..", "..", "quant", "results", "quants"))

  # process design
  if (!is.factor(dd.de.design$Assay)) dd.de.design$Assay <- factor(dd.de.design$Assay, levels = unique(dd.de.design$Assay))
  if (!is.factor(dd.de.design$Condition)) dd.de.design$Condition <- factor(dd.de.design$Condition, levels = unique(dd.de.design$Condition))
  dd.de.design <- as.data.table(merge(dd.assays, dd.de.design))

  # read MCMC samps
  mcmc.quants.all <- array(NA, c(nsamp, nP, nA))
  for (a in 1:nA) {
    print(paste0("[", Sys.time(), " Reading ", a, ".", chain, ".rds ...]"))

    mcmc.quants <- readRDS(file.path(prefix, paste0(a, ".", chain, ".rds")))

    # todo: check baseline
    mcmc.quants.all[, as.integer(sub("\\.[0-9]+$", "", colnames(mcmc.quants))), a] <- mcmc.quants
  }

  # set up parallel processing, seed and go
  doParallel::registerDoParallel(params$nthread)
  set.seed(params$seed * params$quant.nchain + chain - 1)
  cts <- levels(dd.de.design$Condition)[2:length(levels(dd.de.design$Condition))]
  output <- foreach::foreach(ct = rep(1:length(cts), each = nsamp), s = rep(1:nsamp, length(cts)), .packages = "data.table", .options.multicore = list(preschedule = F, silent = T)) %dorng% {
    # process samp s
    dd.0 <- as.data.table(mcmc.quants.all[s,, dd.de.design[Condition == levels(dd.de.design$Condition)[1], AssayID]])
    colnames(dd.0) <- rep("0", ncol(dd.0))

    dd.1 <- as.data.table(mcmc.quants.all[s,, dd.de.design[Condition == cts[ct], AssayID]])
    colnames(dd.1) <- rep("1", ncol(dd.1))

    dd.qprot <- cbind(dd.proteins$ProteinID, dd.0, dd.1)
    colnames(dd.qprot)[1] <- "Protein"

    # remove NAs (check qprot requirements)
    dd.qprot <- dd.qprot[complete.cases(dd.qprot[, 2:ncol(dd.qprot)]),]

    # exponent as qprot needs intensities, not log ratios
    for (j in 2:ncol(dd.qprot)) dd.qprot[[j]] <- 2^dd.qprot[[j]]

    # because we can't pass seed to qprot, randomise the rows to get the desired effect
    dd.qprot <- dd.qprot[sample(1:nrow(dd.qprot), nrow(dd.qprot)),]

    # run qprot
    filename.qprot <- paste0("_", ct, ".", chain, ".", s, ".tsv")
    fwrite(dd.qprot, filename.qprot, sep = "\t")
    if (params$de.paired) {
      system2(paste0(params$qprot.path, "qprot-paired"), args = c(filename.qprot, format(params$qprot.burnin, scientific = F), format(params$qprot.nitt - params$qprot.burnin, scientific = F), "0"), stdout = NULL, stderr = NULL)
    } else {
      system2(paste0(params$qprot.path, "qprot-param"), args = c(filename.qprot, format(params$qprot.burnin, scientific = F), format(params$qprot.nitt - params$qprot.burnin, scientific = F), "0"), stdout = NULL, stderr = NULL)
    }
    system2(paste0(params$qprot.path, "getfdr"), arg = c(paste0(filename.qprot, "_qprot")))
    dd.fdr <- fread(paste0(filename.qprot, "_qprot_fdr"))

    file.remove(filename.qprot)
    file.remove(paste0(filename.qprot, "_qprot"))
    file.remove(paste0(filename.qprot, "_qprot_density"))
    file.remove(paste0(filename.qprot, "_qprot_fdr"))

    setnames(dd.fdr, c("Protein", "fdr", "FDRup", "FDRdown"), c("ProteinID", "PEP", "PEPup", "PEPdown"))
    dd.fdr[, Condition := cts[ct]]
    dd.fdr[, samp := s]
    dd.fdr
  }
  doParallel::stopImplicitCluster()

  # save
  saveRDS(rbindlist(output), file = paste0(chain, ".rds"))

  message(paste0("[", Sys.time(), "] QPROT finished]"))
}


