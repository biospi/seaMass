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
  nsamp <- params$quant.nsample / params$quant.thin

  nP <- length(levels(dd.proteins$ProteinID))
  nA <- length(levels(dd.assays$AssayID))

  # read MCMC samps
  prefix <- ifelse(file.exists(paste0("protein.quants.", chain, ".fst")), ".", file.path("..", "..", "model2", "results"))
  dd.protein.quants <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)

  # generate all combinations of conditions with at least three samples
  cts <- combn(dd.assays[, .N >= 3, by = Condition][V1 == T, Condition], 2)

  # for each contrast
  for (ct in 1:ncol(cts)) {
    message(paste0("[", Sys.time(), "]  ", cts[1, ct], " vs ", cts[2, ct], "..."))

    assays0 <- as.character(dd.assays[Condition == cts[1, ct], Assay])
    assayIDs0 <- as.integer(dd.assays[Condition == cts[1, ct], AssayID])
    assays1 <- as.character(dd.assays[Condition == cts[2, ct], Assay])
    assayIDs1 <- as.integer(dd.assays[Condition == cts[2, ct], AssayID])

    # parallel processing - modelComparisonBatch
    suppressPackageStartupMessages(require(doRNG))

    cl <- parallel::makeCluster(params$nthread)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(levels(dd.proteins$ProteinID)), style = 3)
    setTxtProgressBar(pb, 0)
    dd.output <- foreach::foreach(s = 1:nsamp, .combine = rbind, .options.multicore = list(preschedule = F),
                                  .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
      dd <- dd.protein.quants[samp == s & AssayID %in% c(assayIDs0, assayIDs1), ]
      dd <- merge(dd.assays, dd, by = "AssayID")
      dd <- dcast(dd, ProteinID ~ Assay, value.var = "value")
      dd <- dd[complete.cases(dd),]
      dd.bmc <- as.data.table(bayesmodelquant::modelComparisonBatch(dd, list(assays0, assays1)))
      dd.bmc <- cbind(dd[, .(ProteinID)], dd.bmc[, .(log2fc.lower = lower, log2fc.mean = mean, log2fc.upper = upper, PEP, samp = s)])
      setorder(dd.bmc, PEP)
      dd.bmc[, Discoveries := 1:nrow(dd.bmc)]
      dd.bmc[, FDR := cumsum(PEP) / Discoveries]
      dd.bmc
    }
    close(pb)
    dd.output[, ProteinID := factor(ProteinID)]
    dd.output[, samp := factor(samp)]
    fst::write.fst(dd.output, paste0(ct, ".", chain, ".bmc0.fst"))

    if (params$qprot) {
      # set up parallel processing, seed and go
      message(paste0("[", Sys.time(), "]   Qprot..."))
      pb <- txtProgressBar(max = length(levels(dd.proteins$ProteinID)), style = 3)
      setTxtProgressBar(pb, 0)
      set.seed(params$quant.seed * params$quant.nchain + chain - 1)
      dd.output <- foreach::foreach(s = 1:nsamp, .combine = rbind, .options.multicore = list(preschedule = F),
                                    .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
        # process samp s
        dd.0 <- dd.protein.quants[samp == s & AssayID %in% assayIDs0, ]
        dd.0 <- dcast(dd.0, ProteinID ~ AssayID, value.var = "value")
        colnames(dd.0)[2:ncol(dd.0)] <- paste("0", colnames(dd.0)[2:ncol(dd.0)])
        dd.1 <- dd.protein.quants[samp == s & AssayID %in% assayIDs1, ]
        dd.1 <- dcast(dd.1, ProteinID ~ AssayID, value.var = "value")
        colnames(dd.1)[2:ncol(dd.1)] <- paste("1", colnames(dd.1)[2:ncol(dd.1)])
        dd.qprot <- merge(dd.0, dd.1, by = "ProteinID")
        colnames(dd.qprot)[2:ncol(dd.qprot)] <- substr(colnames(dd.qprot)[2:ncol(dd.qprot)], start = 1, stop = 1)

        # exponent as qprot needs intensities, not log ratios
        for (j in 2:ncol(dd.qprot)) dd.qprot[[j]] <- 2^dd.qprot[[j]]

        # missing data needs to be set as zeros, as in qprot vignette!
        for (j in 2:ncol(dd.qprot)) dd.qprot[[j]][is.na(dd.qprot[[j]])] <- 0

        # remove rows with less than 6 non-zeros
        #dd.qprot$nnz <- apply(dd.qprot, 1, function(x) (length(x) - sum(x == 0)))
        #dd.qprot <- cbind(dd.proteins$ProteinID, dd.qprot)[nnz >= 6, -"nnz"]
        #colnames(dd.qprot)[1] <- "Protein"

        # because we can't pass seed to qprot, randomise the rows to get the desired effect
        set.seed(params$qprot.seed)
        dd.qprot <- dd.qprot[sample(1:nrow(dd.qprot), nrow(dd.qprot)),]

        # run qprot
        filename.qprot <- paste0("_", ct, ".", chain, ".", s, ".tsv")
        fwrite(dd.qprot, filename.qprot, sep = "\t")
        ret <- 1
        if (params$de.paired) {
          while (ret != 0) {
            ret <- system2(ifelse(params$qprot.path == "", "qprot-paired", file.path(params$qprot.path, "qprot-paired")),
                           args = c(filename.qprot, format(params$qprot.burnin, scientific = F), format(params$qprot.nitt - params$qprot.burnin, scientific = F), "0"),
                           stdout = NULL, stderr = NULL)
          }
        } else {
          while (ret != 0) {
            ret <- system2(ifelse(params$qprot.path == "", "qprot-param", file.path(params$qprot.path, "qprot-param")),
                           args = c(filename.qprot, format(params$qprot.burnin, scientific = F), format(params$qprot.nitt - params$qprot.burnin, scientific = F), "0"),
                           stdout = NULL, stderr = NULL)
          }
        }
        ret <- 1
        while (ret != 0) {
          ret <- system2(ifelse(params$qprot.path == "", "getfdr", file.path(params$qprot.path, "getfdr")), arg = c(paste0(filename.qprot, "_qprot")), stdout = NULL, stderr = NULL)
        }
        dd.qprot <- fread(paste0(filename.qprot, "_qprot_fdr"))[, .(ProteinID = Protein, log2fc.mean = LogFoldChange / log(2), Z = Zstatistic, PEP = fdr, samp = s)]
        setorder(dd.qprot, PEP)
        dd.qprot[, Discoveries := 1:nrow(dd.qprot)]
        dd.qprot[, FDR := cumsum(PEP) / 1:nrow(dd.qprot)]

        file.remove(filename.qprot)
        file.remove(paste0(filename.qprot, "_qprot"))
        file.remove(paste0(filename.qprot, "_qprot_density"))
        file.remove(paste0(filename.qprot, "_qprot_fdr"))

        dd.qprot
      }
      close(pb)
      dd.output[, ProteinID := factor(ProteinID)]
      dd.output[, samp := factor(samp)]
      fst::write.fst(dd.output, paste0(ct, ".", chain, ".qprot.fst"))
    }
    parallel::stopCluster(cl)
  }

  message(paste0("[", Sys.time(), "] BMC finished]"))
}


