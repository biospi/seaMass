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
  chain <- formatC(chain, width = ceiling(log10(params$quant.nchain + 1)) + 1, format = "d", flag = "0")
  nA <- length(levels(dd.assays$AssayID))

  # read MCMC samps
  prefix <- ifelse(file.exists(paste0("protein.quants.", chain, ".fst")), ".", file.path("..", "..", "model2", "results"))
  dd.protein.quants <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)

  # generate all combinations of conditions with at least three samples
  cts <- combn(dd.assays[, .N >= 3, by = Condition][V1 == T, Condition], 2)

  # for each contrast
  for (ct in 1:ncol(cts)) {
    message(paste0("[", Sys.time(), "]  ", cts[1, ct], " vs ", cts[2, ct], "..."))

    suppressPackageStartupMessages(require(doRNG))

    assayIDs0 <- dd.assays[Condition == cts[1, ct], AssayID]
    assayIDs1 <- dd.assays[Condition == cts[2, ct], AssayID]
    assayIDs <- dd.assays[Condition == cts[1, ct] | Condition == cts[2, ct], AssayID]
    protein.quants.ct <- split(dd.protein.quants[AssayID %in% assayIDs,], by = "mcmcID")

    # parallel processing - modelComparisonBatch
    message(paste0("[", Sys.time(), "]   BMC0..."))

    cl <- parallel::makeCluster(params$nthread)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(protein.quants.ct), style = 3)
    dd.output <- foreach::foreach(
      dd = iterators::iter(protein.quants.ct), .packages = "data.table", .combine = rbind,
      .options.multicore = list(preschedule = F), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
        s <- dd[1, mcmcID]
        bmc <- function(value, AssayID) {
          out <- bayesmodelquant::modelComparison(value[AssayID %in% assayIDs0], value[AssayID %in% assayIDs1])
          data.table(log2fc.lower = out$HPDI$lower[1], log2fc.mean = out$mean, log2fc.upper = out$HPDI$upper[1], PEP = out$PEP[1])
        }
        dd <- dd[, as.list(bmc(value, AssayID)), by = ProteinID]
        dd[, mcmcID := s]
        setorder(dd, PEP)
        dd[, Discoveries := 1:nrow(dd)]
        dd[, FDR := cumsum(PEP) / Discoveries]
        dd
    }
    close(pb)
    fst::write.fst(dd.output, paste0(ct, ".", chain, ".bmc0.fst"))

    if (params$qprot) {
      # set up parallel processing, seed and go
      message(paste0("[", Sys.time(), "]   Qprot..."))
      pb <- txtProgressBar(max = length(protein.quants.ct), style = 3)
      set.seed(params$qprot.seed)
      dd.output <- foreach::foreach(
        dd = iterators::iter(protein.quants.ct), .packages = "data.table", .combine = rbind,
        .options.multicore = list(preschedule = F), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
          s <- dd[1, mcmcID]
          dd <- dcast(dd, ProteinID ~ AssayID, value.var = "value")
          colnames(dd)[colnames(dd) %in% assayIDs0] <- "0"
          colnames(dd)[colnames(dd) %in% assayIDs1] <- "1"

          # exponent as qprot needs intensities, not log ratios
          for (j in 2:ncol(dd)) dd[[j]] <- 2^dd[[j]]

          # missing data needs to be set as zeros, as in qprot vignette!
          for (j in 2:ncol(dd)) dd[[j]][is.na(dd[[j]])] <- 0

          # because we can't pass seed to qprot, randomise the rows to get the desired effect
          dd <- dd[sample(1:nrow(dd), nrow(dd)),]

          # run qprot
          filename <- paste0("_", ct, ".", chain, ".", s, ".tsv")
          fwrite(dd, filename, sep = "\t")
          ret <- 1
          if (params$de.paired) {
            for (i in 1:3) {
             if (0 == system2(ifelse(params$qprot.path == "", "qprot-paired", file.path(params$qprot.path, "qprot-paired")),
                              args = c(filename, format(params$qprot.nwarmup, scientific = F), format(params$qprot.nsample, scientific = F), "0"),
                              stdout = NULL, stderr = NULL)) break
            }
          } else {
            for (i in 1:3) {
              if (0 == system2(ifelse(params$qprot.path == "", "qprot-param", file.path(params$qprot.path, "qprot-param")),
                               args = c(filename, format(params$qprot.nwarmup, scientific = F), format(params$qprot.nsample, scientific = F), "0"),
                               stdout = NULL, stderr = NULL)) break
            }
          }
          ret <- 1
          for (i in 1:3) {
            if (0 == system2(ifelse(params$qprot.path == "", "getfdr", file.path(params$qprot.path, "getfdr")),
                             args = c(paste0(filename, "_qprot")), stdout = NULL, stderr = NULL)) break
          }
          dd <- fread(paste0(filename, "_qprot_fdr"))[, .(ProteinID = Protein, log2fc.mean = LogFoldChange / log(2), Z = Zstatistic, PEP = fdr, mcmcID = s)]
          setorder(dd, PEP)
          dd[, Discoveries := 1:nrow(dd)]
          dd[, FDR := cumsum(PEP) / 1:nrow(dd)]

          file.remove(filename)
          file.remove(paste0(filename, "_qprot"))
          file.remove(paste0(filename, "_qprot_density"))
          file.remove(paste0(filename, "_qprot_fdr"))

          dd
      }
      close(pb)
      fst::write.fst(dd.output, paste0(ct, ".", chain, ".qprot.fst"))
    }
    parallel::stopCluster(cl)
  }

  message(paste0("[", Sys.time(), "] BMC finished]"))
}


