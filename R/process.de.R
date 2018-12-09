#' process.de (BayesProt internal function)
#'
#' @return .
#' @import data.table
#' @export

process.de <- function() {
  message(paste0("[", Sys.time(), "] DE started"))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
  nsamp <- (params$quant.nitt - params$quant.burnin) / params$quant.thin

  stats.dir <- paste0(params$id, ".de")
  dir.create(stats.dir, showWarnings = F)

  nP <- length(levels(dd.proteins$ProteinID))
  nA <- length(levels(dd.assays$AssayID))

  cts <- combn(levels(dd.assays$Condition), 2)

  prefix <- ifelse(file.exists("protein_quants.csv"), ".", file.path("..", "..", "quant", "results", paste0(params$id, ".quant")))
  dd <- fread(file.path(prefix, "protein_quants.csv"))
  prefix <- ifelse(file.exists("1.fst"), ".", file.path("..", "..", "bmc", "results"))
  for (ct in 1:ncol(cts)) {

    # BMC ON POSTERIOR MEANS ONLY, FOR FUN
    dd.bmc <- as.data.table(bayesmodelquant::modelComparisonBatch(dd, list(paste("x", dd.assays[Condition == cts[1, ct], Assay]), paste("x", dd.assays[Condition == cts[2, ct], Assay]))))
    dd.bmc <- cbind(dd.proteins, dd.bmc[, .(log2fc.lower = lower, log2fc.mean = mean, log2fc.upper = upper, PEP)])
    setorder(dd.bmc, PEP)
    dd.bmc[, FDR := cumsum(PEP) / 1:nrow(dd.bmc)]
    fwrite(dd.bmc, file.path(stats.dir, paste0("protein_de__", cts[1, ct], "_vs_", cts[2, ct], "__bmc__point_est.csv")))

    bmc3 <- bayesmodelquant::populationLevel(dd, list(paste("x", dd.assays[Condition == cts[1, ct], Assay]), paste("x", dd.assays[Condition == cts[2, ct], Assay])))
    dd.bmc3 <- cbind(dd.proteins, data.table(log2fc.mean = bmc3$mean, PEP = bmc3$PEP))
    setorder(dd.bmc3, PEP)
    dd.bmc3[, FDR := cumsum(PEP) / 1:nrow(dd.bmc3)]
    fwrite(dd.bmc3, file.path(stats.dir, paste0("protein_de__", cts[1, ct], "_vs_", cts[2, ct], "__bmc3__point_est.csv")))

    bmc11 <- bayesmodelquant::populationLevel(dd, list(paste("x", dd.assays[Condition == cts[1, ct], Assay]), paste("x", dd.assays[Condition == cts[2, ct], Assay])), K = 11)
    dd.bmc11 <- cbind(dd.proteins, data.table(log2fc.mean = bmc11$mean, PEP = bmc11$PEP))
    setorder(dd.bmc11, PEP)
    dd.bmc11[, FDR := cumsum(PEP) / 1:nrow(dd.bmc11)]
    fwrite(dd.bmc11, file.path(stats.dir, paste0("protein_de__", cts[1, ct], "_vs_", cts[2, ct], "__bmc11__point_est.csv")))

    dds <- list("BayesProt/BMC" = dd.bmc, "BayesProt/BMC3" = dd.bmc3, "BayesProt/BMC11" = dd.bmc11)
    g <- fdr.plot(dds)
    ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__", cts[1, ct], "_vs_", cts[2, ct], "__bmc__point_est.pdf")), g, width=8, height=8)

    if (!is.null(params$truth)) {
      g <- pr.plot(dds, params$truth, 0.21)
      ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__", cts[1, ct], "_vs_", cts[2, ct], "__bmc__point_est__fdp.pdf")) , g, width=8, height=8)
    }

    # GATHER PEPs FROM BMC ON MCMC SAMPLES
    bmc.mcmc <- function(method) {
      dd.bmc.mcmc <- vector("list", params$quant.nchain)
      for (j in 1:params$quant.nchain) {
        dd.bmc.mcmc[[j]] <- fst::read.fst(file.path(prefix, paste0(j, method, ".fst")), as.data.table = T)[Baseline == cts[1, ct] & Condition == cts[2, ct],]
        dd.bmc.mcmc[[j]][, chain := j]
      }
      dd.bmc.mcmc <- rbindlist(dd.bmc.mcmc)

      # Discoveries based on mean FDRs
      dd.bmc.fdr <- merge(dd.proteins, dd.bmc.mcmc[, .(
        log2fc.mean = mean(log2fc.mean),
        PEP.lower = coda::HPDinterval(coda::as.mcmc(PEP))[, "lower"], PEP.mean = mean(PEP), PEP.upper = coda::HPDinterval(coda::as.mcmc(PEP))[, "upper"],
        FDR.mean = mean(FDR)
      ), by = ProteinID], sort = F)
      setorder(dd.bmc.fdr, FDR.mean)
      dd.bmc.fdr[, FDR.mean := NULL]
      dd.bmc.fdr[, Discoveries := 1:nrow(dd.bmc.fdr)]

      # with this Discoveries, recompute FDR for each samp to derive credible interval
      dd.bmc.mcmc.fdr <- merge(dd.bmc.fdr[, .(ProteinID, Discoveries)], dd.bmc.mcmc[, .(ProteinID, PEP, samp, chain)], sort = F)
      setorder(dd.bmc.mcmc.fdr, Discoveries, samp, chain)
      dd.bmc.mcmc.fdr <- dd.bmc.mcmc.fdr[, .(Discoveries, FDR = cumsum(PEP) / Discoveries), by = .(samp, chain)]
      dd.bmc.fdr <- merge(dd.bmc.fdr, dd.bmc.mcmc.fdr[, .(
        FDR.lower = coda::HPDinterval(coda::as.mcmc(FDR))[, "lower"], FDR.mean = mean(FDR), FDR.upper = coda::HPDinterval(coda::as.mcmc(FDR))[, "upper"]
      ), by = Discoveries])
      dd.bmc.fdr[, Discoveries := NULL]
      dd.bmc.fdr[, FDR := FDR.mean]

      fwrite(dd.bmc.fdr, file.path(stats.dir, paste0("protein_de__", cts[1, ct], "_vs_", cts[2, ct], "__bmc", method, ".csv")))

      dd.bmc.fdr
    }

    dd.bmc.fdr <- bmc.mcmc("")
    dd.bmc3.fdr <- bmc.mcmc(".3")
    dd.bmc11.fdr <- bmc.mcmc(".11")
    dds <- list("BayesProtMCMC/BMC" = dd.bmc.fdr, "BayesProtMCMC/BMC3" = dd.bmc3.fdr, "BayesProtMCMC/BMC11" = dd.bmc11.fdr)
    g <- fdr.plot(dds)
    ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__", cts[1, ct], "_vs_", cts[2, ct], "__bmc.pdf")), g, width=8, height=8)

    if (!is.null(params$truth)) {
      g <- pr.plot(dds, params$truth, 0.21)
      ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__", cts[1, ct], "_vs_", cts[2, ct], "__bmc__fdp.pdf")) , g, width=8, height=8)
    }
  }

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), ": DE finished]"))
}
