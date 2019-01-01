#' process.model (BayesProt internal function)
#'
#' @param chain .
#' @return .
#' @import data.table
#' @export

process.model <- function(chain) {
  # load priors as process.output0 has been run
  prefix <- ifelse(file.exists("priors.rds"), ".", file.path("..", "..", "output0", "results"))
  priors <- readRDS(file.path(prefix, "priors.rds"))
  # and exposures, shifting so first assay is zero
  dd.assay.exposures <- fst::read.fst(file.path(prefix, "assay.exposures.fst"), as.data.table = T)
  dd.assay.exposures <- merge(dd.assay.exposures, dd.assay.exposures[AssayID == levels(AssayID)[1], .(chainID, mcmcID, shift = exposure)])
  dd.assay.exposures[, exposure := exposure - shift]
  dd.assay.exposures <- dd.assay.exposures[, .(mean = mean(exposure), var = var(exposure)), by = AssayID]

  # load metadata
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)[!is.na(model.row0),]
  nitt <- params$model.nwarmup + (params$model.nsample * params$model.thin) / params$model.nchain
  message(paste0("[", Sys.time(), "] MODEL started, chain=", chain, "/", params$model.nchain, " nitt=", nitt))

  # our combine function
  rbindlistlist <- function(...) {
    input <- list(...)
    for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
    input[[1]]
  }

  # set up parallel processing, seed and go
  suppressPackageStartupMessages(require(doRNG))
  cl <- parallel::makeCluster(params$nthread)
  doSNOW::registerDoSNOW(cl)
  options(max.print = 99999)
  pb <- txtProgressBar(max = sum(!is.na(dd.proteins$model.row0)), style = 3)
  setTxtProgressBar(pb, 0)
  set.seed(params$model.seed * params$model.nchain + chain - 1)

  output <- foreach::foreach(i = which(!is.na(dd.proteins$model.row0)), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.multicore = list(preschedule = F), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    # prepare dd for MCMCglmm
    dd <- fst::read.fst(file.path(prefix, "data.fst"), as.data.table = T, from = dd.proteins[i, model.row0], to = dd.proteins[i, model.row1])
    dd <- droplevels(dd)
    dd[, Count := round(Count)]
    if (!is.null(dd$MaxCount)) dd[, MaxCount := round(MaxCount)]

    # set prior and run model
    nA <- length(levels(dd$AssayID))
    nT <- length(levels(dd$PeptideID))
    nF <- length(levels(dd$FeatureID))

    prior <- list(
      B = list(mu = matrix(0, nF + nA - 1, 1), V = diag(nF + nA - 1) * 1e+6),
      G = list(
        G1 = list(V = priors$protein.V, nu = priors$protein.nu),
        #G2 = list(V = priors$peptide.V * diag(nT), nu = priors$peptide.nu)
        G2 = list(V = diag(nT), nu = 0.02)
      ),
      #R = list(V = priors$feature.V * diag(nF), nu = priors$feature.nu)
      R = list(V = diag(nF), nu = 0.02)
    )
    prior$B$mu[(nF + 1):(nF + nA - 1)] <- dd.assay.exposures$mean[2:nrow(dd.assay.exposures)]
    diag(prior$B$V)[(nF + 1):(nF + nA - 1)] <- dd.assay.exposures$var[2:nrow(dd.assay.exposures)]

    #run model
    output <- list()
    output$dd.summary <- as.character(Sys.time())
    output$dd.timing <- system.time(model <- (MCMCglmm::MCMCglmm(
      as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ ", ifelse(nF==1, "", "FeatureID-1 +"), " AssayID")),
      random = as.formula(paste0("~ SampleID +", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":DigestID")),
      rcov = as.formula(paste0("~ ", ifelse(nF==1, "AssayID", "idh(FeatureID):AssayID"))),
      family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"), data = dd, prior = prior, singular.ok = T,
      nitt = nitt, burnin = params$model.nwarmup, thin = params$model.thin, pr = T, verbose = F
    )))
    output$dd.timing <- data.table(ProteinID = dd[1, ProteinID], as.data.table(t(as.matrix(output$dd.timing))))
    output$dd.summary <- data.table(ProteinID = dd[1, ProteinID], Summary = paste(c(output$dd.summary, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

    # extract quants, shifting so that denominator is mean of reference assays
    output$dd.protein.quants <- as.data.table(model$Sol[, grep("^SampleID\\.[0-9]+$", colnames(model$Sol)), drop = F])
    output$dd.protein.quants[, mcmcID := factor(formatC(1:nrow(output$dd.protein.quants), width = ceiling(log10(nrow(output$dd.protein.quants))) + 1, format = "d", flag = "0"))]
    output$dd.protein.quants <- melt(output$dd.protein.quants, variable.name = "SampleID", id.vars = "mcmcID")
    output$dd.protein.quants[, ProteinID := dd[1, ProteinID]]
    output$dd.protein.quants[, SampleID := factor(sub("^SampleID\\.", "", SampleID))]
    setcolorder(output$dd.protein.quant, c("ProteinID", "SampleID"))

    # extract peptide deviations from protein quants
    if (nT==1) {
      output$dd.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID:DigestID\\.[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
      output$dd.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$dd.peptide.deviations), width = ceiling(log10(nrow(output$dd.peptide.deviations))) + 1, format = "d", flag = "0"))]
      output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
      output$dd.peptide.deviations[, DigestID := factor(sub("^PeptideID:DigestID\\.[0-9]+\\.([0-9]+)$", "\\1", PeptideID))]
      output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID:DigestID\\.([0-9]+)\\.[0-9]+$", "\\1", PeptideID))]
    } else {
      output$dd.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.DigestID\\.[0-9]+$", colnames(model$Sol)), drop = F])
      output$dd.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$dd.peptide.deviations), width = ceiling(log10(nrow(output$dd.peptide.deviations))) + 1, format = "d", flag = "0"))]
      output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
      output$dd.peptide.deviations[, DigestID := factor(sub("^PeptideID[0-9]+\\.DigestID\\.([0-9]+)$", "\\1", PeptideID))]
      output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.DigestID\\.([0-9]+)$", "\\1", PeptideID))]
    }
    output$dd.peptide.deviations[, ProteinID := dd[1, ProteinID]]
    setcolorder(output$dd.peptide.deviations, c("ProteinID", "PeptideID", "DigestID"))

    model$Sol <- NULL

    # extract peptide variances
    if (nT == 1) {
      output$dd.peptide.vars <- as.data.table(model$VCV[, "PeptideID:DigestID", drop = F])
      output$dd.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$dd.peptide.vars), width = ceiling(log10(nrow(output$dd.peptide.vars))) + 1, format = "d", flag = "0"))]
      output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
      output$dd.peptide.vars[, PeptideID := factor(levels(dd$PeptideID))]
    } else {
      output$dd.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.DigestID$", colnames(model$VCV)), drop = F])
      output$dd.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$dd.peptide.vars), width = ceiling(log10(nrow(output$dd.peptide.vars))) + 1, format = "d", flag = "0"))]
      output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
      output$dd.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.DigestID$", "\\1", PeptideID))]
    }
    output$dd.peptide.vars[, ProteinID := dd[1, ProteinID]]
    setcolorder(output$dd.peptide.vars, c("ProteinID", "PeptideID"))

    # extract feature variances
    if (nF == 1) {
      output$dd.feature.vars <- as.data.table(model$VCV[, "AssayID", drop = F])
      output$dd.feature.vars[, mcmcID := factor(formatC(1:nrow(output$dd.feature.vars), width = ceiling(log10(nrow(output$dd.feature.vars))) + 1, format = "d", flag = "0"))]
      output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
      output$dd.feature.vars[, FeatureID := factor(levels(dd$FeatureID))]
    } else {
      output$dd.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
      output$dd.feature.vars[, mcmcID := factor(formatC(1:nrow(output$dd.feature.vars), width = ceiling(log10(nrow(output$dd.feature.vars))) + 1, format = "d", flag = "0"))]
      output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
      output$dd.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.AssayID$", "\\1", FeatureID))]
    }
    output$dd.feature.vars[, ProteinID := dd[1, ProteinID]]
    setcolorder(output$dd.feature.vars, c("ProteinID", "FeatureID"))

    output
  }
  close(pb)

  # write output
  chain <- formatC(chain, width = ceiling(log10(params$model.nchain + 1)) + 1, format = "d", flag = "0")
  fst::write.fst(output$dd.summary, paste0("summary.", chain, ".fst"))
  output$dd.summary <- NULL
  fst::write.fst(output$dd.timing, paste0("timing.", chain, ".fst"))
  output$dd.timing <- NULL
  fst::write.fst(output$dd.peptide.deviations, paste0("peptide.deviations.", chain, ".fst"))
  output$dd.peptide.deviations <- NULL
  fst::write.fst(output$dd.peptide.vars, paste0("peptide.vars.", chain, ".fst"))
  output$dd.peptide.vars <- NULL
  fst::write.fst(output$dd.feature.vars, paste0("feature.vars.", chain, ".fst"))
  output$dd.feature.vars <- NULL
  fst::write.fst(output$dd.protein.quants, paste0("protein.quants.", chain, ".fst"))

  # DE
  if (!is.null(dd.assays$ConditionID) & params$de.mcmc) {
    bmc <- function(SampleID, value) {
      out <- bayesmodelquant::modelComparison(value[SampleID %in% sampleIDs0], value[SampleID %in% sampleIDs1])
      data.table(
        stdev = out$sd, Z = out$Zstat,
        fc.lower = out$HPDI$lower[1], fc.mean = out$mean, fc.upper = out$HPDI$upper[1],
        PEP.down = out$bpFDR$down[1], PEP.same = out$bpFDR$same[1], PEP.up = out$bpFDR$up[1], PEP = out$PEP[1]
      )
    }

    # all combinations of conditions with at least two samples
    cts <- combn(dd.assays[, length(unique(SampleID)) >= 2, by = ConditionID][V1 == T, ConditionID], 2)
    for (ct in 1:ncol(cts)) {
      ct0 <- dd.assays[ConditionID == cts[1, ct], Condition][1]
      ct1 <- dd.assays[ConditionID == cts[2, ct], Condition][1]
      sampleIDs0 <- dd.assays[ConditionID == cts[1, ct], SampleID]
      sampleIDs1 <- dd.assays[ConditionID == cts[2, ct], SampleID]
      sampleIDs <- dd.assays[ConditionID == cts[1, ct] | ConditionID == cts[2, ct], SampleID]
      protein.quants.ct <- split(output$dd.protein.quants[SampleID %in% sampleIDs,], by = "mcmcID")
      message(paste0("[", Sys.time(), "]  BMC MCMC diffential analysis for ", ct1, " vs ", ct0, "..."))

      pb <- txtProgressBar(max = length(protein.quants.ct), style = 3)
      dd.output <- foreach::foreach(dd = iterators::iter(protein.quants.ct), .packages = "data.table", .combine = rbind, .options.multicore = list(preschedule = F), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
        s <- dd[1, mcmcID]
        dd <- dd[, as.list(bmc(SampleID, value)), by = ProteinID]

        dd[, abs.Z := abs(Z)]
        dd[, abs.fc.mean := abs(fc.mean)]
        setorder(dd, PEP, -abs.Z, -abs.fc.mean, na.last = T)
        dd[, abs.Z := NULL]
        dd[, abs.fc.mean := NULL]

        dd[, Discoveries := 1:nrow(dd)]
        dd[, FDR := cumsum(PEP) / Discoveries]
        dd[, mcmcID := factor(s)]
        dd
      }
      close(pb)
      fst::write.fst(dd.output, paste0("de.bmc.", chain, ".", formatC(ct, width = ceiling(log10(nrow(cts))) + 1, format = "d", flag = "0"), ".fst"))

      if (params$qprot) {
        message(paste0("[", Sys.time(), "]  Qprot MCMC diffential analysis for ", ct1, " vs ", ct0, "..."))

        pb <- txtProgressBar(max = length(protein.quants.ct), style = 3)
        set.seed(params$qprot.seed)
        dd.output <- foreach::foreach(dd = iterators::iter(protein.quants.ct), .packages = "data.table", .combine = rbind, .options.multicore = list(preschedule = F), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
          s <- dd[1, mcmcID]
          dd <- dcast(dd, ProteinID ~ SampleID, value.var = "value")
          colnames(dd)[colnames(dd) %in% sampleIDs0] <- "0"
          colnames(dd)[colnames(dd) %in% sampleIDs1] <- "1"

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
          dd.qprot <- fread(paste0(filename, "_qprot_fdr"))[, .(ProteinID = Protein, fc.mean = LogFoldChange, Z = Zstatistic, PEP.down = FDRdown, PEP.up = FDRup, PEP = fdr, mcmcID = s)]

          file.remove(filename)
          file.remove(paste0(filename, "_qprot"))
          file.remove(paste0(filename, "_qprot_density"))
          file.remove(paste0(filename, "_qprot_fdr"))

          dd.qprot[, abs.Z := abs(Z)]
          dd.qprot[, abs.fc.mean := abs(fc.mean)]
          setorder(dd.qprot, PEP, -abs.Z, -abs.fc.mean, na.last = T)
          dd.qprot[, abs.Z := NULL]
          dd.qprot[, abs.fc.mean := NULL]
          dd.qprot[, FDR := cumsum(PEP) / 1:nrow(dd.qprot)]
          dd.qprot
        }
        close(pb)
        fst::write.fst(dd.output, paste0("de.qprot.", chain, ".", formatC(ct, width = ceiling(log10(nrow(cts))) + 1, format = "d", flag = "0"), ".fst"))
      }
    }
  }
  parallel::stopCluster(cl)

  message(paste0("[", Sys.time(), "] MODEL finished"))
}

