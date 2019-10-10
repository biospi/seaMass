#' process_model0 (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
#' @export
process_model0 <- function(
  fit,
  chain = 1
) {
  control <- control(fit)
  set.seed(control$model.seed + chain)
  start_parallel(fit)

  # EXECUTE MODEL
  execute_model(fit, chain)

  if (length(list.files(file.path(fit, "model0"), "\\.finished$")) == control$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT stage=eb"))
    priors <- list()

    message("[", paste0(Sys.time(), "]  computing empirical Bayes feature variance prior..."))
    DT.feature.vars <- feature_vars(fit, stage = "0", as.data.table = T)

    # plot feature variance fit
    DT.feature.vars[, nPeptide := .N, by = ProteinID]
    g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = FeatureID, y = rhat)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
    ggplot2::ggsave(file.path(fit, "model0", "feature_vars0_rhat.pdf"), g, width = 8, height = 6, limitsize = F)
    g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = FeatureID, y = v)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
    ggplot2::ggsave(file.path(fit, "model0", "feature_vars0_v.pdf"), g, width = 8, height = 6, limitsize = F)
    g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = FeatureID, y = df)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
    ggplot2::ggsave(file.path(fit, "model0", "feature_vars0_df.pdf"), g, width = 8, height = 6, limitsize = F)

    # compute EB feature prior
    priors$DT.feature <- DT.feature.vars[, squeeze_var_func(fit)$value(v, df)]

    if(!is.null(control$peptide.model)) {
      message("[", paste0(Sys.time(), "]  computing empirical Bayes peptide variance prior..."))
      DT.peptide.vars <- peptide_vars(fit, stage = "0", as.data.table = T)

      # plot peptide variance fit
      DT.peptide.vars[, nPeptide := .N, by = ProteinID]
      g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = PeptideID, y = rhat)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model0", "peptide_vars0_rhat.pdf"), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(DT.peptide.vars, ggplot2::aes(x = PeptideID, y = v)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model0", "peptide_vars0_v.pdf"), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(DT.peptide.vars, ggplot2::aes(x = PeptideID, y = df)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model0", "peptide_vars0_df.pdf"), g, width = 8, height = 6, limitsize = F)

      # compute EB peptide prior
      priors$DT.peptide <- DT.peptide.vars[, squeeze_var_func(fit)$value(v, df)]
    }

    if(!is.null(control$assay.model)) {
      message("[", paste0(Sys.time(), "]  computing empirical Bayes assay variance prior..."))
      DT.assay.vars <- assay_vars(fit, stage = "0", as.data.table = T)

      #DT <- assay_vars(fit, 80, stage = "0", summary = F, as.data.table = T)
      #DT.fits <- DT[, dist_invchisq_mcmc(chainID, mcmcID, value, control=list(trace=1, REPORT=1)), by = .(ProteinID, AssayID)]
      #plot_fits(DT, DT.fits, by = "AssayID", ci = c(0.001, 0.999), trans = log2, inv.trans = function(x) 2^x)
      #DT <- protein_quants(fit,  1720, stage = "0", summary = F, as.data.table = T, norm.func.key = NULL)
      #DT <- protein_quants(fit,  1, stage = "0", summary = F, as.data.table = T, norm.func.key = NULL)
      #system.time(DT.fits <- DT[, dist_lst_mcmc(chainID, mcmcID, value), by = .(ProteinID, AssayID)])
      #system.time(DT.fits <- DT[, dist_lst_mcmc(chainID, mcmcID, value, control=list(trace=1, REPORT=1)), by = .(ProteinID, AssayID)])
      #plot_fits(DT, DT.fits, by = "AssayID", c(0.01, 0.99))

      # plot assay variance fit
      DT.assay.vars[, nPeptide := .N, by = ProteinID]
      g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = ProteinID, y = rhat)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model0", "assay_vars0_rhat.pdf"), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(DT.assay.vars, ggplot2::aes(x = ProteinID, y = v)) + ggplot2::geom_point(ggplot2::aes(colour = factor(AssayID)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model0", "assay_vars0_v.pdf"), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(DT.assay.vars, ggplot2::aes(x = ProteinID, y = df)) + ggplot2::geom_point(ggplot2::aes(colour = factor(AssayID)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model0", "assay_vars0_df.pdf"), g, width = 8, height = 6, limitsize = F)

      # compute EB assay prior(s)
      if (control$assay.model == "single") {
        priors$DT.assay <- DT.assay.vars[, squeeze_var_func(fit)$value(v, df)]
      } else {
        priors$DT.assay <- DT.assay.vars[, squeeze_var_func(fit)$value(v, df), by = AssayID]
      }
    }

    # save priors
    saveRDS(priors, file = file.path(fit, "model0", "priors.rds"))

    # plot EB

    # prepare data tables
    DT.priors <- priors$DT.feature[, .(Effect = "Feature", v0, df0, v, df)]
    if (!is.null(priors$DT.peptide)) DT.priors <- rbind(DT.priors, priors$DT.peptide[, .(Effect = "Peptide", v0, df0, v, df)])
    if (!is.null(priors$DT.assay)) {
      if (is.null(priors$DT.assay$AssayID)) {
        DT.priors <- rbind(DT.priors, priors$DT.assay[, .(Effect = "Assay", v0, df0, v, df)])
      } else {
        DT.assay <- merge(design(fit, as.data.table = T)[, .(AssayID, Assay, Sample)], priors$DT.assay, by = "AssayID")
        setorder(DT.assay, Assay)
        DT.assay[, Effect := paste0("Assay ", Assay, " (", Sample, ") ")]
        DT.assay[, Assay := NULL]
        DT.assay[, Sample := NULL]
        DT.assay[, AssayID := NULL]
        DT.priors <- rbind(DT.priors, DT.assay)
      }
    }
    DT.priors[, Effect := factor(Effect, levels = unique(Effect))]

    DT.vars <- DT.feature.vars[, .(Effect = "Feature", v, df)]
    if (!is.null(priors$DT.peptide)) DT.vars <- rbind(DT.vars, DT.peptide.vars[, .(Effect = "Peptide", v, df)])
    if (!is.null(priors$DT.assay)) {
      if (is.null(priors$DT.assay$AssayID)) {
        DT.vars <- rbind(DT.vars, DT.assay.vars[, .(Effect = "Assay", v, df)])
      } else {
        DT.assay <- merge(design(fit, as.data.table = T)[, .(AssayID, Assay, Sample)], DT.assay.vars[, .(AssayID, v, df)], by = "AssayID")
        setorder(DT.assay, Assay)
        DT.assay[, Effect := paste0("Assay ", Assay, " (", Sample, ") ")]
        DT.assay[, Assay := NULL]
        DT.assay[, Sample := NULL]
        DT.assay[, AssayID := NULL]
        DT.vars <- rbind(DT.vars, DT.assay)
      }
    }
    DT.vars[, Effect := factor(Effect, levels = unique(Effect))]

    # plot eb prior fits
    g <- plot_fits(DT.vars, DT.priors,  ci = c(0.05, 0.95), by = "Effect")
    suppressWarnings(ggplot2::ggsave(file.path(fit, "output", "protein_quant_qc_vars.pdf"), g, width = 8, height = 0.5 + 1.5 * nrow(DT.priors), limitsize = F))

    # plot eb prior stdevs
    g <- plot_fits(DT.vars, DT.priors, ci = c(0.01, 0.95), by = "Effect", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
    suppressWarnings(ggplot2::ggsave(file.path(fit, "output", "protein_quant_qc_stdevs.pdf"), g, width = 8, height = 0.5 + 1.5 * nrow(DT.priors), limitsize = F))
  }

  stop_parallel()
}


#' process_model (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
#' @import ggfortify
#' @export
process_model <- function(
  fit,
  chain = 1
) {
  control <- control(fit)
  set.seed(control$model.seed + chain)
  start_parallel(fit)

  # EXECUTE MODEL
  priors <- readRDS(file.path(fit, "model0", "priors.rds"))
  execute_model(fit, chain, priors)

  if (length(list.files(file.path(fit, "model"), "\\.finished$")) == control$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT stage=full"))

    # load parameters
    DT.design <- design(fit, as.data.table = T)
    DT.proteins <- proteins(fit, as.data.table = T)
    DT.peptides <- peptides(fit, as.data.table = T)
    DT.features <- features(fit, as.data.table = T)

    # WRITE MODEL OUTPUT

    # timings
    DT.timings <- timings(fit, as.data.table = T)
    DT.timings <- data.table::dcast(DT.timings, ProteinID ~ chainID, value.var = "elapsed")
    DT.timings <- merge(DT.proteins[, .(ProteinID, nPeptide, nFeature, nMeasure, pred = timing)], DT.timings, by = "ProteinID")
    fwrite(DT.timings, file.path(fit, "output", "protein_timings.csv"))
    rm(DT.timings)

    # feature vars
    if (control$feature.vars == T) {
      message("[", paste0(Sys.time(), "]  computing feature variances..."))
      DT.feature.vars <- feature_vars(fit, as.data.table = T)
      if (control$feature.model == "independent") {
        DT.feature.vars <- merge(DT.features, DT.feature.vars, by = "FeatureID")
      }
      setcolorder(DT.feature.vars, c("ProteinID", "PeptideID"))
      fwrite(DT.feature.vars, file.path(fit, "output", "feature_log2vars.csv"))
      rm(DT.feature.vars)
    }

    # peptide vars
    if(control$peptide.vars == T && !is.null(control$peptide.model)) {
      message("[", paste0(Sys.time(), "]  computing peptide variances..."))
      DT.peptide.vars <- peptide_vars(fit, as.data.table = T)
      if (control$peptide.model == "independent") {
        DT.peptide.vars <- merge(DT.peptides, DT.peptide.vars, by = "PeptideID")
      }
      setcolorder(DT.peptide.vars, "ProteinID")
      fwrite(DT.peptide.vars, file.path(fit, "output", "peptide_log2vars.csv"))
      rm(DT.peptide.vars)
    }

    # peptide deviations
    if(control$peptide.deviations == T && !is.null(control$assay.model) && control$assay.model == "independent") {
      message("[", paste0(Sys.time(), "]   peptide deviations..."))
      for (k in 1:length(control$dist.mean.func)) {
        set.seed(control$model.seed + chain)
        DT.peptide.deviations <- peptide_deviations(fit, summary = k, as.data.table = T)

        # save csv summary
        DT.peptide.deviations <- merge(DT.design[, .(AssayID, Sample, Assay)], DT.peptide.deviations, by = "AssayID")
        DT.peptide.deviations[, ProteinIDPeptideID := paste(ProteinID, PeptideID, sep = "_")]
        DT.peptide.deviations[, SampleAssay := paste(Sample, Assay, sep = "_")]
        setcolorder(DT.peptide.deviations, c("ProteinIDPeptideID", "SampleAssay"))
        DT.peptide.deviations <- dcast(DT.peptide.deviations, ProteinIDPeptideID ~ SampleAssay, drop = FALSE, value.var = colnames(DT.peptide.deviations)[8:ncol(DT.peptide.deviations)])
        DT.peptide.deviations[, ProteinID := as.integer(sub("^([0-9]+)_[0-9]+$", "\\1", ProteinIDPeptideID))]
        DT.peptide.deviations[, PeptideID := as.integer(sub("^[0-9]+_([0-9]+)$", "\\1", ProteinIDPeptideID))]
        DT.peptide.deviations[, ProteinIDPeptideID := NULL]
        DT.peptide.deviations <- merge(DT.peptides[, .(PeptideID, Peptide, nFeature, nMeasure)], DT.peptide.deviations, by = "PeptideID")
        setcolorder(DT.peptide.deviations, "ProteinID")
        fwrite(DT.peptide.deviations, file.path(fit, "output", paste0("peptide_log2deviations", ifelse(k == 1, "", dist_mean_func(fit, k)$index), ".csv")))
        rm(DT.peptide.deviations)
      }
    }

    # protein quants
    message("[", paste0(Sys.time(), "]  computing protein quants..."))
    n <- max(length(control$ref.assays), length(control$norm.func), length(control$dist.mean.func))
    for (k in 1:n) {
      set.seed(control$model.seed + chain)
      message("[", paste0(Sys.time(), "]   ref.assays=", ref_assays(fit, k)$key, " norm.func=", norm_func(fit, k)$key, "..."))
      DT.protein.quants <- protein_quants(fit, norm.func.key = k, ref.assays.key = k, summary = k, as.data.table = T)

      # save csv summary
      DT.protein.quants <- merge(DT.design[, .(AssayID, Sample, Assay)], DT.protein.quants, by = "AssayID")
      DT.protein.quants[, SampleAssay := paste(Sample, Assay, sep = "_")]
      setcolorder(DT.protein.quants, "SampleAssay")
      DT.protein.quants <- dcast(DT.protein.quants, ProteinID ~ SampleAssay, drop = FALSE, value.var = colnames(DT.protein.quants)[6:ncol(DT.protein.quants)])
      DT.protein.quants <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], DT.protein.quants, by = "ProteinID")
      fwrite(DT.protein.quants, file.path(fit, "output", paste0("protein_log2quants__", ref_assays(fit, k)$key, "__", norm_func(fit, k)$key, ".csv")))
      rm(DT.protein.quants)

      # plots
      DT.protein.quants <- protein_quants(fit, norm.func.key = k, ref.assays.key = k, summary = F, as.data.table = T)

      g <- plot_exposures(fit, DT.protein.quants)
      ggplot2::ggsave(file.path(fit, "output", paste0("assay_exposures__", ref_assays(fit, k)$key, "__", norm_func(fit, k)$key, ".pdf")),
                      g, width = 8, height = 0.5 + 0.75 * nrow(DT.design), limitsize = F)

      g <- plot_pca(fit, DT.protein.quants)
      ggplot2::ggsave(file.path(fit, "output", paste0("assay_pca__", ref_assays(fit, k)$key, "__", norm_func(fit, k)$key, ".pdf")),
                      g, width = 12, height = 12, limitsize = F)

      rm(DT.protein.quants)
    }

    # differential expression analysis
    message("[", paste0(Sys.time(), "]  differential expression analysis..."))
    dir.create(file.path(fit, "model", "protein.de"), showWarnings = F)
    n <- max(length(control$ref.assays), length(control$norm.func), length(control$dist.mean.func), length(control$dea.func))
    for (k in 1:n) {
      set.seed(control$model.seed + chain)
      if (!is.null(dea_func(fit, k))) {
        message("[", paste0(Sys.time(), "]   ref.assays=", ref_assays(fit, k)$key, " norm.func=", norm_func(fit, k)$key, " dea.func=", dea_func(fit, k)$key, "..."))

        # precompute DE
        protein_de(fit, dea.func.key = k, dist.mean.func.key = k, norm.func.key = k, ref.assays.key = k, as.data.table = T)
      }
    }

    # fdr
    message("[", paste0(Sys.time(), "]  false discovery rate control..."))
    dir.create(file.path(fit, "model", "protein.fdr"), showWarnings = F)
    n <- max(length(control$ref.assays), length(control$norm.func), length(control$dist.mean.func), length(control$dea.func), length(control$fdr.func))
    for (k in 1:n) {
      set.seed(control$model.seed + chain)
      message("[", paste0(Sys.time(), "]   ref.assays=", ref_assays(fit, k)$key, " norm.func=", norm_func(fit, k)$key, " dea.func=", dea_func(fit, k)$key, " fdr.func=", fdr_func(fit, k)$key, "..."))

      # run fdr
      DTs.fdr <- split(protein_fdr(fit, fdr.func.key = k, dea.func.key = k, dist.mean.func.key = k, norm.func.key = k, ref.assays.key = k, as.data.table = T), by = c("Model", "Effect"), drop = T)

      for (i in 1:length(DTs.fdr)) {
        # save pretty version
        DT.fdr <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], DTs.fdr[[i]], by = "ProteinID")
        setorder(DT.fdr, qvalue)
        fwrite(DT.fdr[, !c("Model", "Covariate")], file.path(fit, "output", paste0("protein_log2DE_", dea_func(fit, k)$key, "_", names(DTs.fdr)[i], ".csv")))
        g <- bayesprot::plot_fdr(DT.fdr, 1.0)
        ggplot2::ggsave(file.path(fit, "output", paste0("protein_log2DE_fdr_", dea_func(fit, k)$key, "_", names(DTs.fdr)[i], ".pdf")), g, width = 8, height = 8)
      }
    }
  }

  stop_parallel()
}


#' process_plots (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @import data.table
#' @import doRNG
#' @import ggplot2
#' @export
process_plots <- function(
  fit,
  chain
) {
  control <- control(fit)
  message(paste0("[", Sys.time(), "] PLOTS set=", chain, "/", control$model.nchain))

  # create subdirs
  dir.create(file.path(fit, "plots", "features"), showWarnings = F)
  dir.create(file.path(fit, "plots", "peptides"), showWarnings = F)

  DT.proteins <- proteins(fit, as.data.table = T)
  DT.design <- design(fit, as.data.table = T)
  DT.protein.quants <- protein_quants(fit, as.data.table = T)
  pids <- levels(DT.protein.quants$ProteinID)
  pids <- pids[seq(chain, length(pids), control$model.nchain)]

  # start cluster and reproducible seed
  pb <- txtProgressBar(max = length(pids), style = 3)
  dfll <- foreach(pid = pids, .packages = c("bayesprot", "data.table", "ggplot2"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    plt.features <- plot_features(fit, proteinID = pid)
    ggplot2::ggsave(file.path(fit, "plots", "features", paste0(pid, ".pdf")), plt.features$g, width = plt.features$width, height = plt.features$height, limitsize = F)

    plt.peptides <- plot_peptides(fit, proteinID = pid)
    ggplot2::ggsave(file.path(fit, "plots", "peptides", paste0(pid, ".pdf")), plt.peptides$g, width = plt.peptides$width, height = plt.peptides$height, limitsize = F)
  }
  setTxtProgressBar(pb, length(pids))
}

