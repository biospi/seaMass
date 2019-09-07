#' process_model1 (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
#' @export
process_model1 <- function(
  fit,
  chain
) {
  # start cluster and reproducible seed
  control <- control(fit)
  cl <- parallel::makeCluster(control$nthread)
  doSNOW::registerDoSNOW(cl)
  RNGkind("L'Ecuyer-CMRG")
  parallel::clusterSetRNGStream(cl, control$model.seed * control$model.nchain + chain - 1)

  # EXECUTE MODEL
  execute_model(fit, chain)

  if (length(list.files(file.path(fit, "model1"), "\\.finished$")) == control$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT stage=1/2"))

    priors <- list()

    if(!is.null(control$peptide.model)) {
      message("[", paste0(Sys.time(), "]  computing empirical Bayes peptide variance prior..."))

      DT.peptide.vars <- peptide_vars(fit, as.data.table = T, stage = 1, parallel = T)

      # fit scaled F Distribution - Limma method fails with low DF as need to calculate mean. Hack - add some DF until it works.
      delta <- 0
      repeat{
        peptides.fit <- limma::fitFDist(DT.peptide.vars$tau2, DT.peptide.vars$nu + delta)
        if(!is.infinite(peptides.fit$df2)) {
          break
        } else if (delta > 2.0) {
          stop(paste0("delta 2.0 still failing, something odd is happening"))
        } else {
          delta <- delta + 0.1
          print(paste0("WARNING: fit failed, trying DF + ", delta))
        }
      }

      priors$peptide <- c(tau2 = peptides.fit$scale, nu = peptides.fit$df2)
    }

    # fit scaled chi squared distribution to each feature
    message("[", paste0(Sys.time(), "]  computing empirical Bayes feature variance prior..."))

    DT.feature.vars <- feature_vars(fit, as.data.table = T, stage = 1, parallel = T)

    # fit scaled F Distribution - Limma method fails with low DF as need to calculate mean. Hack - add some DF until it works.
    delta <- 0.0
    repeat{
      features.fit <- limma::fitFDist(DT.feature.vars$tau2, DT.feature.vars$nu + delta)
      if(!is.infinite(features.fit$df2)) {
        break
      } else if (delta > 2.0) {
        stop(paste0("delta 2.0 still failing, something odd is happening"))
      } else {
        delta <- delta + 0.1
        print(paste0("fit failed, trying nu + ", delta))
      }
    }
    priors$feature <- c(tau2 = features.fit$scale, nu = features.fit$df2)

    # save output
    saveRDS(priors, file = file.path(fit, "model1", "priors.rds"))

    # PLOT VARIANCE FITS
    DT.peptide.vars[, nPeptide := .N, by = ProteinID]
    if (control$peptide.model == "independent") {
      g <- ggplot2::ggplot(DT.peptide.vars, ggplot2::aes(x = as.numeric(PeptideID), y = log2(tau2))) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1)
      ggplot2::ggsave(file.path(fit, "model1", "peptide_vars_tau2.pdf"), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(DT.peptide.vars, ggplot2::aes(x = as.numeric(PeptideID), y = log2(nu))) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1)
      ggplot2::ggsave(file.path(fit, "model1", "peptide_vars_nu.pdf"), g, width = 8, height = 6, limitsize = F)
    } else {
      g <- ggplot2::ggplot(DT.peptide.vars, ggplot2::aes(x = as.numeric(ProteinID), y = log2(tau2))) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1)
      ggplot2::ggsave(file.path(fit, "model1", "peptide_vars_tau2.pdf"), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(DT.peptide.vars, ggplot2::aes(x = as.numeric(ProteinID), y = log2(nu))) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1)
      ggplot2::ggsave(file.path(fit, "model1", "peptide_vars_nu.pdf"), g, width = 8, height = 6, limitsize = F)
    }

    DT.feature.vars[, nPeptide := .N, by = ProteinID]
    if (control$feature.model == "independent") {
      g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = as.numeric(FeatureID), y = log2(tau2))) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1)
      ggplot2::ggsave(file.path(fit, "model1", "feature_vars_tau2.pdf"), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = as.numeric(FeatureID), y = log2(nu))) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1)
      ggplot2::ggsave(file.path(fit, "model1", "feature_vars_nu.pdf"), g, width = 8, height = 6, limitsize = F)
    } else {
      g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = as.numeric(ProteinID), y = log2(tau2))) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1)
      ggplot2::ggsave(file.path(fit, "model1", "feature_vars_tau2.pdf"), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(DT.feature.vars, ggplot2::aes(x = as.numeric(ProteinID), y = log2(nu))) + ggplot2::geom_point(ggplot2::aes(colour = log2(nPeptide)), size = 0.1)
      ggplot2::ggsave(file.path(fit, "model1", "feature_vars_nu.pdf"), g, width = 8, height = 6, limitsize = F)
    }

    # PLOT PRIOR FIT
    plot_priors <- function(priors, DT.peptide.vars, DT.feature.vars) {
      prior.vars.meta <- function(x) {
        m = median(x)
        data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
      }

      if(!is.null(control$peptide.model)) {
        DT.prior.vars.meta <- rbind(
          data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.meta(tau2))]),
          data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(tau2))])
        )
      } else {
        DT.prior.vars.meta <- data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(tau2))])
      }
      DT.prior.vars.meta[, Type := factor(Type, levels = unique(Type))]

      if(!is.null(control$peptide.model)) {
        DT.prior.fit.meta <- rbind(
          data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.meta(priors$peptide["tau2"]))]),
          data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(priors$feature["tau2"]))])
        )
      } else {
        DT.prior.fit.meta <- data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(priors$feature["tau2"]))])
      }
      DT.prior.fit.meta[, Type := factor(Type, levels = unique(Type))]

      prior.vars.density <- function(x) {
        DT <- as.data.table(density(log(x), n = 4096)[c("x","y")])
        DT[, x := exp(x)]
        DT
      }

      if(!is.null(control$peptide.model)) {
        DT.prior.vars.density <- rbind(
          data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.density(tau2))]),
          data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.density(tau2))])
        )
      } else {
        DT.prior.vars.density <- data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.density(tau2))])
      }
      DT.prior.vars.density[, Type := factor(Type, levels = unique(Type))]

      if(!is.null(control$peptide.model)) {
        DT.prior.fit.density <- rbind(
          data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.density(MCMCglmm::rIW(priors$peptide["tau2"] * diag(1), priors$peptide["nu"], n = 1000000)))]),
          data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.density(MCMCglmm::rIW(priors$feature["tau2"] * diag(1), priors$feature["nu"], n = 1000000)))])
        )
      } else {
        DT.prior.fit.density <- data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.density(MCMCglmm::rIW(priors$feature["tau2"] * diag(1), priors$feature["nu"], n = 1000000)))])
      }
      DT.prior.fit.density[, Type := factor(Type, levels = unique(Type))]

      fmt_signif <- function(signif = 2) {
        function(x) formatC(signif(x, digits = signif))
      }

      g <- ggplot2::ggplot(DT.prior.vars.density, ggplot2::aes(x = x, y = y))
      g <- g + ggplot2::theme_bw()
      g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                              panel.grid.major = ggplot2::element_line(size = 0.5),
                              axis.ticks = ggplot2::element_blank(),
                              axis.text.y = ggplot2::element_blank(),
                              plot.title = ggplot2::element_text(size = 10),
                              strip.background = ggplot2::element_blank(),
                              strip.text.y = ggplot2::element_text(angle = 0))
      #g <- g + ggplot2::coord_cartesian(xlim = c(min(DT.prior.vars.density$x) / 1.1, max(DT.prior.vars.density$x) * 1.1), ylim = c(0, max(DT.prior.fit.density$y) * 1.35))
      g <- g + ggplot2::coord_cartesian(xlim = c(0.00001, 100), ylim = c(0, max(DT.prior.fit.density$y) * 1.35))
      g <- g + ggplot2::xlab(expression('Log2 Variance'))
      g <- g + ggplot2::ylab("Probability Density")
      g <- g + ggplot2::facet_grid(Type ~ .)
      g <- g + ggplot2::scale_x_log10(labels = fmt_signif(1), expand = c(0, 0))
      g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
      g <- g + ggplot2::geom_ribbon(data = DT.prior.vars.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
      g <- g + ggplot2::geom_line(data = DT.prior.vars.density, ggplot2::aes(x = x,y = y), size = 1/2)
      g <- g + ggplot2::geom_line(data = DT.prior.fit.density, ggplot2::aes(x = x,y = y), size = 1/2, colour = "red")
      g <- g + ggplot2::geom_vline(data = DT.prior.fit.meta, ggplot2::aes(xintercept = median), size = 1/2, colour = "red")
      g <- g + ggplot2::geom_text(data = DT.prior.fit.meta, ggplot2::aes(x = median, label = fc), y = max(DT.prior.fit.density$y) * 1.25, hjust = 0, vjust = 1, size = 3, colour = "red")
      g
    }
    g <- plot_priors(priors, DT.peptide.vars, DT.feature.vars)
    ggplot2::ggsave(file.path(fit, "output", "peptide_feature_priors.pdf"), g, width = 8, height = 6, limitsize = F)
  }

  # stop cluster
  parallel::stopCluster(cl)
}


#' process_model2 (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
#' @import ggfortify
#' @export
process_model2 <- function(
  fit,
  chain
) {
  # start cluster and reproducible seed
  control <- control(fit)
  cl <- parallel::makeCluster(control$nthread)
  doSNOW::registerDoSNOW(cl)
  RNGkind("L'Ecuyer-CMRG")
  parallel::clusterSetRNGStream(cl, control$model.seed * control$model.nchain + chain - 1)

  # EXECUTE MODEL
  execute_model(fit, chain, readRDS(file.path(fit, "model1", "priors.rds")))

  if (length(list.files(file.path(fit, "model2"), "\\.finished$")) == control$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT stage=2/2"))

    # load parameters
    DT.design <- design(fit, as.data.table = T)
    DT.proteins <- proteins(fit, as.data.table = T)
    DT.peptides <- peptides(fit, as.data.table = T)
    DT.features <- features(fit, as.data.table = T)
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

    # WRITE MODEL OUTPUT
    message("[", paste0(Sys.time(), "]  writing model output to csv..."))

    # timings
    DT.timings <- timings(fit, as.data.table = T)
    DT.timings <- data.table::dcast(DT.timings, ProteinID ~ chainID, value.var = "elapsed")
    DT.timings <- merge(DT.proteins[, .(ProteinID, nPeptide, nFeature, nMeasure, pred = timing)], DT.timings, by = "ProteinID")
    fwrite(DT.timings, file.path(fit, "output", "protein_timings.csv"))
    rm(DT.timings)

    # feature vars
    message("[", paste0(Sys.time(), "]   feature vars..."))
    DT.feature.vars <- feature_vars(fit, as.data.table = T, parallel = T)
    if (control$feature.model == "independent") {
      DT.feature.vars <- merge(DT.features, DT.feature.vars, by = "FeatureID")
      DT.feature.vars <- merge(DT.peptides[, .(PeptideID, Peptide)], DT.feature.vars, by = "PeptideID")
    }
    DT.feature.vars <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.feature.vars, by = "ProteinID")
    fwrite(DT.feature.vars, file.path(fit, "output", "feature_log2SDs.csv"))
    rm(DT.feature.vars)

    if(!is.null(control$peptide.model)) {
      # peptide vars
      message("[", paste0(Sys.time(), "]   peptide vars..."))
      DT.peptide.vars <- peptide_vars(fit, as.data.table = T, parallel = T)
      if (control$peptide.model == "independent") {
        DT.peptide.vars <- merge(DT.peptides, DT.peptide.vars, by = "PeptideID")
      }
      DT.peptide.vars <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.peptide.vars, by = "ProteinID")
      fwrite(DT.peptide.vars, file.path(fit, "output", "peptide_log2SDs.csv"))
      rm(DT.peptide.vars)

      # peptide deviations
      message("[", paste0(Sys.time(), "]   peptide deviations..."))
      DT.peptide.deviations <- peptide_deviations(fit, as.data.table = T, parallel = T)
      # DT.peptide.deviations <- merge(DT.peptide.deviations, DT.design[, .(SampleID, Sample)], by = "SampleID")
      # DT.peptide.deviations[, SampleSE := paste0("SE:", Sample)]
      # DT.peptide.deviations[, Sample := paste0("est:", Sample)]
      # DT.peptide.deviations.ses <- data.table::dcast(DT.peptide.deviations, ProteinID + PeptideID ~ SampleSE, value.var = "SE")
      # DT.peptide.deviations <- data.table::dcast(DT.peptide.deviations, ProteinID + PeptideID ~ Sample, value.var = "est")
      # DT.peptide.deviations <- cbind(DT.peptide.deviations, DT.peptide.deviations.ses[, 3:ncol(DT.peptide.deviations.ses), with = F])
      # setcolorder(DT.peptide.deviations, c("ProteinID", "PeptideID", paste0(c("est:", "SE:"), rep(levels(DT.design$Sample), each = 2))))
      # DT.peptide.deviations <- merge(DT.peptides, DT.peptide.deviations, by = "PeptideID")
      # DT.peptide.deviations <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.peptide.deviations, by = "ProteinID")
      # fwrite(DT.peptide.deviations, file.path(fit, "output", "peptide_log2deviations.csv"))
      rm(DT.peptide.deviations)
    }

    # raw protein quants
    message("[", paste0(Sys.time(), "]   raw protein quants..."))
    # save_summary <- function(DT.protein.quants.summary, name = "") {
    #   DT.protein.quants.summary.out <- merge(DT.design[, .(AssayID, Sample, Assay)], DT.protein.quants.summary, by = "AssayID")
    #   headings <- interaction(DT.protein.quants.summary.out$Sample, DT.protein.quants.summary.out$Assay, drop = T, sep = ":")
    #   DT.protein.quants.summary.out[, Assay := headings]
    #   DT.protein.quants.summary.out[, AssaySE := headings]
    #   DT.protein.quants.summary.out[, AssayNPeptide := headings]
    #   DT.protein.quants.summary.out[, AssayNFeature := headings]
    #   DT.protein.quants.summary.out[, Sample := NULL]
    #   levels(DT.protein.quants.summary.out$Assay) <- paste0("est:", levels(DT.protein.quants.summary.out$Assay))
    #   levels(DT.protein.quants.summary.out$AssaySE) <- paste0("SE:", levels(DT.protein.quants.summary.out$AssaySE))
    #   levels(DT.protein.quants.summary.out$AssayNPeptide) <- paste0("nPeptide:", levels(DT.protein.quants.summary.out$AssayNPeptide))
    #   levels(DT.protein.quants.summary.out$AssayNFeature) <- paste0("nFeature:", levels(DT.protein.quants.summary.out$AssayNFeature))
    #   DT.protein.quants.summary.out.SEs <- data.table::dcast(DT.protein.quants.summary.out, ProteinID ~ AssaySE, value.var = "SE")
    #   DT.protein.quants.summary.out.nPeptides <- data.table::dcast(DT.protein.quants.summary.out, ProteinID ~ AssayNPeptide, value.var = "nPeptide")
    #   DT.protein.quants.summary.out.nFeatures <- data.table::dcast(DT.protein.quants.summary.out, ProteinID ~ AssayNFeature, value.var = "nFeature")
    #   DT.protein.quants.summary.out <- data.table::dcast(DT.protein.quants.summary.out, ProteinID ~ Assay, value.var = "est")
    #   DT.protein.quants.summary.out <- cbind(DT.protein.quants.summary.out, DT.protein.quants.summary.out.SEs[, 2:ncol(DT.protein.quants.summary.out.SEs), with = F])
    #   DT.protein.quants.summary.out <- cbind(DT.protein.quants.summary.out, DT.protein.quants.summary.out.nPeptides[, 2:ncol(DT.protein.quants.summary.out.nPeptides), with = F])
    #   DT.protein.quants.summary.out <- cbind(DT.protein.quants.summary.out, DT.protein.quants.summary.out.nFeatures[, 2:ncol(DT.protein.quants.summary.out.nFeatures), with = F])
    #   setcolorder(DT.protein.quants.summary.out, c("ProteinID", paste0(c("est:", "SE:", "nPeptide:", "nFeature:"), rep(levels(headings), each = 4))))
    #   DT.protein.quants.summary.out <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], DT.protein.quants.summary.out, by = "ProteinID")
    #   fwrite(DT.protein.quants.summary.out, file.path(fit, "output", paste0("protein_log2quants", name, ".csv")))
    # }
    DT.protein.quants.summary <- protein_quants(fit, data.exposures = NULL, as.data.table = T)
    #save_summary(DT.protein.quants.summary, "_unnormalised")

    # normalisation
    if (!is.null(control$norm.func)) {
        message("[", paste0(Sys.time(), "]   normalised protein quants..."))

        # calculate exposures and save normalised protein quants
        DT.assay.exposures <- control$norm.func(fit)
        setDT(DT.assay.exposures)
        fst::write.fst(DT.assay.exposures, file.path(fit, "model2", "assay.exposures.fst"))
        DT.protein.quants.summary <- protein_quants(fit, data.exposures = DT.assay.exposures, as.data.table = T)
        save_summary(DT.protein.quants.summary, "_normalised")

        # plot
        assay.exposures.meta <- function(x) {
          m = median(x)
          data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
        }
        DT.assay.exposures.meta <- merge(DT.design[, .(AssayID, Assay)], DT.assay.exposures[, as.list(assay.exposures.meta(value)), by = AssayID], by = "AssayID")

        assay.exposures.density <- function(x) {
          as.data.table(density(x, n = 4096)[c("x","y")])
        }
        DT.assay.exposures.density <- merge(DT.design[, .(AssayID, Assay)], DT.assay.exposures[, as.list(assay.exposures.density(value)), by = AssayID], by = "AssayID")

        x.max <- max(0.5, max(abs(DT.assay.exposures.density$x)))
        g <- ggplot2::ggplot(DT.assay.exposures.density, ggplot2::aes(x = x, y = y))
        g <- g + ggplot2::theme_bw()
        g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                                panel.grid.major = ggplot2::element_line(size = 0.5),
                                axis.ticks = ggplot2::element_blank(),
                                axis.text.y = ggplot2::element_blank(),
                                plot.title = ggplot2::element_text(size = 10),
                                strip.background = ggplot2::element_blank(),
                                strip.text.y = ggplot2::element_text(angle = 0))
        g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
        g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
        g <- g + ggplot2::coord_cartesian(xlim = c(-x.max, x.max) * 1.1, ylim = c(0, max(DT.assay.exposures.density$y) * 1.35))
        g <- g + ggplot2::facet_grid(Assay ~ .)
        g <- g + ggplot2::xlab(expression('Log2 Ratio'))
        g <- g + ggplot2::ylab("Probability Density")
        g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
        g <- g + ggplot2::geom_ribbon(data = DT.assay.exposures.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
        g <- g + ggplot2::geom_line(data = DT.assay.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
        g <- g + ggplot2::geom_vline(data = DT.assay.exposures.meta, ggplot2::aes(xintercept = median), size = 1/2)
        g <- g + ggplot2::geom_text(data = DT.assay.exposures.meta, ggplot2::aes(x = median, label = fc), y = max(DT.assay.exposures.density$y) * 1.25, hjust = 0, vjust = 1, size = 3)
        ggplot2::ggsave(file.path(fit, "output", paste0("assay_exposures.pdf")), g, width = 8, height = 0.5 + 0.75 * length(levels(DT.assay.exposures.density$Assay)), limitsize = F)
    }

    # compute and write out Rhat
    if (control$model.nchain > 1) {
      message("[", paste0(Sys.time(), "]  calculating Rhats..."))

      rhat <- function(DT) {
        chains <- split(DT[, .(chainID, value)], by = "chainID", keep.by = F, drop = T)
        chains <- coda::as.mcmc.list(lapply(names(chains), function(name) coda::as.mcmc(chains[[name]])))
        coda::gelman.diag(chains, autoburnin = F)$psrf[1]
      }
      DT.protein.quants.rhats <- protein_quants(fit, summary = F, as.data.table = T)[, .(rhat = rhat(.SD)), by = .(AssayID, ProteinID)]
      DT.protein.quants.rhats <- merge(DT.design[, .(AssayID, Assay)], DT.protein.quants.rhats, by = "AssayID")
      DT.protein.quants.rhats <- dcast(DT.protein.quants.rhats, ProteinID ~ Assay, value.var = "rhat")
      colnames(DT.protein.quants.rhats)[2:ncol(DT.protein.quants.rhats)] <- paste0("rhat:", colnames(DT.protein.quants.rhats)[2:ncol(DT.protein.quants.rhats)])
      DT.protein.quants.rhats <- merge(DT.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], DT.protein.quants.rhats, by = "ProteinID")
      fwrite(DT.protein.quants.rhats, file.path(fit, "output", "protein_log2quants_rhats.csv"))
      rm(DT.protein.quants.rhats)
    }

    # differential expression analysis
    dir.create(file.path(fit, "model2", "protein.de"), showWarnings = F)
    if (!is.null(control$dea.func)) {
      for (i in 1:length(control$dea.func)) {
        message("[", paste0(Sys.time(), "]  performing differential expression analysis ", i, "/", length(control$dea.func), "..."))

        # run de
        DTs.de <- control$dea.func[[i]](fit, output = names(control$dea.func)[i])
        fst::write.fst(DTs.de, file.path(fit, "model2", "protein.de", paste0(names(control$dea.func)[i], ".fst")))

        # fdr for pretty version
        DTs.de <- fdr_ash(DTs.de, as.data.table = T, nthread = control$nthread)

        DTs.de <- split(DTs.de, by = c("Model", "Covariate"), drop = T)
        for (j in 1:length(DTs.de)) {
          # save pretty version
          DT.de <- merge(DTs.de[[j]], DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], by = "ProteinID", sort = F)
          setcolorder(DT.de, c("ProteinID", "Protein", "ProteinInfo", "nPeptide", "nFeature", "nMeasure"))
          fwrite(DT.de[, !c("Model", "Covariate")], file.path(fit, "output", paste0("protein_log2DE_", names(control$dea.func)[i], "_", names(DTs.de)[j], ".csv")))
          g <- bayesprot::plot_fdr(DT.de, 1.0)
          ggplot2::ggsave(file.path(fit, "output", paste0("protein_log2DE_fdr_", names(control$dea.func)[i], "_", names(DTs.de)[j], ".pdf")), g, width = 8, height = 8)
        }
      }
    }

    message("[", paste0(Sys.time(), "]  computing PCA plot..."))

    # PCA PLOT
    g <- plot_pca(fit, DT.protein.quants.summary)
    ggplot2::ggsave(file.path(fit, "output", "assay_pca.pdf"), g, width = 12, height = 12, limitsize = F)
  }

  # stop cluster
  parallel::stopCluster(cl)
}


#' process_plots (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
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
  cl <- parallel::makeCluster(control$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(pids), style = 3)
  null <- foreach(pid = pids, .packages = c("bayesprot", "data.table", "ggplot2"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    plt.features <- plot_features(fit, proteinID = pid)
    ggplot2::ggsave(file.path(fit, "plots", "features", paste0(pid, ".pdf")), plt.features$g, width = plt.features$width, height = plt.features$height, limitsize = F)

    plt.peptides <- plot_peptides(fit, proteinID = pid)
    ggplot2::ggsave(file.path(fit, "plots", "peptides", paste0(pid, ".pdf")), plt.peptides$g, width = plt.peptides$width, height = plt.peptides$height, limitsize = F)
  }
  setTxtProgressBar(pb, length(pids))
  close(pb)
  parallel::stopCluster(cl)
}

