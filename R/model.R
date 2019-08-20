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
  # EXECUTE MODEL
  execute_model(fit, chain)

  control <- control(fit)
  if (length(list.files(file.path(fit, "model1"), "\\.finished$")) == control$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT stage=1/2"))

    # load parameters
    DT.proteins <- proteins(fit, as.data.table = T)
    DT.design <- design(fit, as.data.table = T)

    # FIT INVERSE GAMMA DISTRIBUTIONS TO VARIANCES
    message("[", paste0(Sys.time(), "]  computing peptide & feature priors..."))
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

    if(!is.null(control$peptide.model)) {
      # fit peptide posterior medians
      DT.peptide.vars <- peptide_vars1_ln(fit)
      peptide.fit <- fitdistrplus::fitdist(1.0 / DT.peptide.vars$est, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
      peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
      peptide.V <- as.numeric((2.0 * 1.0 / peptide.fit$estimate["scale"]) / peptide.nu)
    } else {
      peptide.nu <- NULL
      peptide.V <- NULL
    }

    # fit feature posterior medians
    DT.feature.vars <- feature_vars1_ln(fit)
    feature.fit <- fitdistrplus::fitdist(1.0 / DT.feature.vars$est, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
    feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
    feature.V <- as.numeric((2.0 * 1.0 / feature.fit$estimate["scale"]) / feature.nu)

    # save output
    saveRDS(list(
      peptide.V = peptide.V, peptide.nu = peptide.nu,
      feature.V = feature.V, feature.nu = feature.nu
    ), file = file.path(fit, "model1", "priors.rds"))


    # PLOT
    prior.vars.meta <- function(x) {
      m = median(x)
      data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
    }

    if(!is.null(control$peptide.model)) {
      DT.prior.vars.meta <- rbind(
        data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.meta(sqrt(est)))]),
        data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(sqrt(est)))])
      )
    } else {
      DT.prior.vars.meta <- data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(sqrt(est)))])
    }
    DT.prior.vars.meta[, Type := factor(Type, levels = unique(Type))]

    if(!is.null(control$peptide.model)) {
      DT.prior.fit.meta <- rbind(
        data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.meta(sqrt(peptide.V)))]),
        data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(sqrt(feature.V)))])
      )
    } else {
      DT.prior.fit.meta <- data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(sqrt(feature.V)))])
    }
    DT.prior.fit.meta[, Type := factor(Type, levels = unique(Type))]

    prior.vars.density <- function(x) {
      DT <- as.data.table(density(log(x), n = 4096)[c("x","y")])
      DT[, x := exp(x)]
      DT
    }

    if(!is.null(control$peptide.model)) {
      DT.prior.vars.density <- rbind(
        data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.density(sqrt(est)))]),
        data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.density(sqrt(est)))])
      )
    } else {
      DT.prior.vars.density <- data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.density(sqrt(est)))])
    }
    DT.prior.vars.density[, Type := factor(Type, levels = unique(Type))]

    if(!is.null(control$peptide.model)) {
      DT.prior.fit.density <- rbind(
        data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.density(sqrt(MCMCglmm::rIW(peptide.V * diag(1), peptide.nu, n = 1000000))))]),
        data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.density(sqrt(MCMCglmm::rIW(feature.V * diag(1), feature.nu, n = 1000000))))])
      )
    } else {
      DT.prior.fit.density <- data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.density(sqrt(MCMCglmm::rIW(feature.V * diag(1), feature.nu, n = 1000000))))])
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
    g <- g + ggplot2::coord_cartesian(xlim = c(min(DT.prior.vars.density$x) / 1.1, max(DT.prior.vars.density$x) * 1.1), ylim = c(0, max(DT.prior.fit.density$y) * 1.35))
    g <- g + ggplot2::xlab(expression('Ln Standard Deviation'))
    g <- g + ggplot2::ylab("Probability Density")
    g <- g + ggplot2::facet_grid(Type ~ .)
    g <- g + ggplot2::scale_x_log10(labels = fmt_signif(1), expand = c(0, 0))
    g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
    g <- g + ggplot2::geom_ribbon(data = DT.prior.vars.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
    g <- g + ggplot2::geom_line(data = DT.prior.vars.density, ggplot2::aes(x = x,y = y), size = 1/2)
    g <- g + ggplot2::geom_line(data = DT.prior.fit.density, ggplot2::aes(x = x,y = y), size = 1/2, colour = "red")
    g <- g + ggplot2::geom_vline(data = DT.prior.fit.meta, ggplot2::aes(xintercept = median), size = 1/2, colour = "red")
    g <- g + ggplot2::geom_text(data = DT.prior.fit.meta, ggplot2::aes(x = median, label = fc), y = max(DT.prior.fit.density$y) * 1.25, hjust = 0, vjust = 1, size = 3, colour = "red")
    ggplot2::ggsave(file.path(fit, "output", "peptide_feature_priors.pdf"), g, width = 8, height = 0.5 + 2 * length(levels(DT.prior.vars.density$Type)), limitsize = F)
  }
}


# for model1 only
peptide_vars1_ln <- function(fit) {
  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")

  DT.peptide.vars <- rbindlist(lapply(c(
    list.files(file.path(fit, "model1", "peptide.stdevs.1"), paste0("^", chains[1], "\\..*fst$"), full.names = T)
  ), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
    DT <- DT[, value := (value * log(2))^2]
    if(control(fit)$peptide.model != "single") {
      DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID, prior)]
    } else {
      DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, prior)]
    }
  }))

  return(DT.peptide.vars)
}


# for model1 only
feature_vars1_ln <- function(fit) {
  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")

  DT.feature.vars <- rbindlist(lapply(c(
    list.files(file.path(fit, "model1", "feature.stdevs.1"), paste0("^", chains[1], "\\..*fst$"), full.names = T)
  ), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
    DT <- DT[, value := (value * log(2))^2]
    if(control(fit)$feature.model != "single") {
      DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, FeatureID, prior)]
    } else {
      DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, prior)]
    }
  }))

  return(DT.feature.vars)
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
  # EXECUTE MODEL
  execute_model(fit, chain, readRDS(file.path(fit, "model1", "priors.rds")))

  control <- control(fit)
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

    # feature stdevs
    DT.feature.stdevs <- feature_stdevs(fit, as.data.table = T)
    if (control$feature.model != "single") {
      DT.feature.stdevs <- merge(DT.features, DT.feature.stdevs, by = "FeatureID")
    }
    DT.feature.stdevs <- merge(DT.peptides[, .(PeptideID, Peptide)], DT.feature.stdevs, by = "PeptideID")
    DT.feature.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.feature.stdevs, by = "ProteinID")
    fwrite(DT.feature.stdevs, file.path(fit, "output", "feature_log2SDs.csv"))
    rm(DT.feature.stdevs)

    if(!is.null(control$peptide.model)) {
      # peptide stdevs
      DT.peptide.stdevs <- peptide_stdevs(fit, as.data.table = T)
      if (control$peptide.model != "single") {
        DT.peptide.stdevs <- merge(DT.peptides, DT.peptide.stdevs, by = "PeptideID")
      }
      DT.peptide.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.peptide.stdevs, by = "ProteinID")
      fwrite(DT.peptide.stdevs, file.path(fit, "output", "peptide_log2SDs.csv"))
      rm(DT.peptide.stdevs)

            # peptide deviations
      DT.peptide.deviations <- peptide_deviations(fit, as.data.table = T)
      DT.peptide.deviations <- merge(DT.peptide.deviations, DT.design[, .(SampleID, Sample)], by = "SampleID")
      DT.peptide.deviations[, SampleSE := paste0("SE:", Sample)]
      DT.peptide.deviations[, Sample := paste0("est:", Sample)]
      DT.peptide.deviations.ses <- data.table::dcast(DT.peptide.deviations, ProteinID + PeptideID + prior ~ SampleSE, value.var = "SE")
      DT.peptide.deviations <- data.table::dcast(DT.peptide.deviations, ProteinID + PeptideID + prior ~ Sample, value.var = "est")
      DT.peptide.deviations <- cbind(DT.peptide.deviations, DT.peptide.deviations.ses[, 4:ncol(DT.peptide.deviations.ses), with = F])
      setcolorder(DT.peptide.deviations, c("ProteinID", "PeptideID", "prior", paste0(c("est:", "SE:"), rep(levels(DT.design$Sample), each = 2))))
      DT.peptide.deviations <- merge(DT.peptides, DT.peptide.deviations, by = "PeptideID")
      DT.peptide.deviations <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.peptide.deviations, by = "ProteinID")
      fwrite(DT.peptide.deviations, file.path(fit, "output", "peptide_log2deviations.csv"))
      rm(DT.peptide.deviations)
    }

    # (unnormalised) protein quants
    save_summary <- function(DT.protein.quants.summary, name = "") {
      DT.protein.quants.summary.out <- merge(DT.design[, .(AssayID, Sample, Assay)], DT.protein.quants.summary, by = "AssayID")
      headings <- interaction(DT.protein.quants.summary.out$Sample, DT.protein.quants.summary.out$Assay, drop = T, sep = ":")
      DT.protein.quants.summary.out[, Assay := headings]
      DT.protein.quants.summary.out[, AssaySE := headings]
      DT.protein.quants.summary.out[, AssayNPeptide := headings]
      DT.protein.quants.summary.out[, AssayNFeature := headings]
      DT.protein.quants.summary.out[, Sample := NULL]
      levels(DT.protein.quants.summary.out$Assay) <- paste0("est:", levels(DT.protein.quants.summary.out$Assay))
      levels(DT.protein.quants.summary.out$AssaySE) <- paste0("SE:", levels(DT.protein.quants.summary.out$AssaySE))
      levels(DT.protein.quants.summary.out$AssayNPeptide) <- paste0("nPeptide:", levels(DT.protein.quants.summary.out$AssayNPeptide))
      levels(DT.protein.quants.summary.out$AssayNFeature) <- paste0("nFeature:", levels(DT.protein.quants.summary.out$AssayNFeature))
      DT.protein.quants.summary.out.SEs <- data.table::dcast(DT.protein.quants.summary.out, ProteinID ~ AssaySE, value.var = "SE")
      DT.protein.quants.summary.out.nPeptides <- data.table::dcast(DT.protein.quants.summary.out, ProteinID ~ AssayNPeptide, value.var = "nPeptide")
      DT.protein.quants.summary.out.nFeatures <- data.table::dcast(DT.protein.quants.summary.out, ProteinID ~ AssayNFeature, value.var = "nFeature")
      DT.protein.quants.summary.out <- data.table::dcast(DT.protein.quants.summary.out, ProteinID ~ Assay, value.var = "est")
      DT.protein.quants.summary.out <- cbind(DT.protein.quants.summary.out, DT.protein.quants.summary.out.SEs[, 2:ncol(DT.protein.quants.summary.out.SEs), with = F])
      DT.protein.quants.summary.out <- cbind(DT.protein.quants.summary.out, DT.protein.quants.summary.out.nPeptides[, 2:ncol(DT.protein.quants.summary.out.nPeptides), with = F])
      DT.protein.quants.summary.out <- cbind(DT.protein.quants.summary.out, DT.protein.quants.summary.out.nFeatures[, 2:ncol(DT.protein.quants.summary.out.nFeatures), with = F])
      setcolorder(DT.protein.quants.summary.out, c("ProteinID", paste0(c("est:", "SE:", "nPeptide:", "nFeature:"), rep(levels(headings), each = 4))))
      DT.protein.quants.summary.out <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], DT.protein.quants.summary.out, by = "ProteinID")
      fwrite(DT.protein.quants.summary.out, file.path(fit, "output", paste0("protein_log2quants", name, ".csv")))
    }
    DT.protein.quants.summary <- protein_quants(fit, data.exposures = NULL, as.data.table = T)
    save_summary(DT.protein.quants.summary, "_unnormalised")

    # normalisation
    if (!is.null(control$norm.func)) {
        message("[", paste0(Sys.time(), "]  calculating exposures and normalising quants..."))

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
          DT.de <- merge(DTs.de[[j]], DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure, prior)], by = "ProteinID", sort = F)
          setcolorder(DT.de, c("ProteinID", "Protein", "ProteinInfo", "nPeptide", "nFeature", "nMeasure", "prior"))
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
}


execute_model <- function(
  fit,
  chain,
  priors = NULL
) {
  stage = ifelse(is.null(priors), "1", "2")
  control <- control(fit)
  message(paste0("[", Sys.time(), "] MODEL stage=", stage, "/2 chain=", chain, "/", control$model.nchain))

  # load metadata
  path.output = file.path(fit, paste0("model", stage))
  DT.design <- design(fit, as.data.table = T)
  DT.proteins <- proteins(fit, as.data.table = T)
  nitt <- control$model.nwarmup + (control$model.nsample * control$model.thin) / control$model.nchain
  chainID <- formatC(chain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

  # create subdirs
  dir.create(file.path(path.output, paste0("protein.quants.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("peptide.deviations.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("peptide.stdevs.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("feature.stdevs.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("summaries.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("timings.", stage)), showWarnings = F)

  # stage 1: nPeptide > control$peptide.prior && nFeature > control$feature.prior, stage 2: others
  DT.proteins <- DT.proteins[prior == (stage == 2)]

  if (nrow(DT.proteins) > 0) {
    # start cluster and reproducible seed
    cl <- parallel::makeCluster(control$nthread)
    doSNOW::registerDoSNOW(cl)
    RNGkind("L'Ecuyer-CMRG")
    parallel::clusterSetRNGStream(cl, control$model.seed * control$model.nchain + chain - 1)

    # go...
    rbindlistlist <- function(...) {
      input <- list(...)
      for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
      input[[1]]
    }
    message(paste0("[", Sys.time(), "]  modelling nprotein=", nrow(DT.proteins), "/", nlevels(DT.proteins$ProteinID), " nitt=", nitt, "/", nitt * control$model.nchain, "..."))
    pb <- txtProgressBar(max = sum(DT.proteins$timing), style = 3)
    progress <- function(n, tag) setTxtProgressBar(pb, getTxtProgressBar(pb) + DT.proteins$timing[tag])
    output <- foreach(i = 1:nrow(DT.proteins), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.snow = list(progress = progress)) %dopar% {
      # prepare DT for MCMCglmm
      DT <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T, from = DT.proteins[i, from], to = DT.proteins[i, to])
      DT <- droplevels(DT)
      DT[, Count := round(Count)]
      if (!is.null(DT$Count1)) DT[, Count1 := round(Count1)]

      # calculate how many real (non-imputed) peptides and features back each Assay and Sample
      DT.assay.n <- DT[, .(nFeature = sum(!is.na(RawCount))), by = .(AssayID, PeptideID)]
      DT.assay.n <- DT.assay.n[, .(nPeptide = sum(nFeature > 0), nFeature = sum(nFeature)), by = AssayID]
      DT.sample.n <- DT[, .(nFeature = sum(!is.na(RawCount))), by = .(SampleID, PeptideID)]
      DT.sample.n <- DT.sample.n[, .(nPeptide = sum(nFeature > 0), nFeature = sum(nFeature)), by = SampleID]
      DT[, RawCount := NULL]

      # create co-occurence matrix of which assays are present in each feature
      DT[, BaselineID := AssayID]
      mat.tmp <- merge(DT, DT, by = "FeatureID", allow.cartesian = T)
      mat.tmp <- table(mat.tmp[, list(AssayID.x, AssayID.y)])
      # matrix multiplication distributes assay relationships
      mat.tmp <- mat.tmp %*% mat.tmp
      # IS THIS NECESSARY? # ignore columns that are not in ref.assays
      # mat.tmp[!DT.design[AssayID %in% colnames(mat.tmp), ref],] <- NA
      # baseline is first non-zero occurence for each assay
      DT[, BaselineID := colnames(mat.tmp)[apply(mat.tmp != 0, 2, which.max)][AssayID]]
      rm(mat.tmp)
      DT[, QuantID := as.character(interaction(DT$AssayID, DT$BaselineID, lex.order = T, drop = T))]
      # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
      DT[AssayID == BaselineID, QuantID := "."]
      DT[, QuantID := factor(QuantID)]

      nQ <- length(levels(DT$QuantID))

      output <- list()
      if (nQ > 1) {
        setcolorder(DT, c("PeptideID", "FeatureID", "AssayID", "SampleID", "QuantID"))
        nT <- length(levels(DT$PeptideID))
        nF <- length(levels(DT$FeatureID))

        # fixed effects
        fixed <- as.formula(paste(ifelse(is.null(DT$Count1), "Count", "c(Count, Count1)"), "~ ", ifelse(nF == 1, "", "FeatureID-1 +"), " QuantID"))

        # random effect
        if (is.null(control$peptide.model)) {
          random <- NULL
          prior.random <- NULL
        } else if (control$peptide.model == "single") {
          random <- as.formula("~PeptideID:SampleID")
          prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
          if (!is.null(priors) && nlevels(DT$PeptideID) <= control$peptide.prior) {
            prior.random <- list(V = priors$peptide.V, nu = priors$peptide.nu)
          } else {
            prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
          }
        } else {
          random <- as.formula(paste0("~", ifelse(nT == 1, "PeptideID", "idh(PeptideID)"), ":SampleID"))
          if (!is.null(priors) && nlevels(DT$PeptideID) <= control$peptide.prior) {
            prior.random <- list(V = priors$peptide.V * diag(nT), nu = priors$peptide.nu)
          } else {
            prior.random <- list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))
          }
        }

        # residual
        if (control$feature.model == "single") {
          rcov <- as.formula("~FeatureID:AssayID")
          if (!is.null(priors) && nlevels(DT$FeatureID) <= control$feature.prior) {
            prior.rcov <- list(V = priors$feature.V, nu = priors$feature.nu)
          } else {
            prior.rcov <- list(V = 1, nu = 0.02)
          }
        } else {
          rcov <- as.formula(paste0("~", ifelse(nF == 1, "FeatureID:AssayID", "idh(FeatureID):AssayID")))
          if (!is.null(priors) && nlevels(DT$FeatureID) <= control$feature.prior) {
            # need to figure out here which should have prior and which not
            prior.rcov <- list(V = priors$feature.V * diag(nF), nu = priors$feature.nu)
          } else {
            prior.rcov <- list(V = diag(nF), nu = 0.02)
          }
        }

        # family
        if (control$error.model == "lognormal") {
          DT$Count <- log(DT$Count)
          if(is.null(DT$Count1)) {
            family <- "gaussian"
          } else {
            DT[, Count1 := log(Count1)]
            family <- "cengaussian"
          }
        } else {
          if(is.null(DT$Count1)) {
            family <- "poisson"
          } else {
            family <- "cenpoisson"
          }
        }

        # prior
        prior <- list(R = prior.rcov)
        if (!is.null(prior.random)) {
          prior$G = list(G1 = prior.random)
        }

        # run model
        output$DT.summaries <- as.character(Sys.time())
        output$DT.timings <- system.time(model <- (MCMCglmm::MCMCglmm(
          fixed, random, rcov, family, data = DT, prior = prior,
          nitt = nitt, burnin = control$model.nwarmup, thin = control$model.thin, pr = T, verbose = F
        )))
        output$DT.timings <- data.table(ProteinID = DT[1, ProteinID], chainID = factor(chainID), as.data.table(t(as.matrix(output$DT.timings))))
        options(max.print = 99999)
        output$DT.summaries <- data.table(ProteinID = DT[1, ProteinID], chainID = factor(chainID), Summary = paste(c(output$DT.summaries, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

        if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
          stop("Some contrasts were dropped unexpectedly")
        }

        # extract protein quants
        output$DT.protein.quants <- as.data.table(model$Sol[, grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
        output$DT.protein.quants[, mcmcID := factor(formatC(1:nrow(output$DT.protein.quants), width = ceiling(log10(nrow(output$DT.protein.quants))) + 1, format = "d", flag = "0"))]
        output$DT.protein.quants <- melt(output$DT.protein.quants, variable.name = "BaselineID", id.vars = "mcmcID")
        output$DT.protein.quants[, ProteinID := DT[1, ProteinID]]
        output$DT.protein.quants[, AssayID := sub("^QuantID([0-9]+)\\.[0-9]+$", "\\1", BaselineID)]
        output$DT.protein.quants[, BaselineID := sub("^QuantID[0-9]+\\.([0-9]+)$", "\\1", BaselineID)]

        # add zeros for baseline assays
        output$DT.protein.quants <- rbind(output$DT.protein.quants, output$DT.protein.quants[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
        output$DT.protein.quants[, chainID := factor(chainID)]
        output$DT.protein.quants[, AssayID := factor(AssayID, levels = levels(DT.design$AssayID))]
        output$DT.protein.quants[, BaselineID := factor(BaselineID, levels = levels(DT.design$AssayID))]
        output$DT.protein.quants[, value := value / log(2)]
        output$DT.protein.quants[, prior := stage == 2]

        # merge with DT.assay.n
        output$DT.protein.quants <- merge(output$DT.protein.quants, DT.assay.n, by = "AssayID")
        setcolorder(output$DT.protein.quants, c("ProteinID", "prior", "AssayID", "BaselineID", "nPeptide", "nFeature", "chainID", "mcmcID"))

        # extract peptide deviations
        if (!is.null(control$peptide.model)) {
          if (nT == 1 || control$peptide.model == "single") {
            output$DT.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID:SampleID\\.[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
            output$DT.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.deviations), width = ceiling(log10(nrow(output$DT.peptide.deviations))) + 1, format = "d", flag = "0"))]
            output$DT.peptide.deviations <- melt(output$DT.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
            output$DT.peptide.deviations[, SampleID := factor(sub("^PeptideID:SampleID\\.[0-9]+\\.([0-9]+)$", "\\1", PeptideID))]
            output$DT.peptide.deviations[, PeptideID := factor(sub("^PeptideID:SampleID\\.([0-9]+)\\.[0-9]+$", "\\1", PeptideID))]
          } else {
            output$DT.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.SampleID\\.[0-9]+$", colnames(model$Sol)), drop = F])
            output$DT.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.deviations), width = ceiling(log10(nrow(output$DT.peptide.deviations))) + 1, format = "d", flag = "0"))]
            output$DT.peptide.deviations <- melt(output$DT.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
            output$DT.peptide.deviations[, SampleID := factor(sub("^PeptideID[0-9]+\\.SampleID\\.([0-9]+)$", "\\1", PeptideID))]
            output$DT.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.SampleID\\.([0-9]+)$", "\\1", PeptideID))]
          }
          output$DT.peptide.deviations[, ProteinID := DT[1, ProteinID]]
          output$DT.peptide.deviations[, chainID := factor(chainID)]
          output$DT.peptide.deviations[, prior := nlevels(DT$PeptideID) <= control$peptide.prior]
          output$DT.peptide.deviations[, value := value / log(2)]

          # merge with DT.n.real
          output$DT.peptide.deviations <- merge(output$DT.peptide.deviations, DT.sample.n, by = "SampleID")
          setcolorder(output$DT.peptide.deviations, c("ProteinID", "prior", "SampleID", "nPeptide", "nFeature", "PeptideID", "chainID", "mcmcID"))
        }

        model$Sol <- NULL

        # extract peptide stdevs
        if (!is.null(control$peptide.model)) {
          if (control$peptide.model == "single" || nT == 1) {
            output$DT.peptide.stdevs <- as.data.table(model$VCV[, "PeptideID:SampleID", drop = F])
            setnames(output$DT.peptide.stdevs, "PeptideID:SampleID", "value")
            output$DT.peptide.stdevs[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.stdevs), width = ceiling(log10(nrow(output$DT.peptide.stdevs))) + 1, format = "d", flag = "0"))]
            if (control$peptide.model != "single") {
              output$DT.peptide.stdevs[, PeptideID := factor(levels(DT$PeptideID))]
            }
          } else {
            output$DT.peptide.stdevs <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.SampleID$", colnames(model$VCV)), drop = F])
            output$DT.peptide.stdevs[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.stdevs), width = ceiling(log10(nrow(output$DT.peptide.stdevs))) + 1, format = "d", flag = "0"))]
            output$DT.peptide.stdevs <- melt(output$DT.peptide.stdevs, variable.name = "PeptideID", id.vars = "mcmcID")
            output$DT.peptide.stdevs[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.SampleID$", "\\1", PeptideID))]
            setcolorder(output$DT.peptide.stdevs, c("value", "mcmcID"))
          }
          output$DT.peptide.stdevs[, ProteinID := DT[1, ProteinID]]
          output$DT.peptide.stdevs[, chainID := factor(chainID)]
          output$DT.peptide.stdevs[, prior := nlevels(DT$PeptideID) <= control$peptide.prior]
          output$DT.peptide.stdevs[, value := sqrt(value) / log(2)]

          if(control$peptide.model != "single") {
            setcolorder(output$DT.peptide.stdevs, c("ProteinID", "prior", "PeptideID", "chainID", "mcmcID"))
          } else {
            setcolorder(output$DT.peptide.stdevs, c("ProteinID", "prior", "chainID", "mcmcID"))
          }
        }

        # extract feature variances
        if (control$feature.model == "single" || nF == 1) {
          output$DT.feature.stdevs <- as.data.table(model$VCV[, "FeatureID:AssayID", drop = F])
          setnames(output$DT.feature.stdevs, "FeatureID:AssayID", "value")
          output$DT.feature.stdevs[, mcmcID := factor(formatC(1:nrow(output$DT.feature.stdevs), width = ceiling(log10(nrow(output$DT.feature.stdevs))) + 1, format = "d", flag = "0"))]
          if (control$feature.model != "single") {
            output$DT.feature.stdevs[, FeatureID := factor(levels(DT$FeatureID))]
            output$DT.feature.stdevs <- merge(output$DT.feature.stdevs, unique(DT[, .(FeatureID, PeptideID, ProteinID)]), by = "FeatureID")
          } else {
            output$DT.feature.stdevs[, ProteinID := DT[1, ProteinID]]
          }
        } else {
          output$DT.feature.stdevs <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
          output$DT.feature.stdevs[, mcmcID := factor(formatC(1:nrow(output$DT.feature.stdevs), width = ceiling(log10(nrow(output$DT.feature.stdevs))) + 1, format = "d", flag = "0"))]
          output$DT.feature.stdevs <- melt(output$DT.feature.stdevs, variable.name = "FeatureID", id.vars = "mcmcID")
          output$DT.feature.stdevs[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.AssayID$", "\\1", FeatureID))]
          setcolorder(output$DT.feature.stdevs, c("value", "mcmcID"))
          output$DT.feature.stdevs <- merge(output$DT.feature.stdevs, unique(DT[, .(FeatureID, PeptideID, ProteinID)]), by = "FeatureID")
        }
        output$DT.feature.stdevs[, chainID := factor(chainID)]
        output$DT.feature.stdevs[, prior := nlevels(DT$FeatureID) <= control$feature.prior]
        output$DT.feature.stdevs[, value := sqrt(value) / log(2)]

        if (control$feature.model != "single") {
          setcolorder(output$DT.feature.stdevs, c("ProteinID", "prior", "PeptideID", "FeatureID", "chainID", "mcmcID"))
        } else {
          setcolorder(output$DT.feature.stdevs, c("ProteinID", "prior", "chainID", "mcmcID"))
        }

        # write out if large enough
        if (object.size(output$DT.protein.quants) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("protein.quants.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))
          fst::write.fst(output$DT.protein.quants, file.path(fit, filename))
          if (chain == 1) output$DT.protein.quants.index <- data.table(ProteinID = DT.proteins[i, ProteinID], file = filename, from = 1, to = nrow(output$DT.protein.quants))
          output$DT.protein.quants <- data.table()
        }

        if (!is.null(output$DT.peptide.deviations) && object.size(output$DT.peptide.deviations) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("peptide.deviations.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))
          fst::write.fst(output$DT.peptide.deviations, file.path(fit, filename))
          if (chain == 1) output$DT.peptide.deviations.index <- data.table(ProteinID = DT.proteins[i, ProteinID], file = filename, from = 1, to = nrow(output$DT.peptide.deviations))
          output$DT.peptide.deviations <- data.table()
        }

        if (!is.null(output$DT.peptide.stdevs) && object.size(output$DT.peptide.stdevs) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("peptide.stdevs.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))
          fst::write.fst(output$DT.peptide.stdevs, file.path(fit, filename))
          if (chain == 1) output$DT.peptide.stdevs.index <- data.table(ProteinID = DT.proteins[i, ProteinID], file = filename, from = 1, to = nrow(output$DT.peptide.stdevs))
          output$DT.peptide.stdevs <- data.table()
        }

        if (object.size(output$DT.feature.stdevs) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("feature.stdevs.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))
          fst::write.fst(output$DT.feature.stdevs, file.path(fit, filename))
          if (chain == 1) output$DT.feature.stdevs.index <- data.table(ProteinID = DT.proteins[i, ProteinID], file = filename, from = 1, to = nrow(output$DT.feature.stdevs))
          output$DT.feature.stdevs <- data.table()
        }
      }

      output
    }
    setTxtProgressBar(pb, sum(DT.proteins$timing))
    close(pb)


    # write out concatenation of smaller output
    message(paste0("[", Sys.time(), "]  writing output..."))
    output$DT.summaries[, ProteinID := factor(as.character(ProteinID))]
    setorder(output$DT.summaries, ProteinID, chainID)
    fst::write.fst(output$DT.summaries, file.path(path.output, file.path(paste0("summaries.", stage), paste0(chainID, ".fst"))))
    output$DT.summaries <- NULL

    setorder(output$DT.timings, ProteinID, chainID)
    output$DT.timings[, ProteinID := factor(as.character(ProteinID))]
    fst::write.fst(output$DT.timings, file.path(path.output, file.path(paste0("timings.", stage), paste0(chainID, ".fst"))))
    output$DT.timings <- NULL

    # write out peptide deviations
    if (!is.null(output$DT.peptide.deviations) && nrow(output$DT.peptide.deviations) > 0) {
      output$DT.peptide.deviations[, ProteinID := factor(as.character(ProteinID))]
      output$DT.peptide.deviations[, PeptideID := factor(as.character(PeptideID))]
      output$DT.peptide.deviations[, SampleID := factor(as.character(SampleID))]
      setorder(output$DT.peptide.deviations, ProteinID, PeptideID, SampleID, chainID, mcmcID)
      filename <- file.path(paste0("model", stage), paste0("peptide.deviations.", stage), paste0(chainID, ".fst"))
      fst::write.fst(output$DT.peptide.deviations, file.path(fit, filename))

      if (chain == 1) {
        output$DT.peptide.deviations.index <- rbind(output$DT.peptide.deviations.index, output$DT.peptide.deviations[, .(
          ProteinID = unique(ProteinID),
          file = filename,
          from = .I[!duplicated(ProteinID)],
          to = .I[rev(!duplicated(rev(ProteinID)))]
        )])
      }
      output$DT.peptide.deviations <- NULL
    }
    # write out peptide deviations index
    if (chain == 1) {
      output$DT.peptide.deviations.index <- merge(proteins(fit, as.data.table = T)[, .(ProteinID, Protein)], output$DT.peptide.deviations.index, by = "ProteinID")
      fst::write.fst(output$DT.peptide.deviations.index, file.path(path.output, paste0("peptide.deviations.", stage, ".index.fst")))
    }

    # write out peptide stdevs
    if (!is.null(output$DT.peptide.stdevs) && nrow(output$DT.peptide.stdevs) > 0) {
      output$DT.peptide.stdevs[, ProteinID := factor(as.character(ProteinID))]

      if(control$peptide.model != "single") {
        output$DT.peptide.stdevs[, PeptideID := factor(as.character(PeptideID))]
        setorder(output$DT.peptide.stdevs, ProteinID, PeptideID, chainID, mcmcID)
      } else {
        setorder(output$DT.peptide.stdevs, ProteinID, chainID, mcmcID)
      }

      filename <- file.path(paste0("model", stage), paste0("peptide.stdevs.", stage), paste0(chainID, ".fst"))
      fst::write.fst(output$DT.peptide.stdevs, file.path(fit, filename))

      if (chain == 1) {
        output$DT.peptide.stdevs.index <- rbind(output$DT.peptide.stdevs.index, output$DT.peptide.stdevs[, .(
          ProteinID = unique(ProteinID),
          file = filename,
          from = .I[!duplicated(ProteinID)],
          to = .I[rev(!duplicated(rev(ProteinID)))]
        )])
      }

      output$DT.peptide.stdevs <- NULL
    }
    # write out peptide stdevs index
    if (chain == 1) {
      output$DT.peptide.stdevs.index <- merge(proteins(fit, as.data.table = T)[, .(ProteinID, Protein)], output$DT.peptide.stdevs.index, by = "ProteinID")
      fst::write.fst(output$DT.peptide.stdevs.index, file.path(path.output, paste0("peptide.stdevs.", stage, ".index.fst")))
    }

    # write out feature stdevs
    if (!is.null(output$DT.feature.stdevs) && nrow(output$DT.feature.stdevs) > 0) {
      output$DT.feature.stdevs[, ProteinID := factor(as.character(ProteinID))]

      if(control$feature.model != "single") {
        output$DT.feature.stdevs[, PeptideID := factor(as.character(PeptideID))]
        output$DT.feature.stdevs[, FeatureID := factor(as.character(FeatureID))]
        setorder(output$DT.feature.stdevs, ProteinID, PeptideID, FeatureID, chainID, mcmcID)
      } else {
        setorder(output$DT.feature.stdevs, ProteinID, chainID, mcmcID)
      }

      filename <- file.path(paste0("model", stage), paste0("feature.stdevs.", stage), paste0(chainID, ".fst"))
      fst::write.fst(output$DT.feature.stdevs, file.path(fit, filename))

      if (chain == 1) {
        output$DT.feature.stdevs.index <- rbind(output$DT.feature.stdevs.index, output$DT.feature.stdevs[, .(
          ProteinID = unique(ProteinID),
          file = filename,
          from = .I[!duplicated(ProteinID)],
          to = .I[rev(!duplicated(rev(ProteinID)))]
        )])
      }

      output$DT.feature.stdevs <- NULL
    }
    # write out feature stdevs index
    if (chain == 1) {
      output$DT.feature.stdevs.index <- merge(proteins(fit, as.data.table = T)[, .(ProteinID, Protein)], output$DT.feature.stdevs.index, by = "ProteinID")
      fst::write.fst(output$DT.feature.stdevs.index, file.path(path.output, paste0("feature.stdevs.", stage, ".index.fst")))
    }

    # write out protein quants
    if (!is.null(output$DT.protein.quants) && nrow(output$DT.protein.quants) > 0) {
      output$DT.protein.quants[, ProteinID := factor(as.character(ProteinID))]
      output$DT.protein.quants[, AssayID := factor(as.character(AssayID))]
      setorder(output$DT.protein.quants, ProteinID, AssayID, chainID, mcmcID)
      filename <- file.path(paste0("model", stage), paste0("protein.quants.", stage), paste0(chainID, ".fst"))
      fst::write.fst(output$DT.protein.quants, file.path(fit, filename))

      if (chain == 1) {
        output$DT.protein.quants.index <- rbind(output$DT.protein.quants.index, output$DT.protein.quants[, .(
          ProteinID = unique(ProteinID),
          file = filename,
          from = .I[!duplicated(ProteinID)],
          to = .I[rev(!duplicated(rev(ProteinID)))]
        )])
      }

      output$DT.protein.quants <- NULL
    }
    # write out protein quants index
    if (chain == 1) {
      output$DT.protein.quants.index <- merge(proteins(fit, as.data.table = T)[, .(ProteinID, Protein)], output$DT.protein.quants.index, by = "ProteinID")
      fst::write.fst(output$DT.protein.quants.index, file.path(path.output, paste0("protein.quants.", stage, ".index.fst")))
    }

    # stop cluster
    parallel::stopCluster(cl)
  }

  write.table(data.frame(), file.path(path.output, paste0(chainID, ".finished")), col.names = F)
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
