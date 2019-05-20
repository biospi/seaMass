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
    message(paste0("[", Sys.time(), "] OUTPUT stage=1"))

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
      DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID, priors)]
    } else {
      DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, priors)]
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
      DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, FeatureID, priors)]
    } else {
      DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, priors)]
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
  execute_model(fit, chain, priors(fit))

  control <- control(fit)
  if (length(list.files(file.path(fit, "model2"), "\\.finished$")) == control$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT stage=2"))

    # load parameters
    DT.design <- design(fit, as.data.table = T)
    DT.proteins <- proteins(fit, as.data.table = T)
    DT.peptides <- peptides(fit, as.data.table = T)
    DT.features <- features(fit, as.data.table = T)
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

    # create nicer assay label for plots
    DT.design[, AssayInfo := paste0(Sample, ":", seq_len(.N)), by = Sample]
    DT.design[, AssayInfo := factor(sub(":1$", "", Sample))]

    # EXPOSURE-CORRECTED PROTEIN QUANTS
    message("[", paste0(Sys.time(), "]  computing exposure-corrected protein quants..."))

    # read in protein quants
    DT.protein.quants <- protein_quants(fit, normalise = F, summary = F, as.data.table = T)

    # assay exposures
    if (!all(!DT.proteins$norm)) {
      DT.assay.exposures <- DT.protein.quants[, .(value = median(value[ProteinID %in% DT.proteins[norm == T, ProteinID]])), by = .(AssayID, chainID, mcmcID)]
      fst::write.fst(DT.assay.exposures, file.path(fit, "model2", "assay.exposures.fst"))

      # plot
      assay.exposures.meta <- function(x) {
        m = median(x)
        data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
      }
      DT.assay.exposures.meta <- merge(DT.design[, .(AssayID, AssayInfo)], DT.assay.exposures[, as.list(assay.exposures.meta(value)), by = AssayID], by = "AssayID")

      assay.exposures.density <- function(x) {
        as.data.table(density(x, n = 4096)[c("x","y")])
      }
      DT.assay.exposures.density <- merge(DT.design[, .(AssayID, AssayInfo)], DT.assay.exposures[, as.list(assay.exposures.density(value)), by = AssayID], by = "AssayID")

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
      g <- g + ggplot2::facet_grid(AssayInfo ~ .)
      g <- g + ggplot2::xlab(expression('Log2 Ratio'))
      g <- g + ggplot2::ylab("Probability Density")
      g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
      g <- g + ggplot2::geom_ribbon(data = DT.assay.exposures.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
      g <- g + ggplot2::geom_line(data = DT.assay.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
      g <- g + ggplot2::geom_vline(data = DT.assay.exposures.meta, ggplot2::aes(xintercept = median), size = 1/2)
      g <- g + ggplot2::geom_text(data = DT.assay.exposures.meta, ggplot2::aes(x = median, label = fc), y = max(DT.assay.exposures.density$y) * 1.25, hjust = 0, vjust = 1, size = 3)
      ggplot2::ggsave(file.path(fit, "output", "assay_exposures.pdf"), g, width = 8, height = 0.5 + 0.75 * length(levels(DT.assay.exposures.density$AssayInfo)), limitsize = F)

      # apply exposures
      DT.protein.quants <- merge(DT.protein.quants, DT.assay.exposures[, .(AssayID, chainID, mcmcID, exposure = value)], by = c("AssayID", "chainID", "mcmcID"))
      DT.protein.quants[, value := value - exposure]
      DT.protein.quants[, exposure := NULL]
    }

    # compute and write out Rhat
    if (control$model.nchain > 1) {
      message("[", paste0(Sys.time(), "]  calculating Rhats..."))

      rhat <- function(DT) {
        chains <- split(DT[, .(chainID, value)], by = "chainID", keep.by = F, drop = T)
        chains <- coda::as.mcmc.list(lapply(names(chains), function(name) coda::as.mcmc(chains[[name]])))
        coda::gelman.diag(chains, autoburnin = F)$psrf[1]
      }
      DT.protein.quants.rhats <- DT.protein.quants[, .(rhat = rhat(.SD)), by = .(AssayID, ProteinID)]
      DT.protein.quants.rhats <- merge(DT.design[, .(AssayID, AssayInfo)], DT.protein.quants.rhats, by = "AssayID")
      DT.protein.quants.rhats <- dcast(DT.protein.quants.rhats, ProteinID ~ AssayInfo, value.var = "rhat")
      colnames(DT.protein.quants.rhats)[2:ncol(DT.protein.quants.rhats)] <- paste0("rhat:", colnames(DT.protein.quants.rhats)[2:ncol(DT.protein.quants.rhats)])
      DT.protein.quants.rhats <- merge(DT.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], DT.protein.quants.rhats, by = "ProteinID")
      fwrite(DT.protein.quants.rhats, file.path(fit, "output", "protein_log2quants_rhats.csv"))
      rm(DT.protein.quants.rhats)
    }

    # summarise MCMC samples
    DT.protein.quants <- DT.protein.quants[, .(priors = any(priors), est = median(value), SE = mad(value)), by = .(ProteinID, AssayID)]
    fst::write.fst(DT.protein.quants, file.path(fit, "model2", "protein.quants.fst"))

    # differential expression analysis
    if (!is.null(control$de.functions)) {
      dir.create(file.path(fit, "model2", "de"), showWarnings = F)
      deID <- formatC(1:length(control$de.functions), width = ceiling(log10(length(control$de.functions) + 1)) + 1, format = "d", flag = "0")

      for (i in 1:length(control$de.functions)) {
        message("[", paste0(Sys.time(), "]  performing differential expression analysis ", i, "/", length(control$de.functions), "..."))

        # run de
        out <- control$de.functions[[i]](fit)

        # save output
        if (class(out) == "bayesprot_de_metafor") {
          saveRDS(out,  file.path(fit, "model2", "de", paste0(deID[i], ".rds")))
          out <- protein_de(out, as.data.table = T)
        }

        if (is.data.frame(out)) {
          out <- list(setDT(out))
        }

        bID <- formatC(1:length(out), width = ceiling(log10(length(out) + 1)) + 1, format = "d", flag = "0")
        for (j in 1:length(out)) {
          if (nrow(out[[j]]) > 0) {
            fst::write.fst(out[[j]],  file.path(fit, "model2", "de", paste0(deID[i], ".", bID[j], ".fst")))

            # pretty version
            out[[j]] <- merge(out[[j]], DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], by = "ProteinID", sort = F)
            out[[j]] <- merge(out[[j]], DT.protein.quants[, .(priors = any(priors)), by = ProteinID], by = "ProteinID", sort = F)
            #out[[j]] <- merge(out[[j]], DT.real.ct, by = "ProteinID", sort = F)
            setcolorder(out[[j]], c("ProteinID", "Protein", "ProteinInfo", "nPeptide", "nFeature", "nMeasure", "priors"))
            fwrite(out[[j]], file.path(fit, "output", paste0("protein_log2DE_", names(control$de.functions)[i], "_", names(out)[j], ".csv")))
            g <- bayesprot::plot_fdr(out[[j]], 1.0)
            ggplot2::ggsave(file.path(fit, "output", paste0("protein_log2DE_fdr_", names(control$de.functions)[i], "_", names(out)[j], ".pdf")), g, width = 8, height = 8)
          }
        }
      }
    }

    message("[", paste0(Sys.time(), "]  writing model output..."))

    # write out
    DT.protein.quants <- merge(DT.design[, .(AssayID, Sample, Assay)], DT.protein.quants, by = "AssayID")

    headings <- interaction(DT.protein.quants$Sample, DT.protein.quants$Assay, drop = T, sep = ":")
    DT.protein.quants[, Assay := headings]
    DT.protein.quants[, AssaySE := headings]
    DT.protein.quants[, Sample := NULL]
    levels(DT.protein.quants$AssaySE) <- paste0("SE:", levels(DT.protein.quants$AssaySE))
    levels(DT.protein.quants$Assay) <- paste0("est:", levels(DT.protein.quants$Assay))
    DT.protein.quants.ses <- dcast(DT.protein.quants, ProteinID + priors ~ AssaySE, value.var = "SE")
    protein.quants.ses <- as.matrix(DT.protein.quants.ses[, 3:ncol(DT.protein.quants.ses), with = F]) # for pca
    DT.protein.quants <- dcast(DT.protein.quants, ProteinID + priors ~ Assay, value.var = "est")
    protein.quants <- as.matrix(DT.protein.quants[, 3:ncol(DT.protein.quants), with = F])  # for pca
    DT.protein.quants <- cbind(DT.protein.quants, DT.protein.quants.ses[, 3:ncol(DT.protein.quants.ses), with = F])
    setcolorder(DT.protein.quants, c("ProteinID", "priors", paste0(c("est:", "SE:"), rep(levels(headings), each = 2))))
    DT.protein.quants <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], DT.protein.quants, by = "ProteinID")
    fwrite(DT.protein.quants, file.path(fit, "output", "protein_log2quants.csv"))

    # write out pca plot
    pca.assays <- prcomp(
      t(protein.quants[complete.cases(protein.quants),]),
      center = T,
      scale = rowMeans(protein.quants.ses[complete.cases(protein.quants.ses),]^2)
    )
    rm(protein.quants.ses)
    rm(protein.quants)
    DT.pca.assays <- ggplot2::fortify(pca.assays)
    DT.pca.assays <- cbind(DT.pca.assays, DT.design)
    #if (any(as.character(DT.pca.assays$Assay) != as.character(DT.pca.assays$Sample))) {
    #  DT.pca.assays$Assay <- paste0(DT.pca.assays$Assay, "; ", DT.pca.assays$Sample)
    #}

    g <- ggplot2::autoplot(pca.assays, data = DT.pca.assays)
    g <- g + ggplot2::theme_bw()
    g <- g + ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", size = 1),
      panel.grid.major = ggplot2::element_line(size = 0.5),
      strip.background = ggplot2::element_blank(),
      aspect.ratio = 1.0
    )
    g <- g + ggplot2::coord_equal()
    if (is.null(DT.pca.assays$Condition)) {
      g <- g + ggplot2::geom_point()
      g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = AssayInfo), size = 3.0)
    } else {
      g <- g + ggplot2::geom_point(ggplot2::aes(colour = Condition))
      g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = AssayInfo, colour = Condition), size = 3.0)
    }
    if (!all(DT.design$isRef)) g <- g + ggtitle(paste("ref.assays =", paste(DT.design[isRef == T, AssayInfo], collapse = "; ")))
    ggplot2::ggsave(file.path(fit, "output", "assays_pca.pdf"), g, width=12, height=12, limitsize = F)


    # THE REST

    if(!is.null(control$peptide.model)) {
      # peptide deviations write out
      DT.peptide.deviations <- peptide_deviations(fit, as.data.table = T)
      fst::write.fst(DT.peptide.deviations, file.path(fit, "model2", "peptide.deviations.fst"))

      DT.peptide.deviations <- merge(DT.design[, .(Sample, SampleID)], DT.peptide.deviations, by = "SampleID")
      DT.peptide.deviations[, SampleSE := paste0("SE:", Sample)]
      DT.peptide.deviations[, Sample := paste0("est:", Sample)]
      DT.peptide.deviations.ses <- dcast(DT.peptide.deviations, ProteinID + PeptideID + priors ~ SampleSE, value.var = "SE")
      DT.peptide.deviations <- dcast(DT.peptide.deviations, ProteinID + PeptideID + priors ~ Sample, value.var = "est")
      DT.peptide.deviations <- cbind(DT.peptide.deviations, DT.peptide.deviations.ses[, 4:ncol(DT.peptide.deviations.ses), with = F])
      setcolorder(DT.peptide.deviations, c("ProteinID", "PeptideID", "priors", paste0(c("est:", "SE:"), rep(levels(DT.design$Sample), each = 2))))
      DT.peptide.deviations <- merge(DT.peptides, DT.peptide.deviations, by = "PeptideID")
      DT.peptide.deviations <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.peptide.deviations, by = "ProteinID")
      fwrite(DT.peptide.deviations, file.path(fit, "output", "peptide_log2deviations.csv"))
      rm(DT.peptide.deviations)

      # peptide stdevs in base 2
      DT.peptide.stdevs <- peptide_stdevs(fit, as.data.table = T)
      fst::write.fst(DT.peptide.stdevs, file.path(fit, "model2", "peptide.stdevs.fst"))

      if (control$peptide.model != "single") {
        DT.peptide.stdevs <- merge(DT.peptides, DT.peptide.stdevs, by = "PeptideID")
      }
      DT.peptide.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.peptide.stdevs, by = "ProteinID")
      fwrite(DT.peptide.stdevs, file.path(fit, "output", "peptide_log2SDs.csv"))
      rm(DT.peptide.stdevs)
    }

    # feature stdevs in base 2
    DT.feature.stdevs <- feature_stdevs(fit, as.data.table = T)
    fst::write.fst(DT.feature.stdevs, file.path(fit, "model2", "feature.stdevs.fst"))

    if (control$feature.model != "single") {
      DT.feature.stdevs <- merge(DT.features, DT.feature.stdevs, by = "FeatureID")
    }
    DT.feature.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.feature.stdevs, by = "ProteinID")
    fwrite(DT.feature.stdevs, file.path(fit, "output", "feature_log2SDs.csv"))
    rm(DT.feature.stdevs)

    # timings
    DT.timings <- timings(fit, as.data.table = T)
    DT.timings <- dcast(DT.timings, ProteinID ~ chainID, value.var = "elapsed")
    DT.timings <- merge(DT.proteins[, .(ProteinID, nPeptide, nFeature, nMeasure, pred = timing)], DT.timings, by = "ProteinID")
    fwrite(DT.timings, file.path(fit, "output", "protein_timings.csv"))
    rm(DT.timings)
  }
}


execute_model <- function(
  fit,
  chain,
  priors = NULL
) {
  stage = ifelse(is.null(priors), "1", "2")
  message(paste0("[", Sys.time(), "] MODEL stage=", stage, " chain=", chain))

  # load metadata
  path.output = file.path(fit, paste0("model", stage))
  control <- control(fit)
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

  # stage 1: nPeptide >= control$model.npeptide, stage 2: nPeptide < control$model.npeptide
  if (stage == 1) {
    DT.proteins <- DT.proteins[nPeptide >= control$model.npeptide]
  } else {
    DT.proteins <- DT.proteins[nPeptide < control$model.npeptide]
  }

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
    message(paste0("[", Sys.time(), "]  modelling nprotein=", nrow(DT.proteins), " nitt=", nitt, "..."))
    pb <- txtProgressBar(max = sum(DT.proteins$timing), style = 3)
    progress <- function(n, tag) setTxtProgressBar(pb, getTxtProgressBar(pb) + DT.proteins$timing[tag])
    output <- foreach(i = 1:nrow(DT.proteins), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.snow = list(progress = progress)) %dopar% {
      # prepare DT for MCMCglmm
      DT <- input(fit, as.data.table = T, from = DT.proteins[i, row], to = DT.proteins[i, row1])
      DT <- droplevels(DT)
      DT[, Count := round(Count)]
      if (!is.null(DT$Count1)) DT[, Count1 := round(Count1)]

      # create co-occurence matrix of which assays are present in each feature
      DT[, BaselineID := AssayID]
      mat.tmp <- merge(DT, DT, by = "FeatureID", allow.cartesian = T)
      mat.tmp <- table(mat.tmp[, list(AssayID.x, AssayID.y)])
      # matrix multiplication distributes assay relationships
      mat.tmp <- mat.tmp %*% mat.tmp
      # ignore columns that are not in ref.assays
      mat.tmp[!DT.design[AssayID %in% colnames(mat.tmp), ref],] <- NA
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
          if (!is.null(priors) && !is.null(control$peptide.prior)) {
            prior.random <- list(V = priors$peptide.V, nu = priors$peptide.nu)
          } else {
            prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
          }
        } else {
          random <- as.formula(paste0("~", ifelse(nT == 1, "PeptideID", "idh(PeptideID)"), ":SampleID"))
          if (!is.null(priors) && !is.null(control$peptide.prior)) {
            prior.random <- list(V = priors$peptide.V * diag(nT), nu = priors$peptide.nu)
          } else {
            prior.random <- list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))
          }
        }

        # residual
        if (control$feature.model == "single") {
          rcov <- as.formula("~FeatureID:AssayID")
          if (!is.null(priors) && !is.null(control$feature.prior)) {
            prior.rcov <- list(V = priors$feature.V, nu = priors$feature.nu)
          } else {
            prior.rcov <- list(V = 1, nu = 0.02)
          }
        } else {
          rcov <- as.formula(paste0("~", ifelse(nF == 1, "FeatureID:AssayID", "idh(FeatureID):AssayID")))
          if (!is.null(priors) && !is.null(control$feature.prior)) {
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
        output$DT.protein.quants[, BaselineID := NULL]

        # mean centre so that denominator is mean of reference assays, and log(2)
        refs <- DT.design[ref == T, AssayID]
        mean.refs <- function(AssayID, value) mean(value[AssayID %in% refs])
        output$DT.protein.quants <- output$DT.protein.quants[, .(AssayID, value = (value - mean.refs(AssayID, value)) / log(2)), by = .(ProteinID, chainID, mcmcID)]
        output$DT.protein.quants[, priors := !is.null(priors)]
        setcolorder(output$DT.protein.quants, c("ProteinID", "AssayID", "priors", "chainID", "mcmcID"))

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
          output$DT.peptide.deviations[, priors := !is.null(priors)]
          output$DT.peptide.deviations[, value := value / log(2)]
          setcolorder(output$DT.peptide.deviations, c("ProteinID", "PeptideID", "SampleID", "priors", "chainID", "mcmcID"))
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
          output$DT.peptide.stdevs[, priors := !is.null(priors)]
          output$DT.peptide.stdevs[, value := sqrt(value) / log(2)]
          if(control$peptide.model != "single") {
            setcolorder(output$DT.peptide.stdevs, c("ProteinID", "PeptideID", "priors", "chainID", "mcmcID"))
          } else {
            setcolorder(output$DT.peptide.stdevs, c("ProteinID", "priors", "chainID", "mcmcID"))
          }
        }

        # extract feature variances
        if (control$feature.model == "single" || nF == 1) {
          output$DT.feature.stdevs <- as.data.table(model$VCV[, "FeatureID:AssayID", drop = F])
          setnames(output$DT.feature.stdevs, "FeatureID:AssayID", "value")
          output$DT.feature.stdevs[, mcmcID := factor(formatC(1:nrow(output$DT.feature.stdevs), width = ceiling(log10(nrow(output$DT.feature.stdevs))) + 1, format = "d", flag = "0"))]
          if (control$feature.model != "single") {
            output$DT.feature.stdevs[, FeatureID := factor(levels(DT$FeatureID))]
          }
        } else {
          output$DT.feature.stdevs <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
          output$DT.feature.stdevs[, mcmcID := factor(formatC(1:nrow(output$DT.feature.stdevs), width = ceiling(log10(nrow(output$DT.feature.stdevs))) + 1, format = "d", flag = "0"))]
          output$DT.feature.stdevs <- melt(output$DT.feature.stdevs, variable.name = "FeatureID", id.vars = "mcmcID")
          output$DT.feature.stdevs[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.AssayID$", "\\1", FeatureID))]
          setcolorder(output$DT.feature.stdevs, c("value", "mcmcID"))
        }
        output$DT.feature.stdevs[, ProteinID := DT[1, ProteinID]]
        output$DT.feature.stdevs[, chainID := factor(chainID)]
        output$DT.feature.stdevs[, priors := !is.null(priors)]
        output$DT.feature.stdevs[, value := sqrt(value) / log(2)]
        if (control$feature.model != "single") {
          setcolorder(output$DT.feature.stdevs, c("ProteinID", "FeatureID", "priors", "chainID", "mcmcID"))
        } else {
          setcolorder(output$DT.feature.stdevs, c("ProteinID", "priors", "chainID", "mcmcID"))
        }

        # write out if large enough
        if (object.size(output$DT.protein.quants) > 2^18) {
          fst::write.fst(output$DT.protein.quants, file.path(path.output, file.path(paste0("protein.quants.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
          output$DT.protein.quants <- data.table()
        }

        if (!is.null(output$DT.peptide.deviations) && object.size(output$DT.peptide.deviations) > 2^18) {
          fst::write.fst(output$DT.peptide.deviations, file.path(path.output, file.path(paste0("peptide.deviations.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
          output$DT.peptide.deviations <- data.table()
        }

        if (!is.null(output$DT.peptide.stdevs) && object.size(output$DT.peptide.stdevs) > 2^18) {
          fst::write.fst(output$DT.peptide.stdevs, file.path(path.output, file.path(paste0("peptide.stdevs.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
          output$DT.peptide.stdevs <- data.table()
        }

        if (object.size(output$DT.feature.stdevs) > 2^18) {
          fst::write.fst(output$DT.feature.stdevs, file.path(path.output, file.path(paste0("feature.stdevs.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
          output$DT.feature.stdevs <- data.table()
        }
      }

      output
    }
    setTxtProgressBar(pb, sum(DT.proteins$timing))
    close(pb)


    # write out concatenation of smaller output
    message(paste0("[", Sys.time(), "]  writing output..."))

    fst::write.fst(output$DT.summaries, file.path(path.output, file.path(paste0("summaries.", stage), paste0(chainID, ".fst"))))
    output$DT.summaries <- NULL

    fst::write.fst(output$DT.timings, file.path(path.output, file.path(paste0("timings.", stage), paste0(chainID, ".fst"))))
    output$DT.timings <- NULL

    if (!is.null(output$DT.peptide.deviations) && nrow(output$DT.peptide.deviations) > 0) {
      fst::write.fst(output$DT.peptide.deviations, file.path(path.output, file.path(paste0("peptide.deviations.", stage), paste0(chainID, ".fst"))))
      output$DT.peptide.deviations <- NULL
    }

    if (!is.null(output$DT.peptide.stdevs) && nrow(output$DT.peptide.stdevs) > 0) {
      fst::write.fst(output$DT.peptide.stdevs, file.path(path.output, file.path(paste0("peptide.stdevs.", stage), paste0(chainID, ".fst"))))
      output$DT.peptide.stdevs <- NULL
    }

    if (!is.null(output$DT.feature.stdevs) && nrow(output$DT.feature.stdevs) > 0) {
      fst::write.fst(output$DT.feature.stdevs, file.path(path.output, file.path(paste0("feature.stdevs.", stage), paste0(chainID, ".fst"))))
      output$DT.feature.stdevs <- NULL
    }

    if (!is.null(output$DT.protein.quants) && nrow(output$DT.protein.quants) > 0) {
      fst::write.fst(output$DT.protein.quants, file.path(path.output, file.path(paste0("protein.quants.", stage), paste0(chainID, ".fst"))))
      output$DT.protein.quants <- NULL
    }

    # stop cluster
    parallel::stopCluster(cl)
  }

  write.table(data.frame(), file.path(path.output, paste0(chainID, ".finished")), col.names = F)
}


#' process_plots (internal)
#'
#' @param chain Number of chain to process.
#' @param fit bayesprot fit object.
#' @import data.table
#' @import foreach
#' @export
process_plots <- function(
  chain,
  fit = NULL
) {
}
