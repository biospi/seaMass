#' Mixed-effects univariate differential expression analysis with 'metafor'
#'
#' Default is to a one-way ANOVA on the column 'Condition' in 'data.design'.
#'
#' @import data.table
#' @import foreach
#' @export
dea_metafor <- function(fit, data.design = design(fit), mods = ~ Condition, random = ~ 1 | Sample, ...) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
  if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")
  if (!("control" %in% names(arguments))) control = list()
  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")

  # work out real n for each level and each protein
  DT.n <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T)[, .(ProteinID, AssayID, RawCount)]
  DT.n <- DT.n[complete.cases(DT.n)]
  DT.n <- unique(merge(DT.n, data.design, by = "AssayID")[, .(ProteinID, Sample)])
  DT.n <- DT.n[, .(N = length(Sample)), by = ProteinID]

  # prepare quants
  DTs <- protein_quants(fit, as.data.table = T)
  DTs <- merge(DTs, data.design, by = "AssayID")
  DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE)])])
  DTs[, SE := ifelse(SE < 0.001, 0.001, SE)] # won't work if SE is 0
  DTs[, BatchID := ProteinID]
  DTs[, ProteinID := as.character(ProteinID)]
  levels(DTs$BatchID) <- substr(levels(DTs$BatchID), 1, nchar(levels(DTs$BatchID)) - 1)
  DTs <- split(DTs, by = "BatchID", keep.by = F)

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(DTs), style = 3)
  output.all <- foreach(DT.chunk = iterators::iter(DTs), .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    lapply(split(DT.chunk, by = "ProteinID"), function(DT) {
      # input
      output <- list(input = as.data.frame(DT))

      # test
      output$fits <- vector("list", 1)
      for (i in 0:9) {
        control$sigma2.init = 0.025 + 0.1 * i
        output$log <- paste0(output$log, "[", Sys.time(), "]  attempt ", i + 1, "\n")
        try({
          output$log <- paste0(output$log, capture.output({
            output$fits[[1]] <- metafor::rma.mv(yi = est, V = SE^2, data = DT, control = control, test = "t", mods = mods, random = random, ...)
          }))
          output$log <- paste0(output$log, "[", Sys.time(), "]   succeeded\n")
          break
        })
      }

      # output
      if (!is.null(output$fits[[1]]) && !is.null(output$fits[[1]]$b) && length(output$fits[[1]]$b) > 0) {
        n.real = DT.n[ProteinID == DT$ProteinID[1], N]

        output$output <- data.table(
          Effect = rownames(output$fits[[1]]$b),
          n.test = length(unique(DT$Sample)),
          n.real = ifelse(length(n.real) > 0, n.real, 0),
          log2FC.lower = output$fits[[1]]$ci.lb,
          log2FC = output$fits[[1]]$b[, 1],
          log2FC.upper = output$fits[[1]]$ci.ub,
          t.value = output$fits[[1]]$zval,
          p.value = output$fits[[1]]$pval
        )
      } else {
        output$output <- NULL
      }

      output
    })
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)
  parallel::stopCluster(cl)

  output.all <- unlist(output.all, recursive = F)
  class(output.all) <- "bayesprot_de_metafor"
  attr(output.all, "bayesprot_fit") <- fit
  return(output.all)
}


#' Pair-wise mixed-effects univariate differential expression analysis with 'metafor'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import foreach
#' @export
dea_metafor_pairwise <- function(fit, data.design = design(fit), output = "dea_metafor_pairwise", as.data.table = F, mods = ~ Condition, random = ~ 1 | Sample, ...) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
  if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")
  if (!("control" %in% names(arguments))) control = list()
  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit))

  # prepare design
  DT.design <- as.data.table(data.design)
  if (!is.null(DT.design$nProtein)) DT.design[, nProtein := NULL]
  if (!is.null(DT.design$nPeptide)) DT.design[, nPeptide := NULL]
  if (!is.null(DT.design$nFeature)) DT.design[, nFeature := NULL]
  if (!is.null(DT.design$nMeasure)) DT.design[, nMeasure := NULL]
  if (!is.null(DT.design$SampleID)) DT.design[, SampleID := NULL]
  if (!is.null(DT.design$ref)) DT.design[, ref := NULL]

  # prepare quants
  DTs <- protein_quants(fit, as.data.table = T)
  DTs <- merge(DTs, DT.design, by = "AssayID")
  DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE, Condition)])])
  cts <- combn(levels(DTs$Condition), 2)

  # batch
  DTs[, BatchID := ProteinID]
  nbatch <- ceiling(nlevels(DTs$ProteinID) / 16)
  levels(DTs$BatchID) <- rep(formatC(1:nbatch, width = ceiling(log10(nbatch)) + 1, format = "d", flag = "0"), each = 16)[1:nlevels(DTs$ProteinID)]
  DTs <- split(DTs, by = "BatchID")

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(DTs), style = 3)
  DT.de <- foreach(DT.chunk = iterators::iter(DTs), .final = rbindlist, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {

    # process chunk
    output.chunk <- lapply(split(DT.chunk, by = "ProteinID", drop = T), function(DT) {
      DT[, ProteinID := NULL]
      DT[, BatchID := NULL]

      output.protein <- lapply(1:ncol(cts), function (j) {
        output.contrast <- list()

        # input
        output.contrast$DT.input <- DT[as.character(Condition) %in% cts[, j]]
        output.contrast$DT.input[, Condition := factor(as.character(Condition), levels = cts[, j])]
        output.contrast$DT.input <- droplevels(output.contrast$DT.input)
        setcolorder(output.contrast$DT.input, c("AssayID", "Assay", "Run", "Channel", "Sample", "est", "SE"))

        # metafor will moan if SE is 0
        output.contrast$DT.input[, SE := ifelse(SE < 0.001, 0.001, SE)]
        # metafor will moan if the ratio between largest and smallest sampling variance is >= 1e7
        output.contrast$DT.input[, SE := ifelse(SE^2 / min(SE^2) >= 1e7, sqrt(min(SE^2) * (1e7 - 1)), SE)]

        # fit
        if (length(unique(output.contrast$DT.input[Condition == cts[1, j], Sample])) >= 2 && length(unique(output.contrast$DT.input[Condition == cts[2, j], Sample])) >= 2) {
          for (i in 0:9) {
            control$sigma2.init = 0.025 + 0.1 * i
            try( {
              output.contrast$fit <- metafor::rma.mv(yi = est, V = SE^2, data = output.contrast$DT.input, control = control, test = "t", mods = mods, random = random)
              output.contrast$log <- paste0("[", Sys.time(), "] attempt ", i + 1, " succeeded\n")
              break
            })
          }
        } else {
          output.contrast$log <- paste0("[", Sys.time(), "] ignored as n < 2 for one or both conditions\n")
        }

        output.contrast
      })
      names(output.protein) <- sapply(1:length(output.protein), function(j) paste(cts[,j], collapse = "_v_"))

      output.protein
    })

    # save chunk
    saveRDS(output.chunk, file.path(fit, "model2", "de", paste0(output, ".", DT.chunk[1, BatchID], ".rds")))

    rbindlist(lapply(output.chunk, function (output.protein) {
      rbindlist(lapply(output.protein, function (output.contrast) {
        DT.out <- output.contrast$DT.input[, .(
          n1.test = length(unique(Sample[Condition == levels(Condition)[1]])),
          n2.test = length(unique(Sample[Condition == levels(Condition)[2]])),
          n1.real = sum(nPeptide[Condition == levels(Condition)[1]] > 0),
          n2.real = sum(nPeptide[Condition == levels(Condition)[2]] > 0)
        )]

        if (!is.null(output.contrast$fit) && !is.null(output.contrast$fit$b) && length(output.contrast$fit$b) > 0) {
          DT.res <- data.table(
            Covariate = rownames(output.contrast$fit$b),
            log2FC.lower = output.contrast$fit$ci.lb,
            log2FC = output.contrast$fit$b[, 1],
            log2FC.upper = output.contrast$fit$ci.ub,
            t.value = output.contrast$fit$zval,
            p.value = output.contrast$fit$pval
          )
          DT.out <- cbind(DT.out[rep(1, nrow(DT.res))], DT.res)
        } else {
          DT.out[, Covariate := NA]
          DT.out[, log2FC.lower := NA]
          DT.out[, log2FC := NA]
          DT.out[, log2FC.upper := NA]
          DT.out[, t.value := NA]
          DT.out[, p.value := NA]
        }
      }), idcol = "Model")
    }), idcol = "ProteinID")
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)
  parallel::stopCluster(cl)

  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Covariate := factor(Covariate, levels = unique(Covariate))]

  # highlights if any proteins are missing
  DT.proteins <- proteins(fit, as.data.table = T)
  DT.de[, ProteinID := factor(ProteinID, levels = levels(DT.proteins$ProteinID))]
  DT.de <- merge(DT.proteins[, .(ProteinID)], DT.de, by = "ProteinID", all.x = T)

  # ensure every ProteinID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(ProteinID, Covariate, unique = T), .SD, by = c("ProteinID", "Covariate"), keyby = "Covariate", all = T), by = Model]
  DT.de <- DT.de[, .(
    Covariate = Covariate,
    n1.test = first(n1.test[order(n1.test)]),
    n2.test = first(n2.test[order(n2.test)]),
    n1.real = first(n1.real[order(n1.real)]),
    n2.real = first(n2.real[order(n2.real)]),
    log2FC.lower,
    log2FC,
    log2FC.upper,
    t.value,
    p.value
  ), by = .(Model, ProteinID)]
  DT.de <- DT.de[!is.na(Covariate)]

  # sort for FDR
  setcolorder(DT.de, c("Model", "Covariate"))
  setorder(DT.de, Model, Covariate, p.value, na.last = T)
  DT.de[, FDR := p.adjust(p.value, method = "BH"), by = c("Model", "Covariate")]

  if (!as.data.table) setDF(DT.de)
  attr(DT.de, "bayesprot_dea_fits") <- file.path(fit, "model2", "de", output)
  return(DT.de)
}


#' Mixed-effects univariate differential expression analysis with 'MCMCglmm'
#'
#' Default is to a one-way ANOVA on the column 'Condition' in 'data.design'.
#'
#' @import data.table
#' @import foreach
#' @export
dea_MCMCglmm <- function(fit, data.design = design(fit), fixed = ~ Condition, prior = list(R = list(V = 1, nu = 0.02)), ...) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("mev" %in% names(arguments)) stop("do not pass a 'mev' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("verbose" %in% names(arguments)) stop("do not pass a 'verbose' argument to metafor")
  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")

  control = control(fit)
  fixed <- as.formula(sub("^.*~", "est ~", deparse(fixed)))

  # work out real n for each level and each protein
  DT.n <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T)[, .(ProteinID, AssayID, RawCount)]
  DT.n <- DT.n[complete.cases(DT.n)]
  DT.n <- unique(merge(DT.n, data.design, by = "AssayID")[, .(ProteinID, Sample)])
  DT.n <- DT.n[, .(N = length(Sample)), by = ProteinID]

  # prepare quants
  DTs <- protein_quants(fit, as.data.table = T)
  DTs <- merge(DTs, data.design, by = "AssayID")
  DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE)])])
  DTs[, SE := ifelse(SE < 0.001, 0.001, SE)] # won't work if SE is 0
  DTs[, BatchID := ProteinID]
  DTs[, ProteinID := as.character(ProteinID)]
  levels(DTs$BatchID) <- substr(levels(DTs$BatchID), 1, nchar(levels(DTs$BatchID)) - 1)
  DTs <- split(DTs, by = "BatchID", keep.by = F)

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  RNGkind("L'Ecuyer-CMRG")
  parallel::clusterSetRNGStream(cl, control$model.seed)

  pb <- txtProgressBar(max = length(DTs), style = 3)
  output.all <- foreach(DT.chunk = iterators::iter(DTs), .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    lapply(split(DT.chunk, by = "ProteinID"), function(DT) {
      # input
      output <- list(input = as.data.frame(DT))

      # test
      output$fits <- vector("list", 1)
      output$log <- capture.output({
        output$fits[[1]] <- MCMCglmm::MCMCglmm(fixed = fixed, mev = DT$SE^2, data = DT, prior = prior, verbose = F, ...)
      })

      # output
      if (!is.null(output$fits[[1]])) {
        sol <- summary(output$fits[[1]])$solutions

        if (nrow(sol) > 0) {
          n.real = DT.n[ProteinID == DT$ProteinID[1], N]

          data.table(
            Effect = rownames(sum$solutions),
            n.test = length(unique(DT$Sample)),
            n.real = ifelse(length(n.real) > 0, n.real, 0),
            log2FC.lower = sol[, "l-95% CI"],
            log2FC = sol[, "post.mean"],
            log2FC.upper = sol[, "u-95% CI"],
            pMCMC = sol[, "pMCMC"]
          )
        } else {
          NULL
        }
      }

      output
    })
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)
  parallel::stopCluster(cl)

  output.all <- unlist(output.all, recursive = F)
  class(output.all) <- "bayesprot_de_MCMCglmm"
  attr(output.all, "bayesprot_fit") <- fit
  return(output.all)
}


#' Pair-wise mixed-effects univariate differential expression analysis with 'MCMCglmm'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import foreach
#' @export
dea_MCMCglmm_pairwise <- function(fit, data.design = design(fit), fixed = ~ Condition, prior = list(R = list(V = 1, nu = 0.02)), ...) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("mev" %in% names(arguments)) stop("do not pass a 'mev' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("verbose" %in% names(arguments)) stop("do not pass a 'verbose' argument to metafor")
  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")

  control = control(fit)
  fixed <- as.formula(sub("^.*~", "est ~", deparse(fixed)))

  # work out real n for each condition of each protein
  DT.n <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T)[, .(ProteinID, AssayID, RawCount)]
  DT.n <- DT.n[complete.cases(DT.n)]
  DT.n <- unique(merge(DT.n, data.design, by = "AssayID")[, .(ProteinID, Sample, Condition)])
  setkey(DT.n, ProteinID, Condition)
  DT.n <- DT.n[CJ(ProteinID, Condition, unique = T), .N, by = .EACHI, allow.cartesian = T]

  # prepare quants
  DTs <- protein_quants(fit, as.data.table = T)
  DTs <- merge(DTs, data.design, by = "AssayID")
  DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE)])])
  cts <- combn(levels(DTs$Condition), 2)
  DTs[, SE := ifelse(SE < 0.001, 0.001, SE)] # won't work if SE is 0
  DTs[, BatchID := ProteinID]
  DTs[, ProteinID := as.character(ProteinID)]
  levels(DTs$BatchID) <- substr(levels(DTs$BatchID), 1, nchar(levels(DTs$BatchID)) - 1)
  DTs <- split(DTs, by = "BatchID", keep.by = F)

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  RNGkind("L'Ecuyer-CMRG")
  parallel::clusterSetRNGStream(cl, control$model.seed)

  pb <- txtProgressBar(max = length(DTs), style = 3)
  output.all <- foreach(DT.chunk = iterators::iter(DTs), .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    lapply(split(DT.chunk, by = "ProteinID"), function(DT) {
      # input
      output <- list(input = as.data.frame(DT))

      # test
      output$fits <- vector("list", ncol(cts))
      for (j in 1:length(output$fit)) {
        DT.input <- DT[as.character(Condition) %in% cts[, j]]
        DT.input[, Condition := factor(as.character(Condition), levels = cts[, j])]

        if (length(unique(DT.input$Sample[DT.input$Condition == cts[1, j]])) >= 2 && length(unique(DT.input$Sample[DT.input$Condition == cts[2, j]])) >= 2) {
          DT.input <- droplevels(DT.input)

          output$log <- paste0(output$log, "[", Sys.time(), "]  ", paste(cts[,j], collapse = "_v_"))
          output$log <- paste0(output$log, capture.output({
            output$fits[[j]] <- MCMCglmm::MCMCglmm(fixed = fixed, mev = DT.input$SE^2, data = DT.input, prior = prior, verbose = F)
          }))
        } else {
          output$log <- paste0(output$log, "[", Sys.time(), "]  ", paste(cts[,j], collapse = "_v_"), " ignored as n < 2 for one or both conditions\n")
        }
       }
      names(output$fits) <- sapply(1:ncol(cts), function(j) paste(cts[,j], collapse = "v"))

      # output
      output$output <- rbindlist(lapply(1:length(output$fits), function(j) {
        if (!is.null(output$fits[[j]])) {
          sol <- summary(output$fits[[j]])$solutions

          if (nrow(sol) > 0) {
            n1.real = DT.n[ProteinID == DT$ProteinID[1] & as.character(Condition) == cts[1, j], N]
            n2.real = DT.n[ProteinID == DT$ProteinID[1] & as.character(Condition) == cts[2, j], N]

            data.table(
              Effect = paste0(paste(cts[,j], collapse = "v"), "_", rownames(sol)),
              n1.test = length(unique(DT[as.character(Condition) == cts[1, j], Sample])),
              n2.test = length(unique(DT[as.character(Condition) == cts[2, j], Sample])),
              n1.real = ifelse(length(n1.real) > 0, n1.real, 0),
              n2.real = ifelse(length(n2.real) > 0, n2.real, 0),
              log2FC.lower = sol[, "l-95% CI"],
              log2FC = sol[, "post.mean"],
              log2FC.upper = sol[, "u-95% CI"],
              pMCMC = sol[, "pMCMC"]
            )
          } else {
            NULL
          }
        }
      }))

      output
    })
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)
  parallel::stopCluster(cl)

  output.all <- unlist(output.all, recursive = F)
  class(output.all) <- "bayesprot_de_MCMCglmm"
  attr(output.all, "bayesprot_fit") <- fit
  return(output.all)
}


#' Return differential expression as a list of FDR-controlled data tables
#'
#' @import data.table
#' @import metafor
#' @export
protein_de <- function(fit, key = 1, as.data.table = F) {
  if (class(fit) == "bayesprot_de_metafor") {

    DTs.de <- rbindlist(lapply(1:length(fit), function(i) {
      if (nrow(fit[[i]]$output) > 0) {
        data.table(ProteinID = factor(names(fit[i])), fit[[i]]$output)
      } else {
        NULL
      }
    }))

    DTs.de <- split(DTs.de, by = "Effect", keep.by = F)
    for (DT in DTs.de) {
      setorder(DT, p.value, na.last = T)
      DT[, FDR := p.adjust(p.value, method = "BH")]
      if (!as.data.table) setDF(DT)
    }
  } else if (class(fit) == "bayesprot_de_MCMCglmm") {
    DTs.de <- rbindlist(lapply(1:length(fit), function(i) {
      if (nrow(fit[[i]]$output) > 0) {
        data.table(ProteinID = factor(names(fit[i])), fit[[i]]$output)
      } else {
        NULL
      }
    }))
    DTs.de <- split(DTs.de, by = "Effect", keep.by = F)
    for (DT in DTs.de) {
      DT[, log2FC.delta := log2FC.upper - log2FC.lower]
      setorder(DT, pMCMC, log2FC.delta, na.last = T)
      DT[, log2FC.delta := NULL]
      DT[, FDR := cumsum(pMCMC) / .I]
      if (!as.data.table) setDF(DT)
    }
  }
  else {
    dea.func <- control(fit)$dea.func
    if (is.character(key)) key = match(key, names(dea.func))
    deID <- formatC(key, width = ceiling(log10(length(dea.func) + 1)) + 1, format = "d", flag = "0")
    DTs.de <- fst::read.fst(file.path(fit, "model2", "de", paste0(deID, ".fst")), as.data.table = T)
    DTs.de <- split(DTs.de, by = "Effect", keep.by = F)
  }

  return(DTs.de)
}


#' Mixed-effects univariate differential expression analysis with 'metafor' (test version)
#'
#' @import data.table
#' @import foreach
#' @export
# dea_metafor2 <- function(fit, data.design = design(fit), ...) {
#   arguments <- eval(substitute(alist(...)))
#   if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
#   if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
#   if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
#   if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
#   if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")
#   if (!("control" %in% names(arguments))) control = list()
#
#   DT.design <- as.data.table(data.design)
#   DT.assays <- design(fit, as.data.table = T)[, .(AssayID, Assay)]
#   DTs <- protein_quants(fit, summary = F, as.data.table = T)
#   DTs <- droplevels(DTs[complete.cases(DTs)])
#   DTs <- split(DTs, by = "ProteinID")
#
#   # start cluster and reproducible seed
#   cl <- parallel::makeCluster(control(fit)$nthread)
#   doSNOW::registerDoSNOW(cl)
#   pb <- txtProgressBar(max = length(DTs), style = 3)
#   fits.out <- foreach(DT = iterators::iter(DTs), .packages = "data.table", .multicombine = T, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
#     mat.cov <- dcast(DT, chainID + mcmcID ~ AssayID)
#     mat.cov[, chainID := NULL]
#     mat.cov[, mcmcID := NULL]
#     mat.cov <- cov(mat.cov)
#
#     #DT.in <- DT[, .(est = median(value), var = mad(value)^2), by = .(AssayID)]
#     DT.in <- DT[, .(est = median(value), var = var(value)), by = .(AssayID)]
#     DT.in <- merge(DT.in, DT.assays, by = "AssayID", sort = F)
#     DT.in <- merge(DT.in, DT.design, by = "Assay", sort = F)
#
#     fit.out <- NULL
#     for (i in 0:9) {
#       control$sigma2.init = 0.025 + 0.1 * i
#       try( {
#         fit.out <- metafor::rma.mv(yi = est, V = var, data = DT.in, control = control, test = "t", ...)
#         break
#       })
#     }
#     fit.out
#   }
#   setTxtProgressBar(pb, length(DTs))
#   close(pb)
#   parallel::stopCluster(cl)
#
#   names(fits.out) <- names(DTs)
#   class(fits.out) <- "bayesprot_de_metafor"
#   return(fits.out)
# }





