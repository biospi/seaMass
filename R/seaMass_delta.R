#' seaMass-Δ
#'
#' Perform seaMass-Δ normalisation, differential expression analysis and/or false discovery rate correction on raw group-level
#' quants output by seaMass-Σ
#' @include seaMass.R
#' @include seaMass_sigma.R
setClass("seaMass_delta", contains = "seaMass", slots = c(
  fit = "seaMass_sigma",
  name = "character"
))


#' @describeIn seaMass_delta-class Runs seaMass-Σ.
#' @param fit A \link{seaMass_sigma} object as returned by seaMass-Σ.
#' @param data.design Optionally, a \link{data.frame} created by \link{new_design} and then customised, which specifies
#'   sample names, RefWeight channels, and any covariates specified by the experimental design.
#' @param name Name of folder on disk where all intermediate and output data will be stored; default is \code{"output"}.
#' @param control A control object created with \link{new_delta_control} specifying control parameters for the differential analysis.
#' @return A \code{seaMass_delta} object that can be interrogated for various metadata and results.
#' @import data.table
#' @export seaMass_delta
seaMass_delta <- function(
  fit,
  data.design = assay_design(fit),
  name = "fit",
  control = delta_control(),
  ...
) {
  # check for finished output and return that
  object <- open_delta(fit, name, quiet = T)
  if (!is.null(object)) {
    message(paste0("returning completed seaMass-delta object - if this wasn't your intention, supply a different 'name' or delete the folder for the returned object with 'del(object)'"))
    return(object)
  }

  ### INIT
  cat(paste0("[", Sys.time(), "] seaMass-delta v", control@version, "\n"))
  control@model.nchain <- control(fit)@model.nchain
  control@model.nsample <- control(fit)@model.nsample
  control@nthread <- control(fit)@nthread
  control@ellipsis <- list(...)
  validObject(control)

  data.table::setDTthreads(control@nthread)
  fst::threads_fst(control@nthread)

  # create fit and output directories
  object <- new("seaMass_delta", fit = fit, name = paste0("delta.", name))
  path <- filepath(object)
  if (file.exists(path)) unlink(path, recursive = T)
  dir.create(path)
  if (file.exists(file.path(dirname(path), "output", basename(path)))) unlink(file.path(dirname(path), "output", basename(path)), recursive = T)
  dir.create(file.path(dirname(path), "output", basename(path)))

  # check and save control
  saveRDS(control, file.path(path, "control.rds"))

  # get design into the format we need
  DT.design <- as.data.table(data.design)[!is.na(Assay)]
  DT.sigma.design <- assay_design(fit, as.data.table = T)[, .(Assay, Block)]
  if ("Block" %in% colnames(DT.design)) {
    DT.design <- merge(DT.design, DT.sigma.design, by = c("Assay", "Block"), sort = F)
  } else {
    DT.design <- merge(DT.design, DT.sigma.design, by = "Assay", sort = F, allow.cartesian = T)
    setcolorder(DT.design, c("Assay", "Block"))
  }
  if (!is.factor(DT.design$Assay)) DT.design[, Assay := factor(Assay, levels = levels(DT.sigma.design$Assay))]
  if (!is.factor(DT.design$Block)) DT.design[, Block := factor(Block, levels = levels(DT.sigma.design$Block))]
  blocks <- grep("^Block\\.", colnames(DT.design))
  if (length(blocks) > 0) set(DT.design, j = blocks, value = NULL)
  cols <- which(colnames(DT.design) %in% c("qG", "uG", "nG", "qC", "uC", "nC", "qM", "uM", "nM", "qD", "uD", "nD"))
  if (length(cols) > 0) set(DT.design, j = cols, value = NULL)
  fst::write.fst(DT.design, file.path(path, "design.fst"))

  # run
  prepare_delta(control(fit)@schedule, object)

  if (file.exists(file.path(filepath(fit), "sigma", "complete.rds"))) {
    run(object)
    cat(paste0("[", Sys.time(), "] finished!\n"))
  } else {
    cat(paste0("[", Sys.time(), "] queued\n"))
  }

  return(invisible(object))
}


#' @describeIn seaMass_delta-class Open a complete \code{seaMass_delta} run from the supplied \link{seaMass_sigma} fir object and \code{name}.
#' @export
open_delta <- function(fit, name = "fit", quiet = FALSE, force = FALSE) {
  path <- file.path(filepath(fit), paste0("delta.", name))

  object <- new("seaMass_delta", fit = fit, name = paste0("delta.", name))
  if (!force && read_completed(file.path(filepath(object))) == 0) {
    if (quiet) {
      return(NULL)
    } else {
      stop("'", path, "' is not complete")
    }
  }

  return(object)
}


#' @describeIn seaMass_delta-class Run.
#' @export
#' @include generics.R
setMethod("run", "seaMass_delta", function(object) {
  job.id <- uuid::UUIDgenerate()
  for (chain in 1:control(object)@model.nchain) process(object, chain, job.id)

  return(invisible(object))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("finish", "seaMass_delta", function(object, job.id) {
  ctrl <- control(object)
  if (!("de.standardised.group.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "standardised.group.deviations*"), recursive = T)
  if (!("de.component.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "component.deviations*"), recursive = T)

  return(invisible(NULL))
})


#' @describeIn seaMass_delta-class Get name.
#' @export
#' @include generics.R
setMethod("name", "seaMass_delta", function(object) {
  return(sub("^delta\\.", "", object@name))
})


#' @describeIn seaMass_delta-class Get path.
#' @include generics.R
#' @export
setMethod("filepath", "seaMass_delta", function(object) {
  return(file.path(filepath(object@fit), object@name))
})


#' @describeIn seaMass_delta-class Get the \link{seaMass_sigma} object.
#' @export
#' @include generics.R
setMethod("parent", "seaMass_delta", function(object) {
  return(object@fit)
})


#' @include generics.R
setMethod("blocks", "seaMass_delta", function(object) {
  return(NULL)
})


#' @describeIn seaMass_delta-class Get the \link{fit_control}.
#' @export
#' @include generics.R
setMethod("control", "seaMass_delta", function(object) {
  if (!file.exists(file.path(filepath(object), "control.rds")))
    stop(paste0("seaMass-delta output '", name(object), "' is missing"))

  return(readRDS(file.path(filepath(object), "control.rds")))
})


#' @describeIn seaMass_delta-class Get the study design as a \code{data.frame}.
#' @export
#' @include generics.R
setMethod("assay_design", "seaMass_delta", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "design.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("groups", "seaMass_delta", function(object, as.data.table = FALSE) {
  return(groups(parent(object), as.data.table))
})


#' @describeIn seaMass_delta-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("components", "seaMass_delta", function(object, as.data.table = FALSE) {
  return(components(parent(object), as.data.table))
})


#' @describeIn seaMass_delta-class Get the assay groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_groups", "seaMass_delta", function(object, as.data.table = FALSE) {
  return(assay_groups(parent(object), as.data.table))
})


#' @describeIn seaMass_delta-class Get the assay groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_components", "seaMass_delta", function(object, as.data.table = FALSE) {
  return(assay_components(parent(object), as.data.table))
})


#' @describeIn seaMass_delta-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurements", "seaMass_delta", function(object, as.data.table = FALSE) {
  return(measurements(parent(object), as.data.table))
})


#' @describeIn seaMass_delta-class Delete the \code{seaMass_delta} run from disk.
#' @import data.table
#' @export
#' @include generics.R
setMethod("del", "seaMass_delta", function(object) {
  return(unlink(filepath(object), recursive = T))
})


#' @describeIn seaMass_delta-class Get the model differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("de_standardised_group_deviations", "seaMass_delta", function(object, groups = NULL, summary = FALSE, type = "standardised.group.deviations", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, ".", type, groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("de_component_deviations", "seaMass_delta", function(object, type = "component.deviations", as.data.table = FALSE) {
  DT <- read_samples(object, ".", type, components, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("fdr_standardised_group_deviations", "seaMass_delta", function(object, groups = NULL, summary = TRUE, type = "standardised.group.deviations", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), paste0("fdr.", type, ".fst")))) {
    DT <- fst::read.fst(file.path(filepath(object), paste0("fdr.", type, ".fst")), as.data.table = as.data.table)
    if (!is.null(groups)) DT <- DT[DT$Group %in% groups,]
    return(DT)
  } else {
    return(NULL)
  }
})


#' @describeIn seaMass_delta-class Get the component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("fdr_component_deviations", "seaMass_delta", function(object, groups = NULL, summary = TRUE, type = "component.deviations", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), paste0("fdr.", type, ".fst")))) {
    DT <- fst::read.fst(file.path(filepath(object), paste0("fdr.", type, ".fst")), as.data.table = as.data.table)
    if (!is.null(groups)) DT <- DT[DT$Group %in% groups,]
    return(DT)
  } else {
    return(NULL)
  }
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_de_standardised_group_deviations", "seaMass_delta", function(
  object, data, limits = NULL, alpha = 1,
  facets = ~ interaction(Effect, interaction(Contrast, Baseline, drop = T, sep = " - ", lex.order = T), drop = T, sep = " : ", lex.order = T),
  sort.cols = "Group", label.cols = "Group", title = NULL, horizontal = TRUE, colour = "G.qC", fill = "G.qC", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, title, value.label = "deviation", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_de_component_deviations", "seaMass_delta", function(
  object, data, limits = NULL, alpha = 1,
  facets = ~ interaction(Effect, interaction(Contrast, Baseline, drop = T, sep = " - ", lex.order = T), drop = T, sep = " : ", lex.order = T),
  sort.cols = "Group", label.cols = "Group", title = NULL, horizontal = TRUE, colour = "G.qC", fill = "G.qC", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, title, value.label = "deviation", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_fdr_standardised_group_deviations", "seaMass_delta", function(
  object, data, limits = NULL, alpha = 1,
  facets = ~ interaction(Effect, interaction(Contrast, Baseline, drop = T, sep = " - ", lex.order = T), drop = T, sep = " : ", lex.order = T),
  sort.cols = NULL, label.cols = c("Group", "qvalue"), title = NULL, horizontal = TRUE, colour = "G.qC", fill = "G.qC", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, title, value.label = "deviation", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_fdr_component_deviations", "seaMass_delta", function(
  object, data, limits = NULL, alpha = 1,
  facets = ~ interaction(Effect, interaction(Contrast, Baseline, drop = T, sep = " - ", lex.order = T), drop = T, sep = " : ", lex.order = T),
  sort.cols = NULL, label.cols = "Group", title = NULL, horizontal = TRUE, colour = "G.qC", fill = "G.qC", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, title, value.label = "deviation", horizontal, colour, fill, file, value.length, level.length))
})


#' @describeIn seaMass_delta-class Add our spikein ground truth to the object of \code{group_fdr} pr \code{component_deviations_fdr}.
#' @import data.table
#' @export
add_seaMass_spikein_truth <- function(data.fdr) {
  # ground truth
  data.truth <- rbind(
    data.frame(set = "Rat [1:1]",      truth =           0, grep = "_RAT$"),
    data.frame(set = "E.coli [10:16]", truth = log2(16/10), grep = "_ECOL[I|X]$"),
    data.frame(set = "[3:1]",          truth =   log2(1/3), grep = "^sp\\|P00330\\|ADH1_YEAST$"),
    data.frame(set = "[5:3]",          truth =   log2(3/5), grep = "^sp\\|P08603\\|CFAH_HUMAN$"),
    data.frame(set = "[3:2]",          truth =   log2(2/3), grep = "^sp\\|P02769\\|ALBU_BOVIN$|^sp\\|P00698\\|LYSC_CHICK$|^sp\\|P00711\\|LALBA_BOVIN$|^sp\\|P00915\\|CAH1_HUMAN$"),
    data.frame(set = "[4:3]",          truth =   log2(3/4), grep = "^sp\\|P46406\\|G3P_RABIT$|^sp\\|P00489\\|PYGM_RABIT$"),
    data.frame(set = "[5:4]",          truth =   log2(4/5), grep = "^sp\\|P00004\\|CYC_HORSE$|^sp\\|P02754\\|LACB_BOVIN$"),
    data.frame(set = "[3:4]",          truth =   log2(4/3), grep = "^sp\\|P02666\\|CASB_BOVIN$|^sp\\|P01012\\|OVAL_CHICK$"),
    data.frame(set = "[15:28]",        truth = log2(28/15), grep = "^sp\\|P00432\\|CATA_BOVIN$"),
    data.frame(set = "[1:2]",          truth =   log2(2/1), grep = "^sp\\|P68082\\|MYG_HORSE$|^sp\\|P06278\\|AMY_BACLI$|^sp\\|Q29443\\|TRFE_BOVIN$")
  )

  # unlist proteins from protein groups
  s <- strsplit(as.character(data.fdr$Group), split = ";")
  data <- data.frame(Group = rep(data.fdr$Group, sapply(s, length)), Protein = unlist(s))
  # initialize new varible with NAs
  data$truth <- NA
  # fill in matching indices
  for (i in 1:nrow(data.truth)) data$truth[grepl(data.truth$grep[i], data$Protein)] <- data.truth$truth[i]
  # remove duplicates
  data <- data[!duplicated(data[, c("Group", "truth")]),]
  # remove groups that have multiple truths
  data <- data[!duplicated(data$Group) & !duplicated(data$Group, fromLast = T),]
  # merge
  cols <- colnames(data.fdr)
  data <- merge(data.fdr, data[, c("Group", "truth")], by = "Group", sort = F, all.x = T)
  setcolorder(data, cols)

  return(data)
}


#' @describeIn seaMass_delta-class Add Navarro spikein ground truth to the object of \code{group_fdr} pr \code{component_deviations_fdr}.
#' @import data.table
#' @export
add_Navarro_spikein_truth <- function(data.fdr) {
  # ground truth
  data.truth <- rbind(
    data.frame(set = "Human [1:1]",  truth = 0, grep = "_HUMAN$"),
    data.frame(set = "E.coli [1:4]", truth = 2, grep = "_ECOLI$"),
    data.frame(set = "Yeast [2:1]",  truth = -1, grep = "_YEAS8$")
  )

  # unlist proteins from protein groups
  s <- strsplit(as.character(data.fdr$Group), split = ";")
  data <- data.frame(Group = rep(data.fdr$Group, sapply(s, length)), Protein = unlist(s))
  # initialize new varible with NAs
  data$truth <- NA
  # fill in matching indices
  for (i in 1:nrow(data.truth)) data$truth[grepl(data.truth$grep[i], data$Protein)] <- data.truth$truth[i]
  # remove duplicates
  data <- data[!duplicated(data[, c("Group", "truth")]),]
  # remove groups that have multiple truths
  data <- data[!duplicated(data$Group) & !duplicated(data$Group, fromLast = T),]
  # merge
  cols <- colnames(data.fdr)
  data <- merge(data.fdr, data[, c("Group", "truth")], by = "Group", sort = F, all.x = T)
  setcolorder(data, cols)

  return(data)
}
