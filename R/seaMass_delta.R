#' seaMass-delta
#'
#' Perform seaMass-delta differential expression analysis and false discovery rate correction on normalised group-level
#' quants output by seaMass-theta
#' @include seaMass.R
#' @include seaMass_sigma.R
setClass("seaMass_delta", contains = "seaMass", slots = c(
  filepath = "character"
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
  name = "default",
  control = delta_control(),
  ...
) {
  # create delta object, checking for finished output
  if (!is.null(open_deltas(root(fit))[[name]])) stop(paste0("ERROR: completed seaMass-delta found at ", name))
  object <- new("seaMass_delta", filepath = file.path(dirname(filepath(fit)), paste0("delta.", name)))

  ### INIT
  control@version <- as.character(packageVersion("seaMass"))
  cat(paste0("[", Sys.time(), "] seaMass-delta v", control@version, "\n"))
  control@fit <- fit
  control@nchain <- control(fit)@nchain
  control@nsample <- control(fit)@nsample
  control@nthread <- control(fit)@nthread
  control@ellipsis <- list(...)
  validObject(control)

  data.table::setDTthreads(control(fit)@nthread)
  fst::threads_fst(control(fit)@nthread)

  # create fit and output directories
  if (file.exists(filepath(object))) unlink(filepath(object), recursive = T)
  if (!dir.create(filepath(object))) stop("ERROR: problem creating folder")

  # check and save control
  saveRDS(control, file.path(filepath(object), "control.rds"))

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
  fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))

  # run
  prepare_delta(control(root(fit))@schedule, object)

  if (file.exists(file.path(filepath(fit), "complete.rds"))) {
    run(object)
  } else {
    cat(paste0("[", Sys.time(), "] queued\n"))
  }

  return(invisible(object))
}


#' @describeIn seaMass_delta-class Run.
#' @export
#' @include generics.R
setMethod("run", "seaMass_delta", function(object) {
  job.id <- uuid::UUIDgenerate()

  cat(paste0("[", Sys.time(), "] processing...\n"))
  for (chain in 1:control(object)@nchain) process(object, chain, job.id)
  cat(paste0("[", Sys.time(), "] finished processing!\n"))

  cat(paste0("[", Sys.time(), "] reporting...\n"))
  for (batch in 1:control(root(object))@plot.nbatch) plots(object, batch, job.id)
  report(object)
  cat(paste0("[", Sys.time(), "] finished reporting!\n"))

  return(invisible(object))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("report", "seaMass_delta", function(object, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  cat(paste0("[", Sys.time(), "]  REPORT\n"))

  if (!("group.quants.de" %in% ctrl@keep)) unlink(file.path(filepath(object), "group.quants.de*"), recursive = T)
  if (!("component.deviations.de" %in% ctrl@keep)) unlink(file.path(filepath(object), "component.deviations.de*"), recursive = T)

  # assemble report
  cat(paste0("[", Sys.time(), "]   assembling html report...\n"))
  assemble_report(root(object))

  return(invisible(NULL))
})


#' @describeIn seaMass_delta-class Get name.
#' @export
#' @include generics.R
setMethod("name", "seaMass_delta", function(object) {
  return(sub("^delta\\.", "", basename(object@filepath)))
})


#' @describeIn seaMass_theta-class Get the \code{seaMass_sigma} object.
#' @export
#' @include generics.R
setMethod("root", "seaMass_delta", function(object) {
  return(root(control(object)@fit))
})


#' @describeIn seaMass_delta-class Get path.
#' @include generics.R
#' @export
setMethod("filepath", "seaMass_delta", function(object) {
  return(object@filepath)
})


#' @describeIn seaMass_delta-class Get the \link{seaMass_sigma} object.
#' @export
#' @include generics.R
setMethod("parent", "seaMass_delta", function(object) {
  return(control(object@fit))
})


#' @describeIn seaMass_delta-class Get the \link{fit_control}.
#' @export
#' @include generics.R
setMethod("control", "seaMass_delta", function(object) {
  if (!file.exists(file.path(filepath(object), "control.rds")))
    stop(paste0("seaMass-delta output '", name(object), "' is missing"))

  return(readRDS(file.path(filepath(object), "control.rds")))
})


#' @include generics.R
setMethod("blocks", "seaMass_delta", function(object) {
  return(NULL)
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
setMethod("group_quants_de", "seaMass_delta", function(object, groups = NULL, summary = TRUE, type = "group.quants.de", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  DT <- read(object, ".", type, groups, chains, summary, summary.func = "lst_ash", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations_de", "seaMass_delta", function(object, components = NULL, summary = TRUE, type = "component.deviations.de", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  DT <- read(object, ".", type, components, chains, summary, summary.func = "lst_ash", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_quants_fdr", "seaMass_delta", function(object, groups = NULL, type = "group.quants.fdr", as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), paste0(type, ".fst")))) {
    DT <- fst::read.fst(file.path(filepath(object), paste0(type, ".fst")), as.data.table = as.data.table)
    if (!is.null(groups)) {
      if (is.data.frame(groups)) {
        DT <- merge(DT, groups, by = colnames(groups), sort = F)
      } else {
        DT <- DT[DT$Group %in% groups,]
      }
    }
    return(DT)
  } else {
    return(NULL)
  }
})


#' @describeIn seaMass_delta-class Get the component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations_fdr", "seaMass_delta", function(object, groups = NULL, type = "component.deviations.fdr", as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), paste0(type, ".fst")))) {
    DT <- fst::read.fst(file.path(filepath(object), paste0(type, ".fst")), as.data.table = as.data.table)
    if (!is.null(groups)) DT <- DT[DT$Group %in% groups,]
    return(DT)
  } else {
    return(NULL)
  }
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_quants_de", "seaMass_delta", function(
  object,
  groups = NULL,
  summary = TRUE,
  colour = "black",
  variable.summary.cols = c("Group", "Effect", "Contrast", "Baseline", "Cont.uS", "Base.uS", "Cont.qS", "Base.qS", "Cont.qC", "Base.qC", "Cont.qM", "Base.qM"),
  variable.label.cols = c("Group", "Contrast", "Baseline"),
  value.label = "fold change",
  ...
) {
  return(plot_dists(
    object,
    data = group_quants_de(object, groups, summary = summary, as.data.table = T),
    colour = colour,
    variable.summary.cols = variable.summary.cols,
    variable.label.cols = variable.label.cols,
    value.label = value.label,
    ...
  ))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_quants_fdr", "seaMass_delta", function(
  object,
  groups = NULL,
  summary = TRUE,
  colour = list("lfdr", "grey"),
  variable.summary.cols = c("qvalue", "Batch", "Effect", "Contrast", "Baseline", "Group", "Cont.uS", "Base.uS", "Cont.qS", "Base.qS",
                            "Cont.qC", "Base.qC", "Cont.qM", "Base.qM", "lfdr", "lfsr", "svalue", "NegativeProb", "PositiveProb"),
  variable.label.cols = c("Group", "qvalue"),
  value.label = "fold change",
  ...
) {
  de.groups <- groups
  if (is.data.frame(de.groups) && "Batch" %in% colnames(de.groups)) {
    de.groups$Batch <- NULL
    de.groups <- unique(de.groups)
    if (nrow(de.groups) == 0) de.groups <- NULL
  }

  return(plot_dists(
    object,
    data = list(
      group_quants_fdr(object, groups, as.data.table = T),
      group_quants_de(object, de.groups, summary = summary, as.data.table = T)
    ),
    colour = colour,
    variable.summary.cols = variable.summary.cols,
    variable.label.cols = variable.label.cols,
    value.label = value.label,
    ...
  ))
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
