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
  set.seed(control@random.seed)

  # create fit and output directories
  object <- new("seaMass_delta", fit = fit, name = paste0("delta.", name))
  path <- filepath(object)
  if (file.exists(path)) unlink(path, recursive = T)
  dir.create(file.path(path, "meta"), recursive = T)
  if (file.exists(file.path(dirname(path), "output", basename(path)))) unlink(file.path(dirname(path), "output", basename(path)), recursive = T)
      dir.create(file.path(dirname(path), "output", basename(path)), recursive = T)

  # check and save control
  saveRDS(control, file.path(path, "meta", "control.rds"))

  # get design into the format we need
  DT.design <- as.data.table(data.design)[!is.na(Assay)]
  blocks <- grep("^Block\\.", colnames(DT.design))
  if (length(blocks) > 0) set(DT.design, j = blocks, value = NULL)
  cols <- which(colnames(DT.design) %in% c("qG", "uG", "nG", "qC", "uC", "nC", "qM", "uM", "nM", "qD", "uD", "nD"))
  if (length(cols) > 0) set(DT.design, j = cols, value = NULL)
  fst::write.fst(DT.design, file.path(path, "meta", "design.fst"))

  # run
  prepare_delta(control(fit)@schedule, object)

  if (completed(fit)) {
    run(object)
  } else {
    cat(paste0("[", Sys.time(), "] queued\n"))
  }

  return(invisible(object))
}


#' @describeIn seaMass_delta-class Open a complete \code{seaMass_delta} run from the supplied \link{seaMass_sigma} fir object and \code{name}.
#' @export
open_delta <- function(fit, name = "fit", quiet = FALSE, force = FALSE) {
  path <- file.path(filepath(fit), paste0("delta.", name))

  if(force || file.exists(file.path(path, "complete"))) {
    object <- new("seaMass_delta", fit = fit, name = paste0("delta.", name))
    return(object)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      if (force) stop("'", path, "' does not contain a completed seaMass-delta run")
      else stop("'", path, "' does not contain a seaMass-delta run")
    }
  }
}


#' @describeIn seaMass_delta-class Run.
#' @export
#' @include generics.R
setMethod("run", "seaMass_delta", function(object) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", name(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  ellipsis <- ctrl@ellipsis
  ellipsis$object <- object

  cat(paste0("[", Sys.time(), "]  running delta with name=", name(object), " for sigma fit=", name(object@fit), "...\n"))

  # group dea
  if (ctrl@dea.model != "" && !all(is.na(assay_design(object, as.data.table = T)$Condition))) do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)

  # summarise group de and perform fdr correction
  if (file.exists(file.path(filepath(object), "group.de.index.fst"))) {
    if(ctrl@fdr.model != "") {
      do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
    } else {
      group_de(object, summary = T, as.data.table = T)
    }
  }
  if (!("group.de" %in% ctrl@keep)) unlink(file.path(filepath(object), "group.de*"), recursive = T)

  # write out group fdr
  if (file.exists(file.path(filepath(object), "group.fdr.fst"))) {
    DTs.fdr <- split(group_fdr(object, as.data.table = T), drop = T, by = "Batch")
    for (name in names(DTs.fdr)) {
      # save
      fwrite(DTs.fdr[[name]], file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_group_de", gsub("\\s", "_", name), "csv", sep = ".")))
      # plot fdr
      plot_fdr(DTs.fdr[[name]])
      ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_group_de", gsub("\\s", "_", name), "pdf", sep = ".")), width = 8, height = 8)
    }
  }
  if (!("group.fdr" %in% ctrl@keep)) unlink(file.path(filepath(object), "group.fdr*"), recursive = T)

  # component deviations
  if (ctrl@component.deviations == T && control(object@fit)@component.model == "independent") {
    # component deviation dea
    ellipsis$type <- "component.deviations.de"
    do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)

    # summarise component deviation de and perform fdr correction
    if (file.exists(file.path(filepath(object), "component.deviations.de.index.fst"))) {
      if(ctrl@fdr.model != "") {
        ellipsis$type <- "component.deviations.fdr"
        do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
      } else {
        component_deviations_de(object, summary = T, as.data.table = T)
      }
    }
    if (!("component.deviations.de" %in% ctrl@keep)) unlink(file.path(filepath(object), "component.deviations.de*"), recursive = T)

    # write out group fdr
    if (file.exists(file.path(filepath(object), "component.deviations.fdr.fst"))) {
      DTs.fdr <- split(component_deviations_fdr(object, as.data.table = T), drop = T, by = "Batch")
      for (name in names(DTs.fdr)) {
        # save
        fwrite(DTs.fdr[[name]], file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_component_deviations_de", gsub("\\s", "_", name), "csv", sep = ".")))
        # plot fdr
        plot_fdr(DTs.fdr[[name]])
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_component_deviations_de", gsub("\\s", "_", name), "pdf", sep = ".")), width = 8, height = 8)
      }
    }
    if (!("component.deviations.fdr" %in% ctrl@keep)) unlink(file.path(filepath(object), "component.deviations.fdr*"), recursive = T)
  }


  # set complete
  write.table(data.frame(), file.path(filepath(object), "complete"), col.names = F)

  return(object@name)
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


#' @include generics.R
setMethod("blocks", "seaMass_delta", function(object) {
  return(NULL)
})


#' @describeIn seaMass_delta-class Get the \link{fit_control}.
#' @export
#' @include generics.R
setMethod("control", "seaMass_delta", function(object) {
  if (!file.exists(file.path(filepath(object), "meta", "control.rds")))
    stop(paste0("seaMass-delta output '", name(object), "' is missing"))

  return(readRDS(file.path(filepath(object), "meta", "control.rds")))
})


#' @describeIn seaMass_delta-class Get the study design as a \code{data.frame}.
#' @export
#' @include generics.R
setMethod("assay_design", "seaMass_delta", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "design.fst"), as.data.table = T)

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
setMethod("group_de", "seaMass_delta", function(object, groups = NULL, summary = FALSE, chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, ".", "group.de", groups, chains, summary, summary.func = "dist_lst_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations_de", "seaMass_delta", function(object, components = NULL, summary = FALSE, chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, ".", "component.deviations.de", components, chains, summary, summary.func = "dist_lst_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the false discovery rate corrected group differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_fdr", "seaMass_delta", function(object, as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), "group.fdr.fst"))) {
    return(fst::read.fst(file.path(filepath(object), "group.fdr.fst"), as.data.table = as.data.table))
  } else {
    return(NULL)
  }
})


#' @describeIn seaMass_delta-class Get the false discovery rate corrected component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations_fdr", "seaMass_delta", function(object, as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), "component.deviations.fdr.fst"))) {
    return(fst::read.fst(file.path(filepath(object), "component.deviations.fdr.fst"), as.data.table = as.data.table))
  } else {
    return(NULL)
  }
})


#' @describeIn seaMass_delta-class Add our spikein ground truth to the object of \code{group_fdr} pr \code{component_deviations_fdr}.
#' @import data.table
#' @export
add_seaMass_spikein_truth <- function(data.fdr) {
  # ground truth
  data.truth <- rbind(
    data.frame(set = "Rat [1:1]",      truth =           0, grep = "_RAT$"),
    data.frame(set = "E.coli [10:16]", truth = log2(16/10), grep = "_ECOLX$"),
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
  data <- merge(data.fdr, data[, c("Group", "truth")], by = "Group", sort = F, all.x = T)

  return(data)
}
