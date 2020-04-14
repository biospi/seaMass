setGeneric("components", function(object, ...) standardGeneric("components"))
setGeneric("groups", function(object, ...) standardGeneric("groups"))
setGeneric("unnormalised_group_quants", function(object, ...) standardGeneric("unnormalised_group_quants"))
setGeneric("normalised_group_quants", function(object, ...) standardGeneric("normalised_group_quants"))
setGeneric("component_deviations", function(object, ...) standardGeneric("component_deviations"))
setGeneric("group_variances", function(object, ...) standardGeneric("group_variances"))
setGeneric("group_de", function(object, ...) standardGeneric("group_de"))
setGeneric("component_deviations_de", function(object, ...) standardGeneric("component_deviations_de"))
setGeneric("group_fdr", function(object, ...) standardGeneric("group_fdr"))
setGeneric("component_deviations_fdr", function(object, ...) standardGeneric("component_deviations_fdr"))


#' seaMass-Δ
#'
#' Perform seaMass-Δ normalisation, differential expression analysis and/or false discovery rate correction on unnormalised group-level
#' quants output by seaMass-Σ
#'
#' @include seaMass_sigma.R
setClass("seaMass_delta", slots = c(
  sigma = "seaMass_sigma",
  name = "character"
))


#' @describeIn seaMass_delta-class Runs seaMass-Σ.
#' @param sigma A \link{seaMass_sigma} object as returned by seaMass-Σ.
#' @param data.design Optionally, a \link{data.frame} created by \link{new_design} and then customised, which specifies
#'   sample names, RefWeight channels, and any covariates specified by the experimental design.
#' @param name Name of folder on disk where all intermediate and output data will be stored; default is \code{"output"}.
#' @param control A control object created with \link{new_delta_control} specifying control parameters for the differential analysis.
#' @return A \code{seaMass_delta} object that can be interrogated for various metadata and results.
#' @import data.table
#' @export seaMass_delta
seaMass_delta <- function(
  sigma,
  data.design = assay_design(sigma),
  name = "fit",
  control = delta_control(),
  ...
) {
  # check for finished output and return that
  object <- open_seaMass_delta(sigma, name, quiet = T)
  if (!is.null(object)) {
    cat(paste0(" returning completed seaMass-delta object - if this wasn't your intention, supply a different 'name' or delete the folder for the returned object with 'del(object)'"))
    return(object)
  }

  ### INIT
  cat(paste0("[", Sys.time(), "] seaMass-delta v", control@version))
  control@model.nchain <- control(sigma)@model.nchain
  control@model.nsample <- control(sigma)@model.nsample
  control@nthread <- control(sigma)@nthread
  control@ellipsis <- list(...)
  validObject(control)

  data.table::setDTthreads(control@nthread)
  fst::threads_fst(control@nthread)
  set.seed(control@random.seed)

  # create fit and output directories
  object <- new("seaMass_delta", sigma = sigma, name = paste0("delta.", name))
  path <- path(object)
  if (file.exists(path)) unlink(path, recursive = T)
  dir.create(path)
  dir.create(file.path(path, "meta"))
  dir.create(file.path(path, "output"))

  # check and save control
  saveRDS(control, file.path(path, "meta", "control.rds"))

  # merged design
  DT.design <- unique(as.data.table(data.design)[!is.na(Assay)], by = "Assay")

  # get design into the format we need
  DT.design <- as.data.table(data.design)[!is.na(Assay)]
  if (!is.factor(DT.design$Assay)) DT.design[, Assay := factor(as.character(Assay), levels = levels(assay_design(sigma)$Assay))]
  if ("nGroup" %in% colnames(DT.design)) DT.design[, nGroup := NULL]
  if ("nComponent" %in% colnames(DT.design)) DT.design[, nComponent := NULL]
  if ("nMeasurement" %in% colnames(DT.design)) DT.design[, nMeasurement := NULL]
  if ("nDatapoint" %in% colnames(DT.design)) DT.design[, nDatapoint := NULL]
  if ("Block" %in% colnames(DT.design)) DT.design[, Block := NULL]
  DT.design <- merge(
    DT.design,
    assay_design(sigma, as.data.table = T)[, .(nGroup = max(nGroup), nComponent = max(nComponent), nMeasurement = max(nMeasurement), nDatapoint = max(nDatapoint)), by = Assay],
    by = "Assay"
  )
  if (!is.factor(DT.design$Assay)) DT.design[, Assay := factor(Assay, levels = unique(Assay))]
  fst::write.fst(DT.design, file.path(path, "meta", "design.fst"))

  # merged groups
  DT.groups <- rbindlist(lapply(fits(sigma), function(fit) groups(fit, as.data.table = T)))
  DT.groups <- DT.groups[, .(GroupInfo = first(GroupInfo), nComponent = sum(nComponent), nMeasurement = sum(nMeasurement), nDatapoint = sum(nDatapoint)), by = Group]
  setorder(DT.groups, -nComponent, -nMeasurement, -nDatapoint, Group)
  DT.groups[, Group := factor(Group, levels = unique(Group))]
  fst::write.fst(DT.groups, file.path(path, "meta", "groups.fst"))

  # merged components
  DT.components <- rbindlist(lapply(fits(sigma), function(fit) components(fit, as.data.table = T)))
  DT.components <- DT.components[, .(nMeasurement = sum(nMeasurement), nDatapoint = sum(nDatapoint)), by = Component]
  setorder(DT.components, -nMeasurement, -nDatapoint, Component)
  DT.components[, Component := factor(Component, levels = unique(Component))]
  fst::write.fst(DT.components, file.path(path, "meta", "components.fst"))

  # run
  prepare_delta(control(sigma)@schedule, object)

  if (completed(sigma)) {
    run(object)
  } else {
    cat(" seaMass-delta will be run when 'run(sigma)' is run")
  }

  return(invisible(object))
}


#' @describeIn seaMass_delta-class Open a complete \code{seaMass_delta} run from the supplied \link{seaMass_sigma} and \code{name}.
#' @export
open_seaMass_delta <- function(sigma, name = "fit", quiet = FALSE, force = FALSE) {
  path <- file.path(path(sigma), paste0("delta.", name))

  if(force || file.exists(file.path(path, "complete"))) {
    object <- new("seaMass_delta", sigma = sigma, name = paste0("delta.", name))
    return(object)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      if (force) stop("'", path, "' does not contain a completed seaMass-Δ run")
      else stop("'", path, "' does not contain a seaMass-Δ run")
    }
  }
}


#' @describeIn seaMass_delta-class Run.
#' @export
setMethod("run", "seaMass_delta", function(object) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
      stop(paste0("version mismatch - 'object' was created with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  DT.design <- assay_design(object, as.data.table = T)
  ellipsis <- ctrl@ellipsis
  ellipsis$object <- object

  cat(paste0("[", Sys.time(), "]  running delta for name=", name(object), "..."))

  # standardise quants using reference weights
  standardise_group_quants(object)
  component.model <- control(object@sigma)@component.model
  if (ctrl@component.deviations == T && component.model == "independent") standardise_component_deviations(object)

  # normalise quants by norm.groups
  if (ctrl@norm.model != "") do.call(paste("norm", ctrl@norm.model, sep = "_"), ellipsis)

  # group quants
  cat("[", paste0(Sys.time(), "]   summarising normalised group quants..."))
  set.seed(ctrl@random.seed)
  DT.groups <- groups(object, as.data.table = T)
  DT.group.quants <- normalised_group_quants(object, summary = T, as.data.table = T)
  DT.group.quants <- dcast(DT.group.quants, Group ~ Assay, drop = F, value.var = colnames(DT.group.quants)[5:ncol(DT.group.quants)])
  DT.group.quants <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants, by = "Group")
  fwrite(DT.group.quants, file.path(path(object), "output", "log2_group_normalised_quants.csv"))
  rm(DT.group.quants)

  # component deviations
  if (ctrl@component.deviations == T && component.model == "independent") {
    cat("[", paste0(Sys.time(), "]   component deviations..."))
    set.seed(ctrl@random.seed)
    DT.component.deviations <- component_deviations(object, summary = T, as.data.table = T)
    DT.component.deviations[, GroupComponent := paste(Group, Component, sep = "_seaMass_")]
    setcolorder(DT.component.deviations, "GroupComponent")
    DT.components <- components(object, as.data.table = T)
    DT.component.deviations <- dcast(DT.component.deviations, GroupComponent ~ Assay, drop = F, value.var = colnames(DT.component.deviations)[7:ncol(DT.component.deviations)])
    DT.component.deviations[, Group := sub("^(.*)_seaMass_.*$", "\\1", GroupComponent)]
    DT.component.deviations[, Component := sub("^.*_seaMass_(.*)$", "\\1", GroupComponent)]
    DT.component.deviations[, GroupComponent := NULL]
    DT.component.deviations <- merge(DT.components[, .(Component, nMeasurement, nDatapoint)], DT.component.deviations, by = "Component")
    setcolorder(DT.component.deviations, c("Group", "Component"))
    fwrite(DT.component.deviations, file.path(path(object), "output", "log2_component_deviations.csv"))
    rm(DT.component.deviations)
  }

  # plot PCA and assay exposures
  cat("[", paste0(Sys.time(), "]   plotting PCA and assay exposures..."))
  if (ctrl@component.deviations == T && component.model == "independent") {
    DT <- component_deviations(object, as.data.table = T)
    DT[, Group := interaction(Group, Component, sep = " : ", lex.order = T, drop = T)]
    DT.summary <- component_deviations(object, summary = T, as.data.table = T)
    DT.summary[, Group := interaction(Group, Component, sep = " : ", lex.order = T, drop = T)]
    g <- plot_pca(object, data = DT, data.summary = DT.summary)
    ggplot2::ggsave(file.path(path(object), "output", "log2_component_deviations_pca.pdf"), width = 12, height = 12, limitsize = F)
  }
  g <- plot_pca(object)
  ggplot2::ggsave(file.path(path(object), "output", "log2_group_quants_pca.pdf"), width = 12, height = 12, limitsize = F)
  g <- plot_assay_exposures(object)
  ggplot2::ggsave(file.path(path(object), "output", "log2_assay_exposures.pdf"), width = 8, height = 0.5 + 1 * nlevels(DT.design$Assay), limitsize = F)

  # differential expression analysis and false discovery rate correction
  if (ctrl@dea.model != "" && !all(is.na(DT.design$Condition))) {
    # group quants
    do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)
    if (file.exists(file.path(path(object), "de.index.fst"))) {
      if (ctrl@fdr.model != "") {
        do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
        DTs.fdr <- split(group_fdr(object, as.data.table = T), drop = T, by = "Batch")
        for (name in names(DTs.fdr)) {
          # save pretty version
          DT.fdr <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DTs.fdr[[name]], by = "Group")
          DT.fdr[, Batch := NULL]
          setcolorder(DT.fdr, c("Effect", "Model"))
          setorder(DT.fdr, qvalue, na.last = T)
          fwrite(DT.fdr, file.path(path(object), "output", paste("log2_group_de", gsub("\\s", "_", name), "csv", sep = ".")))
          # plot fdr
          g <- plot_fdr(DT.fdr, 1.0)
          ggplot2::ggsave(file.path(path(object), "output", paste("log2_group_de", gsub("\\s", "_", name), "pdf", sep = ".")), g, width = 8, height = 8)
        }
      } else {
        group_de(object, summary = T, as.data.table = T)
      }
    }

    # component deviations
    if (ctrl@component.deviations == T && component.model == "independent") {
      ellipsis$input <- "standardised.component.deviations"
      ellipsis$output <- "de.component.deviations"
      ellipsis$type <- "component.deviations"
      do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)
      if (file.exists(file.path(path(object), "de.component.deviations.index.fst"))) {
        if (ctrl@fdr.model != "") {
          ellipsis$input <- "de.component.deviations"
          ellipsis$output <- "fdr.component.deviations"
          do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
          DTs.fdr <- split(component_deviations_fdr(object, as.data.table = T), drop = T, by = "Batch")
          for (name in names(DTs.fdr)) {
            # save pretty version
            DT.fdr <- merge(DT.components, DTs.fdr[[name]], by = "Component")
            DT.fdr[, Batch := NULL]
            setcolorder(DT.fdr, c("Effect", "Model", "Group"))
            setorder(DT.fdr, qvalue, na.last = T)
            fwrite(DT.fdr, file.path(path(object), "output", paste("log2_component_deviations_de", gsub("\\s", "_", name), "csv", sep = ".")))
            # plot fdr
            g <- plot_fdr(DT.fdr, 1.0)
            ggplot2::ggsave(file.path(path(object), "output", paste("log2_component_deviations_de", gsub("\\s", "_", name), "pdf", sep = ".")), g, width = 8, height = 8)
          }
        } else {
          component_deviations_de(object, summary = T, as.data.table = T)
        }
      }
    }
  }

  # set complete
  write.table(data.frame(), file.path(path(object), "complete"), col.names = F)

  return(object@name)
})



#' @describeIn seaMass_delta-class Get name.
#' @export
setMethod("name", "seaMass_delta", function(object) {
  return(sub("^delta\\.", "", object@name))
})


#' @describeIn seaMass_delta-class Get path.
#' @export
setMethod("path", "seaMass_delta", function(object) {
  return(file.path(object@sigma@path, object@name))
})


#' @describeIn seaMass_delta-class Get the \link{sigma_control}.
#' @export
setMethod("control", "seaMass_delta", function(object) {
  if (!file.exists(file.path(path(object), "meta", "control.rds")))
    stop(paste0("seaMass-Δ output '", name(object), "' is missing or zipped"))

  return(readRDS(file.path(path(object), "meta", "control.rds")))
})


#' @describeIn seaMass_delta-class Get the study design as a \code{data.frame}.
#' @export
setMethod("assay_design", "seaMass_delta", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(path(object), "meta", "design.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Delete the \code{seaMass_delta} run from disk.
#' @import data.table
#' @export
setMethod("del", "seaMass_delta", function(object) {
  return(unlink(path(object), recursive = T))
})


#' @describeIn seaMass_delta-class Get the component metadata as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("components", "seaMass_delta", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(path(object), "meta", "components.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the group metadata as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("groups", "seaMass_delta", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(path(object), "meta", "groups.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model unnormalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("unnormalised_group_quants", "seaMass_delta", function(object, groups = NULL, summary = FALSE, input = "standardised", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, input, "Group", "Group", c("Group", "Assay"), groups, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model normalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("normalised_group_quants", "seaMass_delta", function(object, groups = NULL, summary = FALSE, input = "normalised", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if (!dir.exists(file.path(path(object), input))) {
    DT <- unnormalised_group_quants(object, groups, summary, chains = chains, as.data.table = T)
    DT[, exposure := 0]
  } else {
    if(is.null(summary) || summary == F) summary <- NULL
    if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

    DT <- read_mcmc(object, input, "Group", "Group", c("Group", "Assay"), groups, ".", chains, summary)
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the input standardised component deviations as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("component_deviations", "seaMass_delta", function(object, components = NULL, summary = FALSE, input = "standardised.component.deviations", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, input, "Component", c("Group", "Component"), c("Group", "Component", "Assay"), components, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})



#' @describeIn seaMass_delta-class Get the model normalised group variances as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("group_variances", "seaMass_delta", function(object, groups = NULL, summary = FALSE, input = "normalised.group.variances", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_invchisq_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, input, "Group", "Group", "Group", groups, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model differential expression as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("group_de", "seaMass_delta", function(object, groups = NULL, summary = FALSE, input = "de", chains = 1:control(object)@model.nchain, as.data.table = FALSE
) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, input, "Group", "Group", c("Group", "Effect", "Model"), groups, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("component_deviations_de", "seaMass_delta", function(object, components = NULL, summary = FALSE, input = "de.component.deviations", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, input, "Component", c("Group", "Component"), c("Group", "Component", "Model", "Effect"), components, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the false discovery rate corrected group differential expression as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("group_fdr", "seaMass_delta", function(object, input = "fdr", as.data.table = FALSE) {
  if (file.exists(file.path(path(object), paste(input, "fst", sep = ".")))) {
    return(fst::read.fst(file.path(path(object), paste(input, "fst", sep = ".")), as.data.table = as.data.table))
  } else {
    return(group_de(object, summary = T, as.data.table = as.data.table))
  }
})


#' @describeIn seaMass_delta-class Get the false discovery rate corrected component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
setMethod("component_deviations_fdr", "seaMass_delta", function(object, input = "fdr.component.deviations", as.data.table = FALSE) {
  if (file.exists(file.path(path(object), paste(input, "fst", sep = ".")))) {
    return(fst::read.fst(file.path(path(object), paste(input, "fst", sep = ".")), as.data.table = as.data.table))
  } else {
    return(component_deviations_de(object, summary = T, as.data.table = as.data.table))
  }
})


