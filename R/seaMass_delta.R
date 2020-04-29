#' seaMass-Δ
#'
#' Perform seaMass-Δ normalisation, differential expression analysis and/or false discovery rate correction on raw group-level
#' quants output by seaMass-Σ
#' @include seaMass.R
#' @include seaMass_sigma.R
setClass("seaMass_delta", contains = "seaMass", slots = c(
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
    message(paste0("returning completed seaMass-delta object - if this wasn't your intention, supply a different 'name' or delete the folder for the returned object with 'del(object)'"))
    return(object)
  }

  ### INIT
  cat(paste0("[", Sys.time(), "] seaMass-delta v", control@version, "\n"))
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
  path <- filepath(object)
  if (file.exists(path)) unlink(path, recursive = T)
  dir.create(file.path(path, "meta"), recursive = T)
  if (file.exists(file.path(dirname(path), "output", basename(path)))) unlink(file.path(dirname(path), "output", basename(path)), recursive = T)
      dir.create(file.path(dirname(path), "output", basename(path)), recursive = T)

  # check and save control
  saveRDS(control, file.path(path, "meta", "control.rds"))

  # merged design
  DT.design <- unique(as.data.table(data.design)[!is.na(Assay)], by = "Assay")

  # get design into the format we need
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
    cat(paste0("[", Sys.time(), "] queued\n"))
  }

  return(invisible(object))
}


#' @describeIn seaMass_delta-class Open a complete \code{seaMass_delta} run from the supplied \link{seaMass_sigma} and \code{name}.
#' @export
open_seaMass_delta <- function(sigma, name = "fit", quiet = FALSE, force = FALSE) {
  path <- file.path(filepath(sigma), paste0("delta.", name))

  if(force || file.exists(file.path(path, "complete"))) {
    object <- new("seaMass_delta", sigma = sigma, name = paste0("delta.", name))
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

  DT.groups <- groups(object, as.data.table = T)
  DT.design <- assay_design(object, as.data.table = T)
  ellipsis <- ctrl@ellipsis
  ellipsis$object <- object

  cat(paste0("[", Sys.time(), "]  running delta for name=", name(object), " sigma=", name(object@sigma), "...\n"))

  if (!(length(fits(object)) == 1 && identical(assay_design(fits(object)[[1]], as.data.table = T)[, .(Assay, RefWeight)], assay_design(object, as.data.table = T)[, .(Assay, RefWeight)]))) {
    # standardise quants using reference weights
    standardise_group_quants(object)

    # normalise quants by norm.groups
    if (ctrl@norm.model != "") do.call(paste("norm", ctrl@norm.model, sep = "_"), ellipsis)
    if (!("standardised.group.quants" %in% ctrl@keep)) unlink(file.path(filepath(object), "standardised.group.quants*"), recursive = T)

    # summarise group variances
    if(file.exists(file.path(filepath(object), "normalised.group.variances"))) {
      cat(paste0("[", Sys.time(), "]   getting normalised group variance summaries...\n"))
      DT.group.variances <- normalised_group_variances(object, summary = T, as.data.table = T)
      setcolorder(DT.group.variances, "Group")
      fwrite(DT.group.variances, file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_normalised_group_variances.csv"))
      rm(DT.group.variances)
    }

    # summarise group quants
    cat(paste0("[", Sys.time(), "]   getting normalised group quant summaries...\n"))
    DT.group.quants <- normalised_group_quants(object, summary = T, as.data.table = T)
    DT.group.quants <- dcast(DT.group.quants, Group ~ Assay, drop = F, value.var = colnames(DT.group.quants)[5:ncol(DT.group.quants)])
    DT.group.quants <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants, by = "Group")
    fwrite(DT.group.quants, file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_normalised_group_quants.csv"))
    rm(DT.group.quants)

    # plot PCA
    do.call("plot_pca", ellipsis)
    ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_normalised_group_quants_pca.pdf"), width = 300, height = 300 * 9/16, units = "mm")

    # plot assay exposures
    #g <- plot_assay_exposures(object)
    #ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_assay_exposures.pdf"), width = 8, height = 0.5 + 1 * nlevels(DT.design$Assay), limitsize = F)
  }

  # group dea
  if (ctrl@dea.model != "" && !all(is.na(DT.design$Condition))) do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)
  if (!("normalised.group.quants" %in% ctrl@keep)) unlink(file.path(filepath(object), "normalised.group.quants*"), recursive = T)
  if (!("normalised.group.variances" %in% ctrl@keep)) {
    if(file.exists(file.path(filepath(object), "normalised.group.variances"))) unlink(file.path(filepath(object), "normalised.group.variances*"), recursive = T)
  }

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
      # save pretty version
      DT.fdr <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DTs.fdr[[name]], by = "Group")
      DT.fdr[, Batch := NULL]
      setcolorder(DT.fdr, c("Effect", "Model"))
      setorder(DT.fdr, qvalue, na.last = T)
      fwrite(DT.fdr, file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_group_de", gsub("\\s", "_", name), "csv", sep = ".")))
      # plot fdr
      g <- plot_fdr(DT.fdr)
      ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_group_de", gsub("\\s", "_", name), "pdf", sep = ".")), g, width = 8, height = 8)
    }
  }
  if (!("group.fdr" %in% ctrl@keep)) unlink(file.path(filepath(object), "group.fdr*"), recursive = T)


  # component deviations
  if (ctrl@component.deviations == T && control(object@sigma)@component.model == "independent") {

    # standardise component deviations using reference weights
    standardise_component_deviations(object)

    # summarise component deviations
    cat(paste0("[", Sys.time(), "]   getting component deviation summaries...\n"))
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
    fwrite(DT.component.deviations, file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_component_deviations.csv"))
    rm(DT.component.deviations)

    # plot PCA
    # if (ctrl@component.deviations == T && component.model == "independent") {
    #   DT <- component_deviations(object, as.data.table = T)
    #   DT[, Group := interaction(Group, Component, sep = " : ", lex.order = T, drop = T)]
    #   DT.summary <- component_deviations(object, summary = T, as.data.table = T)
    #   DT.summary[, Group := interaction(Group, Component, sep = " : ", lex.order = T, drop = T)]
    #   g <- plot_pca(object, data = DT, data.summary = DT.summary)
    #   ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_component_deviations_pca.pdf"), width = 12, height = 12, limitsize = F)
    # }

    # component deviation dea
    ellipsis$input <- "component.deviations"
    ellipsis$output <- "component.deviations.de"
    ellipsis$type <- "component.deviations"
    do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)
    if (!("component.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "component.deviations*"), recursive = T)

    # summarise component deviation de and perform fdr correction
    if (file.exists(file.path(filepath(object), "component.deviations.de.index.fst"))) {
      if(ctrl@fdr.model != "") {
        ellipsis$input <- "component.deviations.de"
        ellipsis$output <- "component.deviations.fdr"
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
        # save pretty version
        DT.fdr <- merge(DT.components, DTs.fdr[[name]], by = "Component")
        DT.fdr[, Batch := NULL]
        setcolorder(DT.fdr, c("Effect", "Model", "Group"))
        setorder(DT.fdr, qvalue, na.last = T)
        fwrite(DT.fdr, file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_component_deviations_de", gsub("\\s", "_", name), "csv", sep = ".")))
        # plot fdr
        g <- plot_fdr(DT.fdr, 1.0)
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_component_deviations_de", gsub("\\s", "_", name), "pdf", sep = ".")), g, width = 8, height = 8)
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
  return(file.path(filepath(object@sigma), object@name))
})


#' @describeIn seaMass_delta-class Get the \link{sigma_control}.
#' @export
#' @include generics.R
setMethod("control", "seaMass_delta", function(object) {
  if (!file.exists(file.path(filepath(object), "meta", "control.rds")))
    stop(paste0("seaMass-delta output '", name(object), "' is missing or zipped"))

  return(readRDS(file.path(filepath(object), "meta", "control.rds")))
})


#' @describeIn seaMass_sigma-class Get the block name and \code{sigma_fit} object as a named list.
#' @export
#' @include generics.R
setMethod("fits", "seaMass_delta", function(object) {
  return(fits(object@sigma))
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


#' @describeIn seaMass_delta-class Get the component metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("components", "seaMass_delta", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "components.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the group metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("groups", "seaMass_delta", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the input standardised component deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations", "seaMass_delta", function(object, components = NULL, summary = FALSE, input = "component.deviations", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, input, "Component", c("Group", "Component"), c("Group", "Component", "Assay"), components, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_de", "seaMass_delta", function(object, groups = NULL, summary = FALSE, input = "group.de", chains = 1:control(object)@model.nchain, as.data.table = FALSE
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
#' @include generics.R
setMethod("component_deviations_de", "seaMass_delta", function(object, components = NULL, summary = FALSE, input = "component.deviations.de", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
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
#' @include generics.R
setMethod("group_fdr", "seaMass_delta", function(object, input = "group.fdr", as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), paste(input, "fst", sep = ".")))) {
    return(fst::read.fst(file.path(filepath(object), paste(input, "fst", sep = ".")), as.data.table = as.data.table))
  } else {
    return(group_de(object, summary = T, as.data.table = as.data.table))
  }
})


#' @describeIn seaMass_delta-class Get the false discovery rate corrected component deviations differential expression as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations_fdr", "seaMass_delta", function(object, input = "component.deviations.fdr", as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), paste(input, "fst", sep = ".")))) {
    return(fst::read.fst(file.path(filepath(object), paste(input, "fst", sep = ".")), as.data.table = as.data.table))
  } else {
    return(component_deviations_de(object, summary = T, as.data.table = as.data.table))
  }
})


