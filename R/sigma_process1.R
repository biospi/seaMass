setGeneric("process1", function(object, ...) standardGeneric("process1"))


#' @import data.table
#' @include generics.R
setMethod("process1", "sigma_fit", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", path(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # EXECUTE MODEL
  model(object, "model1", chain)

  if (all(sapply(1:ctrl@model.nchain, function(chain) file.exists(file.path(object@path, "model1", paste(".complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    cat(paste0("[", Sys.time(), "]   OUTPUT1 block=", sub("^.*sigma\\.(.*)$", "\\1", object@path), "\n"))

    # load parameters
    DT.design <- fst::read.fst(file.path(object@path, "meta", "design.fst"), as.data.table = T)
    DT.groups <- fst::read.fst(file.path(object@path, "meta", "groups.fst"), as.data.table = T)
    DT.components <- fst::read.fst(file.path(object@path, "meta", "components.fst"), as.data.table = T)
    DT.measurements <- fst::read.fst(file.path(object@path, "meta", "measurements.fst"), as.data.table = T)

    # timings
    DT.timings <- timings(object, as.data.table = T)
    DT.timings <- data.table::dcast(DT.timings, GroupID ~ chainID, value.var = "elapsed")
    DT.timings <- merge(DT.groups[, .(GroupID, Group, nComponent, nMeasurement, nDatapoint, pred = timing)], DT.timings, by = "GroupID")
    DT.timings[, GroupID := NULL]
    fwrite(DT.timings, file.path(dirname(object@path), "output", basename(object@path), "group_timings.csv"))
    rm(DT.timings)

    # measurement var summary
    if ("measurement.variances" %in% ctrl@summarise) {
      # measurement variances
      cat(paste0("[", Sys.time(), "]    getting measurement variance summaries...\n"))
      DT.measurement.variances <- measurement_variances(object, summary = T, as.data.table = T)
    }
    # delete if not in 'keep'
    if (!("measurement.variances" %in% ctrl@keep)) unlink(file.path(object@path, "model1", "measurement.variances*"), recursive = T)

    # component variances summary
    if("component.variances" %in% ctrl@summarise && !is.null(ctrl@component.model)) {
      cat(paste0("[", Sys.time(), "]    getting component variance summaries...\n"))
      DT.component.variances <- component_variances(object, summary = T, as.data.table = T)
    }
    # delete if not in 'keep'
    if (!("component.variances" %in% ctrl@keep)) unlink(file.path(object@path, "model1", "component.variances*"), recursive = T)

    # component deviations summary
    if (!is.null(ctrl@component.model) && ctrl@component.model == "independent") {
      if ("component.deviations" %in% ctrl@summarise || "component.deviations.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    getting component deviation summaries...\n"))
        DT.component.deviations <- component_deviations(object, summary = T, as.data.table = T)
      }

      if("component.deviations" %in% ctrl@summarise) {
        DT.component.deviations <- merge(DT.groups[, .(GroupID, Group)], DT.component.deviations, by = "Group")
        DT.component.deviations <- merge(DT.components[, .(ComponentID, Component)], DT.component.deviations, by = "Component")
        DT.component.deviations[, GroupIDComponentID := paste(GroupID, ComponentID, sep = "_")]
        setcolorder(DT.component.deviations, "GroupIDComponentID")
        DT.component.deviations <- dcast(DT.component.deviations, GroupIDComponentID ~ Assay, drop = F, value.var = colnames(DT.component.deviations)[7:ncol(DT.component.deviations)])
        DT.component.deviations[, GroupID := as.integer(sub("^([0-9]+)_[0-9]+$", "\\1", GroupIDComponentID))]
        DT.component.deviations[, ComponentID := as.integer(sub("^[0-9]+_([0-9]+)$", "\\1", GroupIDComponentID))]
        DT.component.deviations[, GroupIDComponentID := NULL]
        DT.component.deviations <- merge(DT.components[, .(ComponentID, Component, nMeasurement, nDatapoint)], DT.component.deviations, by = "ComponentID")
        DT.component.deviations <- merge(DT.groups[, .(GroupID, Group)], DT.component.deviations, by = "GroupID")
        setcolorder(DT.component.deviations, c("GroupID", "ComponentID"))
        DT.component.deviations[, ComponentID := NULL]
        DT.component.deviations[, GroupID := NULL]
        fwrite(DT.component.deviations, file.path(dirname(object@path), "output", basename(object@path), "log2_component_deviations.csv"))
      }

      if ("component.deviations.pca" %in% ctrl@plot) {
        ellipsis <- ctrl@ellipsis
        ellipsis$object <- object
        ellipsis$type <- "component.deviations"
        do.call("plot_pca", ellipsis)
        ellipsis$type <- NULL
        ggplot2::ggsave(file.path(dirname(object@path), "output", basename(object@path), "log2_component_deviations_pca.pdf"), width = 300, height = 300 * 9/16, units = "mm")
      }

      rm(DT.component.deviations)
    }
    # delete if not in 'keep'
    if (!("component.deviations" %in% ctrl@keep)) unlink(file.path(object@path, "model1", "component.deviations*"), recursive = T)

    # assay deviations summary
    if("assay.deviations" %in% ctrl@summarise && !is.null(ctrl@assay.model) && ctrl@assay.model == "component") {
      cat(paste0("[", Sys.time(), "]    getting assay deviation summaries...\n"))
      DT.assay.deviations <- assay_deviations(object, summary = T, as.data.table = T)
      DT.assay.deviations <- merge(DT.groups[, .(GroupID, Group)], DT.assay.deviations, by = "Group")
      DT.assay.deviations <- merge(DT.components[, .(ComponentID, Component)], DT.assay.deviations, by = "Component")
      DT.assay.deviations[, GroupIDComponentID := paste(GroupID, ComponentID, sep = "_")]
      setcolorder(DT.assay.deviations, "GroupIDComponentID")
      DT.assay.deviations <- dcast(DT.assay.deviations, GroupIDComponentID ~ Assay, drop = F, value.var = colnames(DT.assay.deviations)[7:ncol(DT.assay.deviations)])
      DT.assay.deviations[, GroupID := as.integer(sub("^([0-9]+)_[0-9]+$", "\\1", GroupIDComponentID))]
      DT.assay.deviations[, ComponentID := as.integer(sub("^[0-9]+_([0-9]+)$", "\\1", GroupIDComponentID))]
      DT.assay.deviations[, GroupIDComponentID := NULL]
      DT.assay.deviations <- merge(DT.components[, .(ComponentID, Component, nMeasurement, nDatapoint)], DT.assay.deviations, by = "ComponentID")
      DT.assay.deviations <- merge(DT.groups[, .(GroupID, Group)], DT.assay.deviations, by = "GroupID")
      setcolorder(DT.assay.deviations, c("GroupID", "ComponentID"))
      DT.assay.deviations[, ComponentID := NULL]
      DT.assay.deviations[, GroupID := NULL]
      fwrite(DT.assay.deviations, file.path(dirname(object@path), "output", basename(object@path), "log2_assay_deviations.csv"))
      rm(DT.assay.deviations)
    }
    # delete if not in 'keep'
    if (!("assay.deviations" %in% ctrl@keep)) unlink(file.path(object@path, "model1", "assay.deviations*"), recursive = T)

    # unnormalised group quants summary
    if ("unnormalised.group.quants" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]    getting unnormalised group quant summaries...\n"))
      DT.group.quants <- unnormalised_group_quants(object, summary = T, as.data.table = T)
      DT.group.quants.out <- dcast(DT.group.quants, Group ~ Assay, drop = F, value.var = colnames(DT.group.quants)[6:ncol(DT.group.quants)])
      DT.group.quants.out <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants.out, by = "Group")
      fwrite(DT.group.quants.out, file.path(dirname(object@path), "output", basename(object@path), "log2_unnormalised_group_quants.csv"))
      rm(DT.group.quants.out)
    }

    # normalised group quants
    if ("normalised.group.variances" %in% ctrl@summarise || "normalised.group.variances" %in% ctrl@keep || "normalised.group.quants" %in% ctrl@summarise || "normalised.group.quants" %in% ctrl@keep || "normalised.group.quants.pca" %in% ctrl@plot) {
      # normalise group quants
      ellipsis <- ctrl@ellipsis
      ellipsis$object <- object
      ellipsis$input <- "model1"
      do.call(paste("norm", ctrl@norm.model, sep = "_"), ellipsis)
      ellipsis$input <- NULL

      if ("normalised.group.variances" %in% ctrl@summarise) {
        cat(paste0("[", Sys.time(), "]   getting normalised group variance summaries...\n"))
        DT.groups <- groups(object, as.data.table = T)
        DT.group.variances <- normalised_group_variances(object, summary = T, as.data.table = T)
        setcolorder(DT.group.variances, "Group")
        fwrite(DT.group.variances, file.path(dirname(object@path), "output", basename(object@path), "log2_normalised_group_variances.csv"))
        rm(DT.group.variances)
      }

      if ("normalised.group.quants" %in% ctrl@summarise || "normalised.group.quants.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]   getting group quant summaries...\n"))
        DT.groups <- groups(object, as.data.table = T)
        DT.group.quants <- normalised_group_quants(object, summary = T, as.data.table = T)

        if ("normalised.group.quants" %in% ctrl@summarise) {
          DT.group.quants <- dcast(DT.group.quants, Group ~ Assay, drop = F, value.var = colnames(DT.group.quants)[6:ncol(DT.group.quants)])
          DT.group.quants <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants, by = "Group")
          fwrite(DT.group.quants, file.path(dirname(object@path), "output", basename(object@path), "log2_normalised_group_quants.csv"))
          rm(DT.group.quants)
        }

        if ("normalised.group.quants.pca" %in% ctrl@plot) {
          do.call("plot_pca", ellipsis)
          ggplot2::ggsave(file.path(dirname(object@path), "output", basename(object@path), "log2_normalised_group_quants_pca.pdf"), width = 300, height = 300 * 9/16, units = "mm")
        }
      }
    }

    # delete group quants if not in 'keep'
    if (!("unnormalised.group.quants" %in% ctrl@keep)) unlink(file.path(object@path, "model1", "group.quants*"), recursive = T)
    if (!("normalised.group.quants" %in% ctrl@keep)) unlink(file.path(object@path, "normalised.group.quants*"), recursive = T)
    if (!("normalised.group.variances" %in% ctrl@keep)) {
      if(file.exists(file.path(object@path, "normalised.group.variances"))) unlink(file.path(object@path, "normalised.group.variances*"), recursive = T)
    }

    # set complete
    write.table(data.frame(), file.path(object@path, "complete"), col.names = F)
  }

  return(invisible(NULL))
})

