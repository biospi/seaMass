#' @import data.table
#' @include generics.R
#' @include seaMass_delta.R
setMethod("process", "seaMass_delta", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", name(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  ellipsis <- ctrl@ellipsis
  ellipsis$object <- object
  ellipsis$chains <- chain

  # group dea
  if (ctrl@dea.model != "" && !all(is.na(assay_design(object, as.data.table = T)$Condition))) {
    if (is.null(normalised_group_quants(object@fit, groups(object@fit, as.data.table = T)[1, Group]))) {
      type <- "standardised.group.quants"
    } else {
      type <- "normalised.group.quants"
    }
    ellipsis$type <- type
    do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)
  }

  # component deviation dea
  if (ctrl@component.deviations == T && control(object@fit)@component.model == "independent") {
    ellipsis$type <- "component.deviations"
    do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)
  }

  write.table(data.frame(), file.path(filepath(object), paste("complete", chain, sep = ".")), col.names = F)

  if (all(sapply(1:ctrl@model.nchain, function(chain) file.exists(file.path(filepath(object), paste("complete", chain, sep = ".")))))) {
    # summarise group de and perform fdr correction
    if (file.exists(file.path(filepath(object), paste0(type, ".index.fst")))) {
      if(ctrl@fdr.model != "") {
        ellipsis$type <- type
        do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
      } else {
        read_samples(object, ".", type, summary = T, as.data.table = T)
      }
    }
    if (!(type %in% ctrl@keep)) unlink(file.path(filepath(object), paste0(type, "*")), recursive = T)

    # write out group fdr
    if (file.exists(file.path(filepath(object), paste0(type, ".fst")))) {
      DTs.fdr <- fst::read.fst(file.path(filepath(object), paste0(type, ".fst")), as.data.table = T)
      DTs.fdr <- split(DTs.fdr, drop = T, by = "Batch")
      for (name in names(DTs.fdr)) {
        # save
        fwrite(DTs.fdr[[name]], file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_", type, "_fdr", gsub("\\s", "_", name), "csv", sep = ".")))
        # plot fdr
        plot_fdr(DTs.fdr[[name]])
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2", type, "_fdr", gsub("\\s", "_", name), "pdf", sep = ".")), width = 8, height = 8)
        # plot volcano
        plot_volcano(DTs.fdr[[name]])
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2", type, "_fdr_volcano", gsub("\\s", "_", name), "pdf", sep = ".")), width = 8, height = 8)
        # plot fc_vs_stdev
        plot_volcano(DTs.fdr[[name]], y.col = "s")
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2", type, "_fdr_fc", gsub("\\s", "_", name), "pdf", sep = ".")), width = 8, height = 8)
      }
    }

    if (ctrl@component.deviations == T && control(object@fit)@component.model == "independent") {
      # summarise component deviation de and perform fdr correction
      if (file.exists(file.path(filepath(object), "component.deviations.index.fst"))) {
        if(ctrl@fdr.model != "") {
          ellipsis$type <- "component.deviations"
          do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
        } else {
          component_deviations_de(object, summary = T, as.data.table = T)
        }
      }
      if (!("component.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "component.deviations*"), recursive = T)

      # write out group fdr
      if (file.exists(file.path(filepath(object), "component.deviations.fst"))) {
        DTs.fdr <- split(component_deviations_fdr(object, as.data.table = T), drop = T, by = "Batch")
        for (name in names(DTs.fdr)) {
          # save
          fwrite(DTs.fdr[[name]], file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_component_deviations_fdr", gsub("\\s", "_", name), "csv", sep = ".")))
          # plot fdr
          plot_fdr(DTs.fdr[[name]])
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_component_deviations_fdr", gsub("\\s", "_", name), "pdf", sep = ".")), width = 8, height = 8)
          # plot volcano
          plot_volcano(DTs.fdr[[name]])
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_component_deviations_fdr_volcano", gsub("\\s", "_", name), "pdf", sep = ".")), width = 8, height = 8)
          # plot fc_vs_stdev
          plot_volcano(DTs.fdr[[name]], y.col = "s")
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste("log2_component_deviations_fdr_fc", gsub("\\s", "_", name), "pdf", sep = ".")), width = 8, height = 8)
        }
      }
    }

    # set complete
    write.table(data.frame(), file.path(filepath(object), "complete"), col.names = F)
    cat(paste0("[", Sys.time(), "] seaMass-delta finished!\n"))
  }

  return(invisible(NULL))
})


