setGeneric("plots", function(object, ...) standardGeneric("plots"))


#' @import data.table
#' @include generics.R
setMethod("plots", "sigma_fit", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", name(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # ctrl <- ctrl(object)
  # cat(paste0("[", Sys.time(), "] PLOTS set=", i, "/", ctrl@assay.nblock * ctrl@model.nchain, "\n"))
  #
  # # create subdirs
  # dir.create(file.path(object, "plots", "measurements"), showWarnings = F)
  # dir.create(file.path(object, "plots", "components"), showWarnings = F)
  #
  # DT.groups <- groups(object, as.data.table = T)
  # DT.design <- assay_design(object, as.data.table = T)
  # DT.group.quants <- group_quants(object, as.data.table = T)
  # pids <- levels(DT.group.quants$GroupID)
  # pids <- pids[seq(chain, length(pids), ctrl@model.nchain)]
  #
  # # start cluster and reproducible seed
  # pb <- txtProgressBar(max = length(pids), style = 3)
  # dfll <- foreach(pid = pids, .packages = c("seamassdelta", "data.table", "ggplot2"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
  #   plt.measurements <- plot_measurements(object, groupID = pid)
  #   ggplot2::ggsave(file.path(object, "plots", "measurements", paste0(pid, ".pdf")), plt.measurements$g, width = plt.measurements$width, height = plt.measurements$height, limitsize = F)
  #
  #   plt.components <- plot_components(object, groupID = pid)
  #   ggplot2::ggsave(file.path(object, "plots", "components", paste0(pid, ".pdf")), plt.components$g, width = plt.components$width, height = plt.components$height, limitsize = F)
  # }
  # setTxtProgressBar(pb, length(pids))

  return(invisible(NULL))
})

