.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("seaMass v", packageVersion("seaMass"), "  |  © 2019-2020  BIOSP", utf8::utf8_encode("\U0001f441"), "  Laboratory"))
  packageStartupMessage("This program comes with ABSOLUTELY NO WARRANTY.")
  packageStartupMessage("This is free software, and you are welcome to redistribute it under certain conditions.")
}


#' seaMass object
#'
#' Methods shared between \link{sigma_fit} and \link{seaMass_delta}
setClass("seaMass", contains = "VIRTUAL")


#' @describeIn seaMass_delta-class Get the model standardised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("standardised_group_quants", "seaMass", function(object, groups = NULL, summary = FALSE, chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  if (dir.exists(file.path(filepath(object), "standardised.group.quants"))) {
    DT <- read_mcmc(object, "standardised.group.quants", "Group", "Group", c("Group", "Assay"), groups, ".", chains, summary)
  } else {
    DT <- read_mcmc(fits(object)[[1]], "standardised.group.quants", "Group", "Group", c("Group", "Assay"), groups, ".", chains, summary)
  }
  if (is.null(DT)) stop(paste("standardised group quants were not", ifelse(is.null(summary), "kept", "summarised")))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Get the model normalised group variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_variances", "seaMass", function(object, groups = NULL, summary = FALSE, chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_invchisq_mcmc", paste("dist", summary, sep = "_"))

  if (dir.exists(file.path(filepath(object), "normalised.group.variances"))) {
    DT <- read_mcmc(object, "normalised.group.variances", "Group", "Group", "Group", groups, ".", chains, summary)
  } else {
    DT <- read_mcmc(fits(object)[[1]], "normalised.group.variances", "Group", "Group", "Group", groups, ".", chains, summary)
  }
  if (is.null(DT)) stop(paste("normalised group variances were not", ifelse(is.null(summary), "kept", "summarised")))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass-class Get the model normalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_quants", "seaMass", function(object, groups = NULL, summary = FALSE, chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if (dir.exists(file.path(filepath(object), "normalised.group.quants"))) {
    if(is.null(summary) || summary == F) summary <- NULL
    if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))
    DT <- read_mcmc(object, "normalised.group.quants", "Group", "Group", c("Group", "Assay"), groups, ".", chains, summary)
  } else {
    if (dir.exists(file.path(filepath(object), "standardised_group_quants"))) {
      DT <- standardised_group_quants(object, groups, summary, chains = chains, as.data.table = T)
      DT[, exposure := 0]
    } else {
      if (dir.exists(file.path(filepath(fits(object)[[1]]), "normalised.group.quants"))) {
        if(is.null(summary) || summary == F) summary <- NULL
        if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))
        DT <- read_mcmc(fits(object)[[1]], "normalised.group.quants", "Group", "Group", c("Group", "Assay"), groups, ".", chains, summary)
      } else {
        DT <- standardised_group_quants(fits(object)[[1]], groups, summary, chains = chains, as.data.table = T)
        DT[, exposure := 0]
      }
    }
  }
  if (is.null(DT)) stop(paste("normalised group variances were not", ifelse(is.null(summary), "kept", "summarised")))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @include generics.R
setMethod("read_mcmc", "seaMass", function(object, effect.name, columnID, batchIDs, summaryIDs, itemIDs, input, chains, summary) {
  if (!is.null(summary)) filename <- file.path(file.path(filepath(object), input, paste(summary, effect.name, "fst", sep = ".")))
  if (!is.null(summary) && file.exists(filename)) {
    # load and filter from cache
    DT <- fst::read.fst(filename, as.data.table = T)
    if (!is.null(itemIDs)) DT <- DT[get(columnID) %in% itemIDs]
  } else {
    # load and filter index
    filename.index <- file.path(filepath(object), input, paste(effect.name, "index.fst", sep = "."))
    if (!file.exists(filename.index)) return(NULL)
    DT.index <- fst::read.fst(filename.index, as.data.table = T)
    if (!is.null(itemIDs)) DT.index <- DT.index[get(columnID) %in% itemIDs]
    DT.index <- DT.index[complete.cases(DT.index)]
    if (nrow(DT.index) == 0) return(NULL)
    setorder(DT.index, file, from)

    # read
    ctrl <- control(object)
    if (is.null(summary)) {
      items <- list(DT.index)
    } else {
      items <- batch_split(DT.index, batchIDs, 16)
    }
    DT <- rbindlist(parallel_lapply(items, function(item, object, input, chains, summary, summaryIDs) {
      # minimise file access
      item[, file.prev := shift(file, fill = "")]
      item[, to.prev := shift(to + 1, fill = 0)]
      item[, file.next := shift(file, fill = "", -1)]
      item[, from.next := shift(from - 1, fill = 0, -1)]
      item <- cbind(
        item[!(file == file.prev & from == to.prev), .(file, from)],
        item[!(file == file.next & to == from.next), .(to)]
      )

      # read
      DT <- rbindlist(lapply(1:nrow(item), function(i) {
        rbindlist(lapply(chains, function(chain) {
          DT1 <- NULL
          try({
            DT1 <- fst::read.fst(
              file.path(filepath(object), input, dirname(item[i, file]), sub("^([0-9]+)", chain, basename(item[i, file]))),
              from = item[i, from],
              to = item[i, to],
              as.data.table = T
            )}, silent = T)
          return(DT1)
        }))
      }))

      # optional summarise
      if (!is.null(summary) && nrow(DT) > 0) {
        # average samples if assay run in multiple blocks
        if (length(unique(DT$BlockID)) > 1) DT <- DT[, .(value = mean(value)), by = c(summaryIDs, "chainID", "mcmcID")]
        DT.summary <- DT[, do.call(summary, list(chainID = chainID, mcmcID = mcmcID, value = value)), by = summaryIDs]
        # add back metadata
        DT <- DT[, c(summaryIDs, colnames(DT)[!(colnames(DT) %in% c(colnames(DT.summary), "chainID", "mcmcID", "value"))]), with = F]
        DT <- merge(DT[!duplicated(DT[, summaryIDs, with = F])], DT.summary, by = summaryIDs)
      }

      setcolorder(DT, summaryIDs)
      return(DT)
    }, nthread = ifelse(length(items) == 1, 1, ctrl@nthread)))

    # cache results
    if (!is.null(summary) && is.null(itemIDs) && identical(chains, 1:ctrl@model.nchain)) {
      fst::write.fst(DT, filename)
    }
  }

  return(DT)
})

