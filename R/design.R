#' Imported dataset run information
#'
#' Get or set run and block information from a \link{data.frame} returned by an \code{import} routine.
#' Used to manually add run and block information to imported datasets.
#'
#' @param data \link{data.frame} returned by an \code{import} routine.
#' @import data.table
#' @export
runs <- function(data) {
  data.is.data.table <- is.data.table(data)
  DT <- setDT(data)

  DT <- DT[, .(Run = first(Run)), keyby = Injection]
  DT[, Run := as.character(Run)]
  DT[, Injection := as.character(Injection)]

  if (!data.is.data.table) setDF(data)
  setDF(DT)
  return(DT[])
}


#' @rdname runs
#' @param value Edited \link{data.frame} created by \code{runs}.
#' @import data.table
#' @export
`runs<-` <- function(data, value) {
  data.is.data.table <- is.data.table(data)
  DT <- setDT(data)

  DT[, Run := NULL]
  DT.runs <- as.data.table(value)
  DT.runs <- DT.runs[complete.cases(DT.runs)]
  if (!is.factor(DT.runs$Injection)) DT.runs[, Injection := factor(Injection)]
  if (!is.factor(DT.runs$Run)) DT.runs[, Run := factor(Run)]
  DT <- merge(DT, DT.runs, by = "Injection")

  if (!data.is.data.table) setDF(data)
  setDF(DT)
  return(DT[])
}


#' Specify study design for imported dataset
#'
#' Returns a skeleton \link{data.frame} for customising the study design of an imported dataset.
#'
#' @param data \link{data.frame} returned by  \link{import_ProteinPilot}, \link{import_ProteinPilot},
#' \link{import_ProteomeDiscovery}, \link{import_Progenesis} or \link{import_OpenSWATH}.
#' @return \link{data.frame} that can be edited and passed as parameter \code{data.design} of \link{seaMass_sigma}.
#' @import data.table
#' @export
new_assay_design <- function(data) {
  data.is.data.table <- is.data.table(data)
  setDT(data)

  DT.design <- data[, .(Assay = paste(Run, Channel, sep = ",")), keyby = .(Run, Channel)]
  DT.design[, Assay := sub("^0,", "", Assay)]
  DT.design[, Assay := sub(",0$", "", Assay)]

  # autodetect blocks
  Block. <- merge(DT.design[, .(Run, Channel, Assay)], data[, .(Run, Channel, Measurement)], by = c("Run", "Channel"))
  Block.[, Run := NULL]
  Block.[, Channel := NULL]
  Block.[, N := 1]
  Block. <- dcast(Block., Measurement ~ Assay, sum, value.var = "N")
  Block.[, Measurement := NULL]
  Block. <- as.matrix(Block.)
  # matrix multiplication distributes assay relationships
  Block. <- t(Block.) %*% Block.
  # Block. is recoded first non-zero occurence for each assay
  Block. <- colnames(Block.)[apply(Block. != 0, 2, which.max)]
  Block. <- factor(as.integer(factor(Block., levels = unique(Block.))))
  if (nlevels(Block.) == 1) {
    DT.design[, paste("Block", levels(Block.), sep = ".")] <- as.logical(as.integer(Block.))
  } else if (nlevels(Block.) > 1) {
    Block. <- data.table(model.matrix(~ Block. - 1))
    DT.design <- cbind(DT.design, mapply(Block., FUN = as.logical))
  }

  # For seaMass-Δ
  if ("RefWeight" %in% colnames(data)) {
    DT.design <- merge(DT.design, unique(data[, .(Run, Channel, RefWeight)]), by = c("Run", "Channel"), sort = F)
  } else {
    DT.design[, RefWeight := 1]
  }

  if ("Sample" %in% colnames(data)) {
    DT.design <- merge(DT.design, unique(data[, .(Run, Channel, Sample)]), by = c("Run", "Channel"), sort = F)
  } else {
    DT.design[, Sample := Assay]
  }

  if ("Condition" %in% colnames(data)) {
    DT.design <- merge(DT.design, unique(data[, .(Run, Channel, Condition)]), by = c("Run", "Channel"), sort = F)
  } else {
    DT.design[, Condition := NA]
  }

  if (!data.is.data.table) setDF(data)
  setDF(DT.design)
  return(DT.design[])
}
