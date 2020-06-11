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
  DT <- merge(DT, DT.runs, by = "Injection", sort = F)

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

  DT.design <- data[, .(Run_ = Run, Channel_ = Channel), keyby = .(Run, Channel)]
  if (length(unique(DT.design$Run_)) == 1) DT.design[, Run_ := ""]
  if (length(unique(DT.design$Channel_)) == 1) DT.design[, Channel_ := ""]
  DT.design[, Assay := paste(Run_, Channel_, sep = ",")]
  DT.design[, Assay := sub("^,", "", Assay)]
  DT.design[, Assay := sub(",$", "", Assay)]
  DT.design[, Run_ := NULL]
  DT.design[, Channel_ := NULL]

  # autodetect blocks
  DT.blocks <- merge(DT.design[, .(Run, Channel, Assay)], data[, .(Run, Channel, Measurement)], by = c("Run", "Channel"))
  DT.blocks[, Run := NULL]
  DT.blocks[, Channel := NULL]
  DT.blocks[, N := 1]
  DT.blocks <- dcast(DT.blocks, Measurement ~ Assay, sum, value.var = "N")
  setcolorder(DT.blocks, DT.design$Assay)
  setDF(DT.blocks)
  rownames(DT.blocks) <- paste0("_seaMass_.Measurement.", DT.blocks$Measurement)
  DT.blocks$Measurement <- NULL
  colnames(DT.blocks) <- paste0("_seaMass_.Assay.", colnames(DT.blocks))
  DT.blocks <- igraph::components(igraph::graph.incidence(DT.blocks))
  DT.blocks <- DT.blocks$membership[grep("^_seaMass_\\.Assay\\.", names(DT.blocks$membership))]
  # igraph using lexi order, get blocks back into assay order
  DT.blocks <- data.table(Assay = sub("^_seaMass_\\.Assay\\.", "", names(DT.blocks)), Block. = factor(DT.blocks, levels = unique(DT.blocks), labels = 1:uniqueN(DT.blocks)))
  DT.design <- merge(DT.design, DT.blocks, by = "Assay", sort = F)
  if (nlevels(DT.design$Block.) == 1) {
    DT.design[,paste0("Block.", levels(DT.design$Block.)) := Block. == levels(DT.design$Block.)]
  } else {
    DT.design <- cbind(DT.design, apply(model.matrix(~ Block. - 1, data = DT.design), 2, as.logical))
  }
  DT.design[, Block. := NULL]

  if ("RefWeight" %in% colnames(data)) {
    DT.design <- merge(DT.design, unique(data[, .(Run, Channel, RefWeight)]), by = c("Run", "Channel"), sort = F)
  } else {
    DT.design[, RefWeight := 1]
  }

  # For seaMass-Δ
  if ("Sample" %in% colnames(data)) {
    DT.design <- merge(DT.design, unique(data[, .(Run, Channel, Sample)]), by = c("Run", "Channel"), sort = F)
  } else {
    DT.design[, Sample := Assay]
  }

  if ("Condition" %in% colnames(data)) {
    DT.design <- merge(DT.design, unique(data[, .(Run, Channel, Condition)]), by = c("Run", "Channel"), sort = F)
  } else {
    DT.design[, Condition := NA_character_]
  }

  setcolorder(DT.design, c("Run", "Channel", "Assay"))

  if (!data.is.data.table) setDF(data)
  setDF(DT.design)
  return(DT.design[])
}
