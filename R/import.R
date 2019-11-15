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

  if (data.is.data.table) setDF(data)
  return(setDF(DT))
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

  if (data.is.data.table) setDF(data)
  return(setDF(DT))
}


#' Specify study design for imported dataset
#'
#' Returns a skeleton \link{data.frame} for customising the study design of a dataset imported with  \link{import_ProteinPilot} or
#' \link{import_ProteomeDiscoverer}.
#'
#' @param data \link{data.frame} returned by \link{import_ProteinPilot} or \link{import_ProteomeDiscoverer}.
#' @return \link{data.frame} that can be edited and passed as parameter \code{data.design} of \link{seamassdelta}.
#' @import data.table
#' @export
new_design <- function(data) {
  data.is.data.table <- is.data.table(data)
  DT <- setDT(data)

  DT.design <- DT[, .(Assay = paste(Run, Channel, sep = ",")), keyby = .(Run, Channel)]
  DT.design[, Assay := sub("^,", "", Assay)]
  DT.design[, Assay := sub(",$", "", Assay)]

  # autodetect blocks
  block <- merge(DT.design[, .(Run, Channel, Assay)], DT[, .(Run, Channel, Measurement)], by = c("Run", "Channel"))
  block[, Run := NULL]
  block[, Channel := NULL]
  block[, N := 1]
  block <- dcast(block, Measurement ~ Assay, sum, value.var = "N")
  block[, Measurement := NULL]
  block <- as.matrix(block)
  # matrix multiplication distributes assay relationships
  block <- t(block) %*% block
  # Block is recoded first non-zero occurence for each assay
  DT.design[, Block := as.integer(factor(colnames(block)[apply(block != 0, 2, which.max)]))]
  # default is all assays are reference assays for each block
  DT.design[, BlockRef := T]

  # DEA
  DT.design[, Sample := Assay]
  DT.design[, Condition := NA]

  if (data.is.data.table) setDF(data)
  return(setDF(DT.design))
}


#' Import SCIEX ProteinPilot data
#'
#' Reads in a SCIEX ProteinPilot \code{ComponentSummary.txt} file for processing with \link{seamassdelta}.
#'
#' @param file Location of the \code{ComponentSummary.txt} file.
#' @param shared Include shared components?
#' @param min.conf Measurements with component ID confidence less than \code{min.conf} (between 0 - 100) are filtered out. The default
#'   \code{"auto"} uses the ProteinPilot default threshold.
#' @param filter Other filters to use, which can include \code{"discordant component type"}, \code{"no iTRAQ" and
#'   \code{"weak signal"}}
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seamassdelta}.
#' @import data.table
#' @export

import_ProteinPilot <- function(
  file = NULL,
  shared = F,
  min.conf = "auto",
  filter = c("discordant component type", "no iTRAQ", "weak signal"),
  data = NULL
) {
  suppressWarnings(suppressMessages(library(R.oo)))

  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file = file, showProgress = T)
  }

  # filtering
  if (min.conf == "auto") min.conf <- DT.raw[Annotation == "auto", min(Conf)]
  DT.raw <- DT.raw[Conf >= min.conf]
  if ("discordant component type" %in% filter) DT.raw <- DT.raw[Annotation != "auto - discordant component type"]
  if ("no iTRAQ" %in% filter) DT.raw <- DT.raw[Annotation != "auto - no iTRAQ"]
  if ("weak signal" %in% filter) DT.raw <- DT.raw[Annotation != "no quant - weak signal"]
  DT.raw <- DT.raw[!grepl("^RRRRR.*", DT.raw$Accessions),] # decoys
  if(!("GroupModifications" %in% colnames(DT.raw))) DT.raw[, GroupModifications := ""]

  # deal with shared components
  if (!shared) DT.raw <- DT.raw[Annotation != "auto - shared MS/MS"]

  # create wide data table
  DT <- DT.raw[, .(
    GroupInfo = paste0("[", N, "] ", Names),
    Group = gsub(";", "", Accessions),
    Component = gsub(" ", "", paste0(Sequence, ",", Modifications, ",", ProteinModifications, ",", Cleavages), fixed = T),
    Measurement = Spectrum,
    Injection = as.integer(matrix(unlist(strsplit(as.character(DT.raw$Spectrum), ".", fixed = T)), ncol = 5, byrow = T)[, 1])
  )]
  if("Area 113" %in% colnames(DT.raw)) DT$Channel.113 <- DT.raw$`Area 113`
  if("Area 114" %in% colnames(DT.raw)) DT$Channel.114 <- DT.raw$`Area 114`
  if("Area 115" %in% colnames(DT.raw)) DT$Channel.115 <- DT.raw$`Area 115`
  if("Area 116" %in% colnames(DT.raw)) DT$Channel.116 <- DT.raw$`Area 116`
  if("Area 117" %in% colnames(DT.raw)) DT$Channel.117 <- DT.raw$`Area 117`
  if("Area 118" %in% colnames(DT.raw)) DT$Channel.118 <- DT.raw$`Area 118`
  if("Area 119" %in% colnames(DT.raw)) DT$Channel.119 <- DT.raw$`Area 119`
  if("Area 121" %in% colnames(DT.raw)) DT$Channel.121 <- DT.raw$`Area 121`

  # group ambiguous PSMs so seaMass-Delta treats them as a single component per group
  DT[, Component := paste(sort(as.character(Component)), collapse = " "), by = .(Group, Measurement)]
  DT <- unique(DT)

  # melt
  DT[, GroupInfo := factor(GroupInfo)]
  DT[, Group := factor(Group)]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT[, Run := factor("")]
  DT[, Injection := factor(Injection)]
  DT <- melt(DT, variable.name = "Channel", value.name = "Count", measure.vars = colnames(DT)[grep("^Channel\\.", colnames(DT))])
  levels(DT$Channel) <- sub("^Channel\\.", "", levels(DT$Channel))

  setDF(DT)
  return(DT)
}


#' Import Thermo ProteomeDiscoverer data
#'
#' Reads in a Thermo ProteomeDiscoverer \code{PSMs.txt} file for processing with \link{seamassdelta}.
#'
#' @param file Location of the \code{PSMs.txt} file.
#' @param shared Include shared components?
#' @param used Include only measurements marked as used by ProteomeDiscoverer?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seamassdelta}.
#' @import data.table
#' @export
import_ProteomeDiscoverer <- function(
  file = NULL,
  shared = F,
  used = T,
  data = NULL
) {
  suppressWarnings(suppressMessages(library(R.oo)))

  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file = file, showProgress = T)
  }

  # only use rows that ProteomeDiscoverer uses for quant
  if (!shared) DT.raw <- DT.raw[`Quan Info` == "Unique",]
  if (!used) DT.raw <- DT.raw[`Component Quan Usage` == "Use",]

  # merge fraction info
  #setnames(DT, "Spectrum File", "File")
  #if (is.null(data.runs)) {
  #  warning("No 'data.runs' parameter supplied, assuming no fractionation!")
  #  DT[, Run := File]
  #  DT[, Fraction := ""]
  #} else {
  #  DT <- merge(DT, setDT(data.runs), by = "File")
  #}

  # only retain the most confidenct and most intense spectrum for each measurement (TODO: make option)
  #DT[, Measurement := factor(paste0(DT$Sequence, " : ", DT$Modifications, " : ", DT$Charge, "+ : ", DT$File))]
  #setorder(DT, Measurement, Confidence, -Intensity)
  #DT <- unique(DT, by = "Measurement")

  # create wide data table
  DT <- DT.raw[ , list(
    GroupInfo = `Protein Descriptions`,
    Group = `Master Protein Accessions`,
    Component = gsub(" ", "", paste0(Sequence, ",", Modifications)),
    Measurement = paste0(`Spectrum File`, ",", `First Scan`),
    Injection = `Spectrum File`
  )]
  if ("Light" %in% colnames(DT.raw)) DT$Channel.Light <- DT.raw$Light
  if ("Medium" %in% colnames(DT.raw)) DT$Channel.Medium <- DT.raw$Medium
  if ("Heavy" %in% colnames(DT.raw)) DT$Channel.Heavy <- DT.raw$Heavy
  if ("126N" %in% colnames(DT.raw)) DT$Channel.126N <- DT.raw$`126N`
  if ("126C" %in% colnames(DT.raw)) DT$Channel.126C <- DT.raw$`126C`
  if ("126" %in% colnames(DT.raw)) DT$Channel.126 <- DT.raw$`126`
  if ("127N" %in% colnames(DT.raw)) DT$Channel.127N <- DT.raw$`127N`
  if ("127C" %in% colnames(DT.raw)) DT$Channel.127C <- DT.raw$`127C`
  if ("127" %in% colnames(DT.raw)) DT$Channel.127 <- DT.raw$`127`
  if ("128N" %in% colnames(DT.raw)) DT$Channel.128N <- DT.raw$`128N`
  if ("128C" %in% colnames(DT.raw)) DT$Channel.128C <- DT.raw$`128C`
  if ("128" %in% colnames(DT.raw)) DT$Channel.128 <- DT.raw$`128`
  if ("129N" %in% colnames(DT.raw)) DT$Channel.129N <- DT.raw$`129N`
  if ("129C" %in% colnames(DT.raw)) DT$Channel.129C <- DT.raw$`129C`
  if ("129" %in% colnames(DT.raw)) DT$Channel.129 <- DT.raw$`129`
  if ("130N" %in% colnames(DT.raw)) DT$Channel.130N <- DT.raw$`130N`
  if ("130C" %in% colnames(DT.raw)) DT$Channel.130C <- DT.raw$`130C`
  if ("130" %in% colnames(DT.raw)) DT$Channel.130 <- DT.raw$`130`
  if ("131N" %in% colnames(DT.raw)) DT$Channel.131N <- DT.raw$`131N`
  if ("131C" %in% colnames(DT.raw)) DT$Channel.131C <- DT.raw$`131C`
  if ("131" %in% colnames(DT.raw)) DT$Channel.131 <- DT.raw$`131`

  # group ambiguous PSMs so seaMass-Delta treats them as a single component per group
  DT[, Component := paste(sort(as.character(Component)), collapse = " "), by = .(Group, Measurement)]
  DT <- unique(DT)

  # melt label counts
  DT[, GroupInfo := factor(GroupInfo)]
  DT[, Group := factor(Group)]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT[, Injection := factor(Injection)]
  DT[, Run := factor("")]
  DT <- melt(DT, variable.name = "Channel", value.name = "Count", measure.vars = colnames(DT)[grep("^Channel\\.", colnames(DT))])
  levels(DT$Channel) <- sub("^Channel\\.", "", levels(DT$Channel))

  setDF(DT)
  return(DT)
}


#' Import Waters Progenesis data
#'
#' Reads in a Waters Progenesis \code{pep_ion_measurements.csv} file for processing with \link{seamassdelta}.
#'
#' @param file Location of the \code{pep_ion_measurements.csv} file.
#' @param used Include only measurements marked as used by Progenesis?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seamassdelta}.
#' @import data.table
#' @export

import_Progenesis <- function(
  file = NULL,
  used = T,
  data = NULL
) {
  suppressWarnings(suppressMessages(library(R.oo)))

  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file = file, showProgress = T)
  }

  # sort out column names
  for (j in 2:ncol(DT.raw)) {
    if (DT.raw[[j]][1] == "") DT.raw[[j]][1] <- DT.raw[[j-1]][1]
    if (DT.raw[[j]][2] == "") DT.raw[[j]][2] <- DT.raw[[j-1]][2]
  }
  colnames(DT.raw) <- apply(DT.raw[1:3,], 2, function(x) trimws(paste(x[1], x[2], x[3])))
  DT.raw <- DT.raw[4:nrow(DT.raw),]

  # column names for different file types
  if ("Accession" %in% colnames(DT.raw)) {
    colname.group <- "Accession"
    colname.groupinfo <- "Description"
    colname.sequence <- "Sequence"
    colname.modifications <- "Modifications"

    # only use rows that Progenesis uses for quant
    if (used) DT.raw <- DT.raw[`Use in quantitation` == "True",]
  } else {
    colname.group <- colnames(DT.raw)[grep("Best peptide match Group$", colnames(DT.raw))]
    colname.groupinfo <- colnames(DT.raw)[grep("Best peptide match Description$", colnames(DT.raw))]
    colname.sequence <- colnames(DT.raw)[grep("Best peptide match Sequence$", colnames(DT.raw))]
    colname.modifications <- colnames(DT.raw)[grep("Best peptide match Variable modifications \\(\\[position\\] description\\)$", colnames(DT.raw))]
  }

  # remove decoys and strange missing Accession
  DT.raw <- DT.raw[get(colname.group) != "",]
  DT.raw <- DT.raw[!grepl("^#DECOY#", get(colname.group)),]

  # create wide data table
  DT <- cbind(DT.raw[, .(
    GroupInfo = get(colname.groupinfo),
    Group = get(colname.group),
    Component = gsub(" ", "", paste0(get(colname.sequence), ",", get(colname.modifications))),
    Measurement = `#`
  )], DT.raw[, .SD, .SDcols = names(DT.raw) %like% "^Raw abundance "])

  # group ambiguous PSMs so seaMass-Delta treats them as a single component per group
  DT[, Component := paste(sort(as.character(Component)), collapse = " "), by = .(Group, Measurement)]
  DT <- unique(DT)

  # melt
  DT[, GroupInfo := factor(GroupInfo)]
  DT[, Group := factor(Group)]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT <- melt(DT, variable.name = "Run", value.name = "Count", measure.vars = colnames(DT)[grep("^Raw abundance ", colnames(DT))])
  levels(DT$Run) <- sub("^Raw abundance ", "", levels(DT$Run))
  DT[, Injection := Run]
  DT[, Count := as.numeric(Count)]

  setDF(DT)
  return(DT)
}


#' Import OpenSWATH data
#'
#' Reads in the output of an OpenSWATH -> PyProphet -> TRIC pipeline for processing with \link{seamassdelta}.
#'
#' @param files A \code{csv} file to import.
#' @param m_score.cutoff Include only measurements with PyProphet m_score >= than this?
#' @param data Advanced: Rather than specifying a \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seamassdelta}.
#' @import data.table
#' @export
import_OpenSwath <- function(
  file = NULL,
  shared = FALSE,
  m_score.cutoff = 0.05,
  #missing.cutoff = 0.5,
  data = NULL
) {
  suppressWarnings(suppressMessages(library(R.oo)))

  if (is.null(file) && is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
  if (!is.null(data)) file <- data

  DT <- rbindlist(lapply(file, function(f) {
    if (is.data.frame(f)) {
      DT.raw <- as.data.file(data)
    } else {
      DT.raw <- fread(file = f, showProgress = T)
    }

    # remove decoys and > m_score.cutoff
    DT.raw <- DT.raw[decoy == 0,]
    DT.raw <- DT.raw[m_score <= m_score.cutoff,]

    # create long data table
    DT <- DT.raw[, .(
      Group = ProteinName,
      Component = FullPeptideName,
      Measurement = gsub(";", ";bayesprot;", aggr_Fragment_Annotation),
      Run = paste(run_id, tools::file_path_sans_ext(basename(filename)), sep = ";"),
      Count = gsub(";", ";bayesprot;", aggr_Peak_Area)
    )]
    DT <- DT[, lapply(.SD, function(x) unlist(tstrsplit(x, ";bayesprot;", fixed = T)))]
    DT[, Count := as.numeric(Count)]
    DT
  }))

  # remove measurements that have more than one identification in any assay
  DT[, N := .N, by = .(Measurement, Run)]
  DT <- DT[N == 1]
  DT[, N := NULL]
  assays <- unique(DT$Run)

  # create wide data table
  DT <- dcast(DT, Group + Component + Measurement ~ Run, value.var = "Count")

  # remove transitions that are present is less than missing.cutoff assays
  #DT <- DT[rowSums(is.na(DT)) < missing.cutoff * length(assays)]

  # remove shared
  if (!shared) DT <- DT[grepl("^1/", DT$Group)]

  # group ambiguous transitions so seaMass-Delta treats them as a single component per group ## UNNECCESARY?
  DT[, Component := paste(sort(as.character(Component)), collapse = " "), by = .(Group, Measurement)]
  DT <- unique(DT)

  # melt
  DT[, GroupInfo := factor(Group)]
  DT[, Group := factor(sub("^1/", "", Group))]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT <- melt(DT, variable.name = "Run", value.name = "Count", measure.vars = assays)
  DT[, Injection := Run]
  DT[, Channel := factor("")]
  setcolorder(DT, c("Group", "GroupInfo", "Component", "Measurement", "Run", "Injection", "Channel"))

  setDF(DT)
  return(DT)
}


#' Import data outputed by an MSstats import routine
#'
#' Reads in a set of \code{_with_dscore} datasets processed by OpenSWATH and PyProphet for processing with \link{seamassdelta}.
#'
#' @param data MSstats output
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seamassdelta}.
#' @import data.table
#' @export
import_MSstats <- function(data) {
  DT <- data.table(
    Group = factor(data$ProteinName),
    GroupInfo = "",
    Component = factor(data$PeptideSequence),
    Measurement = factor(data$FragmentIon),
    Run = factor(data$Run, levels = unique(data$Run)),
    Injection = factor(data$Run, levels = unique(data$Run)),
    Channel = factor(data$IsotopeLabelType, levels = unique(data$IsotopeLabelType)),
    Count = as.numeric(data$Intensity)
  )

  setDF(DT)
  return(DT)
}
