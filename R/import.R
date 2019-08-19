#' Imported dataset run information
#'
#' Get or set run information from a \link{data.frame} returned by \link{import_ProteinPilot} or \link{import_ProteomeDiscoverer}.
#' Used to manually add run information to imported datasets
#'
#' @param data \link{data.frame} returned by \link{import_ProteinPilot} or \link{import_ProteomeDiscoverer}.
#' @import data.table
#' @export
runs <- function(data) {
  data <- setDT(data)[, .(Run = first(Run)), keyby = Injection]
  data[, Run := as.character(Run)]
  data[, Injection := as.character(Injection)]
  setDF(data)
  return(data)
}


#' @rdname runs
#' @param value Edited \link{data.frame} created by \code{runs}.
#' @import data.table
#' @export
`runs<-` <- function(data, value) {
  DT <- setDT(data)
  DT[, Run := NULL]
  DT.runs <- as.data.table(value)
  DT.runs <- DT.runs[complete.cases(DT.runs)]
  if (!is.factor(DT.runs$Injection)) DT.runs[, Injection := factor(Injection)]
  if (!is.factor(DT.runs$Run)) DT.runs[, Run := factor(Run)]
  DT <- merge(DT, DT.runs, by = "Injection")
  data <- setDF(DT)
  return(data)
}


#' Specify study design for imported dataset
#'
#' Returns a skeleton \link{data.frame} for customising the study design of a dataset imported with  \link{import_ProteinPilot} or
#' \link{import_ProteomeDiscoverer}.
#'
#' @param data \link{data.frame} returned by \link{import_ProteinPilot} or \link{import_ProteomeDiscoverer}.
#' @return \link{data.frame} that can be edited and passed as parameter \code{data.design} of \link{bayesprot}.
#' @import data.table
#' @export
new_design <- function(data) {
  data.design <- setDT(data)[, .(Assay = paste(Run, Channel, sep = ",")), keyby = .(Run, Channel)]
  data.design[, Assay := sub("^,", "", Assay)]
  data.design[, Assay := sub(",$", "", Assay)]
  data.design[, Sample := Assay]
  data.design[, Condition := NA]
  data.design[, ref := T]
  setDF(data.design)
  setDF(data)
  return(data.design)
}


#' Import SCIEX ProteinPilot data
#'
#' Reads in a SCIEX ProteinPilot \code{PeptideSummary.txt} file for processing with \link{bayesprot}.
#'
#' @param file Location of the \code{PeptideSummary.txt} file.
#' @param shared Include shared peptides?
#' @param min.conf Features with peptide ID confidence less than \code{min.conf} (between 0 - 100) are filtered out. The default
#'   \code{"auto"} uses the ProteinPilot default threshold.
#' @param filter Other filters to use, which can include \code{"discordant peptide type"}, \code{"no iTRAQ" and
#'   \code{"weak signal"}}
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export

import_ProteinPilot <- function(
  file = NULL,
  shared = F,
  min.conf = "auto",
  filter = c("discordant peptide type", "no iTRAQ", "weak signal"),
  data = NULL
) {
  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file = file, showProgress = T)
  }

  # filtering
  if (min.conf == "auto") min.conf <- DT.raw[Annotation == "auto", min(Conf)]
  DT.raw <- DT.raw[Conf >= min.conf]
  if ("discordant peptide type" %in% filter) DT.raw <- DT.raw[Annotation != "auto - discordant peptide type"]
  if ("no iTRAQ" %in% filter) DT.raw <- DT.raw[Annotation != "auto - no iTRAQ"]
  if ("weak signal" %in% filter) DT.raw <- DT.raw[Annotation != "no quant - weak signal"]
  DT.raw <- DT.raw[!grepl("^RRRRR.*", DT.raw$Accessions),] # decoys
  if(!("ProteinModifications" %in% colnames(DT.raw))) DT.raw[, ProteinModifications := ""]

  # deal with shared peptides
  if (!shared) DT.raw <- DT.raw[Annotation != "auto - shared MS/MS"]

  # create wide data table
  DT <- DT.raw[, .(
    ProteinInfo = paste0("[", N, "] ", Names),
    Protein = gsub(";", "", Accessions),
    Peptide = gsub(" ", "", paste0(Sequence, ",", Modifications, ",", ProteinModifications, ",", Cleavages), fixed = T),
    Feature = Spectrum,
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

  # group ambiguous PSMs so BayesProt treats them as a single peptide per protein
  DT[, Peptide := paste(sort(as.character(Peptide)), collapse = " "), by = .(Protein, Feature)]
  DT <- unique(DT)

  # melt
  DT[, ProteinInfo := factor(ProteinInfo)]
  DT[, Protein := factor(Protein)]
  DT[, Peptide := factor(Peptide)]
  DT[, Feature := factor(Feature)]
  DT[, Injection := factor(Injection)]
  DT[, Run := factor("")]
  DT <- melt(DT, variable.name = "Channel", value.name = "Count", measure.vars = colnames(DT)[grep("^Channel\\.", colnames(DT))])
  levels(DT$Channel) <- sub("^Channel\\.", "", levels(DT$Channel))

  setDF(DT)
  return(DT)
}


#' Import Thermo ProteomeDiscoverer data
#'
#' Reads in a Thermo ProteomeDiscoverer \code{PSMs.txt} file for processing with \link{bayesprot}.
#'
#' @param file Location of the \code{PSMs.txt} file.
#' @param shared Include shared peptides?
#' @param used Include only features marked as used by ProteomeDiscoverer?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export
import_ProteomeDiscoverer <- function(
  file = NULL,
  shared = F,
  used = T,
  data = NULL
) {
  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file = file, showProgress = T)
  }

  # only use rows that ProteomeDiscoverer uses for quant
  if (!shared) DT.raw <- DT.raw[`Quan Info` == "Unique",]
  if (!used) DT.raw <- DT.raw[`Peptide Quan Usage` == "Use",]

  # merge fraction info
  #setnames(DT, "Spectrum File", "File")
  #if (is.null(data.runs)) {
  #  warning("No 'data.runs' parameter supplied, assuming no fractionation!")
  #  DT[, Run := File]
  #  DT[, Fraction := ""]
  #} else {
  #  DT <- merge(DT, setDT(data.runs), by = "File")
  #}

  # only retain the most confidenct and most intense spectrum for each feature (TODO: make option)
  #DT[, Feature := factor(paste0(DT$Sequence, " : ", DT$Modifications, " : ", DT$Charge, "+ : ", DT$File))]
  #setorder(DT, Feature, Confidence, -Intensity)
  #DT <- unique(DT, by = "Feature")

  # create wide data table
  DT <- DT.raw[ , list(
    ProteinInfo = `Protein Descriptions`,
    Protein = `Master Protein Accessions`,
    Peptide = gsub(" ", "", paste0(Sequence, ",", Modifications)),
    Feature = paste0(`Spectrum File`, ",", `First Scan`),
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
  if ("131N" %in% colnames(DT.raw)) DT$Channel.131N <- DT.raw$`131N`
  if ("131C" %in% colnames(DT.raw)) DT$Channel.131C <- DT.raw$`131C`
  if ("131" %in% colnames(DT.raw)) DT$Channel.131 <- DT.raw$`131`

  # group ambiguous PSMs so BayesProt treats them as a single peptide per protein
  DT[, Peptide := paste(sort(as.character(Peptide)), collapse = " "), by = .(Protein, Feature)]
  DT <- unique(DT)

  # melt label counts
  DT[, ProteinInfo := factor(ProteinInfo)]
  DT[, Protein := factor(Protein)]
  DT[, Peptide := factor(Peptide)]
  DT[, Feature := factor(Feature)]
  DT[, Injection := factor(Injection)]
  DT[, Run := factor("")]
  DT <- melt(DT, variable.name = "Channel", value.name = "Count", measure.vars = colnames(DT)[grep("^Channel\\.", colnames(DT))])
  levels(DT$Channel) <- sub("^Channel\\.", "", levels(DT$Channel))

  setDF(DT)
  return(DT)
}


#' Import Waters Progenesis data
#'
#' Reads in a Waters Progenesis \code{pep_ion_measurements.csv} file for processing with \link{bayesprot}.
#'
#' @param file Location of the \code{pep_ion_measurements.csv} file.
#' @param used Include only features marked as used by Progenesis?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export

import_Progenesis <- function(
  file = NULL,
  used = T,
  data = NULL
) {
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
    colname.protein <- "Accession"
    colname.proteininfo <- "Description"
    colname.sequence <- "Sequence"
    colname.modifications <- "Modifications"

    # only use rows that Progenesis uses for quant
    if (used) DT.raw <- DT.raw[`Use in quantitation` == "True",]
  } else {
    colname.protein <- colnames(DT.raw)[grep("Best peptide match Protein$", colnames(DT.raw))]
    colname.proteininfo <- colnames(DT.raw)[grep("Best peptide match Description$", colnames(DT.raw))]
    colname.sequence <- colnames(DT.raw)[grep("Best peptide match Sequence$", colnames(DT.raw))]
    colname.modifications <- colnames(DT.raw)[grep("Best peptide match Variable modifications \\(\\[position\\] description\\)$", colnames(DT.raw))]
  }

  # remove decoys and strange missing Accession
  DT.raw <- DT.raw[get(colname.protein) != "",]
  DT.raw <- DT.raw[!grepl("^#DECOY#", get(colname.protein)),]

  # create wide data table
  DT <- cbind(DT.raw[, .(
    ProteinInfo = get(colname.proteininfo),
    Protein = get(colname.protein),
    Peptide = gsub(" ", "", paste0(get(colname.sequence), ",", get(colname.modifications))),
    Feature = `#`
  )], DT.raw[, .SD, .SDcols = names(DT.raw) %like% "^Raw abundance "])

  # group ambiguous PSMs so BayesProt treats them as a single peptide per protein
  DT[, Peptide := paste(sort(as.character(Peptide)), collapse = " "), by = .(Protein, Feature)]
  DT <- unique(DT)

  # melt
  DT[, ProteinInfo := factor(ProteinInfo)]
  DT[, Protein := factor(Protein)]
  DT[, Peptide := factor(Peptide)]
  DT[, Feature := factor(Feature)]
  DT <- melt(DT, variable.name = "Run", value.name = "Count", measure.vars = colnames(DT)[grep("^Raw abundance ", colnames(DT))])
  levels(DT$Run) <- sub("^Raw abundance ", "", levels(DT$Run))
  DT[, Count := as.numeric(Count)]

  setDF(DT)
  return(DT)
}


#' Import OpenSWATH - PyProphet data
#'
#' Reads in a set of \code{_with_dscore} datasets processed by OpenSWATH and PyProphet for processing with \link{bayesprot}.
#'
#' @param files One of more \code{_with_dscore} files to import and merge.
#' @param m_score.cutoff Include only features with PyProphet m_score >= than this?
#' @param data Advanced: Rather than specifying \code{files}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export
import_OpenSwath_PyProphet <- function(
  files = NULL,
  shared = F,
  m_score.cutoff = 0.01,
  data = NULL
) {
  if (is.null(file) && is.null(data)) stop("One of 'data' or 'files' needs to be specified.")
  if (!is.null(data)) files <- data

  DT <- rbindlist(lapply(files, function(file) {
    if (is.data.frame(file)) {
      DT.raw <- as.data.file(data)
    } else {
      DT.raw <- fread(file = file, showProgress = T)
    }

    # remove decoys and > m_score.cutoff
    DT.raw <- DT.raw[decoy == 0,]
    DT.raw <- DT.raw[m_score <= m_score.cutoff,]

    # create long data table
    DT <- DT.raw[, .(
      Protein = ProteinName,
      Peptide = FullPeptideName,
      Feature = gsub(";", ";bayesprot;", aggr_Fragment_Annotation),
      Run = filename,
      Count = gsub(";", ";bayesprot;", aggr_Peak_Area)
    )]
    DT <- DT[, lapply(.SD, function(x) unlist(tstrsplit(x, ";bayesprot;", fixed = T)))]
    DT[, Count := as.numeric(Count)]
    DT
  }))

  # remove features that have more than one identification in any assay
  DT[, N := .N, by = .(Feature, Run)]
  DT <- DT[N == 1]
  DT[, N := NULL]
  assays <- unique(DT$Run)

  # create wide data table
  DT <- dcast(DT, Protein + Peptide + Feature ~ Run, value.var = "Count")

  # remove shared
  if (!shared) DT <- DT[grepl("^1/", DT$Protein)]

  # group ambiguous PSMs so BayesProt treats them as a single peptide per protein
  DT[, Peptide := paste(sort(as.character(Peptide)), collapse = " "), by = .(Protein, Feature)]
  DT <- unique(DT)

  # melt
  DT[, ProteinInfo := factor(Protein)]
  DT[, Protein := factor(sub("^1/", "", Protein))]
  DT[, Peptide := factor(Peptide)]
  DT[, Feature := factor(Feature)]
  DT <- melt(DT, variable.name = "Run", value.name = "Count", measure.vars = assays)
  DT[, Channel := factor("")]
  setcolorder(DT, c("Protein", "ProteinInfo", "Peptide", "Feature", "Run", "Channel"))

  setDF(DT)
  return(DT)
}


#' Import data outputed by an MSstats import routine
#'
#' Reads in a set of \code{_with_dscore} datasets processed by OpenSWATH and PyProphet for processing with \link{bayesprot}.
#'
#' @param files One of more \code{_with_dscore} files to import and merge.
#' @param m_score.cutoff Include only features with PyProphet m_score >= than this?
#' @param data Advanced: Rather than specifying \code{files}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export
import_MSstats <- function(data) {
  DT <- data.table(
    Protein = factor(data$ProteinName),
    ProteinInfo = "",
    Peptide = factor(data$PeptideSequence),
    Feature = factor(data$FragmentIon),
    Run = factor(data$Run, levels = unique(data$Run)),
    Channel = factor(data$IsotopeLabelType, levels = unique(data$IsotopeLabelType)),
    Count = as.numeric(data$Intensity)
  )

  setDF(DT)
  return(DT)
}
