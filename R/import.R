check_dups <- function(DT) {
  # if duplicated rows, warn and use only first
  DT[, N := .N, by = .(Measurement, Run)]
  DT.dup <- DT[N != 1]
  if (nrow(DT.dup) > 0) {
    warning(paste0("Some measurements are duplicated, only the first are used:\n", paste(capture.output(print(DT.dup)), collapse = "\n")))
    DT[N != 1, Use := F]
  }
  DT[, N := NULL]
  return(DT)
}


#' Import OpenSWATH data
#'
#' Reads in the output of an OpenSWATH -> PyProphet -> TRIC pipeline.
#'
#' @param files A \code{csv} file to import.
#' @param max.m_score Include only measurements with PyProphet m_score >= than this?
#' @param data Advanced: Rather than specifying a \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seaMass_sigma}.
#' @import data.table
#' @export
import_OpenSWATH <- function(
  file = NULL,
  max.m_score = 0.05,
  use.decoys = FALSE,
  protein.decoy.prefix = "DECOY",
  data = NULL
) {
  if (is.null(file) && is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
  if (!is.null(data)) file <- data

  DT <- rbindlist(lapply(file, function(f) {
    if (is.data.frame(f)) {
      DT.wide <- as.data.table(data)
    } else {
      DT.wide <- fread(file = f, showProgress = T)
    }

    # consolidate OpenSWATH's 'Protein Groups', removing decoys
    DT.wide[, row := seq_len(nrow(DT.wide))]
    DT.groups <- DT.wide[, .(Group = unlist(strsplit(sub("^[^/]+/", "", ProteinName), "/"))), by = row]
    DT.groups[, Group := trimws(Group)]
    if (!is.null(protein.decoy.prefix) && protein.decoy.prefix != "") DT.groups <- DT.groups[!grep(paste0("^", protein.decoy.prefix, "_"), Group)]
    setorder(DT.groups, row, Group)
    # reassemble to our Protein Groups
    DT.groups <- DT.groups[, .(N = .N, Group = paste(Group, collapse = "/")), by = row]
    setorder(DT.groups, Group)
    DT.groups[N > 1, Group := paste(N, Group, sep = "/")]
    DT.groups[, N := NULL]
    DT.groups[, Group := factor(Group, levels = unique(Group))]
    DT.wide <- merge(DT.wide, DT.groups, by = "row", sort = F, all.x = T)

    # filters
    DT.wide[, Use := T]
    DT.wide[m_score > max.m_score, Use := F]
    if (!use.decoys) DT.wide[decoy != 0, Use := F]

    # create long data table
    DT.wide <- DT.wide[, .(
      Group,
      Component = FullPeptideName,
      Measurement = gsub(";", ";seaMass;", aggr_Fragment_Annotation),
      Run = paste(run_id, tools::file_path_sans_ext(basename(filename)), sep = ";"),
      Count = gsub(";", ";seaMass;", aggr_Peak_Area),
      Use
    )]
    DT.long <- DT.wide[, lapply(.SD, function(x) unlist(tstrsplit(x, ";seaMass;", fixed = T)))]
    DT.long[, Group := factor(Group, levels = levels(DT.wide$Group))]
    DT.long[, Component := factor(Component, levels = unique(Component))]
    DT.long[, Measurement := factor(Measurement, levels = unique(Measurement))]
    DT.long[, Run := factor(Run, levels = unique(Run))]
    DT.long[, Count := as.numeric(Count)]
    DT.long[, Use := as.logical(Use)]

    return(DT.long)
  }))

  DT[, GroupInfo := ""]
  DT[, Channel := factor("1")]
  setcolorder(DT, c("Group", "GroupInfo", "Component", "Measurement", "Run", "Channel"))
  setorder(DT, Group, Component, Measurement, Run, Channel)

  setattr(DT, "group", c("Protein", "Proteins"))
  setattr(DT, "component", c("Peptidoform", "Peptidoforms"))
  setattr(DT, "measurement", c("Transition", "Transitions"))

  DT <- check_dups(DT)
  setDF(DT)
  return(DT[])
}


#### BELOW TO CHECK


#' Import Thermo ProteomeDiscoverer data
#'
#' Reads in a Thermo ProteomeDiscoverer \code{PSMs.txt} file.
#'
#' @param file Location of the \code{PSMs.txt} file.
#' @param shared Include shared peptides?
#' @param used Include only measurements marked as used by ProteomeDiscoverer?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seaMass_sigma}.
#' @import data.table
#' @export
import_ProteomeDiscoverer <- function(
  file = NULL,
  use.shared.peptides = FALSE,
  use.not.unique = TRUE,
  use.excluded.by.method = TRUE,
  use.redundant = FALSE,
  data = NULL
) {
  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file = file, showProgress = T)
  }

  # filtering
  DT.raw[, Use := T]
  if (!use.shared.peptides) DT.raw[`Number of Protein Groups` > 1, Use := F]
  if (!use.not.unique) DT.raw[`Quan Info` == "NotUnique", Use := F]
  if (!use.excluded.by.method) DT.raw[`Quan Info` == "ExcludedByMethod", Use := F]
  if (!use.redundant) DT.raw[`Quan Info` == "Redundant", Use := F]

  # create wide data table
  DT <- DT.raw[ , list(
    Groups = `Master Protein Accessions`,
    GroupInfo = factor(""),
    Component = gsub(" ", "", paste0(Sequence, ",", Modifications)),
    Measurement = paste0(`Spectrum File`, ",", `First Scan`),
    Injection = `Spectrum File`,
    Use
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
  rm(DT.raw)

  # fold out shared peptides
  DT[, row := seq_len(nrow(DT))]
  DT <- merge(DT, DT[, list(Group = unlist(strsplit(Groups, ";"))), by = row], by = "row", sort = F)
  DT[, Group := trimws(Group)]
  DT[, Groups := NULL]
  DT[, row := NULL]

  # we can get a spectrum assigned to multiple peptides in a protein - if this occurs, assign the spectrum as a 'new' ambiguous peptide
  # hack: data.table bug crashes R?
  setorder(DT, Component)
  DT[, Component2 := toString(Component), by = .(Group, Measurement, Use)]
  DT[, Component := NULL]
  setnames(DT, "Component2", "Component")
  DT <- unique(DT)

  # melt label counts
  DT[, GroupInfo := factor(GroupInfo)]
  DT[, Group := factor(Group)]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT[, Run := factor("1")]
  DT[, Injection := factor(Injection)]
  DT <- melt(DT, variable.name = "Channel", value.name = "Count", measure.vars = colnames(DT)[grep("^Channel\\.", colnames(DT))])
  levels(DT$Channel) <- sub("^Channel\\.", "", levels(DT$Channel))

  setcolorder(DT, c("Group", "GroupInfo", "Component", "Measurement", "Run", "Channel", "Count"))
  setorder(DT, Group, Component, Measurement, Run, Channel)

  setattr(DT, "group", c("ProteinGroup", "ProteinGroups"))
  setattr(DT, "component", c("Peptidoform", "Peptidoforms"))
  setattr(DT, "measurement", c("Spectrum", "Spectra"))

  setDF(DT)
  return(DT[])
}


#' Import SCIEX ProteinPilot data
#'
#' Reads in a SCIEX ProteinPilot \code{PeptideSummary.txt} file for processing with \link{seaMass_sigma}.
#'
#' @param file Location of the \code{PeptideSummary.txt} file.
#' @param shared Include shared peptides?
#' @param min.conf Measurements with peptide ID confidence less than \code{min.conf} (between 0 - 100) are filtered out. The default
#'   \code{"auto"} uses the ProteinPilot default threshold.
#' @param filter Other filters to use, which can include \code{"discordant peptide type"}, \code{"no iTRAQ" and
#'   \code{"weak signal"}}
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seaMass_sigma}.
#' @import data.table
#' @export
import_ProteinPilot <- function(
  file = NULL,
  min.conf = "auto",
  use.decoys = FALSE,
  use.shared.peptides = FALSE,
  use.discordant.peptide.type = FALSE,
  use.no.itraq = FALSE,
  use.weak.signal = FALSE,
  data = NULL
) {
  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT.raw(data)
  } else {
    DT.raw <- fread(file = file, showProgress = T)
  }

  # filtering
  DT.raw[, Use := T]
  if (min.conf == "auto") min.conf <- DT.raw[Annotation == "auto", min(Conf)]
  DT.raw[Conf < min.conf, Use := F]
  if (!use.decoys) DT.raw[grepl("^RRRRR.*", DT.raw$Accessions), Use := F]
  if (!use.shared.peptides) DT.raw[Annotation == "auto - shared MS/MS", Use := F]
  if (!use.discordant.peptide.type) DT.raw[Annotation == "auto - discordant peptide type", Use := F]
  if (!use.no.itraq) DT.raw[Annotation == "auto - no iTRAQ", Use := F]
  if (!use.weak.signal) DT.raw[Annotation == "no quant - weak signal", Use := F]

  # create wide data table
  if(!("ProteinModifications" %in% colnames(DT.raw))) DT.raw[, ProteinModifications := ""]
  DT <- DT.raw[, .(
    Group = gsub(";", "", Accessions),
    GroupInfo = paste0("[", N, "] ", Names),
    Component = gsub(" ", "", paste0(Sequence, ",", Modifications, ",", ProteinModifications, ",", Cleavages), fixed = T),
    Measurement = Spectrum,
    Injection = as.integer(matrix(unlist(strsplit(as.character(DT.raw$Spectrum), ".", fixed = T)), ncol = 5, byrow = T)[, 1]),
    Use
  )]
  if("Area 113" %in% colnames(DT.raw)) DT$Channel.113 <- DT.raw$`Area 113`
  if("Area 114" %in% colnames(DT.raw)) DT$Channel.114 <- DT.raw$`Area 114`
  if("Area 115" %in% colnames(DT.raw)) DT$Channel.115 <- DT.raw$`Area 115`
  if("Area 116" %in% colnames(DT.raw)) DT$Channel.116 <- DT.raw$`Area 116`
  if("Area 117" %in% colnames(DT.raw)) DT$Channel.117 <- DT.raw$`Area 117`
  if("Area 118" %in% colnames(DT.raw)) DT$Channel.118 <- DT.raw$`Area 118`
  if("Area 119" %in% colnames(DT.raw)) DT$Channel.119 <- DT.raw$`Area 119`
  if("Area 121" %in% colnames(DT.raw)) DT$Channel.121 <- DT.raw$`Area 121`
  rm(DT.raw)

  # we can get a spectrum assigned to multiple peptides in a protein - if this occurs, assign the spectrum as a 'new' ambiguous peptide
  # hack: data.table bug crashes R?
  setorder(DT, Component)
  DT[, Component2 := toString(Component), by = .(Group, Measurement, Use)]
  DT[, Component := NULL]
  setnames(DT, "Component2", "Component")
  DT <- unique(DT)

  # melt
  DT[, GroupInfo := factor(GroupInfo)]
  DT[, Group := factor(Group)]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT[, Run := factor("1")]
  DT[, Injection := factor(Injection)]
  DT <- melt(DT, variable.name = "Channel", value.name = "Count", measure.vars = colnames(DT)[grep("^Channel\\.", colnames(DT))])
  levels(DT$Channel) <- sub("^Channel\\.", "", levels(DT$Channel))

  setcolorder(DT, c("Group", "GroupInfo", "Component", "Measurement", "Injection", "Run", "Channel", "Count", "Use"))
  setorder(DT, Group, Component, Measurement, Injection, Run, Channel)

  setattr(DT, "group", c("ProteinGroup", "ProteinGroups"))
  setattr(DT, "component", c("Peptidoform", "Peptidoforms"))
  setattr(DT, "measurement", c("Spectrum", "Spectra"))

  setDF(DT)
  return(DT[])
}




#' Import MaxQuant LF data
#'
#' Reads in MaxQuant \code{evidence.txt} and \code{proteinGroups.txt}.
#'
#' @param proteinGroups.file Location of the \code{proteinGroups.txt} file.
#' @param evidence.file Location of the \code{evidence.txt} file.
#' @param shared Include shared peptides?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export
import_MaxQuant <- function(
  proteinGroups.file = NULL,
  evidence.file = NULL,
  use.decoys = FALSE,
  use.shared.peptides = FALSE,
  use.only.identified.by.site = FALSE,
  use.potential.contaminant = FALSE,
  proteinGroups.data = NULL,
  evidence.data = NULL
) {
  warning("import_MaxQuant currently only supports labelfree data with no fractionation")

  # load groups
  if (is.null(proteinGroups.file)) {
    if (is.null(proteinGroups.data)) stop("One of 'proteinGroups.data' or 'proteinGroups.file' needs to be specified.")
    DT.groups <- setDT(proteinGroups.data)
  } else {
    DT.groups <- fread(file = proteinGroups.file, showProgress = T)
  }

  # make sure we keep ALL groups, components and measurements for reference
  DT.groups[, Group := factor(`Protein IDs`, levels = unique(`Protein IDs`))]
  DT.groups[, Use := T]

  # filter groups
  if (!use.only.identified.by.site) DT.groups[`Only identified by site` == "+", Use := F]
  if (!use.decoys) DT.groups[`Reverse` == "+", Use := F]
  if (!use.potential.contaminant) DT.groups[`Potential contaminant` == "+", Use := F]

  DT.groups <- DT.groups[, .(
    Group,
    GroupInfo = factor(`Fasta headers`, levels = unique(`Fasta headers`)),
    GroupID = as.integer(id),
    Use
  )]

  # load wide raw
  if (is.null(evidence.file)) {
    if (is.null(evidence.data)) stop("One of 'evidence.data' or 'evidence.file' needs to be specified.")
    DT <- setDT(evidence.data)
  } else {
    DT <- fread(file = evidence.file, showProgress = T)
  }

  # make sure we keep ALL components and measurements for reference
  DT <- DT[, .(
    GroupID = as.character(`Protein group IDs`),
    MeasurementID = as.integer(id),
    Component = factor(`Modified sequence`, levels = unique(`Modified sequence`)),
    Measurement = paste0(`Modified sequence`, ",", Charge, "+"),
    Run = factor(Experiment, levels = unique(Experiment)),
    Channel = factor(1),
    Count = as.double(Intensity)
  )]
  DT[, Measurement := factor(Measurement, levels = unique(Measurement))]

  # fold out rows with shared features
  DT[, Use := ifelse(grepl(";", GroupID), use.shared.peptides, T)]
  DT = merge(DT[, !"GroupID"], DT[, {
    groupID = strsplit(GroupID, ";")
    list(MeasurementID = rep(MeasurementID, sapply(groupID, length)), GroupID = as.integer(unlist(groupID)))
  }], by = "MeasurementID", sort = F)

  # merge groups and raw
  DT <- merge(DT.groups, DT, by = "GroupID", sort = F, suffixes = c(".groups",""))
  DT[, GroupInfo := paste0("[", GroupID, "] ", GroupInfo)]
  DT[, GroupInfo := factor(GroupInfo, levels = unique(GroupInfo))]
  DT[, Use := Use & Use.groups]
  DT[, Use.groups := NULL]
  DT[, GroupID := NULL]
  DT[, MeasurementID := NULL]

  # sum multiple datapoints per measurement
  DT <- DT[, .(GroupInfo = GroupInfo[1], is.na = all(is.na(Count)), Count = sum(Count, na.rm = T)), by = .(Group, Component, Measurement, Run, Channel, Use)]
  DT[is.na == T, Count := NA_real_]
  DT[, is.na := NULL]

  setcolorder(DT, c("Group", "GroupInfo", "Component", "Measurement", "Run", "Channel", "Count"))
  setorder(DT, Group, Component, Measurement, Run, Channel)

  setattr(DT, "group", c("ProteinGroup", "ProteinGroups"))
  setattr(DT, "component", c("Peptidoform", "Peptidoforms"))
  setattr(DT, "measurement", c("Feature", "Features"))

  setDF(DT)
  return(DT)
}



























#' Import DIA-NN data
#'
#' Reads in DIA-NN \code{report.tsv} for processing.
#'
#' @param proteinGroups.file Location of the \code{proteinGroups.txt} file.
#' @param evidence.file Location of the \code{evidence.txt} file.
#' @param shared Include shared peptides?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export
import_DIANN <- function(
  file = NULL,
  use.shared.peptides = FALSE,
  data = NULL
) {
  stop("todo: not implemented yet")

  # load wide raw
  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT <- setDT(data)
  } else {
    DT <- fread(file = file, showProgress = T)
  }

  setDF(DT)
  return(DT)
}


#' Import Waters Progenesis data
#'
#' Reads in a Waters Progenesis \code{pep_ion_measurements.csv} file.
#'
#' @param file Location of the \code{pep_ion_measurements.csv} file.
#' @param used Include only measurements marked as used by Progenesis?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seaMass_sigma}.
#' @import data.table
#' @export

import_Progenesis <- function(
  file = NULL,
  used = T,
  data = NULL
) {
  stop("todo: needs updating")
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
  #if ("Accession" %in% colnames(DT.raw)) {
    colname.group <- "Accession"
    colname.groupinfo <- "Description"
    colname.sequence <- "Sequence"
    colname.modifications <- "Modifications"

    # only use rows that Progenesis uses for quant
    if (used) DT.raw <- DT.raw[`Use in quantitation` == "True",]
  #} else {
  #  colname.group <- colnames(DT.raw)[grep("Best peptide match Group$", colnames(DT.raw))]
  #  colname.groupinfo <- colnames(DT.raw)[grep("Best peptide match Description$", colnames(DT.raw))]
  #  colname.sequence <- colnames(DT.raw)[grep("Best peptide match Sequence$", colnames(DT.raw))]
  #  colname.modifications <- colnames(DT.raw)[grep("Best peptide match Variable modifications \\(\\[position\\] description\\)$", colnames(DT.raw))]
  #}

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
  # hack: data.table bug crashes R?
  setorder(DT, Component)
  DT[, Component2 := toString(Component), by = .(Group, Measurement)]
  DT[, Component := NULL]
  setnames(DT, "Component2", "Component")
  DT <- unique(DT)

  # melt
  DT[, GroupInfo := factor(GroupInfo)]
  DT[, Group := factor(Group)]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT <- melt(DT, variable.name = "Run", value.name = "Count", measure.vars = colnames(DT)[grep("^Raw abundance ", colnames(DT))])
  levels(DT$Run) <- sub("^Raw abundance ", "", levels(DT$Run))
  DT[, Channel := factor("1")]
  DT[, Injection := Run]
  DT[, Count := as.numeric(Count)]

  setDF(DT)
  return(DT[])
}




#' Import MaxQuant LF data
#'
#' Reads in MaxQuant \code{peptides.txt} and \code{proteinGroups.txt}.
#'
#' @param proteinGroups.file Location of the \code{proteinGroups.txt} file.
#' @param peptides.file Location of the \code{peptides.txt} file.
#' @param shared Include shared peptides?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export
import_MaxQuant_peptides <- function(
  proteinGroups.file = NULL,
  peptides.file = NULL,
  shared = FALSE,
  proteinGroups.data = NULL,
  peptides.data = NULL
) {
  stop("todo: needs updating")
  # load groups
  if (is.null(proteinGroups.file)) {
    if (is.null(proteinGroups.data)) stop("One of 'proteinGroups.data' or 'proteinGroups.file' needs to be specified.")
    DT.groups <- setDT(proteinGroups.data)
  } else {
    DT.groups <- fread(file = proteinGroups.file, showProgress = T)
  }

  # filter groups
  #DT.groups <- DT.groups[`Only identified by site` != "+",]
  DT.groups <- DT.groups[`Reverse` != "+",]
  DT.groups <- DT.groups[`Potential contaminant` != "+",]
  DT.groups <- DT.groups[, .(Group = `Protein IDs`, GroupInfo = `Fasta headers`, GroupID = id)]

  # load wide raw
  if (is.null(peptides.file)) {
    if (is.null(peptides.data)) stop("One of 'peptides.data' or 'peptides.file' needs to be specified.")
    DT.raw <- setDT(peptides.data)
  } else {
    DT.raw <- fread(file = peptides.file, showProgress = T)
  }

  # filter raw
  DT <- DT.raw[ , .(
    GroupID = `Protein group IDs`,
    ComponentID = id,
    #Component = paste(id, Sequence, sep = ","),
    Component = `Protein group IDs`,
    Measurement = paste(id, Sequence, sep = ",")
  )]
  DT.raw <- DT.raw[, colnames(DT.raw)[grep("^Intensity ", colnames(DT.raw))], with = F]
  colnames(DT.raw) <- sub("^Intensity ", "", colnames(DT.raw))
  DT <- cbind(DT, DT.raw)

  # remove or expand out rows with shared features
  if (shared == F) {
    DT <- DT[!grepl(";", GroupID),]
    DT[, GroupID := as.integer(GroupID)]
  } else {
    DT = merge(DT[, !"GroupID"], DT[, {
      groupID = strsplit(GroupID, ";")
      list(ComponentID = rep(ComponentID, sapply(groupID, length)), GroupID = as.integer(unlist(groupID)))
    }], by = "ComponentID")
  }

  # merge groups and raw
  DT <- merge(DT.groups, DT, by = "GroupID")
  DT[, GroupInfo := paste0("[", GroupID, "] ", GroupInfo)]
  DT[, GroupID := NULL]
  DT[, ComponentID := NULL]

  # melt to long data table and convert zeros to NA
  DT <- melt(DT, variable.name = "Run", value.name = "Count", id.vars = c("Group", "GroupInfo", "Component", "Measurement"))
  DT$Count[DT$Count == 0] <- NA
  DT[, GroupInfo := factor(GroupInfo)]
  DT[, Group := factor(Group)]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT[, Run := factor(Run)]
  DT[, Channel := factor("1")]
  DT[, Injection := Run]
  DT[, Count := as.double(Count)]
  setcolorder(DT, c("Group", "GroupInfo", "Component", "Measurement", "Run", "Channel", "Injection"))

  setDF(DT)
  return(DT)
}


#' Import MaxQuant LF data
#'
#' Reads in MaxQuant \code{evidence.txt} and \code{proteinGroups.txt}.
#'
#' @param proteinGroups.file Location of the \code{proteinGroups.txt} file.
#' @param evidence.file Location of the \code{evidence.txt} file.
#' @param shared Include shared peptides?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export
import_MaxQuant_evidence0 <- function(
  proteinGroups.file = NULL,
  evidence.file = NULL,
  shared = F,
  rollup = "measurement",
  proteinGroups.data = NULL,
  evidence.data = NULL
) {
  stop("todo: needs updating")
  suppressWarnings(suppressMessages(library(R.oo)))

  if (is.null(evidence.file) && is.null(evidence.data)) stop("One of 'evidence.data' or 'evidence.file' needs to be specified.")
  if (!is.null(evidence.data)) evidence.file <- evidence.data

  if (is.null(proteinGroups.file) && is.null(proteinGroups.data)) stop("One of 'proteinGroups.data' or 'proteinGroups.file' needs to be specified.")
  if (!is.null(proteinGroups.data)) proteinGroups.file <- proteinGroups.data

  #Use proteinGroups.txt and evidence.txt files to match MS/MS IDs to protein groups
  if (is.null(proteinGroups.file)) {
    DT.proteins <- setDT(proteinGroups.data)
  } else {
    DT.proteins <- fread(file = proteinGroups.file, showProgress = T)
  }

  #DT.proteins <- DT.proteins[`Only identified by site` != "+",]
  DT.proteins <- DT.proteins[`Reverse` != "+",]
  DT.proteins <- DT.proteins[`Potential contaminant` != "+",]

  DT.proteins <- unique(DT.proteins[, c("Protein IDs", "Majority protein IDs", "id")])
  DT.proteins[, `Protein group IDs` := id]
  DT.proteins[, id := NULL]

  DT <- rbindlist(lapply(evidence.file, function(f) {
    if (is.data.frame(f)) {
      DT.raw <- as.data.file(evidence.data)
    } else {
      DT.raw <- fread(file = f, showProgress = T)
    }

    # remove decoys and > max.m_score
    #DT.raw <- DT.raw[decoy == 0,]
    #DT.raw <- DT.raw[m_score <= max.m_score,]

    #DT.raw <- DT.raw[`MS/MS IDs` %in% DT.msms[,`MS/MS IDs`],]
    #DT.raw <- DT.raw[`MS/MS IDs` != "",]
    #DT.raw[, `MS/MS IDs` := factor(`MS/MS IDs`)]
    #DT.raw <- merge(DT.raw, DT.msms, by = "MS/MS IDs")
    #DT.msms <- NULL

    #MaxQuant labels unmodified peptides as "Unmodified"
    DT.raw[, Modifications := gsub("Unmodified", "", Modifications)]
    DT.raw[, Measurement := gsub(" ", "", paste0(Sequence, ",", Modifications, ",", Charge))]

    DT.raw[,Shared := grepl(";", `Protein group IDs`)]
    pGroupIDs <- unique(DT.raw[, c("Protein group IDs", "Measurement")])[, .(`Protein group IDs` = tstrsplit(`Protein group IDs`, ";", fixed = T)), by = .(Measurement)]
    DT.raw[, `Protein group IDs` := NULL]
    DT.raw <- merge(DT.raw, pGroupIDs, by = "Measurement")
    DT.raw[, `Protein group IDs` := strtoi(`Protein group IDs`)]

    # only use rows that MaxQuant uses for quant
    if (!shared) DT.raw <- DT.raw[Shared == F,]
    #if (!used) DT.raw <- DT.raw[`Peptide Quan Usage` == "Use",]

    DT.raw <- merge(DT.raw, DT.proteins, by = "Protein group IDs")

    # create long data table
    DT <- DT.raw[, .(
      Group = `Protein IDs`, #From the proteinGroups file
      Component = gsub(" ", "", paste0(Sequence, ",", Modifications)),
      Measurement = Measurement,
      Run = Experiment,
      Count = Intensity
    )]
    #DT <- DT[, lapply(.SD, function(x) unlist(tstrsplit(x, ";bayesprot;", fixed = T)))]
    DT[, Count := as.numeric(Count)]
    DT

  }))

  #Do rollup
  if (rollup == "measurement") {
    DT[, Count := sum(Count, na.rm = !all(is.na(Count))), by = .(Run, Measurement)]
    DT <- unique(DT)
  }

  assays <- unique(DT$Run)

  # create wide data table
  DT <- dcast(DT, Group + Component + Measurement ~ Run, value.var = "Count")

  # group ambiguous transitions so seaMass-Delta treats them as a single component per group ## UNNECCESARY?
  # hack: data.table bug crashes R?
  setorder(DT, Component)
  DT[, Component2 := toString(Component), by = .(Group, Measurement)]
  DT[, Component := NULL]
  setnames(DT, "Component2", "Component")
  DT <- unique(DT)

  # melt
  DT[, GroupInfo := factor(Group)]
  #DT[, Group := factor(sub("^1/", "", Group))]
  DT[, Group := factor(Group)]
  DT[, Component := factor(Component)]
  DT[, Measurement := factor(Measurement)]
  DT <- melt(DT, variable.name = "Run", value.name = "Count", measure.vars = assays)
  DT[, Injection := Run]
  DT[, Channel := factor("1")]
  setcolorder(DT, c("Group", "GroupInfo", "Component", "Measurement", "Run", "Injection", "Channel"))

  setDF(DT)
  return(DT[])

}


#' Import MaxQuant LF data
#'
#' Reads in MaxQuant \code{evidence.txt} and \code{proteinGroups.txt}.
#'
#' @param proteinGroups.file Location of the \code{proteinGroups.txt} file.
#' @param evidence.file Location of the \code{evidence.txt} file.
#' @param shared Include shared peptides?
#' @param data Advanced: Rather than specifying \code{file}, you can enter a \link{data.frame} preloaded with
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{bayesprot}.
#' @import data.table
#' @export
import_MSqRob <- function(
  fData,
  exprs,
  evidence.file = NULL,
  is.log2 = TRUE,
  evidence.data = NULL
) {
  stop("todo: needs updating")
  if (is.null(evidence.file) & is.null(evidence.data)) {
    DT <- setDT(cbind(fData[, c("Proteins", "Sequence")], exprs))
    DT <- melt(DT, variable.name = "Run", value.name = "Count", id.vars = c("Proteins", "Sequence"))
    setnames(DT, c("Proteins", "Sequence"), c("Group", "Measurement"))
    DT[, Group := factor(Group)]
    DT[, GroupInfo := ""]
    DT[, Component := Group]
    DT[, Measurement := factor(Measurement)]
    DT[, Channel := factor("1")]
    if (is.log2) DT[, Count := 2^Count]
  } else {
    # if MaxQuant evidence file supply, convert from two to three stage hierarchy
    DT <- setDT(fData[, c("Proteins", "Sequence","Evidence.IDs")])
    setnames(DT, "Proteins", "Group")
    DT <- DT[, {
      id = strsplit(as.character(Evidence.IDs), ";")
      list(Group = rep(Group, sapply(id, length)), Sequence = rep(Sequence, sapply(id, length)), id = as.integer(unlist(id)))
    }]

    if (!is.null(evidence.file) | is.null(evidence.data)) {
      if (is.null(evidence.file)) {
        DT.evidence <- setDT(evidence.data)
      } else {
        DT.evidence <- fread(file = evidence.file, showProgress = T)
      }
    }

    DT.evidence <- DT.evidence[, .(
      id,
      GroupInfo = "",
      Component = `Modified sequence`,
      Measurement = paste0(`Modified sequence`, ",", Charge, "+"),
      Fraction,
      Run = Experiment,
      Channel = "",
      Count = Intensity
    )]

    DT <- droplevels(merge(DT, DT.evidence, by = "id")[, !"id"])

    # create wide data table (summing up multiple features per measurement)
    DT <- dcast(DT, Group + GroupInfo + Component + Measurement + Fraction ~ Run, fun.aggregate = sum, value.var = "Count")

    # melt to long data table and convert zeros to NA
    DT <- melt(DT, variable.name = "Run", value.name = "Count", id.vars = c("Group", "GroupInfo", "Component", "Measurement", "Fraction"))
    DT$Count[DT$Count == 0] <- NA
    DT[, GroupInfo := factor(GroupInfo)]
    DT[, Group := factor(Group)]
    DT[, Component := factor(Component)]
    DT[, Measurement := factor(Measurement)]
    DT[, Fraction := factor(Fraction)]
    DT[, Run := factor(Run)]
    DT[, Channel := factor("1")]
    DT[, Count := as.double(Count)]
  }

  setcolorder(DT, c("Group", "GroupInfo", "Component", "Measurement", "Run"))
  setDF(DT)
  return(DT)
}


#' Import data outputed by an MSstats import routine
#'
#' Reads in a set of \code{_with_dscore} datasets processed by OpenSWATH and PyProphet.
#'
#' @param data MSstats output
#'   \link[data.table]{fread} default parameters.
#' @return A \link{data.frame} for input into \link{seaMass_sigma}.
#' @import data.table
#' @export
import_MSstats <- function(data) {
  stop("todo: needs updating")
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
  return(DT[])
}
