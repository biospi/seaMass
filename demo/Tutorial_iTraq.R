library(deamass)

#######################################

# Import tutorial iTraq dataset.
file <- system.file(file.path("demo", "Tutorial_PeptideSummary.txt.bz2"), package = "deamass")
data <- import_ProteinPilot(file)

#######################################

# Get skeleton injection-run table from imported data.
data.runs <- runs(data)

# Indicate which injection refers to which run; use 'NA' to indicate injections to ignore.
data.runs$Run[1:26] <- NA
data.runs$Run[27:52] <- "1"
data.runs$Run[53:78] <- "2"

# Update the imported data with the run information and remove any ignored injections.
runs(data) <- data.runs

#######################################

# Get skeleton design matrix
data.design <- new_design(data)
print(data.design)

# You can rename assays, or remove them from the analysis with 'NA'.
data.design$Assay <- factor(c(
  NA, NA, NA, NA, "1.1", "1.2", "1.3", "1.4",
  NA, NA, NA, NA, "2.1", "2.2", "2.3", "2.4"
))

#######################################

# If you have a large amount of assays, you may group them into 'Blocks' so that deaMass will be
#  robust to batch effects or measurement drift etc. This will also help memory consumption, as
#  blocks are processed by deaMass independently. If you are processing iTraq or TMT data, each
#  run will also be treated as a seperate block, which deaMass will autodetect. To compare across
#  blocks correctly, you need to normalise them via reference assays. Here we can specify
#  specific reference assays e.g. pooled reference samples, of if the study design has been
#  blocked appropriately, we can just use all relevant assays:
data.design$BlockRef <- factor(c(
  F, F, F, F, T, T, T, T,
  F, F, F, F, T, T, T, T)
)

#######################################

# specify a list of one of more differential expression analysis functions. deaMass currently
#  implements tests between conditions using the 'dea_limma' and the 'dea_MCMCglmm' functions.
# By default these are standard t.tests for limma and Welch's t.tests for MCMCglmm (limma does
#  not support Welch's t.tests), but you can add covariates.
dea.func <- list(t.tests = dea_limma)

# Assign samples to assays (a 1-to-1 mapping if all your samples are separately digested;
#  pure technical replicates should share the same sample name).
data.design$Sample <- factor(c(
  NA, NA, NA, NA, "A1", "A2", "B1", "B2",
  NA, NA, NA, NA, "A3", "A4", "B3", "B4"
))

# 'dea_limma' expects a column 'Condition' to have been specified in 'data.design'.
#  You can use 'NA' to ignore irrelevant samples.
data.design$Condition <- factor(c(
  NA, NA, NA, NA, "A", "A", "B", "B",
  NA, NA, NA, NA, "A", "A", "B", "B"
))

#######################################

# By default deaMass uses median normalisation. If you want to normalise just to a specific
# set of protein groups, do this (here to all the rat protein groups).
norm_truth <- function(...) {
  norm_median(..., ref.groups = levels(data$Group)[grep("_RAT", levels(data$Group))])
}

#######################################

# run deaMass.
fit <- deamass(
  data,
  data.design = data.design,
  norm.func = list(truth = norm_truth),
  dea.func = dea.func,
  output = "Tutorial.deamass",
  measurement.vars = T,
  component.vars = T,
  component.deviations = T,
  control = new_control(nthread = 4)
)

#######################################

# Output list of groups analysed.
data.groups <- groups(fit)
print(data.groups)

# Output processed design matrix
data.design <- design(fit)
print(data.design)

# Output group quants with accessions and assay/sample names
data.group.quants <- group_quants(fit)
data.group.quants <- merge(data.group.quants, data.design[, c("AssayID", "Assay", "Sample")])
data.group.quants <- merge(data.group.quants, data.groups[, c("GroupID", "Group")])
print(data.group.quants)

# Output fdr-controlled differential expression for the 't.test' analysis, with accessions.
data.de <- group_fdr(fit)
data.de <- merge(data.de, data.groups[, c("GroupID", "Group")], sort = F)
print(data.de)

# Plot precision-recall curve (sensitivity against false discovery rate.
data.de$truth <- ifelse(grepl("_RAT", data.de$Group), 0, 1)
plot_pr(data.de)
