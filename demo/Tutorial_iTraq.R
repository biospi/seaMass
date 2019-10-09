library(bayesprot)

#######################################

# Import tutorial iTraq dataset.
file <- system.file(file.path("demo", "Tutorial_PeptideSummary.txt.bz2"), package = "bayesprot")
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

# You can rename assays, or remove them from the analysis with 'NA'.
data.design$Assay <- factor(c(
  NA, NA, NA, NA, "1.1", "1.2", "1.3", "1.4",
  NA, NA, NA, NA, "2.1", "2.2", "2.3", "2.4"
))

# Assign samples to assays (a 1-to-1 mapping if all your samples are separately digested;
#  pure technical replicates should share the same sample name).
data.design$Sample <- factor(c(
  NA, NA, NA, NA, "A1", "A2", "B1", "B2",
  NA, NA, NA, NA, "A3", "A4", "B3", "B4"
))

#######################################

# specify a list of one of more differential expression analysis functions. Bayesprot currently
#  implements tests between conditions using the 'dea_MCMCglmm' function. By default
#  these are t.tests but you can add covariates, random effects etc using the 'MCMCglmm::MCMCglmm'
#  syntax
dea.func <- list(t.tests = dea_limma)

# 'dea_MCMCglmm' expects a column 'Condition' to have been specified in 'data.design'.
#  You can use 'NA' to ignore irrelevant samples.
data.design$Condition <- factor(c(
  NA, NA, NA, NA, "A", "A", "B", "B",
  NA, NA, NA, NA, "A", "A", "B", "B"
))

#######################################

# iTraq/TMT/SILAC only: Since we have more than one iTraq run we need to normalise across them.
#  Here we can specify specific reference assays e.g. pooled reference samples, of if the study
#  design has been blocked appropriately, we can just use all relevant assays:
data.design$ref <- factor(c(
  F, F, F, F, T, T, T, T,
  F, F, F, F, T, T, T, T)
)

#######################################

# By default BayesProt uses median normalisation. If you want to normalise just to a specific
# set of proteins, do this (here to all the rat proteins).
norm_truth <- function(...) {
  norm_median(..., ref.proteins = levels(data$Protein)[grep("_RAT", levels(data$Protein))])
}

#######################################

# run BayesProt.
fit <- bayesprot(
  data,
  data.design = data.design,
  norm.func = list(truth = norm_truth),
  dea.func = dea.func,
  output = "Tutorial.bayesprot",
  control = new_control(nthread = 4)
)

#######################################

# Output list of proteins analysed.
data.proteins <- proteins(fit)
print(data.proteins)

# Output processed design matrix
data.design <- design(fit)
print(data.design)

# Output protein quants with accessions and assay/sample names
data.protein.quants <- protein_quants(fit)
data.protein.quants <- merge(data.protein.quants, data.design[, c("AssayID", "Assay", "Sample")])
data.protein.quants <- merge(data.protein.quants, data.proteins[, c("ProteinID", "Protein")])
print(data.protein.quants)

# Output fdr-controlled differential expression for the 't.test' analysis, with accessions.
data.de <- protein_fdr(fit)
data.de <- merge(data.de, data.proteins[, c("ProteinID", "Protein")], sort = F)
print(data.de)

# Plot precision-recall curve (sensitivity against false discovery rate.
data.de$truth <- ifelse(grepl("_RAT", data.de$Protein), 0, 1)
plot_pr(data.de)
