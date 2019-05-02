library(bayesprot)

# import tutorial iTraq dataset
file <- system.file(file.path("demo", "Tutorial_PeptideSummary.txt.bz2"), package = "bayesprot")
data <- import_ProteinPilot(file)

# get skeleton injection-run table from imported data
data.runs <- runs(data)
# indicate which injection refers to which run
data.runs$Run <- factor(c(rep_len(1, 10), rep_len(2, 10)))
# update the imported data with this information
runs(data) <- data.runs

# get skeleton design matrix
data.design <- design(data)
# you can assign samples to assays
data.design$Sample <- factor(c(
  "A1", "A2", "B1", "B2", "A3", "A4", "B3", "B4",
  "O1", "O2", "O3", "O4", "A5", "A6", "B5", "B6"
))
# specify the conditions for differential expression analysis (use NA to ignore irrelevant channels)
data.design$Condition <- factor(c(
  "A", "A", "B", "B", "A", "A", "B", "B",
   NA,  NA,  NA,  NA, "A", "A", "B", "B")
)
# specify the reference channels
data.design$ref <- !is.na(data.design$Condition)

# run BayesProt
fit <- bayesprot(
  data,
  data.design = data.design,
  normalisation.proteins = levels(data$Protein)[grep("_RAT$", levels(data$Protein))],
  output = "Tutorial.bayesprot",
  control = control(nthread = 4)
)

# get protein-level quants
dd.quants <- protein_quants(fit)
# get FDR-controlled differential expression
dds.de <- de_metafor(fit)

