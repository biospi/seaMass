# seaMass-delta v1.0-0.0 (alpha1 - 30th October 2019)
Differential expression analysis for mass spectrometry proteomics and metabolomics using a mixed-effects model and uncertainty propagation, achieving sensitive protein and metabolite group-level quantification and differential expression analysis. Currently works with imported data from SCIEX ProteinPilot, Thermo ProteomeDiscoverer, Waters Progenesis and OpenSWATH/PyProphet across iTraq/TMT, SILAC, Label-free and SWATH data types.

## Current citation
Xu et al, Nature Communications Biology 2019, 2:43, [https://doi.org/10.1038/s42003-018-0254-9]

## Requirements

seaMass-delta is an R package that works on Windows, MacOS and Linux. Small studies (n=4 vs n=4) will take about one to two hours to process with default control parameters. Large studies (e.g. 100 samples) could take overnight or longer. Memory requirements are dependent on the amount of data in the largest block of assays. Typically 8-32Gb of memory is required. To reduce memory usage at the expense of speed, reduce the number of CPU threads via the 'nthread' parameter of the seaMass-delta 'control' object (see tutorial below).

## Installation

Please install the R package directly from our Github repository. To do this you first need to install the 'devtools' package: 

```
install.packages("devtools")
```

Then to install seaMass-delta, or upgrade to the latest version, simply run:

```
devtools::install_github("biospi/seamassdelta", dependencies = TRUE)
```

To upgrade to the latest version of seaMass-delta, simply run this line again.

## Usage

Firstly, you need to use an 'import' function to convert from an upstream tool format to seaMass-delta's standardised 'data.frame' format. Then you can assign injections to runs if you've used fractionation, and specify a study design if you'd like to do differential expression analysis. Finally, you use the 'seaMass-delta' function to fit the model and generate the results.

### Tutorial

Load the included ProteinPilot iTraq dataset (note this is a small subset of proteins from our spike-in study and as such is not useful for anything else than this tutorial).

```
library(seamassdelta)

# Import tutorial iTraq dataset.
file <- system.file(file.path("demo", "Tutorial_PeptideSummary.txt.bz2"), package = "seamassdelta")
data <- import_ProteinPilot(file)
```

Unfortunately the input file does not contain information for linking fractions to runs, so seaMass-delta allows you to update the imported data with this information. If your study does not employ fractionation, you can skip this section.

```
# Get skeleton injection-run table from imported data.
data.runs <- runs(data)

# Indicate which injection refers to which run; use 'NA' to indicate injections to ignore.
data.runs$Run[1:26] <- NA
data.runs$Run[27:52] <- "1"
data.runs$Run[53:78] <- "2"

# Update the imported data with the run information and remove any ignored injections.
runs(data) <- data.runs
```

Next, you can give the assays the names you prefer, and assign samples to assays (in case you have pure technical replicates, for example).

```
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
```

If you have a large amount of assays, you may group them into blocks so that seaMass-delta will be robust to batch effects or measurement drift etc. This will also help memory consumption, as blocks are processed by seaMass-delta independently. If you are processing iTraq or TMT data, each run will be treated as a seperate block by default, which seaMass-delta will autodetect. To compare across blocks correctly, you need to normalise them via reference assays. Here we can specify specific reference assays e.g. pooled reference samples, of if the study design has been blocked appropriately, we can just use all relevant assays:

```
# Define reference assays for each block
data.design$BlockRef <- c(
  F, F, F, F, T, T, T, T,
  F, F, F, F, T, T, T, T
)
```

Optionally, you can do differential expression analysis (or any mixed-effects model supported by the 'MCMCglmm' package). You can run a 'dea_' function post-hoc using the 'fit' object returned by 'seaMass-delta', or you can supply one or more 'dea_' functions to run during the 'seaMass-delta' call. The 'dea' functions will perform the model seperately for each pair of conditions specified by a 'Condition' column in your experiment design:

```
# specify a list of one of more differential expression analysis functions.
dea.func <- list(t.tests = dea_MCMCglmm)

# 'dea_' functions expect a column 'Condition' to have been specified in 'data.design'.
#  You can use 'NA' to ignore irrelevant samples.
data.design$Condition <- factor(c(
  NA, NA, NA, NA, "A", "A", "B", "B",
  NA, NA, NA, NA, "A", "A", "B", "B"
))
```

If you have multiple channels per run (e.g. iTraq, TMT, SILAC) you need to specify the reference assay(s) for each run so that quants can be linked between runs. seaMass-delta does not need specific reference samples (e.g. pooled samples) to have been run as long as you have a blocked design with the same proportion of samples in each treatment group in each run. For example, below we have designated 4 A and 4 B samples as reference channels for Run 1, and 2 A and 2 B samples as reference channels for Run 2.

```
# iTraq/TMT/SILAC only: Since we have more than one iTraq run we need to normalise across them.
#  Here we can specify specific reference assays e.g. pooled reference samples, of if the study
#  design has been blocked appropriately, we can just use all relevant assays:
data.design$BlockRef <- factor(c(
  F, F, F, F, T, T, T, T,
  F, F, F, F, T, T, T, T)
)
```

For differential expression analysis, a Bayesian version of median normalisation will be used. By default all proteins are considered in the normalisation, but you can choose a subset if required as illustrated below (here to all the rat proteins):

```
# custom normalisation function
norm_truth <- function(...) {
  norm_median(..., ref.groups = levels(data$Group)[grep("_RAT", levels(data$Group))])
}
```

Finally, run the model. Intermediate and results data is stored in the directory specified by the seaMass-delta 'output' parameter, and any internal control parameters (such as the number of CPU threads to use) can be specified through a 'control' object'. 

```
# run seaMass-delta
fit <- seamassdelta(
  data,
  data.design = data.design,
  norm.func = norm_truth,
  dea.func = dea.func,
  output = "Tutorial.seamassdelta",
  control = new_control(nthread = 4)
)
```

Once run, results tables (in csv format) and diagnostic plots are available in the 'output' subdirectory of the output directory specified above. Or you can use the R package functions to retrieve results and generate plots from the fit object created by the 'seaMass-delta' call:

```
# Output list of protein groups analysed.
data.groups <- groups(fit)
print(data.groups)

# Output processed design matrix
data.design <- design(fit)
print(data.design)

# Output protein group quants with accessions and assay/sample names
data.group.quants <- group_quants(fit)
data.group.quants <- merge(data.group.quants, data.design[, c("AssayID", "Assay", "Sample")])
data.group.quants <- merge(data.group.quants, data.groups[, c("GroupID", "Group")])
print(data.group.quants)

# Output fdr-controlled protein group differential expression for the 't.test' analysis, with accessions.
data.de <- group_fdr(fit)
data.de <- merge(data.de, data.groups[, c("GroupID", "Group")], sort = F)
print(data.de)

# Plot precision-recall curve (sensitivity against false discovery rate.
data.de$truth <- ifelse(grepl("_RAT", data.de$Group), 0, 1)
plot_pr(data.de)
```
