# BayesProt v2.0.0 (beta1 - 29th October 2019)
Bayesian mixed-effects model and uncertainty propagation for mass spectrometry proteomics, achieving sensitive protein-level quantification and differential expression analysis. Currently works with imported data from SCIEX ProteinPilot, Thermo ProteomeDiscoverer, Waters Progenesis and OpenSWATH/PyProphet across iTraq/TMT, SILAC, Label-free and SWATH data types.

## Current citation
Xu et al, Nature Communications Biology 2019, 2:43, [https://doi.org/10.1038/s42003-018-0254-9]

## Requirements

BayesProt is an R package that works on Windows, MacOS and Linux. Small studies (n=4 vs n=4) will take about one to two hours to process with default control parameters. Large studies (e.g. 100 samples) could take overnight or longer. Small studies will run with 16Gb of memory, but larger studies could take 64Gb or more. To reduce memory usage at the expense of speed, reduce the number of CPU threads via the 'nthread' parameter of the BayesProt 'control' object (see tutorial below).

## Installation

Please install the R package directly from our Github repository (use the same command to upgrade to the latest version):

```
install.packages("devtools")
devtools::install_github("biospi/bayesprot", dependencies = TRUE)
```

## Usage

Firstly, you need to use an 'import' function to convert from an upstream tool format to BayesProt's standardised 'data.frame' format. Then you can assign injections to runs if you've used fractionation, and specify a study design if you'd like to do differential expression analysis. Finally, you use the 'bayesprot' function to fit the model and generate the results.

### Tutorial

Load the included ProteinPilot iTraq dataset (note this is a small subset of proteins from our spike-in study and as such is not useful for anything else than this tutorial).

```
library(bayesprot)

# Import tutorial iTraq dataset.
file <- system.file(file.path("demo", "Tutorial_PeptideSummary.txt.bz2"), package = "bayesprot")
data <- import_ProteinPilot(file)
```

Unfortunately the input file does not contain information for linking fractions to runs, so BayesProt allows you to update the imported data with this information. If your study does not employ fractionation, you can skip this section.

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

Optionally, you can do differential expression analysis (or any mixed-effects model supported by the 'metafor' package). You can run a 'dea' function post-hoc using the 'fit' object returned by 'bayesprot', or you can supply one or more 'dea' functions to run during the 'bayesprot' call. If your experiment contains two or more treatment groups, we suggest you use a 'pairwise' function (only function implemented at present), which will perform the model seperately for each pair of conditions specified by a 'Condition' column in your experiment design:

```
# specify a list of one of more differential expression analysis functions. Bayesprot currently
#  implements tests between conditions in a pairwise fashion using 'limma'
dea.func <- list(t.tests = dea_limma)

# 'dea_limma' expects a column 'Condition' to have been specified in 'data.design'.
#  You can use 'NA' to ignore irrelevant samples.
data.design$Condition <- factor(c(
  NA, NA, NA, NA, "A", "A", "B", "B",
  NA, NA, NA, NA, "A", "A", "B", "B"
))
```

If you have multiple channels per run (e.g. iTraq, TMT, SILAC) you need to specify the reference assay(s) for each run so that quants can be linked between runs. BayesProt does not need specific reference samples (e.g. pooled samples) to have been run as long as you have a blocked design with the same proportion of samples in each treatment group in each run. For example, below we have designated 4 A and 4 B samples as reference channels for Run 1, and 2 A and 2 B samples as reference channels for Run 2.

```
# iTraq/TMT/SILAC only: Since we have more than one iTraq run we need to normalise across them.
#  Here we can specify specific reference assays e.g. pooled reference samples, of if the study
#  design has been blocked appropriately, we can just use all relevant assays:
data.design$ref <- factor(c(
  F, F, F, F, T, T, T, T,
  F, F, F, F, T, T, T, T)
)
```

Finally, run the model. Intermediate and results data is stored in the directory specified by the BayesProt 'output' parameter, and any internal control parameters (such as the number of CPU threads to use) can be specified through a 'control' object'. For differential expression analysis, a Bayesian version of median normalisation will be used. By default all proteins are considered in the normalisation, but you can choose a subset if required as illustrated below.

```
# By default BayesProt uses median normalisation. If you want to normalise just to a specific
# set of proteins, do this (here to all the rat proteins).
norm_truth <- function(...) {
  norm_median(..., ref.proteins = levels(data$Protein)[grep("_RAT", levels(data$Protein))])
}

# run BayesProt.
fit <- bayesprot(
  data,
  data.design = data.design,
  norm.func = list(truth = norm_truth),
  dea.func = dea.func,
  output = "Tutorial.bayesprot",
  control = new_control(nthread = 4)
)
```

Once run, results tables (in csv format) and diagnostic plots (PCA, exposures aka 'normalisation coefficients') are available in the 'output' subdirectory of the output directory specified above. Or you can use the R package functions to retrieve results and generate plots from the fit object created by the 'bayesprot' call:

```
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
```
