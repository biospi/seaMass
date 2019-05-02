# BayesProt v1.2.0 (beta)
Bayesian mixed-effects model and uncertainty propagation for mass spectrometry proteomics, achieving sensitive protein-level quantification and differential expression analysis. Currently works with imported data from SCIEX ProteinPilot, Thermo ProteomeDiscoverer and Waters Progenesis across iTraq/TMT, SILAC, Label-free and SWATH data.

## Current citation
Xu et al, Nature Communications Biology 2019, 2:43, [https://doi.org/10.1038/s42003-018-0254-9]

## Installation

```
install.packages("remotes")
remotes::install_github("biospi/bayesprot", dependencies = T)
```

## Usage

Firstly, you need to use an 'import' function to convert from an upstream tool format to BayesProt's standardised data.frame format. Then you can assign runs if you've used fractionation and specify a study design if you'd like to do differential expression analysis (pair-wise t-tests between conditions only at present). Finally, you use the 'bayesprot' function to fit the model and generate the results.

### Tutorial

Load the included ProteinPilot iTraq dataset (note this is a small subset of the fractions from our spike-in study and as such is not useful for anything else than this tutorial).

```
library(bayesprot)

# import tutorial iTraq dataset
file <- system.file(file.path("demo", "Tutorial_PeptideSummary.txt.bz2"), package = "bayesprot")
data <- import_ProteinPilot(file)
```

Unfortunately the input file does not contain information for linking fractions to runs. So BayesProt allows you to update the imported data with this information:

```
# get skeleton injection-run table from imported data
data.runs <- runs(data)
# indicate which injection refers to which run
data.runs$Run <- factor(c(rep_len(1, 10), rep_len(2, 10)))
# update the imported data with this information
runs(data) <- data.runs
```

Next, you can assign samples to conditions for pairwise t-testing, and if you are processing iTraq/TMT/SILAC data, you can define which samples are used as reference channels to link the multiplexes. BayesProt allows you to select multiple samples within each multiplex as reference channels - in this instance, the mean of these samples will be used.

```
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
```

Finally, run the model. Intermediate and results data is stored in the directory specified by 'output', and any internal control parameters (such as the number of CPU threads to use) can be specified through a 'control' object':

```
# run BayesProt
fit <- bayesprot(
  data,
  data.design = data.design,
  normalisation.proteins = levels(data$Protein)[grep("_RAT$", levels(data$Protein))],
  output = "Tutorial.bayesprot",
  control = control(nthread = 4)
)
```

Once run, results tables (in csv format) and diagnostic plots (PCA, exposures aka 'normalisation coefficients') are available in the output subdirectory of the output directory specified above. Or you can use the R package to retrieve results from the fit object created by the 'bayesprot' call:

```
# get protein-level quants
dd.quants <- protein_quants(fit)
# get FDR-controlled differential expression
dds.de <- de_metafor(fit)
```
