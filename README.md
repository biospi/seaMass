# BayesProt v1.2.0 (beta)
Bayesian mixed-effects model and uncertainty propagation for mass spectrometry proteomics, achieving sensitive protein-level quantification and differential expression analysis. Currently works with imported data from SCIEX ProteinPilot, Thermo ProteomeDiscoverer, Waters Progenesis and OpenSWATH/PyProphet across iTraq/TMT, SILAC, Label-free and SWATH data types.

## Current citation
Xu et al, Nature Communications Biology 2019, 2:43, [https://doi.org/10.1038/s42003-018-0254-9]

## Requirements

BayesProt is an R package that works on Windows, MacOS and Linux. Small studies (n=4 vs n=4) will take about one to two hours to process with default control parameters. Large studies (e.g. 100 samples) could take overnight or longer. Small studies will run with 16Gb of memory, but larger studies could take 64Gb or more. To reduce memory usage at the expense of speed, reduce the number of CPU threads via the 'nthread' parameter of the BayesProt 'control' object (see tutorial below).

## Installation

```
install.packages("remotes")
remotes::install_github("biospi/bayesprot", dependencies = T)
```

## Usage

Firstly, you need to use an 'import' function to convert from an upstream tool format to BayesProt's standardised 'data.frame' format. Then you can assign runs if you've used fractionation and specify a study design if you'd like to do differential expression analysis. Finally, you use the 'bayesprot' function to fit the model and generate the results.

### Tutorial

Load the included ProteinPilot iTraq dataset (note this is a small subset of the fractions from our spike-in study and as such is not useful for anything else than this tutorial).

```
library(bayesprot)

# import tutorial iTraq dataset
file <- system.file(file.path("demo", "Tutorial_PeptideSummary.txt.bz2"), package = "bayesprot")
data <- import_ProteinPilot(file)
```

Unfortunately the input file does not contain information for linking fractions to runs, so BayesProt allows you to update the imported data with this information. If your study does not employ fractionation, you can skip this section.

```
# get skeleton injection-run table from imported data
data.runs <- runs(data)

# indicate which injection refers to which run
data.runs$Run <- factor(c(rep_len(1, 10), rep_len(2, 10)))

# update the imported data with this information
runs(data) <- data.runs
```

Next, you can give the assays the names you prefer, and assign samples to assays (in case you have technical replicates, for example).

```
# get skeleton design matrix
data.design <- new_design(data)

# you can give the assays different names
data.design$Assay <- factor(c(
  "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8",
  "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8"
))

# you can assign samples to assays (here we have 6 'A' samples, 6 'B' samples and 4 'O' samples - others)
data.design$Sample <- factor(c(
  "A1", "A2", "B1", "B2", "A3", "A4", "B3", "B4",
  "O1", "O2", "O3", "O4", "A5", "A6", "B5", "B6"
))
```

Optionally, you can do differential expression analysis (or any mixed-effects model supported by the 'metafor' package). You can do it post-hoc using the 'fit' object returned by 'bayesprot', or you can supply one or more functions to run during the 'bayesprot' call. There is a pre-canned function for doing a Student's t-test (as this test uses the uncertainty in the quants, it has much better performance than doing the t-tests yourself!), which needs you to specify a 'Condition' column in your experiment design:

```
# specify the pre-canned t-test function
de.func <- protein_de_ttest

# specify the conditions for differential expression analysis (use NA to ignore irrelevant assays)
data.design$Condition <- factor(c(
  "A", "A", "B", "B", "A", "A", "B", "B",
  NA,  NA,  NA,  NA, "A", "A", "B", "B")
)
```

If you have multiple channels per run (e.g. iTraq, TMT, SILAC) you need to specify the reference assay(s) for each run so that quants can be linked between runs. BayesProt does not need specific reference samples (e.g. pooled samples) to have been run as long as you have a blocked design with the same proportion of samples in each treatment group in each run. For example, below we have designated 4 A and 4 B samples as reference channels for Run 1, and 2 A and 2 B samples as reference channels for Run 2.

```
data.design$ref <- factor(c(
  T, T, T, T, T, T, T, T,
  F, F, F, F, T, T, T, T)
)
```

Finally, run the model. Intermediate and results data is stored in the directory specified by the BayesProt 'output' parameter, and any internal control parameters (such as the number of CPU threads to use) can be specified through a 'control' object'. For differential expression analysis, a Bayesian version of median normalisation will be used. By default all proteins are considered in the normalisation, but you can choose a subset if required with the 'normalisation.proteins' parameter as shown below.

```
# run BayesProt
fit <- bayesprot(
  data,
  data.design = data.design,
  normalisation.proteins = levels(data$Protein)[grep("_RAT$", levels(data$Protein))],
  de.func = de.func,
  output = "Tutorial.bayesprot",
  control = new_control(nthread = 16)
)
```

Once run, results tables (in csv format) and diagnostic plots (PCA, exposures aka 'normalisation coefficients') are available in the 'output' subdirectory of the output directory specified above. Or you can use the R package functions to retrieve results and generate plots from the fit object created by the 'bayesprot' call:

```
# output list of proteins analysed
print(proteins(fit))

# output protein quants
print(protein_quants(fit))

# output fdr-controlled differential expression
print(protein_de(fit))

# view the plot for the top differential expression candidate
plot_peptides(fit, data.de$ProteinID[5])
```
