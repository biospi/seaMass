# seaMass
The seaMass suite of tools for quantification and differential expression analysis in mass spectrometry proteomics. The current status of the individual tools is:

- **seaMass-Σ**: A Bayesian protein group-level quantification technique univerally supporting label-free, SILAC, iTraq/TMT and DIA data. Currently we support proteomics input from Waters Progenesis (label-free), SCIEX ProteinPilot (iTraQ), Thermo ProteomeDiscoverer (SILAC/TMT) and OpenSWATH (SWATH). Other packages supported on request (MaxQuant coming soon). The model provides automatic quality control by downweighting problematic samples and peptides/features, can scale to massive study sizes, and propagates quantification uncertainty downstream to the differential expression analysis stage. 

- **seaMass-Θ**: Bayesian protein group normalisation on the output of seaMass-Σ.

- **seaMass-Δ**: Bayesian differential expression analysis on the output of seaMass-Σ or seaMass-Θ. By harnessing the generic [MCMCglmm](https://cran.r-project.org/web/packages/MCMCglmm) package for Bayesian mixed-effect modelling, the tool allows the user to perform many kinds of univariate analysis on the data, from simple Welch's t-tests and two-way ANOVA to timecourse and multi-level models. Studies analysed with an earlier versions of seaMass-Σ/Θ/Δ which was called BayesProt are published in [[Freeman et al, Diabetes, 2016]](	
https://doi.org/10.2337/db15-0835), [[Xu et al, Nature Comms Biology, 2019]](https://doi.org/10.1038/s42003-018-0254-9), [[Kassab et al, Molecular Metabolism, 2019]](	
https://doi.org/10.1016/j.molmet.2019.08.003) and [[Philbert et al, Biochem Biophys Res Comm, 2021]](https://doi.org/10.1016/j.bbrc.2021.02.107).

- **seaMass-α**: Still in development, this will provide sensitive feature extraction directly from raw mass spectrometry data. Please see [[Liao et al. IEEE ISBI, 2014]](https://doi.org/10.1109/ISBI.2014.6868123) and [[Zhang et al., Proteomics, 2015]](https://doi.org/10.1002/pmic.201400428) for technical details. The C++ codebase resides at [here](https://github.com/biospi/seaMass-alpha) and when ready will have an R interface through this package.

- **mzMLb**: A future-proof binary HDF5 encoding of the [mzML](http://www.psidev.info/mzML) standard interchange format for raw mass spectrometry data flexibly providing best-in-class compression or best-in-class random access speeds. The manuscript is [[Bhamber et al, 2020]](https://doi.org/10.1021/acs.jproteome.0c00192), the code is implemented in ProteoWizard [here](https://github.com/biospi/pwiz) and is moving towards approval by the [Proteomics Standards Intiative](http://www.psidev.info/).  

## Requirements

seaMass is an R package that works on Windows, MacOS and Linux, and supports local processing as well as on SGE, PBS and SLURM High Performance Computing (HPC) clusters. Small studies (n=4 vs n=4) will take about one to two hours to process with default control parameters. Large studies (e.g. 100 samples) could take overnight or longer unless run on HPC. Memory requirements are dependent on the amount of data in the largest block of assays (Each iTraq/TMT run is a separate block, whereas in label-free/SILAC/DIA the blocks are user-defined). Typically 16-32Gb of memory is required. To reduce memory usage at the expense of speed, reduce the number of CPU threads via the 'nthread' parameter of the *sigma_control* object (see tutorial below).

## Installation

Please install the R package directly from our Github repository. To do this you first need to install the *devtools* package: 

```
install.packages("devtools")
```

Then to install seaMass or upgrade to the latest version, simply run:

```
devtools::install_github("biospi/seaMass", dependencies = TRUE)
```

## seaMass-Σ, seaMass-Θ and seaMass-Δ

Firstly, you need to use an *import* function to convert from an upstream tool format to seaMass' standardised *data.frame* format. Then you can assign injections to runs if you've used fractionation, and specify a block structure to your assays if this is a large study (TMT/iTraq studies are blocked automatically). Finally, you use the *seaMass-sigma* function to fit the model and generate raw group-level quants.

### Tutorial 1 - basic analysis

Load the included MaxQuant labelfree Orbitrap dataset (note this is a small subset of proteins from our spike-in study and as such is not useful for anything else than this tutorial).

```
library(seaMass)

# Import MaxQuant tutorial dataset from the seaMass R package.
proteinGroups.file <- system.file("proteinGroups.txt.bz2", package = "seaMass")
evidence.file <- system.file("evidence.txt.bz2", package = "seaMass")
data <- import_MaxQuant(proteinGroups.file, evidence.file)
```

By default, differential expression analysis functions look for columns *Sample* and *Condition* in *data.design*. If you want to perform differential expression analysis, assign samples to assays, and conditions to samples:

```
# Get skeleton design matrix
data.design <- new_assay_design(data)

# Define Sample and Condition mappings
data.design$Sample <- c("A1", "A3", "A5", "A6", "B2", "B3", "B4", "B5")
data.design$Condition <- c("A", "A", "A", "A", "B", "B", "B", "B")
```

Now you can run the seaMass-Σ model for raw protein group-level quantification, followed by the seaMass-Θ normalisation model: 

```
# run seaMass-sigma and seaMass-theta
fit.sigma <- seaMass_sigma(data, data.design)
fit.theta <- seaMass_theta(fit.sigma)
```

A directory has been created called *fit.seaMass* with a *csv* subdirectory that contains numerical results and a *report.zip* file that contains an interactive html report full of quality control and results plots. It is strongly recommend to use a utility such as [Pismo File Mount](https://pismotec.com/pfm/ap/) on Windows, [Archive Mounter](https://github.com/biospi/ArchiveMounter) on MacOS or [fuse-zip](https://linux.die.net/man/1/fuse-zip) on Linux to mount the zip report as a drive for browsing, as the unzipped contents can be more than 20 times the size as the zip file itself. You can also get all these results in R. For example, you can output the raw or normalised protein group quants and generate a PCA plot coloured by 'Assay.SD', which is the unexplained variation in each assay - the higher this is for an assay, the more uncertain the results are (but note, this tutorial study is well run, the highest Assay.SD of 0.06 is still considered very low): 

```
# output raw protein group quant summaries
group_quants(fit.sigma, summary = T)

# output normalised protein group quant summaries
group_quants(fit.theta, summary = T)

# plot PCA with "Assay.SD"
plot_robust_pca(fit.theta)
```

Now you can run seaMass-Δ on the results of seaMass-Θ:

```
# run seaMass-delta
fit.delta <- seaMass_delta(fit.theta)
```

Results and plots of this analysis are added to the *csv* directory and *report.zip*. Again, you can get these results in R and generate plots such as Volcano plots:

```
# get protein group FDR results
group_quants_fdr(fit.delta)

# plot volcano
plot_volcano(fit.delta)
```

### Tutorial 2 - running on a HPC cluster

To run seaMass on a HPC cluster, make sure that R and the seaMass R package is installed on the cluster nodes. Then you need to add a *schedule_slurm*, *schedule_pbs* or *schedule_sge* object to *sigma_control*. We also can specify the path of the output directory, in this case to 'hpc':

```
# prepare seaMass-Σ
fit.sigma <- seaMass_sigma(data, data.design, path = "hpc", run = FALSE, 
  control = sigma_control(schedule = schedule_slurm(cpus_per_task = 14, mem = "64000m",
                                                    mail_user = "username@domain.com")))
```

Note above we also set set *run = FALSE* which prepares but does not run seaMass-Σ. This gives you the opportunity to add one or more seaMass-Θ and seaMass-Δ runs so that everything can be run in one go later.

```
# prepare seaMass-Θ
fit.theta <- seaMass_theta(fit.sigma)

# prepare seaMass-Δ
fit.delta <- seaMass_delta(fit.theta)

# run seaMass-Σ and then seaMass-Δ
run(fit.sigma)
```

Note, if you are not doing this from the HPC submit node, instead of *run(sigma)* you can copy the output folder *Tutorial_SLURM.seaMass* to your HPC submit node and execute *Tutorial_SLURM.seaMass/submit.sh* from the command prompt.

### Tutorial 3 - more advanced analysis

In this tutorial we will analyse a fractionated ProteinPilot iTraq study conducted over three multiplexed runs, removing one of these runs. 

```
# Import tutorial iTraq dataset.
file <- system.file("PeptideSummary.txt.bz2", package = "seaMass")
data <- import_ProteinPilot(file)
```

Unfortunately the input file does not contain information for linking fractions to runs, so seaMass allows you to update the imported data with this information.

```
# Get skeleton injection-run table from imported data
data.runs <- runs(data)

# Indicate which injection refers to which run; use 'NA' to indicate injections to ignore
data.runs$Run[1:26] <- NA
data.runs$Run[27:52] <- "1"
data.runs$Run[53:78] <- "2"

# Update the imported data with the run information and remove any ignored injections
runs(data) <- data.runs
```

If you have a large amount of assays, you may group them into blocks so that seaMass will be robust to batch effects or measurement drift etc. This will also help memory consumption, as blocks are processed by seaMass independently. If you are processing iTraq or TMT data, each run will be treated as a seperate block by default, which seaMass-Σ will autodetect.

```
# Get skeleton design matrix
data.design <- new_assay_design(data)

# You can rename assays, or remove them from the analysis with 'NA'
data.design$Assay <- factor(c(
  NA, NA, NA, NA, "1.1", "1.2", "1.3", "1.4",
  NA, NA, NA, NA, "2.1", "2.2", "2.3", "2.4"
))
```

Next, run the seaMass-Σ model. The *seaMass_sigma* *path* parameter allows you to specify the output directory, and any internal control parameters (such as the number of CPU threads to use) can be specified through a *sigma_control* object. 

```
# run seaMass-Σ
fit.sigma <- seaMass_sigma(
  data,
  data.design,
  path = "tutorial.seaMass",
  control = sigma_control(nthread = 4)
)
```

seaMass-Σ computes raw protein group quants (together with peptide deviations from the protein group quants, and peptide/feature stdevs) for each block. To compare across blocks correctly you need to normalise them via 'reference weights' for each assay. You can add columns to *data.design* before or after seaMass-Σ; here we specify the specific reference assays e.g. pooled reference samples, of if the study design has been blocked appropriately, we can just use all relevant assays:

```
# Define reference weights for each block
data.design$RefWeight <- c(
  0, 0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 0, 1, 1, 1, 1
)

# run seaMass-Θ
fit.theta <- seaMass_theta(fit.sigma, data.design)
```

You can also add columns to *data.design* before or after seaMass-Θ to supply seaMass-Δ with samples and conditions for different differential analyses:

```
# Define Sample and Condition mappings
data.design$Sample <- factor(c(
  NA, NA, NA, NA, "A1", "A2", "B1", "B2",
  NA, NA, NA, NA, "A3", "A4", "B3", "B4"
))
data.design$Condition <- factor(c(
  NA, NA, NA, NA, "A", "A", "B", "B",
  NA, NA, NA, NA, "A", "A", "B", "B"
))
```

Using [MCMCglmm](https://cran.r-project.org/package=MCMCglmm) formula syntax you define all kinds of different models (note that the prior specification is the most difficult part to get right!). For example, to perform a standard t-test rather than a Welch's t-test, do: 

```
# run seaMass-Δ
fit.delta <- seaMass_delta(
  fit.theta, data.design,
  rcov = ~ units,
  prior = list(R = list(V = 1, nu = 2e-4))  
)
```

Since we know the ground truth, lets visualise our performance with a Precision-Recall plot.

```
# get protein group FDR results
group_quants_fdr(fit.delta)

# plot precision-recall curve with ground truth
plot_pr(fit.delta, "truth_seaMass_spikein")
```
