# seaMass
The seaMass suite of tools for quantification and differential expression analysis in mass spectrometry proteomics. The current status of the individual tools is:

- **seaMass-α**: Still in development, this will provide sensitive feature extraction directly from raw mass spectrometry data. Please see [[Liao et al. IEEE ISBI, 2014]](https://doi.org/10.1109/ISBI.2014.6868123) and [[Zhang et al., Proteomics, 2015]](https://doi.org/10.1002/pmic.201400428) for technical details. The C++ codebase resides at (here)[https://github.com/biospi/seaMass-alpha] and when ready will have an R interface through this package.

- **seaMass-Σ**: A Bayesian protein group-level quantification technique univerally supporting label-free, SILAC, iTraq/TMT and DIA data. Currently we support proteomics input from Waters Progenesis (label-free), SCIEX ProteinPilot (iTraQ), Thermo ProteomeDiscoverer (SILAC/TMT) and OpenSWATH (SWATH). Other packages supported on request (MaxQuant coming soon). The model provides automatic quality control by downweighting problematic samples and peptides/features, can scale to massive study sizes, and propagates quantification uncertainty downstream to the differential expression analysis stage. 

- **seaMass-Δ**: Bayesian normalisation and differential expression analysis on the output of seaMass-Σ. By harnessing the generic [MCMCglmm](https://cran.r-project.org/web/packages/MCMCglmm) package for Bayesian mixed-effect modelling, the tool allows the user to perform many kinds of univariate analysis on the data, from simple Welch's t-tests and two-way ANOVA to timecourse and multi-level models. Studies analysed with earlier versions of seaMass-Σ and seaMass-Δ are published in [[Freeman et al, Diabetes, 2016]](	
https://doi.org/10.2337/db15-0835), [[Xu et al, Nature Comms Biology, 2019]](https://doi.org/10.1038/s42003-018-0254-9) and [[Kassab et al, Molecular Metabolism, 2019]](	
https://doi.org/10.1016/j.molmet.2019.08.003).

- **mzMLb**: A future-proof binary HDF5 encoding of the [mzML](http://www.psidev.info/mzML) standard interchange format for raw mass spectrometry data flexibly providing best-in-class compression or best-in-class random access speeds. The code is implemented in ProteoWizard [here](https://github.com/biospi/mzmlb) and is moving towards approval by the [Proteomics Standards Intiative](http://www.psidev.info/).  

## Requirements

seaMass is an R package that works on Windows, MacOS and Linux. Small studies (n=4 vs n=4) will take about one to two hours to process with default control parameters. Large studies (e.g. 100 samples) could take overnight or longer. Memory requirements are dependent on the amount of data in the largest block of assays (Each iTraq/TMT run is a separate block, whereas in label-free/SILAC/DIA the blocks are user-preference). Typically 16-32Gb of memory is required. To reduce memory usage at the expense of speed, reduce the number of CPU threads via the 'nthread' parameter of the seaMass 'control' object (see tutorial below).

## Installation

Please install the R package directly from our Github repository. To do this you first need to install the *devtools* package: 

```
install.packages("devtools")
```

Then to install seaMass, or upgrade to the latest version, simply run:

```
devtools::install_github("biospi/seaMass", dependencies = TRUE)
```

To upgrade to the latest version of seaMass, simply run this line again.
