# BayesProt v1.2.0
Bayesian mixed-effects model and uncertainty propagation for mass spectrometry proteomics, achieving sensitive protein-level quantification and differential expression analysis. Currently works with imported data from SCIEX ProteinPilot, Thermo ProteomeDiscoverer and Waters Progenesis across iTraq/TMT, SILAC, Label-free and SWATH data.

## Installation

```
install.packages("devtools")
library(devtools)
install_github("biospi/bayesprot", dependencies = T)
```

## Usage

Firstly, you need to use an 'import' function to convert from the upstream tool format to BayesProt's standardised data.frame format. Then you use the 'bayesprot' function output a submission zip file for transfer to your HPC cluster.

### Importing from MSstats

BayesProt can import data prepared by MSstats. For example:

```
library(bayesprot)
library(data.table)

dd.input <- fread("E1508100903_6600_32Fixed.tsv", check.names=T)
dd <- importMSstats(dd.input)
bayesprot(dd, id = "SWATHBenchmark")
```

### Importing from ProteomeDiscoverer (TMT)

If you are analysing more than one 6 or 12-plex, please ensure you run all through ProteomeDiscoverer at the same time. BayesProt then takes in the PSM output. For example:

```
todo
```

### Importing from Protein Pilot (iTraq)

If you are analysing more than one 4 or 8-plex, please ensure you run all through Protein Pilot at the same time. Since Protein Pilot doesn't understand datasets merged in this way, you also need to supply a data.frame which associates each fraction with each run. For example

Fraction | Run
---------|-----
1        |  A 
2        |  A 
3        |  A 
4        |  B
5        |  B 
6        |  B  

```
library(readxl)
library(bayesprot)

fractions <- read_excel("fractions.xlsx")
dd <- importProteinPilot("20140910_NR_JXU_GRP 1 SET 1-2-3_combined_PeptideSummary.txt", fractions)
bayesprot(dd, id = "JXU1")
```
