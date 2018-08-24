# BayesProt
Bayesian linear mixed-effects model for protein-level quantification in proteomics. Currently works with PSM output from ProteinPilot, ProteomeDiscoverer and MSstats.

## Installation

BayesProt requires the use of a HPC cluster. We support the popular cluster managers - SLURM, PBS, SGE (todo), and HTCondor (todo). Before BayesProt can run on the cluster, you need to ensure R and all of BayesProt's R library dependencies are installed there. The later can be done by installing the BayesProt package on your cluster's submit node as follows:

```
install.packages("devtools")
library(devtools)
install_github("biospi/bayesprot", ref = "develop", dependencies = T)
```

On the desktop PC you are using to prepare and analyse datasets, you need to do the same, but you don't need to install all the dependencies: 

```
install.packages("devtools")
library(devtools)
install_github("biospi/bayesprot", ref = "develop")
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

### Submitting to HPC cluster

Transfer the submission zip file to your cluster, unzip and then execute the relevant batch script (e.g. slurm.sh or pbs.sh). WhenFor example:

```
scp submit.bayesprot.zip hpc.myuni.ac.uk:

ssh hpc.myuni.ac.uk
unzip submit.bayesprot.zip
submit.bayesprot/slurm.sh
exit
```

When finished you will be notified by email, and a zip file of your results can be downloaded:

```
scp hpc.myuni.ac.uk:submit.bayesprot.output.zip .
unzip submit.bayesprot.output.zip
```
