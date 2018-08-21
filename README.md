# BayesProt
Bayesian linear mixed-effects model for protein-level quantification and
study-level statistical testing in proteomics.

Contact: Al Phillips - A.M.Phillips@liverpool.ac.uk

## Usage

For Protein Pilot iTRAQ data:

```
Rscript bayesprot.R my_iTRAQ_experiment.xlsx PeptideSummary.txt
```

For Proteome Discoverer TMT data:

```
Rscript bayesprotTMT.R my_TMT_experiment.xlsx PSMS.txt
```




The design workbooks are formatted with three spreadsheets as follows:

### Parameters
Key-Value pairs:

Key           | Value
--------------|--------
email         | Kendal.Mintcake@footlights.ac.uk
model_fc      | 1.05
norm_nitt     | 13000
norm_nburnin  | 3000
norm_nchain   | 10
norm_nsamp    | 1000
model_nitt    | 13000
model_nburnin | 3000
model_nchain  | 100
model_nsamp   | 10000
nworker       | 1000

### Design

Experimental design.

Run | Channel | Sample | Volume | Digestion | Population | Condition
----|---------|--------|--------|-----------|------------|-----------
A   | 113     | P      | 1      | P1        | Pool       | 2.pool
A   | 114     | P      | 1      | P2        | Pool       | 2.pool
A   | 115     | S1     | 1      | S1        | All        | 0.A
A   | 116     | S2     | 1      | S2        | All        | 0.A
A   | 117     | S3     | 1      | S3        | All        | 0.A
A   | 118     | S4     | 1      | S4        | All        | 1.B
A   | 119     | S5     | 1      | S5        | All        | 1.B
A   | 121     | S6     | 1      | S6        | All        | 1.B
B   | 113     | P      | 1      | P1        | Pool       | 2.pool
B   | 114     | P      | 1      | P2        | Pool       | 2.pool
B   | 115     | S7     | 1      | S7        | All        | 0.A
B   | 116     | S8     | 1      | S8        | All        | 0.A
B   | 117     | S9     | 1      | S9        | All        | 0.A
B   | 118     | S10    | 1      | S10       | All        | 1.B
B   | 119     | S11    | 1      | S11       | All        | 1.B
B   | 121     | S12    | 1      | S12       | All        | 1.B


Populations with a single sample (such as "Pool" above) will be assigned a
biological variance of zero. This is used for e.g. a pooled reference sample.
Condition 0 is used as the reference for differential expression testing.
Conditions with names starting with a lowercase will not be output in plot results.

### Fractions
The mapping between fractions and iTRAQ/TMT Runs is described here. For Protein
Pilot data, the fraction number should be used. For Proteome Discoverer TMT
data, values from the "Spectrum File" column of *PSMs.txt should be used.


#### Protein Pilot
Fraction | Run
---------|-----
1        |  A 
2        |  A 
3        |  A 
4        |  B
5        |  B 
6        |  B  

#### Proteome Discoverer
Fraction | Run
---------|-----
file1.raw|  A 
file2.raw|  A 
file3.raw|  A 
file3.raw|  B
file5.raw|  B 
file6.raw|  B  
