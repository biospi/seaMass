import ClusterManagers

# import
symlink("../../input/data.Rdata", "import/results/data.Rdata")
symlink("../../input/parameters.Rdata", "import/results/parameters.Rdata")
symlink("../../input/design.Rdata", "import/results/design.Rdata")
ClusterManagers.addprocs_sge(1)
@fetch cd(() -> run(`Rscript ../../import.R HPC`), "import/results")

# norm

# exposures

# model

# plots

# output
