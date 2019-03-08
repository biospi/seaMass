# load params
params <- readRDS(file.path(output, "input", "params.rds"))
wd <- getwd()

# run model0
setwd(file.path(wd, output, "model0", "results"))
sapply(1:params$model0.nchain, function(chain) {
  commandArgs <- function(...) chain
  source(system.file("hpc/model0.R", package = "bayesprot"), local = T)
})

# run output0
setwd(file.path(wd, output, "output0", "results"))
source(system.file("hpc/output0.R", package = "bayesprot"), local = T)

# run model
setwd(file.path(wd, output, "model", "results"))
sapply(1:params$model.nchain, function(chain) {
  commandArgs <- function(...) chain
  source(system.file("hpc/model.R", package = "bayesprot"), local = T)
})

# run output
setwd(file.path(wd, output, "output", "results"))
source(system.file("hpc/output.R", package = "bayesprot"), local = T)

if (file.exists(file.path(output, "plots"))) {
  # run plots
  setwd(file.path(wd, output, "plots", "results"))
  sapply(1:params$model.nchain, function(chain) {
    commandArgs <- function(...) chain
    source(system.file("hpc/plots.R", package = "bayesprot"), local = T)
  })
}

setwd(wd)
