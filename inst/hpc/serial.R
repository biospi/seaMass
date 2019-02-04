# load params
params <- readRDS(file.path(id, "input", "params.rds"))
wd <- getwd()

# run model0
setwd(file.path(wd, id, "model0", "results"))
sapply(1:params$model0.nchain, function(i) bayesprot::process.model0(i))

# run output0
setwd(file.path(wd, id, "output0", "results"))
bayesprot::process.output0()

# run model
setwd(file.path(wd, id, "model", "results"))
sapply(1:params$model.nchain, function(i) bayesprot::process.model(i))

# run output
setwd(file.path(wd, id, "output", "results"))
bayesprot::process.output()

if (file.exists(file.path(id, "plots"))) {
  # run plots
  setwd(file.path(wd, id, "plots", "results"))
  sapply(1:params$model.nchain, bayesprot::process.plots)
}

setwd(wd)
