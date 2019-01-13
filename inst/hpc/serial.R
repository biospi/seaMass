# load params
params <- readRDS(file.path(id, "input", "params.rds"))

# run model0
setwd(file.path(id, "model0", "results"))
sapply(1:params$model0.nchain, function(i) bayesprot::process.model0(i))
setwd(file.path("..", "..", ".."))

# run output0
setwd(file.path(id, "output0", "results"))
bayesprot::process.output0()
setwd(file.path("..", "..", ".."))

# run model
setwd(file.path(id, "model", "results"))
sapply(1:params$model.nchain, function(i) bayesprot::process.model(i))
setwd(file.path("..", "..", ".."))

# run output
setwd(file.path(id, "output", "results"))
bayesprot::process.output()
setwd(file.path("..", "..", ".."))

if (file.exists(file.path(id, "plots"))) {
  # run plots
  setwd(file.path(id, "plots", "results"))
  sapply(1:params$model.nchain, bayesprot::process.plots)
  setwd(file.path("..", "..", ".."))
}
