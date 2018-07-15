# default parameters
id <- sub("[.][^.]*$", "", basename(args[1]), perl=T)

dd.params <- rbind(
  c("seed", "0"),
  c("nbatch", "100"),
  
  c("quant.nitt", "1300"),
  c("quant.burnin", "300"),
  c("quant.thin", "1"),
  c("quant.nchain", "10"),

  c("de.nitt", "13000"),
  c("de.burnin", "3000"),
  c("de.thin", "100"),
  c("de.nchain", "100"),
  
  c("hpc", "SLURM"),
  c("hpc.ncpu", "14"),
  c("hpc.nnode", "1"),
  c("hpc.mem", "3G"),
  c("hpc.himem", "16G"),
  c("hpc.longq", "cpu"),
  c("hpc.shortq", "serial"),
  c("hpc.totaljobs", "1"),
  c("hpc.lowncpu", "6"),
  
  c("internal.nrun", length(unique(dd$Run))),
  c("internal.nproteingroup", length(unique(dd$ProteinGroup))),
  c("internal.npeptidoform", length(unique(dd$Peptidoform))),
  c("internal.nfeature", length(unique(dd$Feature))),
  c("internal.nlabel", length(unique(dd$Label))),
  c("internal.id", id)
)
colnames(dd.params) <- c("Key", "Value")
dd.params <- as.data.table(dd.params)

# build indices
dd.proteingroups <- dd[, .(
  nPeptidoform = length(unique(Peptidoform)),
  nFeature = length(unique(Feature)),
  nCount = .N
), by = ProteinGroup]
dd.proteingroups <- dd.proteingroups[order(-nPeptidoform, -nFeature, -nCount, ProteinGroup)]
dd.proteingroups$ProteinGroup <- factor(dd.proteingroups$ProteinGroup, unique(dd.proteingroups$ProteinGroup))
dd.proteingroups$ProteinGroupID <- as.integer(dd.proteingroups$ProteinGroup)

dd$ProteinGroupID <- dd.proteingroups$ProteinGroupID[match(dd$ProteinGroup, dd.proteingroups$ProteinGroup)]
dd.proteingroups$ProteinGroupID <- factor(dd.proteingroups$ProteinGroupID)

dd.peptidoforms <- dd[, .(
  ProteinGroupID,
  nFeature = length(unique(Feature)),
  nCount = .N
), by = Peptidoform]
dd.peptidoforms  <- dd.peptidoforms[order(-nFeature, -nCount, ProteinGroupID, Peptidoform)]
dd.peptidoforms$Peptidoform <- factor(dd.peptidoforms$Peptidoform, unique(dd.peptidoforms$Peptidoform))
dd.peptidoforms$PeptidoformID <- as.integer(dd.peptidoforms$Peptidoform)

dd$PeptidoformID <- dd.peptidoforms$PeptidoformID[match(dd$Peptidoform, dd.peptidoforms$Peptidoform)]
dd.peptidoforms$ProteinGroupID <- factor(dd.peptidoforms$ProteinGroupID)
dd.peptidoforms$PeptidoformID <- factor(dd.peptidoforms$PeptidoformID)

dd.features <- dd[, .(
  PeptidoformID,
  nCount = .N
), by = Feature]
dd.features  <- dd.features[order(-nCount, PeptidoformID, Feature)]
dd.features$Feature <- factor(dd.features$Feature, unique(dd.features$Feature))
dd.features$PeptidoformID <- factor(dd.features$PeptidoformID)
dd.features$FeatureID <- as.integer(dd.features$Feature)

dd$FeatureID <- dd.features$FeatureID[match(dd$Feature, dd.features$Feature)]
dd.features$FeatureID <- factor(dd.features$FeatureID)

dd.preparations <- dd[, .(
  Preparation = paste(Run, ":", Label)
), by = list(Run, Label)]
dd.preparations$RunID <- as.integer(factor(dd.preparations$Run))
dd.preparations$LabelID <- as.integer(factor(dd.preparations$Label))
dd.preparations$PreparationID <- as.integer(factor(dd.preparations$Preparation))

dd$RunID <- dd.preparations$RunID[match(dd$Run, dd.preparations$Run)]
dd$LabelID <- dd.preparations$LabelID[match(dd$Label, dd.preparations$Label)]
dd$PreparationID <- dd.preparations$PreparationID[match(dd$Preparation, dd.preparations$PreparationID)]


# estimate how long each ProteinGroup will take to process
# intercept, peptidoforms, peptidoforms^2, features, features^2, peptidoforms*features 
a <- c(1.462154e+01, 6.002002e-01, 3.886181e+00, 5.376223e-02, 2.671169e-04, 7.516105e+001)
score <- a[1] + a[2]*dd.index$peptidoforms + a[3]*dd.index$peptidoforms^2 +
  a[4]*dd.index$features +  + a[5]*dd.index$features^2 + a[6]*dd.index$peptidoforms*dd.index$features 

# assign ProteinGroups to batches
message("batching protein groups...")
nbatch <- as.integer(dd.params[Key=="nbatch", Value])
scores <- rep(0, nbatch)
dds = vector("list", nbatch)
for (j in 1:nbatch) dds[[j]] <- vector("list", nrow(dd.index))
for (i in 1:nrow(dd.index)) {
  j <- which.min(scores)
  dds[[j]][[i]] <- dd[ForeignKey == dd.index$ForeignKey[i],]
  scores[j] <- scores[j] + score[i]
}

# build HPC submission zip
id <- sub("[.][^.]*$", "", basename(args[1]), perl=T)
out_zip <- paste0(id, ".bayesprot.zip")
tmp_dir <- tempfile("bayesprot.")
out_dir <- file.path(tmp_dir, paste0(id, ".bayesprot"))
dir.create(out_dir, recursive = T)
for (file in list.files(file.path(script_dir, "submit")))
  file.copy(file.path(script_dir, "submit", file), out_dir, recursive = T)

# save batches
for (j in 1:nbatch) {
  dd <- rbindlist(dds[[j]])
  
  dd$Run <- factor(dd$Run)
  dd$Label <- factor(dd$Label)
  dd$ProteinGroup <- factor(dd$ProteinGroup)
  dd$Peptidoform <- factor(dd$Peptidoform)
  dd$Feature <- factor(dd$Feature)
  
  save(dd, file = file.path(out_dir, "input", paste0(j-1, ".Rdata")))
}

# save metadata
save(dd.index, dd.params, file = file.path(out_dir, "input", "metadata.Rdata"))

# this is where the SLURM/PBS scripts should be generated 
source(file.path(script_dir, "submit", "ScheduleHPC.R"))

# New xlsx reader library so Key = X__1 and Value = X__2
systemHPC <- as.character(dd.params[Key=="hpc", Value])

if (systemHPC == "SLURM") {
  norm_chains <- as.integer(dd.params[Key=="quant.nchain", Value])
  model_chains <- as.integer(dd.params[Key=="de.nchain", Value])

  cpu_num <- as.integer(dd.params[Key=="hpc.ncpu", Value])
  node <- as.integer(dd.params[Key=="hpc.nnode", Value])
  mem <- dd.params[Key=="hpc.mem", Value]
  himem <- dd.params[Key=="hpc.himem", Value]
  long_que <- dd.params[Key=="hpc.longq", Value]
  short_que <- dd.params[Key=="hpc.shortq", Value]
  total_jobs <- as.integer(dd.params[Key=="hpc.totaljobs", Value])
  low_cpu_num <- as.integer(dd.params[Key=="hpc.lowncpu", Value])

  clusterHPC <- new(systemHPC, batch = nbatch, normChain = norm_chains, modelChain = model_chains, path = out_dir,
                    cpuNum = cpu_num, node = node, mem = mem, himem = himem, longQue = long_que,shortQue = short_que,
                    totalJobs = total_jobs,lowCPUNum = low_cpu_num)
} else if (systemHPC == "PBS") {
  stop("PBS HPC system yet to be implemented!...")
} else if (systemHPC == "SGE") {
  stop("SGE HPC system yet to be implemented!...")
} else {
  stop("Error: Unknown HPC system. Possible HPC systems = SLURM, PBS, SGE.")
}

#message(paste("batch",batch,"norm_chain",norm_chains,"model_chains",model_chains,"cpu_num",cpu_num,"node",node,"mem",
#              mem,"himem",himem,"long_que",long_que,"short_que",short_que,"total_jobs",total_jobs,"low_cpu_num",low_cpu_num))

# Norm: 'Rscript norm.R <batch> <norm_chain> <norm_chains>' where <batch> is from 0 to batches-1 and <norm_chain> is from 0 to norm_chaines-1
normHPC(clusterHPC)
# Exposures: 'Rscript exposures.R'
exposuresHPC(clusterHPC)
# Model: 'Rscript model.R <batch> <model_chain> <model_chains>' where <batch> is from 0 to batches-1 and <model_chain> is from 0 to model_chains-1
modelHPC(clusterHPC)
# Plots: 'Rscript plots.R <batch>' where <batch> is from 0 to batches-1
plotsHPC(clusterHPC)
# Output: 'Rscript output.R'
outputHPC(clusterHPC)
# Generate bash script for chained HPC submit job:
genJobFileHPC(clusterHPC)

# create zip file and clean up
message(paste0("writing: ", out_zip, "..."))
wd <- getwd()
setwd(tmp_dir)
zip(file.path(wd, out_zip), ".", flags="-r9Xq")
setwd(wd)
unlink(tmp_dir, recursive = T)
message("")  
