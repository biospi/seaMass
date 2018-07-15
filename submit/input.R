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
  
  c("bayesprot.id", id),
  c("bayesprot.version", "1.2")
)
colnames(dd.params) <- c("Key", "Value")
dd.params <- as.data.table(dd.params)

# ensure factors for building indicies
dd$ProteinGroup <- factor(dd$ProteinGroup)
dd$Peptidoform <- factor(dd$Peptidoform)
dd$Feature <- factor(dd$Feature)
dd$Run <- factor(dd$Run)
dd$Label <- factor(dd$Label)

# build ProteinGroup index
dd.proteingroups <- dd[, .(
  nPeptidoform = length(unique(as.character(Peptidoform))),
  nFeature = length(unique(as.character(Feature))),
  nCount = .N
), by = ProteinGroup]
dd.proteingroups <- dd.proteingroups[order(-nPeptidoform, -nFeature, -nCount, ProteinGroup)]
dd.proteingroups$ProteinGroupID <- factor(as.integer(factor(dd.proteingroups$ProteinGroup, unique(dd.proteingroups$ProteinGroup))))

dd <- merge(dd, dd.proteingroups[, list(ProteinGroup, ProteinGroupID)], by = "ProteinGroup")
dd$ProteinGroup <- NULL

# build Peptidoform index
dd.peptidoforms <- dd[, .(
  nFeature = length(unique(as.character(Feature))),
  nCount = .N
), by = Peptidoform]
dd.peptidoforms  <- dd.peptidoforms[order(-nFeature, -nCount, Peptidoform)]
dd.peptidoforms$PeptidoformID <- factor(as.integer(factor(dd.peptidoforms$Peptidoform, unique(dd.peptidoforms$Peptidoform))))

dd <- merge(dd, dd.peptidoforms[, list(Peptidoform, PeptidoformID)], by = "Peptidoform")
dd$Peptidoform <- NULL

# build Feature index
dd.features <- dd[, .(
  nCount = .N
), by = Feature]
dd.features  <- dd.features[order(-nCount, Feature)]
dd.features$FeatureID <- factor(as.integer(factor(dd.features$Feature, unique(dd.features$Feature))))

dd <- merge(dd, dd.features[, list(Feature, FeatureID)])
dd$Feature <- NULL

# build Preparation index
dd.preparations <- dd[, .(
  Preparation = paste(Run, ":", Label)
), by = list(Run, Label)]
dd.preparations$Preparation <- factor(dd.preparations$Preparation)
dd.preparations$RunID <- factor(as.integer(dd.preparations$Run))
dd.preparations$LabelID <- factor(as.integer(dd.preparations$Label))
dd.preparations$PreparationID <- factor(as.integer(dd.preparations$Preparation))

dd <- merge(dd, dd.preparations[, list(Run, Label, RunID, LabelID, PreparationID)], by = c("Run", "Label"))
dd$Run <- NULL
dd$Label <- NULL

# factors are appaulingly slow to split, so change to strings as we want to drop levels anyway
dd$ProteinGroupID <- as.integer(dd$ProteinGroupID)
dd$PeptidoformID <- as.integer(dd$PeptidoformID)
dd$FeatureID <- as.integer(dd$FeatureID)
dd$RunID <- as.integer(dd$RunID)
dd$LabelID <- as.integer(dd$LabelID)
dd$PreparationID <- as.integer(dd$PreparationID)

# estimate how long each ProteinGroup will take to process
# Intercept, ProteinGroups, Peptidoforms, ProteinGroups^2, Peptidoforms^2, ProteinGroups*Peptidoforms 
a <- c(746.09, 32.76, 7.36, 5017.36, 4975.00, 5038.31)

# assign ProteinGroups to batches
message("batching protein groups...")
nbatch <- as.integer(dd.params[Key=="nbatch", Value])
dd.batches <- data.table(nP = rep(0, nbatch), nT = rep(0, nbatch))
dds = vector("list", nbatch)
for (j in 1:nbatch) dds[[j]] <- vector("list", nrow(dd.proteingroups))
for (i in 1:nrow(dd.proteingroups)) {
  scores <- a[1] +
    a[2]*(dd.batches$nP+1) +
    a[3]*(dd.batches$nT+dd.proteingroups$nPeptidoform[i]) +
    a[4]*(dd.batches$nP+1)*(dd.batches$nP+1) +
    a[5]*(dd.batches$nT+dd.proteingroups$nPeptidoform[i])*(dd.batches$nT+dd.proteingroups$nPeptidoform[i]) +
    a[6]*(dd.batches$nP+1)*(dd.batches$nT+dd.proteingroups$nPeptidoform[i])
    
  j <- which.min(scores)
  dds[[j]][[i]] <- dd[ProteinGroupID == dd.proteingroups$ProteinGroupID[i],]
  dd.batches$nP[j] <- dd.batches$nP[j] + 1
  dd.batches$nT[j] <- dd.batches$nT[j] + dd.proteingroups$nPeptidoform[i]
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
  
  # and back to factors...
  dd$ProteinGroupID <- factor(dd$ProteinGroupID)
  dd$PeptidoformID <- factor(dd$PeptidoformID)
  dd$FeatureID <- factor(dd$FeatureID)
  dd$RunID <- factor(dd$RunID)
  dd$LabelID <- factor(dd$LabelID)
  dd$PreparationID <- factor(dd$PreparationID)
  
  save(dd, file = file.path(out_dir, "input", paste0(j-1, ".Rdata")))
}

# save metadata
save(dd.params, dd.proteingroups, dd.peptidoforms, dd.features, dd.preparations, file = file.path(out_dir, "input", "metadata.Rdata"))

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
