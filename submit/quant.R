invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[", Sys.time(), " Starting]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MCMCglmm))
suppressPackageStartupMessages(library(methods))
options(max.print = 99999)

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

nbatch <- as.integer(dd.params[Key=="nbatch", Value])
nitt <- as.integer(dd.params[Key=="quant.nitt", Value])
burnin <- as.integer(dd.params[Key=="quant.burnin", Value])
thin <- as.integer(dd.params[Key=="quant.thin", Value])

# process arguments
args <- commandArgs(T)
if (length(args) == 0) args <- c("1")
batch <- ((as.integer(args[1]) - 1) %% nbatch) + 1
chain <- ((as.integer(args[1]) - 1) %/% nbatch) + 1

# set seed
set.seed(as.integer(dd.params[Key=="seed", Value]) + chain - 1)

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

## prepare dd for MCMCglmm
nP <- length(levels(dd$ProteinID))
nT <- length(levels(dd$PeptideID))
nF <- length(levels(dd$FeatureID))  
nA <- length(levels(dd$AssayID))
dd$Count = round(dd$Count)

# setup QuantID contrasts
dd[, BaselineID := dd.assays$AssayID[dd.assays$RunID == RunID[1] & dd.assays$LabelID == unique(LabelID)[1]], by = c("ProteinID", "RunID")] # todo: work for no label case
dd$QuantID <- as.character(interaction(dd$ProteinID, dd$BaselineID, dd$AssayID, lex.order = T, drop = T))
dd[AssayID == BaselineID, QuantID := "."]
dd$QuantID <- factor(dd$QuantID)
nQ <- length(levels(dd$QuantID))

# run model!
prior <- list(
  G = list(G1 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(1000, nT))),
  R = list(V = diag(nF), nu = 0.002)
)
time.mcmc <- system.time(model <- (MCMCglmm(
  Count ~ FeatureID + QuantID - 1,
  random = as.formula(paste0("~ ", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
  rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
  family = "poisson", data = dd, prior = prior, nitt = nitt, burnin = burnin, thin = thin, verbose = F
)))
summary(model)
message("")

if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != length(levels(dd$QuantID)) - 1) {
  stop("Some contrasts were dropped unexpectedly")
}

# extract quants, converting to log2
mcmc.quants <- mcmc(matrix(0.0, nrow(model$Sol), nQ-1), start = start(model$Sol), end = end(model$Sol), thin = thin(model$Sol))
colnames(mcmc.quants) <- paste0("QuantID", levels(dd$QuantID)[2:nQ])
for (i in grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$$", colnames(model$Sol))) mcmc.quants[, colnames(model$Sol)[i]] <- model$Sol[, i] / log(2)
model$Sol <- NULL

# extract variances, converting to log2 stdevs
mcmc.peptides <- sqrt(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV))] / log(2))
if (nT==1) {
  colnames(mcmc.peptides) <- levels(dd$PeptideID)
} else {
  colnames(mcmc.peptides) <- gsub("^PeptideID([0-9]+)\\.AssayID$", "\\1", colnames(mcmc.peptides))
}
mcmc.features <- sqrt(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV))] / log(2))
if (nF==1) {
  colnames(mcmc.features) <- levels(dd$FeatureID)
} else {
  colnames(mcmc.features) <- gsub("^FeatureID([0-9]+)\\.units$", "\\1", colnames(mcmc.features))
}
model$VCV <- NULL

# save to files
dir.create("quants", showWarnings = F)
for (i in 1:nA) {
  mcmc.quants.assay <- mcmc.quants[, grep(paste0("^QuantID[0-9]+\\.[0-9]+\\.", i, "$"), colnames(mcmc.quants))]
  if (length(mcmc.quants.assay) > 0) {
    colnames(mcmc.quants.assay) <- sub("^QuantID([0-9]+\\.[0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.quants.assay))
    saveRDS(mcmc.quants.assay, file=file.path("quants", paste0(batch, ".", chain, ".", i, ".rds")))
  }
} 
dir.create("stdevs", showWarnings = F)
save(mcmc.peptides, mcmc.features, time.mcmc, file=file.path("stdevs", paste0(batch, ".", chain, ".Rdata")))

message(paste0("[",Sys.time()," Finished]"))

