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

## set to dd
nP <- length(levels(dd$ProteinID))
nT <- length(levels(dd$PeptideID))
nF <- length(levels(dd$FeatureID))  
nA <- length(levels(dd$AssayID))
dd$Count = round(dd$Count)

# MCMCglmm doesn't keep treatment contrast order
if (nP > 1 & nA > 1) dd$AssayID <- factor(dd$AssayID, rev(levels(dd$AssayID)))

dd$ProteinID <- factor(dd$ProteinID)
dd$PeptideID <- factor(dd$PeptideID)
dd$FeatureID <- factor(dd$FeatureID) 
dd$AssayID <- factor(dd$AssayID)

# run model!
prior <- list(
  G = list(G1 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(1000, nT))),
  R = list(V = diag(nF), nu = 0.002)
)
time.mcmc <- system.time(model <- (MCMCglmm(
  as.formula(paste0("Count ~ FeatureID + ", ifelse(nP==1, "", "ProteinID:"), "AssayID - 1")),
  random = as.formula(paste0("~ ", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
  rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
  family = "poisson", data = dd, prior = prior, nitt = nitt, burnin = burnin, thin = thin, verbose = F
)))
summary(model)
message("")

# refer to dd to figure out what mcmc samples we need to save
dd.ids <- dd[, list(ProteinID, AssayID)]
dd.ids <- dd.ids[!duplicated(dd.ids),]
if (nP>1) dd.ids$variable <- paste0("ProteinID", dd.ids$ProteinID, ":")
dd.ids$variable <- paste0(dd.ids$variable, "AssayID", dd.ids$AssayID)

# save mcmc samples for exposures.R, converting to log2
mcmc.quants <- as.data.table(model$Sol[, colnames(model$Sol) %in% dd.ids$variable] / log(2))
mcmc.peptides <- as.data.table(sqrt(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV))] / log(2)))
mcmc.features <- as.data.table(sqrt(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV))] / log(2)))
model$Sol <- NULL
model$VCV <- NULL

mcmc.quants$sampID <- factor(1:nrow(mcmc.quants))
mcmc.quants <- melt(mcmc.quants, id.vars = "sampID", value.name = "log2quant")
mcmc.quants <- merge(dd.ids[, list(variable, ProteinID, AssayID)], mcmc.quants)
mcmc.quants$variable <- NULL
mcmc.quants$AssayID <- factor(mcmc.quants$AssayID)
dir.create("quants", showWarnings = F)
for (i in levels(mcmc.quants$AssayID)) saveRDS(mcmc.quants[AssayID == i, !"AssayID"], file=file.path("quants", paste0(batch, ".", chain, ".", i, ".rds")))

colnames(mcmc.peptides) <- gsub("^PeptideID([0-9]+)\\.LabelID$", "\\1", colnames(mcmc.peptides))
mcmc.peptides$sampID <- factor(1:nrow(mcmc.peptides))
mcmc.peptides <- melt(mcmc.peptides, id.vars = "sampID", variable.name = "PeptideID", value.name = "log2stdev")
setcolorder(mcmc.peptides, c("PeptideID", "sampID", "log2stdev"))

colnames(mcmc.features) <- gsub("^FeatureID([0-9]+)\\.units$", "\\1", colnames(mcmc.features))
mcmc.features$sampID <- factor(1:nrow(mcmc.features))
mcmc.features <- melt(mcmc.features, id.vars = "sampID", variable.name = "FeatureID", value.name = "log2stdev")
setcolorder(mcmc.features, c("FeatureID", "sampID", "log2stdev"))

dir.create("stdevs", showWarnings = F)
save(mcmc.peptides, mcmc.features, time.mcmc, file=file.path("stdevs", paste0(batch, ".", chain, ".Rdata")))

message(paste0("[",Sys.time()," Finished]"))

