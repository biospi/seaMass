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
nitt <- as.integer(dd.params[Key=="model.nitt", Value])
burnin <- as.integer(dd.params[Key=="model.burnin", Value])
thin <- as.integer(dd.params[Key=="model.thin", Value])

# process arguments
args <- commandArgs(T)
if (length(args) == 0) args <- c("95")
batch <- ((as.integer(args[1]) - 1) %% nbatch) + 1
chain <- ((as.integer(args[1]) - 1) %/% nbatch) + 1

# set seed
set.seed(as.integer(dd.params[Key=="seed", Value]) + chain - 1)

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

# remove single measurements per feature as cannot constuct a ratio
#dd <- merge(dd, dd[, list(nMeasurement = .N), by = list(FeatureID)])
#dd <- dd[nMeasurement > 1,]
#dd$nMeasurement <- NULL

# model needs at least one Assay per Protein to have more than one Feature measurement
dd <- merge(dd, dd[, list(nFeature = length(unique(as.character(FeatureID)))), by = list(ProteinID, AssayID)], by = c("ProteinID", "AssayID"))
dd <- merge(dd, dd[, list(nFeatureMax = max(nFeature)), by = FeatureID], by = "FeatureID")
dd <- dd[nFeatureMax > 1, !c("nFeature", "nFeatureMax")]

# need to drop lost levels
dd$ProteinID <- factor(dd$ProteinID)
dd$PeptideID <- factor(dd$PeptideID)
dd$FeatureID <- factor(dd$FeatureID)
dd$AssayID <- factor(dd$AssayID)

# prepare dd for MCMCglmm
nP <- length(levels(dd$ProteinID))
nT <- length(levels(dd$PeptideID))
nF <- length(levels(dd$FeatureID))
nR <- length(levels(dd$RunID))
nL <- length(levels(dd$LabelID))
nA <- length(levels(dd$AssayID))
dd$Count = round(dd$Count)

# setup QuantID contrasts
if (nL == 1) {
  dd[, BaselineID := dd.assays$AssayID[dd.assays$RunID == unique(RunID)[1]], by = list(ProteinID)]
} else {
  dd[, BaselineID := dd.assays$AssayID[dd.assays$RunID == RunID[1] & dd.assays$LabelID == unique(LabelID)[1]], by = list(ProteinID, RunID)]
}
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
  family = "poisson", data = dd, prior = prior, nitt = nitt, burnin = burnin, thin = thin, pr = T, verbose = F
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
colnames(mcmc.quants) <- sub("^QuantID", "", colnames(mcmc.quants))

# extract peptide level quants, converting to log2
mcmc.peptides.quants <- model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol))] / log(2)
colnames(mcmc.peptides.quants) <- sub("^PeptideID([0-9]+)\\.AssayID(\\.[0-9]+)$", "\\1\\2", colnames(mcmc.peptides.quants))
model$Sol <- NULL

# extract peptide variances, converting to log2 stdevs
mcmc.peptides.sd <- sqrt(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV))] / log(2))
if (nT==1) {
  colnames(mcmc.peptides.sd) <- levels(dd$PeptideID)
} else {
  colnames(mcmc.peptides.sd) <- gsub("^PeptideID([0-9]+)\\.AssayID$", "\\1", colnames(mcmc.peptides.sd))
}

# extract feature variances, converting to log2 stdevs
mcmc.features.sd <- sqrt(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV))] / log(2))
if (nF==1) {
  colnames(mcmc.features.sd) <- levels(dd$FeatureID)
} else {
  colnames(mcmc.features.sd) <- gsub("^FeatureID([0-9]+)\\.units$", "\\1", colnames(mcmc.features.sd))
}
model$VCV <- NULL

# save
save(mcmc.quants, mcmc.peptides.quants, mcmc.peptides.sd, mcmc.features.sd, time.mcmc, file = paste0(batch, ".", chain, ".Rdata"))

message(paste0("[",Sys.time()," Finished]"))

