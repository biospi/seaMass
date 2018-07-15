invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[", Sys.time(), " Starting]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MCMCglmm))
suppressPackageStartupMessages(library(methods))
options(max.print = 99999)

# process arguments
args <- commandArgs(T)
if (length(args) == 0) args <- c("6", "0")
batch <- as.integer(args[1])
chain <- as.integer(args[2])

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

nitt <- as.integer(dd.params[Key=="quant.nitt", Value])
burnin <- as.integer(dd.params[Key=="quant.burnin", Value])
thin <- as.integer(dd.params[Key=="quant.thin", Value])

# set seed
set.seed(as.integer(dd.params[Key=="seed", Value]) + chain)

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

## set to dd
nP <- length(levels(dd$ProteinGroupID))
nT <- length(levels(dd$PeptidoformID))
nF <- length(levels(dd$FeatureID))  
nR <- length(levels(dd$RunID))
nL <- length(levels(dd$LabelID))
dd$Count = round(dd$Count)

# MCMCglmm doesn't keep treatment contrast order
if (nP > 1 & nR == 1) levels(dd$LabelID) = rev(levels(dd$LabelID))
if (nP > 1 & nL == 1) levels(dd$RunID) = rev(levels(dd$RunID))
if (nP == 1 & nR > 1 & nL > 1) levels(dd$LabelID) = rev(levels(dd$LabelID))

# run model!
prior <- list(
  G = list(G1 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(1000, nT))),
  R = list(V = diag(nF), nu = 0.002)
)
time.mcmc <- system.time(model <- (MCMCglmm(
  as.formula(paste0("Count ~ FeatureID + ", ifelse(nP==1, "", "ProteinGroupID:"), ifelse(nR<=1, ifelse(nL<=1, "", "LabelID"), ifelse(nL<=1, "RunID", "RunID:LabelID")), " - 1")),
  random = as.formula(paste0("~ ", ifelse(nT==1, "PeptidoformID", "idh(PeptidoformID)"), ifelse(nR<=1, "", ":RunID"), ifelse(nL<=1, "", ":LabelID"))),
  rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
  family = "poisson", data = dd, prior = prior, nitt = nitt, burnin = burnin, thin = thin, verbose = F
)))
summary(model)
message("")

# refer to dd to figure out what mcmc samples we need to save
dd.ids <- dd[!duplicated(dd[, list(ProteinGroupID, PreparationID)]), list(ProteinGroupID, RunID, LabelID)]
dd.ids$ProteinGroupMCMCglmm <- ""
dd.ids$RunMCMCglmm <- ""
dd.ids$LabelMCMCglmm <- ""
if (nP>1) dd.ids$ProteinGroupMCMCglmm <-paste0("ProteinGroupID", dd.ids$ProteinGroupID, ":")
if (nR>1) dd.ids$RunMCMCglmm <-paste0("RunID", dd.ids$RunID, ":")
if (nL>1) dd.ids$LabelMCMCglmm <-paste0("LabelID", dd.ids$LabelID, ":")
dd.ids$MCMCglmm <- paste0(dd.ids$ProteinGroupMCMCglmm, dd.ids$RunMCMCglmm, dd.ids$LabelMCMCglmm)
dd.ids$MCMCglmm <- substr(dd.ids$MCMCglmm, 1, nchar(dd.ids$MCMCglmm) - 1)

# save mcmc samples for exposures.R
denominators <- dd.ids$MCMCglmm[!dd.ids$MCMCglmm %in%  colnames(model$Sol)] 
mcmc.quants <- as.matrix(model$Sol[, colnames(model$Sol) %in% dd.ids$MCMCglmm])
model$Sol <- NULL

mcmc.quants.demoninator <- matrix(0, nrow(mcmc.quants), length(denominators))
colnames(mcmc.quants.demoninator) <- denominators

mcmc.quants <- cbind(mcmc.quants.demoninator, mcmc.quants)
j <- match(colnames(mcmc.quants), dd.ids$MCMCglmm)
colnames(mcmc.quants) <- paste0(dd.ids$ProteinGroupID[j], ":", dd.ids$RunID[j], ":", dd.ids$LabelID[j])

mcmc.peptidoforms <- as.matrix(model$VCV[, grep("^Peptidoform", colnames(model$VCV))])
model$VCV <- NULL
colnames(mcmc.peptidoforms) <- gsub("^PeptidoformID([0-9]+)\\.LabelID$", "\\1", colnames(mcmc.peptidoforms))

save(mcmc.quants, mcmc.peptidoforms, time.mcmc, file=paste0(batch, ".", chain, ".Rdata"))

message(paste0("[",Sys.time()," Finished]"))

