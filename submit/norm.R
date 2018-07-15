invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[", Sys.time(), " Starting]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MCMCglmm))
suppressPackageStartupMessages(library(methods))
options(max.print = 99999)

# process arguments
args <- commandArgs(T)
if (length(args) == 0) args <- c("1", "0")
batch <- as.integer(args[1])
chain <- as.integer(args[2])

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

# set seed
set.seed(as.integer(dd.params[Key=="seed", Value]) + chain)

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

## create dd.model
nP <- length(levels(dd$ProteinGroup))
nT <- length(levels(dd$Peptidoform))
nF <- length(levels(dd$Feature))  
nL <- length(levels(dd$Label))
nR <- length(levels(dd$Run))

dd.model <- data.table(
  ProteinGroup = factor(as.integer(dd$ProteinGroup)),
  Peptidoform = factor(as.integer(dd$Peptidoform)),
  Feature = factor(as.integer(dd$Feature)),
  Run = factor(as.integer(dd$Run)),
  Label = factor(as.integer(dd$Label)),
  Count = round(dd$Count)
)

# MCMCglmm doesn't keep treatment contrast order
if (nP > 1 & nR == 1) levels(dd.model$Label) = rev(levels(dd.model$Label))
if (nP > 1 & nL == 1) levels(dd.model$Run) = rev(levels(dd.model$Run))
if (nP == 1 & nR > 1 & nL > 1) levels(dd.model$Label) = rev(levels(dd.model$Run))

# run model!
nitt <- as.integer(dd.params[Key=="quant.nitt", Value])
burnin <- as.integer(dd.params[Key=="quant.burnin", Value])
thin <- as.integer(dd.params[Key=="quant.thin", Value])

prior <- list(
  G = list(G1 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(1000, nT))),
  R = list(V = diag(nF), nu = 0.002)
)

time.mcmc <- system.time(model <- (MCMCglmm(
  as.formula(paste0("Count ~ Feature + ", ifelse(nP==1, "", "ProteinGroup:"), ifelse(nR<=1, ifelse(nL<=1, "", "Label"), ifelse(nL<=1, "Run", "Run:Label")), " - 1")),
  random = as.formula(paste0("~ ", ifelse(nT==1, "Peptidoform", "idh(Peptidoform)"), ifelse(nR<=1, "", ":Run"), ifelse(nL<=1, "", ":Label"))),
  rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(Feature):units"))),
  family = "poisson", data = dd.model, prior = prior, nitt = nitt, burnin = burnin, thin = thin, verbose = F
)))

summary(model)
message("")

# figure out what mcmc samples we should save
dd.id <- data.table(
  ID = CJ(paste("ProteinGroup", levels(dd$ProteinGroup)), paste("Run", levels(dd$Run)), paste("Label", levels(dd$Label)))[, paste(V1, V2, V3, sep =" : ")],
  MCMCglmmID = CJ(paste0("ProteinGroup", 1:nP), paste0("Run", 1:nR), paste0("Label", 1:nL))[, paste(V1, V2, V3, sep =":")]
)
if (nP == 1) dd.id$MCMCglmmID <- sub("^ProteinGroup1:", "", dd.id$MCMCglmmID)
if (nR == 1) dd.id$MCMCglmmID <- sub("Run1:", "", dd.id$MCMCglmmID)
if (nR == 1) dd.id$MCMCglmmID <- sub("Run1", "", dd.id$MCMCglmmID)
if (nL == 1) dd.id$MCMCglmmID <- sub("Label1", "", dd.id$MCMCglmmID)
dd.id$MCMCglmmID <- sub("::", ":", dd.id$MCMCglmmID)

# save mcmc samples for exposures.R
model$Sol <- model$Sol[, colnames(model$Sol) %in% dd.id$MCMCglmmID]
colnames(model$Sol) <- dd.id$ID[match(colnames(model$Sol), dd.id$MCMCglmmID)]
mcmc.quants <- model$Sol

model$VCV <- model$VCV[, grep("^Peptidoform", colnames(model$VCV))]  
colnames(model$VCV) <- levels(dd$Peptidoform)
mcmc.peptidoforms <- model$VCV

save(mcmc.quants, mcmc.peptidoforms, time.mcmc, file=paste0(batch, ".", chain, ".Rdata"))

message(paste0("[",Sys.time()," Finished]"))

