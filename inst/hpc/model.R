invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[", Sys.time(), " Starting]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MCMCglmm))
suppressPackageStartupMessages(library(methods))
options(max.print = 99999)

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

# process arguments
args <- commandArgs(T)
if (length(args) == 0) args <- c("2")
chain <- ((as.integer(args[1]) - 1) %% params$model.nchain) + 1
batch <- ((as.integer(args[1]) - 1) %/% params$model.nchain) + 1

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

# setup QuantID contrasts, accounting for complete separability between assays (nb: not efficient!)
# first, figure out baseline for each measurement (might be more than one per feature if separable [e.g. iTraq runs])
baseline.func <- function(dd.protein) {
  # create co-occurence matrix of which assays are present in each feature
  dd.tmp <- copy(dd.protein)
  dd.tmp <- merge(dd.tmp, dd.tmp, by = "FeatureID", allow.cartesian = T)
  mat.tmp <- table(dd.tmp[, list(AssayID.x, AssayID.y)])
  # matrix multiplication distributes assay relationships
  mat.tmp <- mat.tmp %*% mat.tmp
  # ignore columns that are not in ref.assays
  mat.tmp[!dd.assays$isRef,] <- NA
  # baseline is first non-zero occurence for each assay
  apply(mat.tmp != 0, 2, which.max)[dd.protein$AssayID]
}
dd$BaselineID <- dd$AssayID
dd[, BaselineID := baseline.func(.SD), by = ProteinID, .SDcols = c("FeatureID", "AssayID")]
dd$QuantID <- as.character(interaction(dd$ProteinID, dd$BaselineID, dd$AssayID, lex.order = T, drop = T))
# and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
dd[AssayID == BaselineID, QuantID := "."]
dd$QuantID <- factor(dd$QuantID)

# prepare dd for MCMCglmm
nP <- length(levels(dd$ProteinID))
nT <- length(levels(dd$PeptideID))
nF <- length(levels(dd$FeatureID))
nA <- length(levels(dd$AssayID))
nQ <- length(levels(dd$QuantID))
dd$Count <- round(dd$Count)
if (!is.null(dd$MaxCount)) dd$MaxCount <- round(dd$MaxCount)

# set seed
set.seed(params$seed * params$model.nchain + chain - 1)

# run model!
prior <- list(
 G = list(G1 = list(V = diag(nA), nu = nT, alpha.mu = rep(0, nA), alpha.V = diag(1000, nA)),
          G2 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(1000, nT))),
  R = list(V = diag(nF), nu = 0.002)
)
time.mcmc <- system.time(model <- (MCMCglmm(
  as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ FeatureID + QuantID - 1")),
  random = as.formula(paste0("~ idh(AssayID):PeptideID + ", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
  rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
  family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"),
  data = dd, prior = prior, nitt = params$model.nitt, burnin = params$model.burnin, thin = params$model.thin, pr = T, verbose = T
)))
summary(model)
message("")

if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
  stop("Some contrasts were dropped unexpectedly")
}

# extract quants, converting to log2
mcmc.quants <- mcmc(matrix(0.0, nrow(model$Sol), nQ-1), start = start(model$Sol), end = end(model$Sol), thin = thin(model$Sol))
colnames(mcmc.quants) <- paste0("QuantID", levels(dd$QuantID)[2:nQ])
for (i in grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$$", colnames(model$Sol))) mcmc.quants[, colnames(model$Sol)[i]] <- model$Sol[, i] / log(2)
colnames(mcmc.quants) <- sub("^QuantID", "", colnames(mcmc.quants))

# extract peptide deviations from protein quants, converting to log2
mcmc.peptides.quants <- model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol))] / log(2)
colnames(mcmc.peptides.quants) <- sub("^PeptideID([0-9]+)\\.AssayID(\\.[0-9]+)$", "\\1\\2", colnames(mcmc.peptides.quants))

# extract assay deviations from peptide quants, converting to log2
mcmc.assays.quants <- model$Sol[, grep("^AssayID[0-9]+\\.PeptideID\\.[0-9]+$", colnames(model$Sol))] / log(2)
colnames(mcmc.assays.quants) <- sub("^AssayID([0-9]+)\\.PeptideID(\\.[0-9]+)$", "\\1\\2", colnames(mcmc.assays.quants))

model$Sol <- NULL

# extract peptide variances, converting to log2 stdevs
mcmc.peptides.sd <- sqrt(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV))] / log(2))
if (nT==1) {
  colnames(mcmc.peptides.sd) <- levels(dd$PeptideID)
} else {
  colnames(mcmc.peptides.sd) <- gsub("^PeptideID([0-9]+)\\.AssayID$", "\\1", colnames(mcmc.peptides.sd))
}

# extract assay variances, converting to log2 stdevs
mcmc.assays.sd <- sqrt(model$VCV[, grep("^AssayID[0-9]+\\.PeptideID$", colnames(model$VCV))] / log(2))
colnames(mcmc.assays.sd) <- gsub("^AssayID([0-9]+)\\.PeptideID$", "\\1", colnames(mcmc.assays.sd))

# extract feature variances, converting to log2 stdevs
mcmc.features.sd <- sqrt(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV))] / log(2))
if (nF==1) {
  colnames(mcmc.features.sd) <- levels(dd$FeatureID)
} else {
  colnames(mcmc.features.sd) <- gsub("^FeatureID([0-9]+)\\.units$", "\\1", colnames(mcmc.features.sd))
}
model$VCV <- NULL

# save
save(
  mcmc.quants,
  mcmc.peptides.quants,
  mcmc.peptides.sd,
  mcmc.assays.quants,
  mcmc.assays.sd,
  mcmc.features.sd,
  time.mcmc,
  file = paste0(batch, ".", chain, ".Rdata")
)

message(paste0("[", Sys.time(), " Finished]"))

