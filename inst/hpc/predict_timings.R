library(brms)
library(ggplot2)
library(data.table)

dd <- fread("protein_timings.csv")
dd[, Intercept := 1]
dd[, PeptideFeature := nPeptide * nFeature]
dd[, Peptide2 := nPeptide * nPeptide]
dd[, Feature2 := nFeature * nFeature]

fit <- brm(
  total ~ Intercept + nPeptide + nFeature + PeptideFeature + Peptide2 + Feature2 - 1,
  data = dd, family = gaussian(), prior = set_prior("uniform(0, 10000)", class = "b", lb = 0.0, ub = 10000),
  chains = 4, cores = 4
)
summary(fit)

dd$secondsPred <- predict(fit)[,"Estimate"]
dd$secondsDiff <- dd$total - dd$secondsPred

g <- ggplot(aes(x=nPeptide, y=secondsDiff), data=dd)
g <- g + geom_point()
g

g <- ggplot(aes(x=nFeature, y=secondsDiff), data=dd)
g <- g + geom_point()
g


# Intercept Peptides  Features
# 2.67      0.87      3.30
