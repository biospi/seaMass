library(brms)
library(ggplot2)
library(data.table)

dd <- fread("protein_timings.csv")
dd$Intercept <- 1

fit <- brm(
  seconds ~ Intercept + nPeptide + nFeature - 1,
  data = dd, family = gaussian(), prior = set_prior("uniform(0, 10000)", class = "b", lb = 0.0, ub = 10000)
)
summary(fit)

dd$secondsPred <- predict(fit)[,"Estimate"]

g <- ggplot(aes(x=nPeptide, y=seconds), data=dd)
g <- g + geom_point()
g <- g + geom_point(aes(y = secondsPred), colour="blue")
g

g <- ggplot(aes(x=nFeature, y=seconds), data=dd)
g <- g + geom_point()
g <- g + geom_point(aes(y = secondsPred), colour="blue")
g

# Intercept Peptides  Features
# 2.67      0.87      3.30
