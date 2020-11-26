# Power analysis for Taya's FLF app

# load packages
library(GLMMmisc)
library(ggplot2)
library(hrbrthemes)
library(lme4)
library(lmerTest)

# model parameters and design choices

# SD at the observation level
SD <- 0.2
# SD between sites (random intercepts)
SD.site <- 0.2
# number of replicates (no of sites per source type)
n.rep <- 3
# mean relative abundance in each source type
source.effect <- c(Human = 2, Livestock = sqrt(2), Aquaculture = 1)
# relative decline with each additional distance unit
dist.effect <- 0.4
# mean log abundance at distance = 0 and source = human or livestock
intercept <- log(4)

# set up study design
sources <- c("Human", "Livestock", "Aquaculture")
distances <- 0:2
dat <- expand.grid(rep = 1:n.rep, Distance = distances, Source = sources)
dat$site <- paste(dat$Source, dat$rep, sep = "-")
dat <- dat[order(dat$site), ]

# simulate log relative abundance
simdat <-
  sim.glmm(
    design.data = dat, 
    fixed.eff = 
      list(
        intercept = intercept,
        Distance = log(dist.effect), 
        Source = log(source.effect)),
    SD = SD,
    rand.V = c(site = SD.site^2))
simdat$Abundance <- exp(simdat$response)

# plot the simulated data
ggplot(data = simdat, mapping = aes(x = Distance, y = Abundance, group = site, color = Source)) +
  geom_point() +
  geom_line() 

# calculate a p-value for the null hypothesis that all sources have the same mean log abundance
mod1 <- lmer(response ~ Source + Distance + (1 | site), data = simdat, REML = FALSE)
mod0 <- update(mod1, ~ . - Source)
anova(mod0, mod1)[2, "Pr(>Chisq)"]
