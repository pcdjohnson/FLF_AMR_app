# Power analysis for Taya's FLF app
library(GLMMmisc)
library(ggplot2)
library(hrbrthemes)
library(lme4)
library(lmerTest)

# set up study design
sources <- c("Human", "Livestock", "Aquaculture")
distances <- 0:2
n.rep <- 3
dat <- expand.grid(rep = 1:n.rep, Distance = distances, Source = sources)
dat$site <- paste(dat$Source, dat$rep, sep = "-")

# model parameters
SD <- 0.2
SD.site <- 0.2

# mean relative abundance in each source type
source.effect <- c(Human = 1, Livestock = 1, Aquaculture = 2)

# % decline with each distance unit
dist.effect <- 0.4


# simulate log relative abundance
simdat <-
  sim.glmm(
    design.data = dat, 
    fixed.eff = 
      list(
        intercept = log(4),
        Distance = log(dist.effect), 
        Source = log(source.effect)),
    SD = SD,
    rand.V = c(site = SD.site^2))
simdat$Abundance <- exp(simdat$response)

ggplot(data = simdat, mapping = aes(x = Distance, y = Abundance, group = site, color = Source)) +
  geom_point() +
  geom_line() 


mod1 <- lmer(response ~ Source + Distance + (1 | site), data = simdat, REML = FALSE)
mod0 <- update(mod1, ~ . - Source)
anova(mod0, mod1)[2, "Pr(>Chisq)"]
  