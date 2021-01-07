# Power analysis for Taya's FLF app

# load packages
library(GLMMmisc)
library(ggplot2)
library(hrbrthemes)
library(lme4)
library(lmerTest)

# model parameters and design choices


###### triple design across three countries - assume country effect is low enough to be
###### subsumed into between site noise


# SD at the observation level
SD <- 0.2 #### make 1.1x (Rowe et al 2016 found only about 10% changes over time)
# SD between sites (random intercepts)
SD.site <- 0.2 #### this should equate to ~2x diff between sites of same source based on Chu
# number of replicates (no of sites per source type)
n.rep <- 9
# mean relative abundance in each source type
source.effect <- c(Human = 5, Livestock = 3, Aquaculture = 1)
# relative decline with each additional distance unit
dist.effect <- 0.4 ### 60% decline per distance unit (based on Chu)
# mean log abundance at distance = 0 and source = human or livestock
intercept <- log(4)   ##### probably doesn't make a difference - check  

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
#mod1 <- lmer(response ~ Source + Distance + (1 | site), data = simdat, REML = FALSE)
#mod0 <- update(mod1, ~ . - Source)
anova(mod0, mod1)[2, "Pr(>Chisq)"]


simdat0 <- droplevels(simdat[simdat$Distance == 0, ])


###### distance/spatial effect #####
#' for a subset of sites sample at distances 1 & 2
#' for each sampling site, at least 2 samples at dist 1 and at dist 2
#' do this just for one country (TZ), so 1 sample at dist 0, 2 at dist 1 and 2 at dist 2 per site
#' (n = 45 across the whole distance analysis)

