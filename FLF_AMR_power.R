rm(list = ls())
##### Power analysis for Taya Forde's FLF app #####

### Aims

# This code will estimate power for the following objective:
# Obj 1: Determine the abundance and diversity of antimicrobial 
# resistance genes (ARG) in water and sediments in at-risk sites.
# Specifically, the aims of the analyses are to detect:
# 1a: Differences in ARG abundance among sources 
#     (null hypothesis: equal ARG abundance between sources)
# 1b: Changes in ARG abundance in relation to the distance from sources 
#     (null hypothesis: equal ARG abundance across distances)


### General approach

# Because ARG abundance is generally right-skewed, we assume that residual variation will be
# approximately normally distributed on the log scale, which will convert 
# multiplicative effects on the abundance scale to additive effects on the log scale.
# An alternative would be to treat the abundances as counts and assume a count distribution 
# allowing for overdispersion (e.g. Poisson-lognormal, negative binomial, CMP), but given 
# the likely large scale of the counts (10s to 100s) this would likely make no substantial
# difference and would slow down the simulations.

# The power analysis will simulate data from the GLMM as specified below and estimate
# power as the proportion of tests giving P < 0.05.


### load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(ggplot2)
library(lme4)
library(parallel)


### settings, model parameters and design choices

# start timer
start.time <- Sys.time()

# how many simulations to run?
nsim <- 10000

# which countries?
countries <- c("TZ", "UG", "KE")
# we assume the country effect is low enough to be subsumed into between-site noise
#??????????? on what basis ?????????????
#??????????? on what basis ?????????????
#??????????? on what basis ?????????????
#??????????? on what basis ?????????????

# number of replicates (no of sites per source type per country)
n.rep <- 3

# source types and their mean relative abundances of AMR genes.
# (choose effect sizes by commenting out lines that aren't used.)
#source.effect <- c(Human = 5, Livestock = 5, Aquaculture = 1)
source.effect <- c(Human = 5, Livestock = 3, Aquaculture = 1)
#source.effect <- c(Human = 3, Livestock = 3, Aquaculture = 1)
sources <- names(source.effect)

# samples will be taken in al countries at the source (distance = 0),
# and in TZ samples will be taken at distance = 1 (n = 2)
# and distance = 2 (n = 2)
distances <- c(0, 1, 1, 2, 2)

# the n.rep replicates of the three sources will be repeated across each country
# so the total number of sites will be
length(source.effect) * n.rep * length(countries)

# relative decline of abundance with each additional distance unit
# in another study (Chu et al.) a ~60% decline per km has been observed.
# we will power the study to detect effects in this range, and smaller, down to
# a 25% reduction per distance unit
dist.effect <- 0.75

# mean log abundance at distance = 0 and source = human or livestock
intercept <- log(4)   
# seems reasonable based on Chu et al. 
# (in fact the power estimate isn't affected by this assumption)

# SD at the observation level (variation within sites over repeated sampling)
# Rowe et al 2016 found only about 10% changes over time. However we are sampling in
# different directions at the same distance, so we assume 50%. Convert this from a 
# relative rate to a SD following Biometrics, Vol. 56, No. 3 (Sep., 2000), pp. 909-914.
SD <- sqrt(inv.mrr(1.5))

# SD between sites (random intercepts)
SD.site <- sqrt(inv.mrr(2)) 
# this equates to ~2x diff between sites of same source based on Chu et al.

# note that the observation-level noise will contribute to the inter-site variation
# so that for objective 1a, the actual residual SD will be 
sqrt(SD^2 + SD.site^2)
# which equates to a rate ratio of
mrr(SD^2 + SD.site^2)


# set up study design
dat <- expand.grid(Country = countries, rep = 1:n.rep, Distance = distances, Source = sources)
dat$site <- paste(dat$Country, dat$Source, dat$rep, sep = "-")
dat <- dat[order(dat$site), ]

# distance sampling will only be done for TZ, so delete observations with distance > 0
# in the other two countries
dat <- droplevels(dat[!(dat$Country != "TZ" & dat$Distance > 0), ])

# dat encodes the samplng design
dat

# function to simulate log relative abundance
simdat.fn <-
  function() {
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
    simdat
  }

simdat.fn()

# plot an example of the simulated data
ggplot(data = simdat.fn(), mapping = aes(x = Distance, y = Abundance, group = site, 
                                         color = Source, shape = Country)) +
  geom_point() +
  stat_summary(fun = "mean", geom = "line")
  
# looping over multiple simulated data sets, test the two null hypotheses
# (objectives 1a and 1b), outputting the p-values
simres.list <- 
  mclapply(1:nsim, function(i) {
    
    # create subset data sets for each objective:
    simdat <- simdat.fn()
    simdat.1a <- droplevels(simdat[simdat$Distance == 0, ])
    simdat.1b <- droplevels(simdat[simdat$Country == "TZ", ])
    
    # objective 1a:
    # calculate a p-value for the null hypothesis that all sources have the same mean log abundance.
    # fit a Gaussian GLM to log abundance. Note that there is no random effect for site
    # here because, for this objective, there is only one observation per site. noise between sites 
    # and noise within sites are therefore combined in this model (conservatively).
    mod1.1a <- glm(response ~ Source, data = simdat.1a)
    mod0.1a <- update(mod1.1a, ~ . - Source)
    p.1a <- anova(mod0.1a, mod1.1a, test = "Chisq")[2, "Pr(>Chi)"]
    
    # objective 1b:
    # calculate a p-value for the null hypothesis that all distances have the same mean log abundance.
    # fit a random effect between the sites and test for a distance effect.
    # use maximum likelihood, not REML, as is required for comparable likelihoods. 
    mod1.1b <- lmer(response ~ Source + Distance + (1 | site), data = simdat.1b, REML = FALSE)
    mod0.1b <- update(mod1.1b, ~ . - Distance)
    p.1b <- anova(mod0.1b, mod1.1b)[2, "Pr(>Chisq)"]
    
    # output p-values
    c(power.1a = p.1a, power.1b = p.1b, MRR.1a = mrr(sigma(mod1.1a)^2))
  }, mc.cores = detectCores() - 1) # parallise across all but one cores
simres <- do.call("rbind", simres.list)

# calculate power as the proportion of p-values < 0.05
apply(simres[, 1:2] < 0.05, 2, mean)

# look at distribution over all simulated data sets of multiplicative 
# effect of combined within-site and between site SDs
hist(simres[, "MRR.1a"])

# stop timer and show time elapsed
finish.time <- Sys.time()
print(finish.time - start.time)
