rm(list = ls())
##### Power analysis for Taya Forde's FLF app #####

### Aims

# This code will estimate power for the following objective:
# Obj 1: Determine the abundance and diversity of antimicrobial 
# resistance genes (ARG) in water and sediments in at-risk sites.
# Specifically, the aims of the analyses are to detect:
# 1a.i: Differences in ARG abundance among sources 
#     (null hypothesis: equal ARG abundance between sources)
# 1a.ii: Differences in ARG abundance between medium and high intensity sites 
#     (null hypothesis: equal ARG abundance between intensities)
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
# we assume the country effect is low enough to be subsumed into between-site noise.
# differences in AMU tend to be complex, and reflect communities that cross national 
# borders.

# source types and their mean relative abundances of AMR genes.
source.effect.list <-
  list(
    c(Human = 5, Pig = 5, Poultry = 5, Aquaculture = 1),
    c(Human = 5, Pig = 3, Poultry = 3, Aquaculture = 1),
    c(Human = 3, Pig = 3, Poultry = 3, Aquaculture = 1)) # this should give the least power

# choose effect sizes
source.effect <- source.effect.list[[3]]
sources <- names(source.effect)

# samples will be taken at the source (distance = 0) at all replicate sites,
# and one of each replicate, samples will be taken at distance = 1 (n = 3)
# and distance = 2 (n = 3)
distances <- c(0, 1, 1, 1, 2, 2, 2)

# two types of site will be sampled: medium and high intensity 
intensity.fold.diff <- 2
intensity.effect <- c(M = sqrt(intensity.fold.diff)/intensity.fold.diff, H = sqrt(intensity.fold.diff))
intensity.effect["H"]/intensity.effect["M"] # fold effect
exp(mean(log(intensity.effect))) # the geometric mean is 1, so no effect on mean log(ARG count)
intensities <- names(intensity.effect)

# number of replicates (no of sites per source type per country)
n.rep <- 3 * length(intensities)

# the n.rep replicates of the three sources will be repeated across each country
# so the total number of sites will be
length(source.effect) * n.rep * length(countries)

# relative decline of abundance with each additional distance unit
# in another study (Chu et al.) a ~60% decline per km has been observed.
#(this study was done using sediments; similar studies of water column unavailable)
# we will power the study to detect effects in this range, and smaller, down to
# a 25% reduction per distance unit
dist.effect <- 0.75

# mean log abundance at distance = 0 and source = human or livestock
intercept <- log(4)   
# based on Chu et al. (but the power estimate isn't affected by this assumption)

# SD at the observation level (variation within sites over repeated sampling)
# Rowe et al 2016 found only about 10% changes over time (across different years). However we are sampling in
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
dat <- 
  expand.grid(
    Country = countries, 
    rep = 1:n.rep, 
    Distance = distances, 
    Source = sources)
dat$Intensity <- factor(intensities[(dat$rep %% length(intensities)) + 1])


dat$site <- paste(dat$Country, dat$Source, dat$rep, sep = "-")
dat <- dat[order(dat$site), ]
dim(dat)

# distance sampling will only be done for one replicate, so delete observations with distance > 0
# in the other replicates
dat <- droplevels(dat[!(!dat$rep %in% c(1, n.rep) & dat$Distance > 0), ])
dim(dat)

# how many samples per country?
nrow(dat)/length(countries)

table(dat$site)
table(dat$Source, dat$Distance, dat$rep)

# dat encodes the sampling design
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
            Intensity = log(intensity.effect),
            Source = log(source.effect)),
        SD = SD,
        rand.V = c(site = SD.site^2))
    simdat$Abundance <- exp(simdat$response)
    simdat
  }

# test function
simdat.fn()

# plot an example of the simulated data
ggplot(data = simdat.fn(), 
       mapping = aes(x = Distance, y = Abundance, group = site, 
                     color = Intensity, shape = Country)) +
  geom_point() +
  stat_summary(fun = "mean", geom = "line") + 
  facet_wrap(~ Source)

# looping over multiple simulated data sets, test the two null hypotheses
# (objectives 1a and 1b), outputting the p-values
simres.list <- 
  mclapply(1:nsim, function(i) {
    
    # create subset data sets for each objective:
    simdat <- simdat.fn()
    # for objective 1a.i, keep all sites at distance = 0 
    simdat.1a.i <- droplevels(simdat[simdat$Distance == 0, ])
    # for objective 1a.ii, keep all sites at distance = 0 within a randomly selected source
    # (so we are averaging across within-source tests, which is more conservative than
    # assuming all sources have the same intensity effect)
    simdat.1a.ii <- 
      droplevels(simdat[simdat$Distance == 0 & simdat$Source == sample(levels(simdat$Source), 1), ])
    # for objective 1a.i, keep only sites that have spatial data (i.e. sampling at distance > 0) 
    # which is all rep = 1 and rep = n.rep sites
    simdat.1b <- droplevels(simdat[simdat$rep %in% c(1, n.rep), ])
    
    # objective 1a.i:
    # calculate a p-value for the null hypothesis that all sources have the same mean log abundance.
    # fit a Gaussian LMM to log abundance. Note that there is no random effect for site
    # here because, for this objective, there is only one observation per site. noise between sites 
    # and noise within sites are therefore combined in this model (conservatively).
    mod1.1a.i <- glm(response ~ Source + Intensity, data = simdat.1a.i)
    mod0.1a.i <- update(mod1.1a.i, ~ . - Source)
    p.1a.i <- anova(mod0.1a.i, mod1.1a.i, test = "Chisq")[2, "Pr(>Chi)"]
    
    # objective 1a.ii:
    # calculate a p-value for the null hypothesis that both intensities have the same mean log abundance.
    # fit a Gaussian GLM to log abundance. As above, there is no random effect for site
    # here because, for this objective, there is only one observation per site. noise between sites 
    # and noise within sites are therefore combined in this model (conservatively).
    mod1.1a.ii <- glm(response ~ Intensity, data = simdat.1a.ii)
    mod0.1a.ii <- update(mod1.1a.ii, ~ . - Intensity)
    p.1a.ii <- anova(mod0.1a.ii, mod1.1a.ii, test = "Chisq")[2, "Pr(>Chi)"]
    
    # objective 1b:
    # calculate a p-value for the null hypothesis that all distances have the same mean log abundance.
    # fit a random effect between the sites and test for a distance effect.
    # use maximum likelihood, not REML, as required for calculating comparable likelihoods. 
    mod1.1b <- lmer(response ~ Source + Distance + Intensity + (1 | site), data = simdat.1b, REML = FALSE)
    mod0.1b <- update(mod1.1b, ~ . - Distance)
    p.1b <- anova(mod0.1b, mod1.1b)[2, "Pr(>Chisq)"]
    
    # output p-values
    c(power.1a.i = p.1a.i, power.1a.ii = p.1a.ii, power.1b = p.1b, MRR.1a.i = mrr(sigma(mod1.1a.i)^2))
  }, mc.cores = detectCores() - 1) # parallelise across all but one core
simres <- do.call("rbind", simres.list)

# calculate power as the proportion of p-values < 0.05
apply(simres[, 1:3] < 0.05, 2, mean)

# look at distribution over all simulated data sets of multiplicative 
# effect of combined within-site and between site SDs
hist(simres[, "MRR.1a.i"])

# stop timer and show time elapsed
finish.time <- Sys.time()
print(finish.time - start.time)

#Chu et al. Metagenomics Reveals the Impact of Wastewater Treatment Plants on the Dispersal of Microorganisms and Genes in Aquatic Sediments.
#Applied and Environmental Microbiology (2018) 84(5):e02168-17
#Rowe et al. Comparative metagenomics reveals a diverse range of antimicrobial resistance genes in effluents entering a river catchment.
#Water Science and Technology (2016) 73(7):1541-9
