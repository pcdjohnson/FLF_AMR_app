rm(list = ls())
##### Power analysis for Taya Forde's FLF app #####

### General approach

# Because ARG abundance is generally right-skewed, we assume that residual variation will be
# approximately normally distributed on the log scale, which will convert 
# multiplicative effects on the count scale to additive effects on the log scale.
# An alternative would be to treat the abundances as counts and assume a count distribution 
# allowing for overdispersion (e.g. Poisson-lognormal, negative binomial, CMP), but given 
# the likely large scale of the counts (100s to 1000s) this would likely make no substantial
# difference and would slow down the simulations.

# The power analysis will simulate data from the GLMM as specified below and estimate
# power as the proportion of tests giving P < 0.05.

### Aims

# This code will estimate power for the following objective:
# Obj 1: Determine the abundance and diversity of ARG in water and sediments in at-risk sites
# Specifically the aims of the analyses are to detect:
# 1a: Differences in ARG abundance among sources 
#     (null hypothesis: equal ARG abundance between sources)
# 1b: Changes in ARG abundance in relation to the distance from sources 
#     (null hypothesis: equal ARG abundance across distances)


### load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(ggplot2)
#library(hrbrthemes)
library(lme4)
library(lmerTest)


### settings, model parameters and design choices

# start timer
start.time <- Sys.time()

# how many simulations to run?
nsim <- 100

# which countries?
countries <- c("TZ", "UG", "KE")

# number of replicates (no of sites per source type)
n.rep <- 3

# source types and their mean relative abundances of AMR genes
source.effect <- c(Human = 5, Livestock = 3, Aquaculture = 1)
sources <- names(source.effect)

# samples will be taken at the source (distance = 0)
# and in TZ samples will be taken at distance = 1 (n = 2)
# and distance = 2 (n = 2)
distances <- c(0, 1, 1, 2, 2)

# the n.rep replicates of the three sources will be repeated across each country
# so the total number of sites will be
length(source.effect) * n.rep * length(countries)
# we assume the country effect is low enough to be subsumed into between-site noise
#??????????? on what basis ?????????????

# relative decline of abundance with each additional distance unit
dist.effect <- 0.4 ### 60% decline per distance unit (based on Chu)

# mean log abundance at distance = 0 and source = human or livestock
intercept <- log(4)   ##### probably doesn't make a difference - check  

# SD at the observation level (variation within site over repeated sampling)
# Rowe et al 2016 found only about 10% changes over time. Convert this to a 
# relative rate following Biometrics, Vol. 56, No. 3 (Sep., 2000), pp. 909-914.
SD <- inv.mrr(1.1)

# SD between sites (random intercepts)
SD.site <- inv.mrr(2) #### this should equate to ~2x diff between sites of same source based on Chu

# set up study design
dat <- expand.grid(Country = countries, rep = 1:n.rep, Distance = distances, Source = sources)
dat$site <- paste(dat$Country, dat$Source, dat$rep, sep = "-")
dat <- dat[order(dat$site), ]

# distance sampling will only be done for TZ, so delete observations with distance > 0
# in the other two countries
dat <- droplevels(dat[!(dat$Country != "TZ" & dat$Distance > 0), ])

# dat encodes the design of the study
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
  geom_line() 


simres.list <- 
  lapply(1:nsim, function(i) {
    # create subset data sets for each objective:
    simdat <- simdat.fn()
    simdat.1a <- droplevels(simdat[simdat$Distance == 0, ])
    simdat.1b <- droplevels(simdat[simdat$Country == "TZ", ])
    
    # objective 1a:
    # calculate a p-value for the null hypothesis that all sources have the same mean log abundance
    mod1.1a <- glm(response ~ Source, data = simdat.1a)
    mod0.1a <- update(mod1.1a, ~ . - Source)
    p.1a <- anova(mod0.1a, mod1.1a, test = "Chisq")[2, "Pr(>Chi)"]
    
    # objective 1b:
    # calculate a p-value for the null hypothesis that all distances have the same mean log abundance
    mod1.1b <- lmer(response ~ Source + Distance + (1 | site), data = simdat.1b, REML = FALSE)
    mod0.1b <- update(mod1.1b, ~ . - Distance)
    p.1b <- anova(mod0.1b, mod1.1b)[2, "Pr(>Chisq)"]
    
    # output p-values
    c(power.1a = p.1a, power.1b = p.1b)
  })

simres <- do.call("rbind", simres.list)

apply(simres < 0.05, 2, mean)

# stop timer and show time elapsed
finish.time <- Sys.time()
print(finish.time - start.time)
