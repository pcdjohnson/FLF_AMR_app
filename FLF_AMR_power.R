# Power analysis for Taya's FLF app
library(GLMMmisc)
library(ggplot2)
library(hrbrthemes)

# set up study design
sources <- c("Human", "Livestock", "Aquaculture")
distances <- 0:2
n.rep <- 3
dat <- expand.grid(rep = 1:n.rep, Distance = distances, source = sources)
dat$rep <- paste(dat$source, dat$rep)

# model parameters
SD <- 0.2

# mean relative abundance in each source type
source.effect <- c(Human = 1, Livestock = 1, Aquaculture = 2)

# % decline with each distance unit
dist.effect <- 0.4


# siumulate log relative abundance
simdat <-
  sim.glmm(
    design.data = dat, 
    fixed.eff = 
      list(
        intercept = log(4),
        Distance = log(dist.effect), 
        source = log(source.effect)),
    SD = SD)
simdat$Abundance <- exp(simdat$response)

ggplot(data = simdat, mapping = aes(x = Distance, y = Abundance, group = rep, color = source)) +
  geom_point() +
  geom_line() #+
  #theme_ipsum()


