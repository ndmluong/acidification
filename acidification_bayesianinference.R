##### Package #####
library("ggplot2")
library("plotly")
library("grid")
library("lattice")
library("IDPmisc")
library("rjags")
library("emdbook")
library("GGally")
library("dplyr")
library("nlstools")

####################### DATA IMPORT ############################################
## The dataset corresponding to the pH measurement of different meat samples made
## in different batches ("Lot"), with different formulations (e.g. "Lactate"), 
## packed under several modified atmosphere ("Atm"), and measured at different time points ("Time")
## This dataset has to be structured as table with at least the following columns
##  - $SampleCode: the unique id the sample
##  - $Atm: the modified atmosphere
##  - $Lactate (can be renamed in the dataset and in the script if applied for other formulations): lactate concentrations
##  - $Time: sampling time (storage time)
##  - $pH
df <- read.table("Data_pH_mean.txt", h=T, sep="\t")


## Data plotting
gTimePH <- ggplot(data=ddf,
                  mapping=aes(Time, pH)) +
  labs(title="") +
  theme(axis.ticks=element_blank(),
        legend.position = "right",
        axis.title = element_text(size=14, face="bold"),
        strip.background = element_rect(fill="gray20"),
        strip.text = element_text(size = 11, face="bold", colour="white"),
        panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_line(colour="lightgrey"),
        panel.grid.major.x=element_line(colour="lightgrey"),
        panel.grid.minor.y=element_line(colour="lightgrey"),
        panel.grid.minor.x=element_line(colour="lightgrey")) +
  coord_cartesian(xlim=c(0,22),
                  ylim=c(4.8,7)) +
  geom_point(aes(Time,pH), colour="gray40", size=1) +
  geom_line(aes(Time,pH, group=Lot), colour="gray60") +
  scale_x_continuous(breaks=c(2,8,15,22)) +
  stat_summary(mapping = aes(Time, pH), fun = mean, geom="line", na.rm=T, size=1.5, colour="black") +
  facet_grid(Atm ~ Lactate) +
  xlab("Storage time") + 
  ylab("pH")
gTimePH




################################ INFERENCE BAYESIENNE #######################
##### Coding the observed data ##
ddf <- arrange(df, Lot, Atm, Time) ## arrange by batches, atmosphere, sampling times
ddf$Lot <- as.factor(ddf$Lot)
ddf$Atm <- as.factor(ddf$Atm)

Nr <- 10 ## total number of batches
Np <- 3 ## total number of packaging
Nt <- 12 ## total number of sampling points across all lactate conditions

pH <- array(NA, dim=c(Nr,Np,Nt))
Lactate <- array(NA, dim=c(Nr,Np,Nt))
Time <- array(NA, dim=c(Nr,Np,Nt))

for (r in 1:Nr) {
  for (p in 1:Np) {
    pH[r,p,] <- subset(ddf, Lot == levels(ddf$Lot)[r] & Atm == levels(ddf$Atm)[p])$pH
    Time[r,p,] <- subset(ddf, Lot == levels(ddf$Lot)[r] & Atm == levels(ddf$Atm)[p])$Time
    Lactate[r,p,] <- subset(ddf, Lot == levels(ddf$Lot)[r] & Atm == levels(ddf$Atm)[p])$Lactate
  }
}

rm(r,p)

dd <- list(Time=Time, Lactate=Lactate, pH=pH,
           Nr=Nr, Np=Np, Nt=Nt)

## initial values for each MCMC chain
ini <- list(list(deltaAir=0.05, deltaMAP1=0.05, deltaMAP2=0.05,
                 n=0.5,
                 theta=9,
                 mupH0=5.9, sigmapH0 = 0.05, sigmapH = 0.05),
            list(deltaAir=-0.05, deltaMAP1=-0.05, deltaMAP2=-0.05,
                 n=1,
                 theta=15,
                 mupH0=6, sigmapH0 = 0.04, sigmapH = 0.04),
            list(deltaAir=0.01, deltaMAP1=0.01, deltaMAP2=0.01,
                 n=1.5,
                 theta=13,
                 mupH0=6.1, sigmapH0=0.05, sigmapH = 0.03))

##### MCMC inference #####
set.seed(408)

## selection of the model
selected_model <- "acidification_model1.txt"
selected_model <- "acidification_model2.txt"
selected_model <- "acidification_model3.txt"
selected_model <- "acidification_model4.txt"

## initializing the model
model <- jags.model(file=selected_model, data=dd, inits=ini, n.chains=3)
update(model, 5000) ## burn-in phase

## Model 1
mcmc <- coda.samples(model, c("deltaAir","deltaMAP1", "deltaMAP2", "lambda", "mupH0", "sigmapH0", "sigmapH"), 
                     n.iter=60000, thin=6)
## Model 2
mcmc <- coda.samples(model, c("deltaAir","deltaMAP1", "deltaMAP2", "n", "lambda", "mupH0", "sigmapH0", "sigmapH"), 
                     n.iter=60000, thin=6)
## Model 3
mcmc <- coda.samples(model, c("deltaAir","deltaMAP1", "deltaMAP2", "theta", "lambda", "mupH0", "sigmapH0", "sigmapH"), 
                     n.iter=60000, thin=6)
## Model 4
mcmc <- coda.samples(model, c("deltaAir","deltaMAP1", "deltaMAP2", "n", "theta", "lambda", "mupH0", "sigmapH0", "sigmapH"), 
                     n.iter=60000, thin=6)


############### EXPLORING POSTERIOR DISTRIBUTION ########################
# Save full joint distribution 
mcmctot <- as.data.frame(as.matrix(mcmc))
mcmctotsample <- mcmctot[sample.int(nrow(mcmctot), size=500), ]

##### Chains visualization (posterior distribution)
xyplot(mcmc, col=c("navyblue", "darkgreen", "darkorange"), strip.left=T, strip=F, par.strip.text=list(cex=0.8), scales=list(cex=0.5))

##### Posterior distribution
summary(mcmc)

##### Critere Gelman-Rubin
## (ratio in standard deviation between/within chains: convergence if close to 1)
gelman.diag(mcmc)
gelman.plot(mcmc)

##### Critere Geweke
## (convergence if values between -2 et 2)
geweke.diag(mcmc)


##### Autocorrelation within chains
par(mfrow=c(3,ncol(mcmctot)))
autocorr.plot(mcmc, auto.layout = F)


##### Visualization of density function and correlation between parameters
# Graphical functions to be used
panel.hist <- function(x, col.hist="orange", ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr=c(usr[1:2], 0, 1.5))
  h <- hist(x, plot=F)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col=col.hist)
}

panel.dens <- function(x, col.dens="black", lwd.dens=2, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr=c(usr[1:2], 0, 1.5))
  densx <- density(x)
  vx <- densx$x
  vy <- densx$y
  lines(vx, vy/max(vy), col=col.dens, lwd=lwd.dens)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor=1.5, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr=c(0,1,0,1))
  r <- cor(x, y, method="spearman") ## correlation de spearman
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  
  #if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  
  text(0.5, 0.5, txt, cex=cex.cor)
}

panel.xy <- function(x, y, pch.xy=20, col.xy="gray20", cex.xy=0.5, ...) {
  points(x,y,pch=pch.xy, col=col.xy, cex=cex.xy)
}

pairs(mcmctotsample,
      upper.panel=panel.cor,
      diag.panel=panel.dens,
      lower.panel=panel.xy)


##### Comparison between models
# Deviance Information Criterion
DIC <- dic.samples(model, n.iter=10000, type="pD")


##### Save image 
save.image(file = "output_model1.RData")
save.image(file = "output_model2.RData")
save.image(file = "output_model3.RData")
save.image(file = "output_model4.RData")








################# POSTERIOR DISTRIBUTION: PLOT AND PREDICTION ################
##############################################################################
source("acidification_simulation_functions.R")
load(file = "output_model3.RData") ## chosen model
set.seed(408)

## Extraction a point estimate (median value) of each parameter
pt_estim <- apply(mcmctot, 2, quantile, probs=0.5)

## Examples for simulations of pH values from the joint distribution for one lactate concentration and one time point 
lactate_for_sim <- 0.5 ## chosen lactate concentrations for simulation
time_for_sim <- 10  ## chosen lactate concentrations for simulation
estimated_pH <- f_pred_pH(jd = mcmctot, 
                          time_pred = time_for_sim, 
                          lactate_pred = lactate_for_sim)
quantile(estimated_pH$pH_pred, probs=c(0.025, 0.5, 0.975)) ## median and IC 95%


## Examples for simulations of pH values from the joint distribution for a ranges of lactate concentrations and different time points 
lactate_for_sim <- c(0.25, 1.5) ## chosen lactate concentrations for simulation
time_for_sim <- seq(from=0, to=23, by=0.5)  ## chosen lactate concentrations for simulation
predicted_pH <- f_pred_pH(jd = mcmctot, time_pred = time_for_sim, lactate_pred = lactate_for_sim)

## Calculation of the 95% credibility band from the predicted value
## Using the implemented functions f_pred_pH and f_IC_band (cf. "Acidification_simulation_functions.R")
pH_IC_band <- f_IC_Band(predicted_data = predicted_pH)



##########################################################################
############## Simulation of the acidification rates #####################
##########################################################################
par(mfrow=c(1,1))
plot(NULL, xlab="Lactate concentration (% w/w)", ylab="Computed acidification rate", 
     xlim=c(0,2),ylim=c(0.05, 0.1), bty="n")
curve(expr = exp( median(mcmctot$lambda) * x + median(mcmctot$deltaAir)),
      from=0,to=5, col="darkgreen", add=T, lwd=2, lty=6)
curve(expr = exp(median(mcmctot$lambda) * x + median(mcmctot$deltaMAP1)),
      from=0,to=5, col="darkorange", add=T, lwd=2, lty=6)
curve(expr = exp(median(mcmctot$lambda) * x + median(mcmctot$deltaMAP2)),
      from=0,to=5, col="navyblue", add=T, lwd=2, lty=6)
text(0.1, y = 0.09, labels = "Air", col="darkgreen", offset = 0.5, cex=1)
text(0.2, y = 0.078, labels = "MAP1", col="darkorange", offset = 0.5, cex=1)
text(0.2, y = 0.098, labels = "MAP2", col="navyblue", offset = 0.5, cex=1)



##############################################################################
############ COMPARING DATA AND ADJUSTMENT BY THE MODEL ######################
##############################################################################
ddf_pred <- f_point_pred(obs = ddf, parms = pt_estim)

ddf_pred_residues <- data.frame(ddf_pred, Residuals = ddf_pred$pH_pred-ddf$pH)

par(mfrow=c(1,2))
plot(ddf_pred$pH, ddf_pred$pH_pred, 
     pch=20, col="gray60", bty="n", xlim=c(4.9,7), ylim=c(4.9, 7),
     xlab="observed pH", ylab="adjusted pH", bty="n")
abline(coef=c(0,1), col="gray10", lwd=1.5)
plot(ddf_pred$pH_pred, (ddf_pred$pH_pred-ddf_pred$pH), 
     pch=20, ylim=c(-1,1), col="gray60", bty="n", xlim=c(4.9,7), 
     xlab="adjusted pH", ylab="residuals")
abline(lm((ddf_pred$pH_pred-ddf_pred$pH)~ddf_pred$pH_pred), col="gray10", lwd=1.5)


##############################################################################
############ COMPARING PRIOR AND POSTERIOR DISTRIBUTIONS######################
##############################################################################
## Parameters: acidification rates
g_rate_prior_posterior <- ggplot(data=mcmctot,
                                 mapping=aes(x=deltaAir)) +
  theme(axis.ticks=element_blank(),
        legend.position = "right",
        axis.title = element_text(size=14, face="bold"),
        #axis.text.y = element_blank(),
        strip.background = element_rect(fill="gray20"),
        strip.text = element_text(size = 11, face="bold", colour="white"),
        panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_line(colour="lightgrey"),
        panel.grid.major.y=element_line(colour="lightgrey"),
        panel.grid.minor.y=element_line(colour="lightgrey"),
        panel.grid.minor.x=element_line(colour="lightgrey")) +
  coord_cartesian(xlim=c(-2.7, -2.1), ylim=c(0,0.015)) + 
  geom_density(mapping=aes(x=deltaAir, ..count../sum(..count..)), fill="darkgreen", colour="darkgreen", alpha=0.6) + 
  geom_density(mapping=aes(x=deltaMAP1, ..count../sum(..count..)), fill="darkorange", colour="darkorange", alpha=0.6) + 
  geom_density(mapping=aes(x=deltaMAP2, ..count../sum(..count..)), fill="navyblue", colour="navyblue", alpha=0.6) + 
  geom_density(data=data.frame(runif=runif(30000,-2.8,-2)),
               mapping=aes(x=runif, ..count../sum(..count..)), colour="gray20") + 
  annotate("text", x=median(mcmctot$deltaAir), y=0.0135, label="Air", colour="darkgreen") +
  annotate("text", x=median(mcmctot$deltaMAP1), y=0.013, label="MAP1", colour="darkorange") +
  annotate("text", x=median(mcmctot$deltaMAP2), y=0.0145, label="MAP2", colour="darkblue") +
  labs(title="") + xlab(expression(paste("acidification rate ",delta, sep=" "))) + ylab("")
g_rate_prior_posterior

## Parameters: average initial pH across batches
g_muph0_prior_posterior <- ggplot(data=mcmctot, mapping=aes(x=mupH0)) +
  theme(axis.ticks=element_blank(),
        legend.position = "right",
        axis.title = element_text(size=14, face="bold"),
        #axis.text.y = element_blank(),
        strip.background = element_rect(fill="gray20"),
        strip.text = element_text(size = 11, face="bold", colour="white"),
        panel.background=element_rect(fill="white"),
        panel.grid = element_line(colour="lightgrey")) +
  geom_density(data=data.frame(rnorm=rnorm(30000,5.84,0.11)),
               mapping=aes(x=rnorm, ..count../sum(..count..)), colour="gray20", size=1.2) +
  geom_density(mapping=aes(x=mupH0, ..count../sum(..count..)), fill="gray20", colour="gray20", alpha=0.6) +
  annotate("text", x=median(mcmctot$mupH0), y = 0.027, label=expression(mu[pH[0]]), colour="gray20") +
  #geom_density(data=subset(ddf, Time==2), aes(x=pH, ..scaled..)) +
  labs(title="") + xlab("average initial pH") + ylab("")
g_muph0_prior_posterior

## Parameters: standard deviation of pH across samples and batches
g_sigma_prior_posterior <- ggplot(data=mcmctot, mapping=aes(x=sigmapH0)) +
  theme(axis.ticks=element_blank(),
        legend.position = "right",
        axis.title = element_text(size=14, face="bold"),
        #axis.text.y = element_blank(),
        strip.background = element_rect(fill="gray20"),
        strip.text = element_text(size = 11, face="bold", colour="white"),
        panel.background=element_rect(fill="white"),
        panel.grid = element_line(colour="lightgrey")) +
  geom_vline(xintercept = 0) +
  geom_density(mapping=aes(x=sigmapH0, ..count../sum(..count..)), fill="gray20", colour="gray20", alpha=0.8) +
  geom_density(mapping=aes(x=sigmapH, ..count../sum(..count..)), fill="gray70", colour="gray70", alpha=0.8) +
  geom_density(data=data.frame(rnorm=abs(rnorm(30000,0,0.1))),
               mapping=aes(x=rnorm, ..count../sum(..count..)), colour="black", size=1.2) +
  geom_hline(yintercept = 0) +
  annotate("text", x= median(mcmctot$sigmapH0), y =0.024, label=expression(sigma[pH[0]])) +
  annotate("text", x= median(mcmctot$sigmapH), y =0.062, label=expression(sigma[pH])) +
  labs(title="") + xlab("standard deviation parameters") + ylab("")
g_sigma_prior_posterior

## Parameter: stabilisation time
g_theta_prior_posterior <- ggplot(data=mcmctot, mapping=aes(x=theta)) +
  theme(axis.ticks=element_blank(),
        legend.position = "right",
        axis.title = element_text(size=14, face="bold"),
        #axis.text.y = element_blank(),
        strip.background = element_rect(fill="gray20"),
        strip.text = element_text(size = 11, face="bold", colour="white"),
        panel.background=element_rect(fill="white"),
        panel.grid = element_line(colour="lightgrey")) +
  #geom_vline(xintercept = 0) +
  geom_density(mapping=aes(x=theta, ..count../sum(..count..)), fill="gray20", colour="gray20", alpha=0.8) +
  geom_density(data=data.frame(runif=runif(30000,8,22)),
               mapping=aes(x=runif, ..count../sum(..count..)), colour="gray20", size=1.2) + 
  #geom_hline(yintercept = 0) +
  annotate("text", x= median(mcmctot$theta), y =0.029, label=expression(theta)) +
  labs(title="") + xlab("stabilisation time (days)") + ylab("")
g_theta_prior_posterior

## All plots
pushViewport(viewport(layout = grid.layout(1,3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(g_muph0_prior_posterior, vp = vplayout(1, 1))
print(g_sigma_prior_posterior, vp = vplayout(1, 2))
print(g_theta_prior_posterior, vp = vplayout(1, 3))




