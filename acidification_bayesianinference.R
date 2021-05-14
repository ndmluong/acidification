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
# Full joint distribution 
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
# Fonctions necessaires pour les graphiques
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

load(file = "output_model1.RData")
