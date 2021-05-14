## Loading the estimated joing distribution (object: mcmctot)
load("RData/output_model3.RData")
source("acidification_simulation_functions.R")
set.seed(408)

########################################################################################################################################
########################################################################################################################################
## Estimated parameters from rjags
summary(mcmc)

########################################################################################################################################
########################################################################################################################################
## Extraction a point estimate (median value) of each parameter
pt_estim <- apply(mcmctot, 2, quantile, probs=0.5)
pt_estim

########################################################################################################################################
########################################################################################################################################
## Examples for simulations of pH values from the joint distribution for one lactate concentration and one time point 
lactate_for_sim <- 0.5 ## chosen lactate concentrations for simulation
time_for_sim <- 12  ## chosen lactate concentrations for simulation
simulated_pH <- f_pred_pH(jd = mcmctot, 
                          time_pred = time_for_sim, 
                          lactate_pred = lactate_for_sim)
estimated_pH <- as.data.frame(t(apply(simulated_pH$pH_pred, 3, quantile, probs=c(0.025, 0.5, 0.975))))
rownames(estimated_pH) <- c("Air", "MAP1:70%O2-30%CO2", "MAP2:50%CO2-50%N2")
colnames(estimated_pH) <- c("2.5% quantile", "Median", "97.5% quantile")
estimated_pH

########################################################################################################################################
########################################################################################################################################
## Examples for simulations of pH values from the joint distribution for a ranges of lactate concentrations and different time points 
lactate_for_sim <- c(0, 1, 2) ## chosen lactate concentrations for simulation
time_for_sim <- seq(from=0, to=10, by=2)  ## chosen lactate concentrations for simulation
simulated_pH_range <- f_pred_pH(jd = mcmctot, 
                                time_pred = time_for_sim, 
                                lactate_pred = lactate_for_sim)

## Calculation of the 95% credibility band from the predicted value
## Using the implemented functions f_pred_pH and f_IC_band (cf. "Acidification_simulation_functions.R")
pH_IC_band <- f_IC_Band(predicted_data = simulated_pH_range)

## Plotting simulations with point estimation and credible intervals
gTime_predPH <- ggplot(data=pH_IC_band, # change lactate concentrations for showing simulation if necessary
                       mapping=aes(Time, qInf)) +
  labs(title="Simulation example for different lactate concentrations") +
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
  coord_cartesian(xlim=c(1,23),
                  ylim=c(4.7,7)) +
  geom_ribbon(aes(x=Time, ymax=qSup, ymin=qInf, alpha=Lactate, fill=Atm), size=0.5) +
  geom_line(aes(Time, qMed, alpha=Lactate), colour="black", size=1) +
  scale_x_continuous(breaks=c(2,6,10,14,18,22)) +
  scale_alpha_manual(name = "Lactate (% w/w)",
                     values = seq(0.1,0.8,length.out=length(unique(pH_IC_band$Lactate)))) +
  scale_colour_manual(values = c("darkgreen", "darkorange", "navyblue")) +
  scale_fill_manual(values = c("darkgreen", "darkorange", "navyblue")) +
  guides(colour=F, fill=F) + 
  facet_grid(. ~ Atm) +
  xlab("Storage time (day)") + 
  ylab("Simulated pH")
gTime_predPH
