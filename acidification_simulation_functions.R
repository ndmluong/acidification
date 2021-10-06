####################################################################
####################################################################
## Fonction for predicting, from the joint distribution, pH value at one or several given time points
f_pred_pH <- function(jd, time_pred, lactate_pred) {
  ## Fonction enabling to predict, from the posterior joint distribution of parameters,
  ## the pH values for each process conditions and at different given time points.
  ## ## INPUT:
  ##      - jd : data frame. Posterior joint distribution of parameters (MCMC) with the columns (model parameters) below
  ##            + deltaAir: ## effect of packaging under Air on acidification rate
  ##            + deltaMAP1: ## effect of packaging under MAP1:70%O2-30%CO2 on acidification rate
  ##            + deltaMAP2: ## effect of packaging under MAP2:70%O2-30%CO2 on acidification rate
  ##            + lambda: 
  ##            + mupJH0: ## mean initial pH value
  ##            + sigmapH0: ## standard deviation of the initial pH across batches
  ##            + sigmapH: ## standard deviation of the pH across samples
  ##            + theta: stabilisation time
  ##      - time_pred (numeric vector). Time points at which the pH values are predicted for each process conditions
  ##      - lactate_pred (numeric vector). Lactate concentration for prediction
  ## ## OUTPUT:
  ##      - pH_pred: three-way array [i,p,tp] containing the predicted value for the iteration i, 
  ##                packaging p, for the lactate concentration lactate_pred and at time point tp
  
  Niter <- nrow(jd) ## total number of MCMC iterations
  Np <- 3 ## number of packaging conditions
  Nt_pred <- length(time_pred) ## total number of time points at which predictions are computed
  Nl_pred <- length(lactate_pred)
  
  
  pH_pred <- array(NA, dim=c(Niter, Nl_pred, Np, Nt_pred)) ## four-way array : iteration, lactate, packaging and time points
  
  # print("Total number of computed iterations") ## show the computation progress every 500 iterations
  
  for (i in 1:Niter) { ## for the iteration i
    # calculation of an initial pH value (pH0_pred) randomly sampled from a normal law with mean = mupH0[i] and
    # sd = sigmapH0[i] from joint distribution (iteration i)
    pH0_pred <- rnorm(1, mean = jd$mupH0[i], sd = jd$sigmapH0[i])
    
    # extracting the acidification rates corresponding to the iteration i
    delta_vec <- c(jd$deltaAir[i], jd$deltaMAP1[i], jd$deltaMAP2[i])
    
    for (l in 1:Nl_pred) {
      for (p in 1:length(delta_vec)) { ## for each above process condition
        
        for (tp in 1:Nt_pred) {## at each sampling time point
          if (time_pred[tp] < jd$theta[i]) { # if the sampling time point is earlier than the stabilisation time theta
            mean_pH_tp <- -exp(jd$lambda[i] * lactate_pred[l] + delta_vec[p]) * time_pred[tp] + pH0_pred 
          } else { ## if the sampling time point is later than the stabilisation time theta
            mean_pH_tp <- -exp(jd$lambda[i] * lactate_pred[l] + delta_vec[p]) * jd$theta[i] + pH0_pred
          }
          
          # the predicted pH for the iteration i, condition k and time point tp is randomly sampled from the normal law
          # with the mean mean_pH_tp and standard deviation sigmapH corresponding to the iteration i
          pH_pred[i, l, p, tp] <- rnorm(1, mean_pH_tp, jd$sigmapH[i]) 
        }
      }
    }
    
    # if  (i %% 500 == 0) {print(i)} ## show the computation progress every 500 iterations
    
  }
  
  return(list("lactate_pred" = lactate_pred,
              "time_pred" = time_pred,
              "pH_pred" = pH_pred)) ## return the predicted values for e
  
}


####################################################################
####################################################################
## Calculation of the 95% credibility band from the predicted value
f_IC_Band <- function(predicted_data) {
  
  IC_band <- data.frame()
  
  Atm <- c("Air", "MAP1", "MAP2") ## to be renamed conveniently
  
  for (l in 1:length(predicted_data$lactate_pred)) {
    for (p in 1:length(Atm)) {
      for (tp in 1:length(predicted_data$time_pred)) {
        IC_point <- cbind(Atm[p],
                          predicted_data$lactate_pred[l],
                          predicted_data$time_pred[tp], 
                          unname(quantile(predicted_data$pH_pred[,l,p,tp], probs=0.025)), ## quantile 2.5%
                          unname(quantile(predicted_data$pH_pred[,l,p,tp], probs=0.25)), ## quantile 25%
                          unname(quantile(predicted_data$pH_pred[,l,p,tp], probs=0.5)), ## median value 
                          unname(quantile(predicted_data$pH_pred[,l,p,tp], probs=0.75)),  ## quantile 75%
                          unname(quantile(predicted_data$pH_pred[,l,p,tp], probs=0.975)))  ## quantile 97.5%
        IC_band <- rbind(IC_band, IC_point)
      }
    }
  }
  
  colnames(IC_band) <- c("Atm", "Lactate", "Time", "qInf", "qInfMed", "qMed", "qSupMed", "qSup") ## to be renamed if needed
  
  IC_band$Time <- as.numeric(IC_band$Time)
  IC_band$qInf <- as.numeric(IC_band$qInf)
  IC_band$qInfMed <- as.numeric(IC_band$qInfMed)
  IC_band$qMed <- as.numeric(IC_band$qMed)
  IC_band$qSupMed <- as.numeric(IC_band$qSupMed)
  IC_band$qSup <- as.numeric(IC_band$qSup)
  
  colnames(IC_band) <- c("Packaging", "Lactate", "Time", "2.5% quantile", "25% quantile", "Median", "75% quantile", "97.5% quantile")
  
  return(IC_band)
}

## Calculation of the 95% credibility band from the predicted value
f_point_pred <- function(obs, parms) {
  ## INPUT
  ##  - obs: dataframe. Observed data with at least the columns
  ##    + Time
  ##    + Atm (factor) with the following levels (to be adapted if needed)
  ##        + Air
  ##        + MAP1
  ##        + MAP2
  ##  - pt_estim: point estimate of model parameters
  ## OUTPUT
  ##  - pred: data.frame with added column containing predicted data
  
  pH_pred <- rep(NA, nrow(obs))
  
  Nr <- length(levels(obs$Lot)) ## number total of batches
  
  pH0_all <- rnorm(Nr, parms["mupH0"], parms["sigmapH0"]) ## initial value for each batch
  pH0_tab <- data.frame(Lot = as.vector(levels(obs$Lot)),
                        pH0 = pH0_all)
  
  for (i in 1:nrow(obs)) { ## for each row of the observed data
    switch (as.character(obs$Atm[i]), ## extract the process condition and the corresponding acidificiation rate  
            'Air' = {delta <- parms["deltaAir"]},
            'MAP1' = {delta <- parms["deltaMAP1"]},
            'MAP2' = {delta <- parms["deltaMAP2"]}
    )
    
    pH0 <- subset(pH0_tab, Lot == as.character(obs$Lot[i]))$pH0
    
    if (obs$Time[i] <= parms["theta"]) {
      pH_pred[i] <- -exp(parms["lambda"] * obs$Lactate[i] + delta) * obs$Time[i] + pH0
    } else {
      pH_pred[i] <- -exp(parms["lambda"] * obs$Lactate[i] + delta) * parms["theta"] + pH0
    }
  }
  
  pred.data <- data.frame(obs,
                          pH_pred = pH_pred)
  
  return(pred.data)
}
