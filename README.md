[![DOI:10.3390/foods11081114](doi.org/10.3390/foods11081114)](doi.org/10.3390/foods11081114)


## Overview
This repository provides several R scripts enabling to perform numerical simulations for pH kinetics of meat products under different user-defined formulation-atmosphere conditions. These scripts give possibilities to apply a Bayesian modelling procedure to describe pH changes and to estimate acidification rates for experimental data collected on food matrices, with user-defined atmosphere factors and/or formulations (expressed as concentrations). 

## The provided dataset ## 
The dataset corresponding to our study is provided: pH measurement of different meat samples made in different batches ("Lot"), with different formulations ("Lactate"), packed under several modified atmosphere ("Atm"), and measured at different time points ("Time").

The formulation used in our is potassium lactate. The lactate contents, denoted Lactate, is included in the model as a variate. The effect of the different modified atmosphere conditions modelled as multi-level factor (3 levels in our study) are denoted 'deltaAir', 'deltaMAP1', 'deltaMAP2'. For further use of the model, rename these variates in the scripts conveniently as well as in the dataset. 

## To get started
The collected experimental data should be structured in table (Excel spreadsheet or .txt files) with at least the following columns
- SampleCode: the unique ID of the sample
- Atm: the modified atmosphere
- Lactate (can be renamed in the dataset and in the scripts if applied for other formulations): lactate concentrations
- Time: sampling time (storage time)
- pH

## Modelling procedure
Check the R script: "acidification_bayesianinference.R"

## Associated article (in submission)
A Bayesian approach to describe and simulate the pH evolution of fresh meat products depending on the preservation conditions.

Authors: Ngoc-Du Martin Luong, Louis Coroller, Monique Zagorec, Nicolas Moriceau, Valérie Anthoine, Sandrine Guillou, Jeanne-Marie Membré

R/git contributor: Ngoc-Du Martin Luong (ngoc-du.luong@anses.fr)

Corresponding author: Jeanne-Marie Membré (jeanne-marie.membre@oniris-nantes.fr)
