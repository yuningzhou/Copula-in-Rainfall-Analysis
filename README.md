# Copula-in-Rainfall-Analysis

This is a stochastic rainfall data analysis during 2005 - 2017. 
From the Raw Data, Kendall's correlation helps to select 3 independent factor. 
These three variables are assumed to be stochastically independent, and thus were selected to fit Copula models. 
When fitting into Copula models, the joint distribution of each pair was used to derive a multivariate rainfall time series. 
This requires the check of transformation to Uniform Distribution. 
Then the Copula models are built out of the Uniform Distributions, and checked with a p–value, the significance level of a hypothesis test. 
A proper copula family can be selected out of this.
This result allows for more rainfall data to be generated following the same pattern, in order to simulate the effect of rainfall on infrastructures for a longer time scale. 

