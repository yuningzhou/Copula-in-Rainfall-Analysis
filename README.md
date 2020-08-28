# Copula-in-Rainfall-Analysis

This is a stochastic rainfall data analysis during 2005 - 2017. 
From the Raw Data, Kendall's correlation helps to select 3 independent factor. 
These three variables are assumed to be stochastically independent, and thus were selected to fit Copula models. 
When fitting into Copula models, the joint distribution of each pair was used to derive a multivariate rainfall time series. 
This requires the check of transformation to Uniform Distribution. 
Then the Copula models are built out of the Uniform Distributions, and checked with a p–value, the significance level of a hypothesis test. 
A proper copula family can be selected out of this.
This result allows for more rainfall data to be generated following the same pattern, in order to simulate the effect of rainfall on infrastructures for a longer time scale. <br>

I tested several copula families with LSC data in R, including: <br>
(1)	t-copula/Norm Copula (from ellipCopula) <br>

(2)	Archimedean Copula: <br>
{Clayton Copula, Frank Copula, and Gumbel Copula) <br>

(3)	Asymmetric Copula: <br>
PCC <br>

## 1.	Analyse observed data:
### 1.1 Analyzation  
The raw data, LSC, is the rain gauge data from 2005 to 2017 with 19 variables. In these variables, only quantitative measurements are reserved for further calculation. 
And the number of storm, whether it has a new start, and the variables with NaN value are omitted. In order to rule out the influence of the extreme value, top 5% value of the measurements such as major axis, mmaxR, and velocity are trimmed off. 
After data cleaning, correlation between the remaining variables are calculated to figure out which variables will show more strong correlation with each other. Here Kendall’s tau is used with: <br>
corrKendall = cor(LSC, method = c("kendall")). <br>
And pairs plot is used to visualize the results, including the Kendall’s correlation coefficient and the ggplot between each two variables. <br>
 
From the result of the correlations, mmaxR, Axismaj, and durT are selected to fit the copula model. These three variables should be transformed into a uniform distribution before fitting a copula model. Both empirical and theoretical methods are tried in this transformation. The goodness of transformation is checked through resample. <br>

The resampled data can be checked through calculating and comparing their parameters and making scatter plots. The parameters of resample data is close to the original data through both empirical and theoretical method. However, on scatter plots, the shape of distributions shows that empirical method using “to.uniform” function is a better choice, because the shape of distribution using CDF transformation has clear lower limits and leaves some spaces white. <br>

K plots are used to check the goodness of fit through comparing the similarity of the plots of resample data and original data. K plot of each copula family is plotted to make comparison with the original plot below. <br>

### 1.2	Correlation Visualization

Two pairs plots presenting the same result while adopting different format. One uses “ggpairs” in “ggplot2”, the other uses “chart.correlation”. 
A preferable plot can be chosen from them. I would suggest the first one because its linear regression line shows the correlation more clearly. 

## 2.	Copula Fitting
### 2.1	T Copula & Normal Copula

For elliptical copula, t copula and normal copula is tested. <br>
f.t <- fitCopula(tCopula(dim=3, dispstr="un"), U, method="itau") <br>
f.t <- fitCopula(normalCopula(dim=3, dispstr="un"), U, method="itau") <br>

### 2.2	Frank Copula
f.t <- fitCopula(frankCopula(dim=3, use.indepC="FALSE"), U, method="itau") 

### 2.3	 Gumbel Copula
f.t <- fitCopula(gumbelCopula(dim=3, use.indepC="FALSE"), U, method="itau")

### 2.4	Clayton Copula
f.t <- fitCopula(claytonCopula(dim=3, use.indepC="FALSE"), U, method="itau")

### 2.5	 Vine Copula

After the same process of transformation to uniform distribution, dataset U is input into function “RVineStructureSelect” in package “VineCopula”. The input of this command includes type of copula to try for each pair-copula group, where 1,3,4,5,6 represents Gaussian, Clayton, Gumbel, Frank, and Joe respectively. The index of these families of copula is in “RVineMatrix” in the same package. <br>

RVM <- RVineStructureSelect(U, c(1,3,4,5,6), progress = TRUE) <br>



