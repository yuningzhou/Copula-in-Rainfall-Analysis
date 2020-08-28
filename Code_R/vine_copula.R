library("dplyr")
library("copula") # Copula package
library("lcopula")
library("tiger")
library("gsl")
library("PerformanceAnalytics")
library("VineCopula")

LSC <- read.csv("result_2005_2017_LSC.csv")
# SSC <- read.csv("result_2005_2017_SSC.csv")

## Filter
LSC = LSC%>% filter(LSC$durT > 1 | LSC$newstart == 0)
LSC <- LSC[,-c(1,2,3,4,16,17,18,19)]
LSC <- na.omit(LSC)

LSC <- LSC %>% filter(
  mmaxR < quantile(LSC$mmaxR, 0.95),
  Axismaj < quantile(LSC$Axismaj, 0.95),
  vel < quantile(LSC$vel, 0.95)
)

## Compute Kendall tau for all statistics
# 
# corrKendall = cor(LSC, method = c("kendall"))
# pairs(corrKendall)
# chart.Correlation(corrKendall, histogram=TRUE, method = c("kendall"))
# # histfit


## mmaxR, Axismaj, durT
X = LSC$mmaxR - 35
shape_X = 1.81776244814184
scale_X = 3.99936970775806

Y = LSC$Axismaj
mu_Y = 2.19381480694163
sigma_Y = 0.456145826153028

Z = LSC$durT
mu_Z = 0.813023327672115
sigma_Z = 0.829317405252584
# random term
aa = runif(lengths(Z), max = dlnorm(Z,mu_Z,sigma_Z), min = dlnorm(Z+1,mu_Z,sigma_Z))
Z = Z+aa

# transform to uniform
ux = to.uniform(X) #mmaxR
uy = to.uniform(Y) #Axismaj
uz = to.uniform(Z) #durT
U = cbind(ux, uy, uz) # U(0,1)^d
Um = U[1:1000,]
pairs(Um)
# BiCopSelect(uy, uz, familyset = c(1,3,4,5), selectioncrit = "AIC",
#             indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
#             se = FALSE, presel = TRUE, method = "mle")
# ux, uy frank
# uy. uz frank

RVM <- RVineStructureSelect(U, c(1,3,4,5,6), progress = TRUE)

# plot(RVM)
# plot(RVM, tree = "ALL", type = 0,
#      edge.labels = NULL, legend.pos = "bottomleft", interactive = FALSE)
# this line cannot run because:  NULL value passed as symbol address

## see the object's content or a summary
str(RVM)
# 
# List of 20
# $ Matrix     : num [1:3, 1:3] 1 2 3 0 2 3 0 0 3
# $ family     : num [1:3, 1:3] 0 3 5 0 0 5 0 0 0
# $ par        : num [1:3, 1:3] 0 0.409 3.68 0 0 ...
# $ par2       : num [1:3, 1:3] 0 0 0 0 0 0 0 0 0
# $ names      : chr [1:3] "ux" "uy" "uz"
# $ MaxMat     : num [1:3, 1:3] 1 2 3 0 2 3 0 0 3
# $ CondDistr  :List of 2
# ..$ direct  : logi [1:3, 1:3] FALSE TRUE TRUE FALSE TRUE TRUE ...
# ..$ indirect: logi [1:3, 1:3] FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ type       : chr "C-vine"
# $ tau        : num [1:3, 1:3] 0 0.17 0.364 0 0 ...
# $ taildep    :List of 2
# ..$ upper: num [1:3, 1:3] 0 0 0 0 0 0 0 0 0
# ..$ lower: num [1:3, 1:3] 0 0.184 0 0 0 ...
# $ beta       : num [1:3, 1:3] 0 0.166 0.407 0 0 ...
# $ call       : language RVineStructureSelect(data = U, familyset = c(1, 3, 4, 5, 6), progress = TRUE)
# $ nobs       : int 88605
# $ logLik     : num 27636
# $ pair.logLik: num [1:3, 1:3] 0 3599 14115 0 0 ...
# $ AIC        : num -55265
# $ pair.AIC   : num [1:3, 1:3] 0 -7195 -28229 0 0 ...
# $ BIC        : num -55237
# $ pair.BIC   : num [1:3, 1:3] 0 -7186 -28220 0 0 ...
# $ emptau     : num [1:3, 1:3] 0 0.183 0.368 0 0 ...
# - attr(*, "class")= chr "RVineMatrix"

summary(RVM)

# tree   edge |  No. family   par  par2 |  tau   UTD   LTD
# --------------------------------------------------------
#   1    3,1 |   5       F  3.68  0.00 | 0.36     -     -
#   3,2 |   5       F  3.00  0.00 | 0.31     -     -
#   2  2,1;3 |   3       C  0.41  0.00 | 0.17     -  0.18
# type: C-vine    logLik: 27635.71    AIC: -55265.43    BIC: -55237.25
# # inspect the fitted model using plots
# # Not run: plot(RVM) # tree structure

contour(RVM) # contour plots of all pair-copulas


set.seed(123)
simdata <- RVineSim(1000, RVM)
## determine the pair-copula families and parameters
RVM1 <- RVineCopSelect(simdata, familyset = c(1:6), RVM$Matrix)
## see the object's content or a summary
str(RVM1)
# List of 20
# $ Matrix           : num [1:3, 1:3] 1 2 3 0 2 3 0 0 3
# $ family           : num [1:3, 1:3] 0 16 5 0 0 5 0 0 0
# $ par              : num [1:3, 1:3] 0 1.31 3.79 0 0 ...
# $ par2             : num [1:3, 1:3] 0 0 0 0 0 0 0 0 0
# $ names            : chr [1:3] "ux" "uy" "uz"
# $ MaxMat           : num [1:3, 1:3] 1 2 3 0 2 3 0 0 3
# $ CondDistr        :List of 2
# ..$ direct  : logi [1:3, 1:3] FALSE TRUE TRUE FALSE TRUE TRUE ...
# ..$ indirect: logi [1:3, 1:3] FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ type             : chr "C-vine"
# $ tau              : num [1:3, 1:3] 0 0.149 0.372 0 0 ...
# $ taildep          :List of 2
# ..$ upper: num [1:3, 1:3] 0 0 0 0 0 0 0 0 0
# ..$ lower: num [1:3, 1:3] 0 0.302 0 0 0 ...
# $ beta             : num [1:3, 1:3] 0 0.141 0.416 0 0 ...
# $ nobs             : int 1000
# $ logLik           : num 323
# $ pair.logLik      : num [1:3, 1:3] 0 60.5 161.5 0 0 ...
# $ AIC              : num -641
# $ pair.AIC         : num(0) 
# $ BIC              : num -626
# $ pair.BIC         : num(0) 
# $ emptau           : num [1:3, 1:3] 0 0.15 0.368 0 0 ...
# $ p.value.indeptest: num [1:3, 1:3] 0.0 1.4e-12 0.0 0.0 0.0 ...
# - attr(*, "class")= chr "RVineMatrix"

summary(RVM1)

# tree   edge |  No. family   par  par2 |  tau   UTD   LTD 
# -------------------------------------------------------- 
#   1    3,1 |   5       F  3.79  0.00 | 0.37     -     -
#   3,2 |   5       F  2.91  0.00 | 0.30     -     -
#   2  2,1;3 |  16      SJ  1.31  0.00 | 0.15     -  0.30
# ---
#   type: C-vine    logLik: 323.45    AIC: -640.9    BIC: -626.18  

pairs(simdata)

corr1 = cor(U, method = c("kendall"))
corr2 = cor(simdata, method = c("kendall"))
write.csv(corr1,"vine Kendall correlation of the original data.csv")
write.csv(corr2,"vine Kendall correlation of the resample data.csv")