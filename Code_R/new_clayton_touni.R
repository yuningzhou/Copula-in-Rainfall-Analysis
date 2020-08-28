library("dplyr")
library("copula") # Copula package
library("lcopula")
library("tiger")
library("gsl")

LSC <- read.csv("result_2005_2017_LSC.csv")
# SSC <- read.csv("result_2005_2017_SSC.csv")

## Filter
LSC = LSC%>% filter(LSC$durT > 1 | LSC$newstart == 0)
LSC <- LSC[1:10000,-c(1,4,17,18,19)]
LSC <- na.omit(LSC)

LSC <- LSC %>% filter(
  mmaxR < quantile(LSC$mmaxR, 0.95),
  Axismaj < quantile(LSC$Axismaj, 0.95),
  vel < quantile(LSC$vel, 0.95)
)

## Compute Kendall tau for all statistics

corrKendall = cor(LSC, method = c("kendall"))
library(GGally)
corrKendall = as.data.frame(corrKendall)

ggpairs(corrKendall,title = NULL,
        upper = list(continuous = "cor", combo = "box_no_facet", discrete =
                       "facetbar", na = "na"), 
        lower = list(continuous = "smooth", combo ="dot", discrete = "facetbar", na = "na"), diag = list(continuous =
                                                                                                                            "densityDiag", discrete = "barDiag", na = "naDiag"), params = NULL,
        xlab = NULL, ylab = NULL, axisLabels = c("show", "internal", "none"),
        labeller = "label_value",
        switch = NULL, showStrips = TRUE, legend = NULL,
        cardinality_threshold = 15, progress = NULL)
# histfit


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

# x <- round(rgamma(100000,shape_X,1/scale_X),1)
# y <- rlnorm(100000, mu_Y, sigma_Y)
# z <- rlnorm(100000, mu_Z, sigma_Z)  
# # hist(x)
# # hist(y)
# # hist(z)
# xcdf = pgamma(x, shape_X, 1/scale_X)
# ycdf = plnorm(y, mu_Y,sigma_Y)
# zcdf = plnorm(z, mu_Z,sigma_Z)
# # hist(xcdf)
# # hist(ycdf)
# # hist(zcdf)
# U = cbind(xcdf, ycdf, zcdf)


# try also pseudo likelihood

## check uniform
#hist(ux)
#hist(uy)
#hist(uz)

## Copula fitting
# f.t <- fitCopula(tCopula(dim = 3), U, method = c("itau"), start = NULL)
summary(f.t <- fitCopula(claytonCopula(dim=3, use.indepC="FALSE"), U, method="itau"))
# to.uniform
#Estimate Std. Error
# param  0.9878      0.006
# cdf
#Estimate Std. Error
#param 0.9874      0.007
Ut <- cCopula(U, copula = f.t@copula) # conditional copula
splom2(Ut, cex = 0.2) 


## resample

# gofC = gofCopula(claytonCopula(dim=3, use.indepC="FALSE"), U, method="Sn")

# para_gof = matrix(c(0.542, 96.4, 5e-04), nrow = 1, ncol = 3)
# 
# dimnames(para_gof) = list(c(""), c("parameter", "Statistic", "p-value"))
# 
# write.csv(para_gof, "para_gof_clayton_touni.csv")

r = mvdc(copula=claytonCopula(f.t@estimate, dim = 3),margins=c("gamma","lnorm","lnorm"), paramMargins = list(list(shape=shape_X, scale=scale_X),list(meanlog=mu_Y,sdlog=sigma_Y),list(meanlog=mu_Z,sdlog=sigma_Z)))
samp <- rMvdc(length(X), r)
x.samp = samp[, c(1)]
y.samp = samp[, c(2)]
z.samp = samp[, c(3)]
# x.para = fitdistr(x.samp, 'gamma') 
# fitdistr(x.samp, 'gamma') 
# fitdistr(y.samp, 'log-normal') 
# fitdistr(z.samp, 'log-normal') 
# #shape         rate    
# # 1.822372165   0.250048225 
# # (0.007991871) (0.001260807)
# y.para = fitdistr(y.samp, densfun = "log-normal") 
# # meanlog        sdlog   
# # 2.194009011   0.456414069 
# # (0.001533310) (0.001084214)
# z.para = fitdistr(z.samp, densfun = "log-normal") 
# meanlog        sdlog   
# 0.819965322   0.829293820 
#(0.002785988) (0.001969991)
ux.samp = to.uniform(x.samp) #mmaxR
uy.samp = to.uniform(y.samp) #Axismaj
uz.samp = to.uniform(z.samp) #durT
U.samp = cbind(ux, uy, uz) # U(0,1)^d
U.samp = U.samp[1:1000,]
pairs(U.samp)

## K plot
 K.plot(U)
 K.plot(samp)
 
 corr1 = cor(U, method = c("kendall"))
 corr2 = cor(U.samp, method = c("kendall"))
 write.csv(corr1,"clayton Kendall correlation of the original data.csv")
 write.csv(corr2,"clayton Kendall correlation of the resample data.csv")