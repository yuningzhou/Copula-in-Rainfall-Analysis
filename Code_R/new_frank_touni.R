library("dplyr")
library("copula") # Copula package
library("lcopula")
library("tiger")
library("gsl")

LSC <- read.csv("result_2005_2017_LSC.csv")
# SSC <- read.csv("result_2005_2017_SSC.csv")

## Filter
LSC = LSC%>% filter(LSC$durT > 1 | LSC$newstart == 0)
LSC <- LSC[,-c(1,4,17,18,19)]
LSC <- na.omit(LSC)

LSC <- LSC %>% filter(
  mmaxR < quantile(LSC$mmaxR, 0.95),
  Axismaj < quantile(LSC$Axismaj, 0.95),
  vel < quantile(LSC$vel, 0.95)
)

## Compute Kendall tau for all statistics

#corrKendall = cor(LSC, method = c("kendall"))
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


# XYZ = cbind(X, Y, Z)

# b = cor(XYZ, method = c("kendall"))
# xy 0.3155247, xz 0.4054274, yz 0.3560683

# param = matrix(c(shape_X, scale_X, mu_Y, sigma_Y, mu_Z, sigma_Z), nrow = 3, ncol = 2)

#tau <- 0.5

# transform to uniform
# 1. to.uniform
ux = to.uniform(X) #mmaxR
uy = to.uniform(Y) #Axismaj
uz = to.uniform(Z) #durT
U = cbind(ux, uy, uz) # U(0,1)^d

## check uniform
#hist(ux)
#hist(uy)
#hist(uz)

# 2. CDF
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




## Copula fitting
# f.t <- fitCopula(tCopula(dim = 3), U, method = c("itau"), start = NULL)
summary(f.t <- fitCopula(frankCopula(dim=3, use.indepC="FALSE"), U, method="itau"))
# to.uniform
#Estimate Std. Error
# param    3.268      0.004
# cdf
#Estimate Std. Error
#param 3.268      0.004
Ut <- cCopula(U, copula = f.t@copula) # conditional copula
splom2(Ut, cex = 0.2) 

## resample

# gofC = gofCopula(frankCopula(dim=3, use.indepC="FALSE"), U, method="Sn")

# para_gof = matrix(c(3.1, 8.04, 5e-04), nrow = 1, ncol = 3)
# 
# dimnames(para_gof) = list(c(""), c("parameter", "Statistic", "p-value"))
# 
# write.csv(para_gof, "para_gof_frank_touni.csv")

r = mvdc(copula=frankCopula(f.t@estimate, dim = 3),margins=c("gamma","lnorm","lnorm"), paramMargins = list(list(shape=shape_X, scale=scale_X),list(meanlog=mu_Y,sdlog=sigma_Y),list(meanlog=mu_Z,sdlog=sigma_Z)))
samp <- rMvdc(length(X), r)
x.samp = samp[, c(1)]
y.samp = samp[, c(2)]
z.samp = samp[, c(3)]
# x.para = fitdistr(x.samp, 'gamma') 
# #shape         rate    
# # 1.816275520   0.249832783 
# # (0.007963340) (0.001260021)
# y.para = fitdistr(y.samp, densfun = "log-normal") 
# # meanlog        sdlog   
# #2.193539443   0.457107617 
# # (0.001535640) (0.001085861)
# z.para = fitdistr(z.samp, densfun = "log-normal") 
# meanlog        sdlog   
# 0.814829015   0.826762192 
# (0.002777484) (0.001963977)
ux.samp = to.uniform(x.samp) #mmaxR
uy.samp = to.uniform(y.samp) #Axismaj
uz.samp = to.uniform(z.samp) #durT
U.samp = cbind(ux, uy, uz) # U(0,1)^d
# U.samp = U.samp[1:1000,]
pairs(U.samp)

pairs(samp)
K.plot(U)
K.plot(samp)

corr1 = cor(U, method = c("kendall"))
corr2 = cor(U.samp, method = c("kendall"))
write.csv(corr1,"frank Kendall correlation of the original data.csv")
write.csv(corr2,"frank Kendall correlation of the resample data.csv")