library("dplyr")
library("copula") # Copula package
library("lcopula")
library("tiger")
library("gsl")
library("PerformanceAnalytics")

LSC <- read.csv("result_2005_2017_LSC.csv")
# SSC <- read.csv("result_2005_2017_SSC.csv")

## Filter
LSC = LSC%>% filter(LSC$durT > 1 | LSC$newstart == 0)
LSC <- LSC[1:10000,-c(1,2,3,4,16,17,18,19)]
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


# XYZ = cbind(X, Y, Z)

# b = cor(XYZ, method = c("kendall"))
# xy 0.3155247, xz 0.4054274, yz 0.3560683

# param = matrix(c(shape_X, scale_X, mu_Y, sigma_Y, mu_Z, sigma_Z), nrow = 3, ncol = 2)

#tau <- 0.5

# # transform to uniform
# ux = to.uniform(X) #mmaxR
# uy = to.uniform(Y) #Axismaj
# uz = to.uniform(Z) #durT
# U = cbind(ux, uy, uz) # U(0,1)^d


xcdf = pgamma(X, shape_X, 1/scale_X)
ycdf = plnorm(Y, mu_Y,sigma_Y)
zcdf = plnorm(Z, mu_Z,sigma_Z)
# hist(xcdf)
# hist(ycdf)
# hist(zcdf)
U = cbind(xcdf, ycdf, zcdf)


# try also pseudo likelihood

## check uniform
#hist(ux)
#hist(uy)
#hist(uz)

## Copula fitting
# f.t <- fitCopula(tCopula(dim = 3), U, method = c("itau"), start = NULL)
summary(f.t <- fitCopula(tCopula(dim=3, dispstr="un"), U, method="itau")) 
# to.uniform
#Estimate Std. Error
# param   0.9907      0.006
# cdf
#Estimate Std. Error
#param -0.0003286     0.002
Ut <- cCopula(U, copula = f.t@copula) # conditional copula
splom2(Ut, cex = 0.2) 

## K plot
K.plot(Ut)

## resample

gofC = gofCopula(tCopula(dim=3, dispstr="un"), U, method="Sn")

r = mvdc(copula=tCopula(f.t@estimate, dim = 3),margins=c("gamma","lnorm","lnorm"), paramMargins = list(list(shape=shape_X, scale=scale_X),list(meanlog=mu_Y,sdlog=sigma_Y),list(meanlog=mu_Z,sdlog=sigma_Z)))
samp <- rMvdc(length(X), r)
x.samp = samp[, c(1)]
y.samp = samp[, c(2)]
z.samp = samp[, c(3)]
x.para = fitdistr(x.samp, 'gamma')  
#shape         rate    
# 1.816255394   0.248564947 
# (0.007963245) (0.001253627)
y.para = fitdistr(y.samp, densfun = "log-normal") 
# meanlog        sdlog   
# 2.193151835   0.455207252 
# (0.001529256) (0.001081347)
z.para = fitdistr(z.samp, densfun = "log-normal") 
# meanlog        sdlog   
# 0.817512917   0.828075333 
# (0.002781895) (0.001967097)


