setwd("C:/Users/Olympos/Documents/Programs/R/STATS-Final")
library(doParallel)

set.seed(87539319)
CL = detectCores() - 2
registerDoParallel(cores = CL)

ndata = 5000
nboot = 2500
sim = 250
b0 = 1
b1 = 3
x = runif(ndata)

eps1 = rnorm(ndata)
yNorm = b0 + b1 * x+eps1

esp2 = rt(ndata,df=3)/sqrt(3)
yTail = b0 + b1 * x+esp2

esp3 = rgamma(ndata,shape=4,rate=2)-2
ySkew = b0 + b1 * x+esp3

esp4 = rnorm(ndata)*x^2/sqrt(mean(x^4))
esp4 = esp4-mean(esp4)
yWide = b0 + b1 * x+esp4


dataNorm<-data.frame(x,yNorm)
dataTail<-data.frame(x,yTail)
dataSkew<-data.frame(x,ySkew)
dataWide<-data.frame(x,yWide)

write.csv(dataNorm, file = 'dataNorm.csv', row.names = F)
write.csv(dataTail, file = 'dataTail.csv', row.names = F)
write.csv(dataSkew, file = 'dataSkew.csv', row.names = F)
write.csv(dataWide, file = 'dataWide.csv', row.names = F)


bootstrap = function(df, n){
  foreach(i=1:n, .combine=c) %dopar% {
    bootstrap_data<-df[sample(nrow(df),nrow(df),replace=T),]
    unname(lm(bootstrap_data[,2]~bootstrap_data[,1])$coef[2])
  }
}


parametric = function(df, n){
  mod = lm(df[,2] ~ df[,1])
  yhat = mod$coefficients[1] + mod$coefficients[2] * df[,1]
  foreach(i=1:n, .combine=c) %dopar% {
    tmp = rnorm(length(yhat), sd = (summary(mod)$sigma)^2)
    ystar = unlist(yhat) + tmp
    cat(paste0(c(length(tmp), length(ystar))))
    unname(lm(ystar~df[,1])$coef[2])
  }
} 

resampling = function(df, n){
  mod = lm(df[,2] ~ df[,1])
  yhat = mod$coefficients[1] + mod$coefficients[2] * df[,1]
  e = (df[,2] - yhat)/(1 - hatvalues(mod))
  r = e - mean(e)
  foreach(i=1:n, .combine=c) %dopar% {
    boot = sample(length(r), replace = T)
    ystar = yhat[boot] + r[boot]
    xstar = df[,1][boot]
    summary(lm(ystar~ df[,1]))
    unname(lm(ystar~xstar)$coef[2])
  }
}

smoothed = function(df, n){
  mod = lm(df[,2] ~ df[,1])
  yhat = mod$coefficients[1] + mod$coefficients[2] * df[,1]
  e = (df[,2] - yhat)/(1 - hatvalues(mod))
  r = e - mean(e)
  f = ecdf(r)
  foreach(i=1:n, .combine=c) %dopar% {
    u = runif(length(r))
    err = unname(quantile(f,u))
    ystar = yhat + err
    unname(lm(ystar~df[,1])$coef[2])

  }
}

wildstrap = function(df, n){
  mod = lm(df[,2] ~ df[,1])
  yhat = mod$coefficients[1] + mod$coefficients[2] * df[,1]
  e = (df[,2] - yhat)/(1 - hatvalues(mod))
  r = e - mean(e)
  foreach(i=1:n, .combine=c) %dopar% {
    v = rnorm(length(r))
    ystar = yhat + r * v
    unname(lm(ystar~df[,1])$coef[2])

  }
}

simNormBoot = replicate(sim, bootstrap(dataNorm, nboot))
simNormPara = replicate(sim, parametric(dataNorm, nboot))
simNormRsam = replicate(sim, resampling(dataNorm, nboot))
simNormSmth = replicate(sim, smoothed(dataNorm, nboot))
simNormWild = replicate(sim, wildstrap(dataNorm, nboot))

simTailBoot = replicate(sim, bootstrap(dataTail, nboot))
simTailPara = replicate(sim, parametric(dataTail, nboot))
simTailRsam = replicate(sim, resampling(dataTail, nboot))
simTailSmth = replicate(sim, smoothed(dataTail, nboot))
simTailWild = replicate(sim, wildstrap(dataTail, nboot))

simSkewBoot = replicate(sim, bootstrap(dataSkew, nboot))
simSkewPara = replicate(sim, parametric(dataSkew, nboot))
simSkewRsam = replicate(sim, resampling(dataSkew, nboot))
simSkewSmth = replicate(sim, smoothed(dataSkew, nboot))
simSkewWild = replicate(sim, wildstrap(dataSkew, nboot))

simWideBoot = replicate(sim, bootstrap(dataWide, nboot))
simWidePara = replicate(sim, parametric(dataWide, nboot))
simWideRsam = replicate(sim, resampling(dataWide, nboot))
simWideSmth = replicate(sim, smoothed(dataWide, nboot))
simWideWild = replicate(sim, wildstrap(dataWide, nboot))


write.csv(simNormBoot, file = 'simNormBoot.csv', row.names = F)
write.csv(simNormPara, file = 'simNormPara.csv', row.names = F)
write.csv(simNormRsam, file = 'simNormRsam.csv', row.names = F)
write.csv(simNormSmth, file = 'simNormSmth.csv', row.names = F)
write.csv(simNormWild, file = 'simNormWild.csv', row.names = F)

write.csv(simTailBoot, file = 'simTailBoot.csv', row.names = F)
write.csv(simTailPara, file = 'simTailPara.csv', row.names = F)
write.csv(simTailRsam, file = 'simTailRsam.csv', row.names = F)
write.csv(simTailSmth, file = 'simTailSmth.csv', row.names = F)
write.csv(simTailWild, file = 'simTailWild.csv', row.names = F)

write.csv(simSkewBoot, file = 'simSkewBoot.csv', row.names = F)
write.csv(simSkewPara, file = 'simSkewPara.csv', row.names = F)
write.csv(simSkewRsam, file = 'simSkewRsam.csv', row.names = F)
write.csv(simSkewSmth, file = 'simSkewSmth.csv', row.names = F)
write.csv(simSkewWild, file = 'simSkewWild.csv', row.names = F)

write.csv(simWideBoot, file = 'simWideBoot.csv', row.names = F)
write.csv(simWidePara, file = 'simWidePara.csv', row.names = F)
write.csv(simWideRsam, file = 'simWideRsam.csv', row.names = F)
write.csv(simWideSmth, file = 'simWideSmth.csv', row.names = F)
write.csv(simWideWild, file = 'simWideWild.csv', row.names = F)