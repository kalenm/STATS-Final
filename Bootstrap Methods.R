library(doParallel)

set.seed(87539319)
CL = detectCores() - 2
registerDoParallel(cores = cl)

n = 5000
b0 = 1
b1 = 3
x = runif(n)

eps1 = rnorm(n)
yNorm = b0 + b1 * x+eps1

esp2 = rt(n,df=3)/sqrt(3)
yTail = b0 + b1 * x+esp2

esp3 = rgamma(n,shape=4,rate=2)-2
ySkew = b0 + b1 * x+esp3

esp4 = rnorm(n)*x^2/sqrt(mean(x^4))
esp = esp4-mean(esp4)
yWide = b0 + b1 * x+esp4

