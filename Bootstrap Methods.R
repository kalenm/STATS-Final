library(doParallel)

set.seed(87539319)
CL = detectCores() - 2
registerDoParallel(cores = CL)

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
esp4 = esp4-mean(esp4)
yWide = b0 + b1 * x+esp4


data<-data.frame(x,yNorm)
start1<-Sys.time()
boot_b <- foreach(i=1:n, .combine=c) %dopar% {
  bootstrap_data<-data[sample(nrow(data),nrow(data),replace=T),]
  unname(lm(yNorm~x,bootstrap_data)$coef[2])
}
end1<-Sys.time()
boot_b<-numeric()
start2<-Sys.time()
for(i in 1:n){
  bootstrap_data<-data[sample(nrow(data),nrow(data),replace=T),]
  boot_b[i]<-lm(yNorm~x,bootstrap_data)$coef[2]
}
end2<-Sys.time()
start1-end1
start2-end2
as.numeric(start1-end1)/as.numeric(start2-end2)

foreach(i=4:1, .combine='c', .inorder=FALSE) %dopar% {
  Sys.sleep(3 * i)
  i
}
