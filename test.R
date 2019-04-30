#Norm(mu,sigma), theta=(mu,sigma)^tau,mu1~Norm(0,3^2),sigma1~gamma(shape=3,scale=2)
#xi~Norm(c(0,0),diag(2)),mu <- mu1+x1[1];sigma <- sigma1+xi[2]
library(matrixcalc) #vec()
library(mgcv) #rmvn()
n <- 20
m <- 20

mui <- c(0,3)
sigmai <- c(4,0.5)
filename <- paste0("norm","n",n,"m",m,".pdf")
pdf(filename,width = 30,height = 10)
par(mfrow=c(1,2))
xj <- seq(0.5,2.4,by=0.1)
xiij <- matrix(xj,m,2)
normsam(n,m,mui,sigmai,xiij)
xiij <- matrix(1,m,2)
normsam(n,m,mui,sigmai,xiij)
dev.off()
normsam <- function(n,m,mui,sigmai,xiij){
  mu1 <- rnorm(n,mean=mui[1],sd=mui[2])
  sigma1 <- rgamma(n,shape=sigmai[1],rate = sigmai[2])
  xi <- array(NA,dim = c(n,m,2))
  for(j in 1:m)
    xi[,j,] <- rmvn(n,rep(0,2),diag(c(xiij[j,])))
 
  mu1 <- matrix(mu1,n,m)
  sigma1 <- matrix(sigma1,n,m)
  mu <- mu1+xi[,,1]
  sigma <- sigma1+xi[,,2]
  fmu <- vec(mu)
  fsigma <- vec(sigma)
  if(sum(fsigma<0)){
  fmu <- fmu[-which(fsigma<0)]
  fsigma <- fsigma[-which(fsigma<0)]
  }
  N <- length(fsigma)
  mimu <- min(fmu)
  mamu <- max(fmu)
  msigma <- max(fsigma)
  yuper <- max(dnorm(fmu,mean=fmu,sd=sqrt(fsigma)))
  xlow <- mimu-3*sqrt(msigma)
  xup <- mamu+3*sqrt(msigma)
  x <- seq(xlow,xup,by=0.01)
  
  curve(dnorm(x,mean = fmu[1],sd=sqrt(fsigma[1])),xlim=c(xlow,xup),ylim=c(0,yuper),xlab = "x",ylab = "density",main=paste("norm N=",N,"n=",n))
  for (i in 2:N) {
    curve(dnorm(x,mean = fmu[i],sd=sqrt(fsigma[i])),xlim=c(xlow,xup),ylim=c(0,yuper),add=TRUE)
  }
  return(yuper)
}


# mu1 <- rnorm(n,mean=0,sd=3)
# sigma1 <- rgamma(n,shape =3,rate=0.5)
# plot(mu)
# plot(sigma)
# var(sigma1)
# xi <- rmvn(n,rep(0,2),diag(2))
# dim(xi)
# mu <- mu1+xi[,1]
# sigma <- sigma1+xi[,2]
# sum(sigma<0)
# x <- seq(-15,15,by=0.01)
# dev.new()
# 
# curve(dnorm(x,mean = mu[i],sd=sqrt(sigma[i])),xlim=c(-15,15),ylim=c(0,1),xlab = "x",ylab = "density",main=paste("norm n=",n))
# for (i in 2:n) {
#   curve(dnorm(x,mean = mu[i],sd=sqrt(sigma[i])),xlim=c(-15,15),ylim=c(0,1),add=TRUE)
# }

##########################################################################################

n <- 1
m <- 10

mui <- c(6,0.7)
sigmai <- c(5,0.5)
filename <- paste0("gamm","n",n,"m",m,".pdf")
pdf(filename,width = 30,height = 10)
par(mfrow=c(1,2))
xj <- seq(0.5,2.4,by=0.1)
xiij <- matrix(xj,m,2)
gammasam(n,m,mui,sigmai,xiij)
xiij <- matrix(1,m,2)
gammasam(n,m,mui,sigmai,xiij)
dev.off()
gammasam <- function(n,m,mui,sigmai,xiij){
  mu1 <- rgamma(n,shape=mui[1],rate=mui[2])
  sigma1 <- rgamma(n,shape=sigmai[1],rate = sigmai[2])
  xi <- array(NA,dim = c(n,m,2))
  for(j in 1:m)
    xi[,j,] <- rmvn(n,rep(0,2),diag(c(xiij[j,])))
  
  mu1 <- matrix(mu1,n,m)
  sigma1 <- matrix(sigma1,n,m)
  mu <- mu1+xi[,,1]
  sigma <- sigma1+xi[,,2]
  fmu <- vec(mu)
  fsigma <- vec(sigma)
  if(sum(fmu<0)||sum(fsigma<0)){
    tep <- sort(c(which(fmu<0),which(fsigma<0)))
  fmu <- fmu[-tep]
  fsigma <- fsigma[-tep]
  }
  N <- length(fsigma)
  xmode <- (fmu-1)/fsigma
  if(sum(xmode<0)){
    xmode[xmode<0] <- 1
  }
  xdata <- cbind(xmode,fmu,fsigma)
  yuper <- max(dgamma(xmode,shape=fmu,rate=fsigma))
  xuper <- 10#max(xmode)+10#max((1+3/fsigma)*fmu/fsigma)
  x <- seq(0.00001,xuper,by=0.0001)
  curve(dgamma(x,shape = fmu[1],rate=fsigma[1]),xlim=c(0,xuper),ylim=c(0,yuper),xlab = "x",ylab = "density",main=paste("gamma N=",N,"n=",n))
  for (i in 2:N) {
    curve(dgamma(x,shape = fmu[i],rate=fsigma[i]),xlim=c(0,xuper),ylim=c(0,yuper),add=TRUE)
  }
  return(yuper)
}


# ffun <- function(a){
#   x <- dgamma(a[1],shape = a[2],rate=a[3])
#   return(x)
# }
# am <- cbind(xuper,fmu,fsigma)
# yuper <- apply(am,1,  ffun)
# yuper <- which.max(yuper)
# n <- 100
# m <- 10
# mu1 <- rgamma(n,shape = 7,rate=0.7)
# sigma1 <- rgamma(n,shape =5,rate=0.4)
# plot(mu)
# plot(sigma)
# var(sigma1)
# xi <- rmvn(n,rep(0,2),diag(2))
# dim(xi)
# mu <- mu1+xi[,1]
# sigma <- sigma1+xi[,2]
# sum(sigma<0)
# sum(mu<0)
# x <- seq(0,10,by=0.001)
# dev.new()
# curve(dgamma(x,shape = mu[i],rate = (sigma[i])),xlim=c(0,10),ylim=c(0,5),xlab = "x",ylab = "density",main=paste("gamma n=",n))
# for (i in 2:n) {
#   curve(dgamma(x,shape = mu[i],rate=(sigma[i])),xlim=c(0,10),ylim=c(0,5),add=TRUE)
# }

