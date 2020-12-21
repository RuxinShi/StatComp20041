## ----random, echo=TRUE--------------------------------------------------------
A=rnorm(1000,mean=5,sd=2)
hist(A)  

## ----table, echo=TRUE---------------------------------------------------------
T=matrix(0,nrow=2,ncol=2)
colnames(T)=c("mean","sd")
rownames(T)=c("expect","real")
T[1,]=c(5,2)
T[2,]=c(mean(A),sd(A))
knitr::kable(T,digits=c(3,3),align="c")

## -----------------------------------------------------------------------------
n=3000#The number of X we choose to generate.
U=runif(n,min=0,max=pi/3)
X=pi/3*sin(U)#get X
E=mean(X)#the mean of X, which is also the value we want
realE=0.5#the real value
rbind(c("right","estimate","difference"),c(realE,E,abs(realE-E)))#compare

## -----------------------------------------------------------------------------
set.seed(726)
n=3000#The number of X we choose to generate.
U=runif(n,min=0,max=1)
X=exp(U)#get X
E1=mean(X)#the mean of X, which is also the value we want
realE=exp(1)-1#the real value
rbind(c("right","estimate","difference"),c(realE,E1,abs(realE-E1)))#compare

## -----------------------------------------------------------------------------
n=1500#The number of X we choose to generate.
set.seed(112)
U=runif(n,min=0,max=1)
X1=exp(U)#get X
X2=exp(1-U)
X=(X1+X2)/2
E2=mean(X)#the mean of X, which is also the value we want
realE=exp(1)-1#the real value
rbind(c("right","estimate","difference"),c(realE,E2,abs(realE-E2)))#compare

## -----------------------------------------------------------------------------
n=3000
V=NULL
for(i in 1:1000){
  Var2=Var1=0
  U=runif(n,0,1)
  X1=exp(U)#simple Monte Carlo method
  X2=(exp(U[1:(n/2)])+exp(1-U[1:(n/2)]))/2#antithetic variate approach
  Var1=var(X1)
  Var2=var(X2)
  V[i]=Var2/Var1
}
mean(V)

## ----echo=TRUE----------------------------------------------------------------
n=3000#numbers of variables
g=function(x){
  x^0.5*exp(-x/2)/2/(2*pi)^0.5*(x>1)
}
#f1=e^(-x/2)
##cdf F1=1-e^(-(x-1)/2)
u=runif(n)
x=-2*log(1-u)+1
fg=g(x)/exp(-x/2)
E1=mean(fg)
SD1=sd(fg)
#f2=1/x^2
##cdf F2=1-1/x
u=runif(n)
x=1/(1-u)
fg=g(x)*x^2
E2=mean(fg)
SD2=sd(fg)
#results
cbind("E"=c(E1,E2),"SD"=c(SD1,SD2))

## ----echo=FALSE---------------------------------------------------------------
gf2=function(x){
  x^0.5*exp(-x/2)/2/(2*pi)^0.5*x^2
}
gf1=function(x){
  x^0.5/2/(2*pi)^0.5
}
curve(gf1(x),ylim=c(0,1),xlim=c(1,6),col="BLUE")
curve(gf2(x),ylim=c(0,1),xlim=c(1,6),add=TRUE)

## -----------------------------------------------------------------------------
n=3000
k=5
m=n/k#replicates per stratum
N=100#number of times to repeat the estimation
g=function(x){
  exp(-x)/(1+x^2)*(x>0)*(x<1)
}
E1=NULL
E2=NULL
for(i in 1:N){
  U=runif(n)
  x=g(U)
  E1[i]=mean(x)
  E=NULL
  for(j in 1:k){
    u=runif(m,(j-1)/k,j/k)
    x=g(u)
    E[j]=mean(x)
  }
  E2[i]=mean(E)
}
cbind("mean"=c(mean(E1),mean(E2)),"sd"=c(sd(E1),sd(E2)))

## -----------------------------------------------------------------------------
miu=2
s=1
n=10
alpha=0.05
X=exp(rnorm(n,miu,s))
miu_hat=mean(log(X))
S_hat=sum((log(X)-miu_hat)^2)/(n-1)
CI=c(miu_hat-S_hat*qt(1-alpha/2,n-1)/n^0.5,miu_hat+S_hat*qt(1-alpha/2,n-1)/n^0.5)
c((CI[1]),(CI[2]))#95%CI

## -----------------------------------------------------------------------------
m=3000#repeat times
CI=matrix(0,nrow=m,ncol=2)
q=0
for(i in 1:m){
  X=exp(rnorm(n,miu,s))
  miu_hat=mean(log(X))
  S_hat=(sum((log(X)-miu_hat)^2)/(n-1))^0.5
  CI[i,]=c(miu_hat-S_hat*qt(1-alpha/2,n-1)/n^0.5,miu_hat+S_hat*qt(1-alpha/2,n-1)/n^0.5)
  if(CI[i,1]<=miu&&CI[i,2]>=miu)q=q+1
}
q/m

## -----------------------------------------------------------------------------
n=20
alpha=0.05
m=1000
q=p=0
U1=matrix(0,nrow=m,ncol=2)
U2=matrix(0,nrow=m,ncol=2)
for(i in 1:m){
  X1=rchisq(n,df=2)
  S=(sum((X1-mean(X1))^2)/(n-1))^0.5
  U1[i,]=c(mean(X1)-S*qt(1-alpha/2,n-1)/n^0.5,mean(X1)+S*qt(1-alpha/2,n-1)/n^0.5)
  U2[i,]=c(mean(X1)-S*qnorm(1-alpha/2)/n^0.5,mean(X1)+S*qnorm(1-alpha/2)/n^0.5)
  if(U1[i,1]<=2&&U1[i,2]>=2)q=q+1
  if(U2[i,1]<=2&&U2[i,2]>=2)p=p+1
}
q/m#t-interval
p/m#interval for variance

## -----------------------------------------------------------------------------
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    a <- sample(c(1, 10), replace = TRUE,
                    size = n, prob = c(1-e, e))
    x <- rbeta(n, a, a)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    v <- sample(c(5, 10), replace = TRUE,
                    size = n, prob = c(1-e, e))
    x <- rt(n, v)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
Ftest=function(x,y,n){
  X <- x - mean(x)
  Y <- y - mean(y)
  S1=sum(X^2)/(length(x)-1)
  S2=sum(Y^2)/(length(y)-1)
  return(as.integer(S1/S2>=qf(1-0.055,n-1,n-1)))
}
sigma1 <- 1
sigma2 <- 1.5
m=10000
n=c(20,50,200)
power=matrix(0,nrow=3,ncol=2)
for(i in 1:3){
  power[i,] <- rowMeans(replicate(m, expr={
    x <- rnorm(n[i], 0, sigma1)
    y <- rnorm(n[i], 0, sigma2)
    c(count5test(x, y),Ftest(x,y,n[i]))
  }))
}
power

## -----------------------------------------------------------------------------
n=c(5,10,50)
cv=qchisq(0.975,3*4*5/6)
sk <- function(x,n) {
  #computes the sample skewness coeff.
  xbar=c(mean(x[,1]),mean(x[,2]),mean(x[,3]))
  E=solve(cov(x,x))
  m=0
  y=x-xbar
  for(l1 in 1:(n-1)){
    for(l2 in (l1+1):n){
      m=m+2*(t(y[l1,])%*%E%*%y[l2,])^3
    }
    m=m+(t(y[l1,])%*%E%*%y[l1,])^3
  }
  m=m+(t(y[n,])%*%E%*%y[n,])^3
  return(m/n^2)
}
p.reject <- numeric(length(n))#to store sim. results
m <- 10000 #num. repl. each sim.
for (i in 1:length(n)) {
  sktests <- numeric(m) #test decisions
  for (j in 1:m) {
    x <- cbind(rnorm(n[i]),rnorm(n[i]),rnorm(n[i]))
    #test decision is 1 (reject) or 0
    sktests[j] <- as.integer(abs(n[i]*sk(x,n[i])/6) >= cv )
  }
  p.reject[i] <- mean(sktests) #proportion rejected  }
}
p.reject

## -----------------------------------------------------------------------------
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
   sigma <- sample(c(1, 10), replace =TRUE,size = n, prob = c(1-e, e))
   x <- rnorm(n, 0, sigma)
   sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
p1=0.676
p2=0.651
pc=(p1+p2)/2
n=10000
1-pnorm(p1-p2,0,(pc*(1-pc)*(2/n))^0.5)

## -----------------------------------------------------------------------------
library(bootstrap)
n=length(law[,1])#the length of data
theta.hat=cor(law$LSAT,law$GPA)#original data correlation
theta.jack=numeric(n)#set to store statistics from jackknife
for(i in 1:n){
  theta.jack[i] <- cor(law$LSAT[-i],law$GPA[-i])#leave the ith data out
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)#compute the bias
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))#compute the standard error
round(c(original=theta.hat,bias.jack=bias.jack,
        se.jack=se.jack),3)#to show results

## -----------------------------------------------------------------------------
library(boot)
m=1e3
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
boot.mean <- function(x,i){
  mean(x[i])
}
dat=c(3,5,7,18,43,85,91,98,100,130,230,487)
for(i in 1:m){
  de <- boot(dat,statistic=boot.mean, R = 1000)
  ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
  ci.bca[i,]<-ci$bca[4:5]
}
cat('norm =',c(mean(ci.norm[,1]),mean(ci.norm[,2])),
    'basic =',c(mean(ci.basic[,1]),mean(ci.basic[,2])),
    'perc =',c(mean(ci.perc[,1]),mean(ci.perc[,2])),
    'BCa =',c(mean(ci.bca[,1]),mean(ci.bca[,2]))
)

## -----------------------------------------------------------------------------
library(bootstrap)
n=length(scor[,1])
lambda_hat=eigen(cov(scor))$values
theta_hat=lambda_hat[1]/sum(lambda_hat)
theta_jack=rep(0,n)
for (i in 1:n) {
 lambda=eigen(cov(scor[-i,]))$values
 theta_jack[i]=lambda[1]/sum(lambda)
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)#compute the bias
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))#compute the standard error
round(c(original=theta.hat,bias.jack=bias.jack,
        se.jack=se.jack),3)#to show results

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)/2)
p=0
for (i in 1:(n-1))
  for (j in (i+1):n) {
    p=p+1
    y=magnetic[-c(i,j)]
    x=chemical[-c(i,j)]
    
    P1=lm(y~x)
    y1_1=chemical[i]*P1$coef[2] + P1$coef[1]
    y1_2=chemical[j]*P1$coef[2] + P1$coef[1]
    e1[p]=(magnetic[i]-y1_1)^2+(magnetic[j]-y1_2)^2
    
    P2=lm(y~x+I(x^2))
    y2_1=P2$coef[1] + P2$coef[2] * chemical[i] + P2$coef[3] * chemical[i]^2
    y2_2=P2$coef[1] + P2$coef[2] * chemical[j] + P2$coef[3] * chemical[j]^2
    e2[p]=(magnetic[i]-y2_1)^2+(magnetic[j]-y2_2)^2
    
    P3=lm(log(y)~x)
    y3_1=exp(P3$coef[1] + P3$coef[2] * chemical[i])
    y3_2=exp(P3$coef[1] + P3$coef[2] * chemical[j])
    e3[p]=(magnetic[i]-y3_1)^2+(magnetic[j]-y3_2)^2
    
    P4=lm(log(y)~log(x))
    y4_1=exp(P4$coef[1] + P4$coef[2] * log(chemical[i]))
    y4_2=exp(P4$coef[1] + P4$coef[2] * log(chemical[j]))
    e4[p]=(magnetic[i]-y4_1)^2+(magnetic[j]-y4_2)^2
  }
e = c(mean(e1)/2,mean(e2)/2,mean(e3)/2,mean(e4)/2)
matrix(e, nrow=1,dimnames=list("prediction error", c("Linear","Quadratic"," Exponential","Log-Log")))
detach(ironslag)

## -----------------------------------------------------------------------------
count5=function(x1,x2){
  x1=x1-mean(x1)
  x2=x2-mean(x2)
  out1 <- sum(x1 > max(x2)) + sum(x1 < min(x2))
  out2 <- sum(x2 > max(x1)) + sum(x2 < min(x1))
  return(max(c(out1, out2)))
}
count5_per=function(z){
  n=length(z)
  x=z[1:(n/2)]
  y=z[-(1:(n/2))]
  a=count5(x,y)
  return(a)
}
per=function(z,R){
  n=length(z)
  result=numeric(R)
  for (i in 1:R){
    q=sample(1:n,n,replace=FALSE)
    result[i]=count5_per(z[q])
  }
  return(sum(result>5)/R)
}
n1=20
n2=40
mu1=mu2=0
s1=s2=2
R=1e3
test1=test2=numeric(R)
#compare
for(i in 1:R){
  x=rnorm(n1,mu1,s1)
  y=rnorm(n2,mu2,s2)
  test1[i]=(count5(x,y)>5)
  z=c(x,y)
  test2[i]=per(z,500)
}
c(mean(test1),mean(test2<0.05))

## ----warning=FALSE------------------------------------------------------------
library(RANN)
library(boot)
library(Ball)
library(energy)
library(MASS)
#NN-methods
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1)
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
NN=function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R, sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  return(c(statistic=ts[1],p.value=p.value))
}
#energy
# eqdist.etest in package energy
#ball
# bd.test in package Ball
#settings
n1=n2=20
n=n1+n2 
N=c(n1,n2)
k=3
R=500
m=100
alpha=0.05
#compare
com=function(mu1,mu2,s1,s2,alpha){
  p.values=matrix(NA,m,3)
  for(i in 1:m){
    data1=mvrnorm(n1,mu1,sigma1)
    data2=mvrnorm(n2,mu2,sigma2)
    data=rbind(data1,data2)
    p.values[i,1]=NN(data,N,k)[2]
    p.values[i,2]=eqdist.etest(data,sizes=N,R=R)$p.value
    p.values[i,3]=bd.test(x=data1,y=data2,num.permutations=R,seed=i*2846)$p.value
  }
  p=colMeans(p.values<alpha)
  return(p)
}

#Unequal variances and equal expectations  
mu1=mu2=c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(2,0,0,0,3,0,0,0,4),nrow=3,ncol=3)
result1=com(mu1,mu2,s1,s2,alpha)

#Unequal variances and unequal expectations  
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,1)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
result2=com(mu1,mu2,s1,s2,alpha)

#Non-normal distributions
#t distribution with 1 df (heavy-tailed distribution)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  data1=rt(n1,1,2)
  data2=rt(n2,2,5)
  data=c(data1,data2)
  p.values[i,1]=NN(data,N,k)[2]
  p.values[i,2]=eqdist.etest(data,sizes=N,R=R)$p.value
  p.values[i,3]=bd.test(x=data1,y=data2,num.permutations=R,seed=i*2846)$p.value
}
result3=colMeans(p.values<alpha)
#bimodel distribution (mixture of two normal distributions)  
bimodel<-function(n,mu1,mu2,sd1,sd2){
  index=sample(1:2,n,replace=TRUE)
  x=numeric(n)
  index1<-which(index==1)
  x[index1]<-rnorm(length(index1),mu1,sd1)
  index2<-which(index==2)
  x[index2]<-rnorm(length(index2),mu2,sd2)
  return(x)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
  data1=bimodel(n1,0,0,1,2)
  data2=bimodel(n2,2,2,2,4)
  data=c(data1,data2)
  p.values[i,1]=NN(data,N,k)[2]
  p.values[i,2]=eqdist.etest(data,sizes=N,R=R)$p.value
  p.values[i,3]=bd.test(x=data1,y=data2,num.permutations=R,seed=i*2846)$p.value
}
result4=colMeans(p.values<alpha)

#Unbalanced samples (say, 1 case versus 10 controls)
n2=10*n1
n=n1+n2 
N=c(n1,n2)
result5=com(mu1,mu2,s1,s2,alpha)
rbind(result1,result2,result3,result4,result5)

## -----------------------------------------------------------------------------
#standard Laplace dendity distribution 
f=function(x) 0.5*exp(-abs(x))
#random walk
Metro=function(sd,x0,N){
  x=numeric(N)
  x[1]=x0
  u=runif(N)
  p=0
  for (i in 2:N) {
    y=rnorm(1,x[i-1],sd)
    if(u[i]<=(f(y)/f(x[i-1]))) x[i]=y 
    else {
      x[i]=x[i-1]
      p=p+1
    }
  }
  return(c(x,p))
}
#compare the chains generated with different variances
N=2000
sd=c(0.1,1,5,10)
x0=10
ch1=Metro(sd[1],x0,N)
ch2=Metro(sd[2],x0,N)
ch3=Metro(sd[3],x0,N)
ch4=Metro(sd[4],x0,N)
#Compare the chains
par(mfrow=c(1,2))
ch=cbind(ch1[1:N],ch2[1:N],ch3[1:N],ch4[1:N])
for (j in 1:4) {
  plot(ch[,j], type="l")
  boxplot(ch[,j])
}
#acceptance rates
R=c(ch1[N+1],ch2[N+1],ch3[N+1],ch4[N+1])
Acc=(N-R)/N
Acc

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

k=4 #number of chains to generate
b=1000 #burn-in length
N=15000
#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)
sd=c(0.1,1,5,10)
par(mfrow=c(1,1))
for(z in 1:4){
 X=matrix(nrow=k,ncol=N)
 for (i in 1:k){
   X[i,]=Metro(sd[z],x0[i],N)[1:N]
 }
 psi=t(apply(X, 1, cumsum))
 for (i in 1:nrow(psi))
   psi[i,]=psi[i,]/(1:ncol(psi))
 rhat=rep(0, N)
 for (j in (b+1):N)
   rhat[j]=Gelman.Rubin(psi[,1:j])
 plot(rhat[(b+1):N], type="l", xlab=paste("sigma=",sd[z]), ylab="R_hat", ylim=c(1,max(max(rhat[(b+1):N]),1.4)))
 abline(h=1.2,lty=2)
}

## -----------------------------------------------------------------------------
k=c(4:25,100,500,1000)
S=function(a,k){
  z=sqrt(a^2*k/(k+1-a^2))
  p=pt(z,df=k,lower.tail=FALSE)
  return(p)
}
f=function(a,k){
  S(a,k)-S(a,k-1)
}
solve = function(k){
  out=uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
  return(out$root)
}
n=length(k)
A_k=matrix(0,ncol=n,nrow=2)
for (i in 1:n){
  A_k[2,i]=solve(k[i])
}
A_k[1,]=k
A_k

## ----warning=FALSE------------------------------------------------------------
library(nloptr)
eval_f0 = function(x,x1,nA=444,nB=132,nOO=361,nAB=63) {
  r1=1-sum(x1)
  nAA=nA*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB=nB*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r=1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-(nA-nAA)*log(2*x[1]*r)-(nB-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}
eval_g0 = function(x,x1,nA=444,nB=132,nOO=361,nAB=63){
  return(sum(x)-0.99999)
}
mle=NULL
r=matrix(c(0.5,0.5),1,2)#the beginning value of p and q
j=1
while(sum(abs(r[j,]-r[j-1,]))>1e-7||j<=1){
  res = nloptr( x0=c(0.2,0.25),eval_f=eval_f0,lb = c(0,0),ub = c(1,1),eval_g_ineq = eval_g0, 
                opts=list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-7),x1=r[j,],nA=444,nB=132,nOO=361,nAB=63)
  j=j+1
  r=rbind(r,res$solution)
  mle=c(mle,eval_f0(x=r[j,],x1=r[j-1,]))
}
#results
n=length(r[,1])
p=r[n,1]
q=r[n,2]
list("results"=c(p,q),"mle"=-mle,"each step"=r)

## -----------------------------------------------------------------------------
formulas = list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
attach(mtcars)
#loops
r=list()
for (i in length(formulas)){
  r=c(r,list(lm(formulas[[i]],data=mtcars)))
}
r
#lapply
l=lapply(formulas,function(x) lm(formula=x, data=mtcars))
l

## -----------------------------------------------------------------------------
trials = replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
sapply(trials,function(x) x[["p.value"]])

## -----------------------------------------------------------------------------
la=function(X,FUN,FUN.VALUE){
  re=Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  re=unlist(re)
  return(re)
}

## ----eval=FALSE---------------------------------------------------------------
#  
#  #include <Rcpp.h>
#  #include <cmath>
#  using namespace Rcpp;
#  // [[Rcpp::export]]
#  //standard Laplace dendity distribution
#  float f(float x){
#    a=0.5*exp(-abs(x));
#    return a;
#  }
#  //random walk
#  NumericVector Metro(float sd,float x0,int N){
#     int x[N+1];
#      x[0]=x0;
#      float u[]=runif(N);
#      p=0;
#      for (int i=1; i<N;i++) {
#        float y[]=rnorm(1,x[i-1],sd);
#        if(u[i]<=(f(y)/f(x[i-1]))){
#          x[i]=y[0];
#        }
#        else {
#            x[i]=x[i-1];
#            p=p+1;
#          }
#      }
#      x[N]=p;
#      return(x);
#    }

## ----warning=FALSE------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
#C
dir_cpp='C:/Users/Shiru/Desktop/ustc/final/Ruxin/src/'
sourceCpp(paste0(dir_cpp,"CMetro.cpp"))
#R
#standard Laplace dendity distribution 
f=function(x) 0.5*exp(-abs(x))
#random walk
Metro=function(sd,x0,N){
  x=numeric(N)
  x[1]=x0
  u=runif(N)
  for (i in 2:N) {
    y=rnorm(1,x[i-1],sd)
    if(u[i]<=(f(y)/f(x[i-1]))) x[i]=y 
    else {
      x[i]=x[i-1]
    }
  }
  return(x)
}
x0=10
N=3000
sd=1
R=Metro(sd,x0,N)[-(1:500)]
C=CMetro(sd,x0,N)[-(1:500)]
qqplot(R,C)

## ----warning=FALSE------------------------------------------------------------
x0=10
N=3000
sd=1
(time=microbenchmark(R=Metro(sd,x0,N),C=CMetro(sd,x0,N)))

