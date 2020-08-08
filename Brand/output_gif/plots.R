library(mvnfast)
library(mcclust)
library(mcclust.ext)
library(fda)
library(MCMCpack)
set.seed(12345)

# Rcpp::sourceCpp("code/functional_extension_code/OfficialVersion/01_SLICE_functional_TRAIN_STOCHASTIC_v2.cpp")
# source("code/functional_extension_code/OfficialVersion/01_SLICE_functional_TRAIN_STOCHASTIC_v2.R")

t <- 100
x <- seq(0, 2 * pi, length.out = t)
# Define the basis
n_basis <- 100
basis_order <- 5 # spline
basis <-
  create.bspline.basis(rangeval = range(x),
    nbasis = n_basis,
    norder = basis_order)
basis_evaluated <- eval.basis(evalarg = x, basisobj = basis)
plot(basis)


A <- replicate(10, expr = sin(x) + rnorm(t, sd = .25), simplify = T)
B <- replicate(20, expr = cos(x) + rnorm(t, sd = .25), simplify = T)
C <- replicate(30, expr = rnorm(t, sd = .25), simplify = T)
extract_tr_coef_splines <- Vectorize(function(single_func_set, x=x){
  smooth_representation  =  smooth.basis(argvals = x, y = single_func_set, fdParobj = basis)
  rowMeans(smooth_representation$fd$coefs)
}, vectorize.args = "single_func_set")

beta_coefs_tr <- extract_tr_coef_splines(single_func_set = list(A,B,C),x = x)

#Tr set
X <- cbind(A, B, C)
cl_train <- rep(1:3, c(10, 20, 30))
matplot(X,
  type = "l",
  col = cl_train,
  lty = 1)
J <- length(unique(cl_train))

matplot(basis_evaluated%*%beta_coefs_tr, col=4, add = T,lty = 2, type = "l")

# Te set
n <- 20
A_test <- replicate(n, expr = sin(x) + rnorm(t, sd = .25), simplify = T)
B_test <- replicate(n, expr = cos(x) + rnorm(t, sd = .25), simplify = T)
C_test <- replicate(n, expr = rnorm(t, sd = .25), simplify = T)
E_test <- replicate(n, expr = cos(10*x)-1+rnorm(t, sd = .25), simplify = T)
F_test <- replicate(n, expr = tanh(x^2)-3+rnorm(t, sd = .25), simplify = T)

Y <- cbind(A_test, B_test, C_test, E_test, F_test)
cl_test <- rep(1:5, each=n)

matplot(Y,
  type = "l",
  col = cl_test,
  lty = 1)
sigma.trainA0 <- apply(A,1,var)
sigma.trainB0 <- apply(B,1,var)
sigma.trainC0 <- apply(C,1,var)
sigma.train0 <- cbind(sigma.trainA0,sigma.trainB0,sigma.trainC0)

matplot(sigma.train0)

XX1<-basis_evaluated%*%beta_coefs_tr

A1<-A-XX1[,1]
B1<-B-XX1[,2]
C1<-C-XX1[,3]
sigma.trainA <- apply(A1,1,function(x) sum(x^2)/(length(x)-1) )
sigma.trainB <- apply(B1,1,function(x) sum(x^2)/(length(x)-1) )
sigma.trainC <- apply(C1,1,function(x) sum(x^2)/(length(x)-1) )
sigma.train <- cbind(sigma.trainA,sigma.trainB,sigma.trainC)

matplot(sigma.train[,1],type="l")
matplot(sigma.train0[,1],type="l",add=T)
matplot(sigma.train[,2],type="l")
matplot(sigma.train0[,2],type="l",add=T)
matplot(sigma.train[,3],type="l")
matplot(sigma.train0[,3],type="l",add=T)




prior <- list(
  basis       = basis_evaluated,
  beta_train  = beta_coefs_tr,
  sigma_train = sigma.train,
  aDir     = c(1, 1, 1, .1),
  aDP      = 1,
  a_H      = 3,
  # prior varianza gruppi nuovi
  b_H      = 2,
  a_tau    = 3,
  # varianza dei betini
  b_tau    = 2,
  a_alphaDP = 1,
  b_alphaDP = 1,
  aDP = 1,
  s_tau = 10,
  vg0 = .001,
  KappaG = .001
)


RES <- brand::functional_brand(Y = Y,
  prior = prior,
  L=10, # L numero possibili gruppi
  nsim=200,
  thinning=1,
  burn_in=200,
  verbose=1,fixed_alphaDP = T,kappa = .5,
  learning_type = "Indutive")


major.vote <- function(x){
  as.numeric(names(sort(table(x),decreasing = T))[1])
}

CL <- apply(RES$AB[,1,],1,major.vote)

matplot(apply(RES$FTr,c(1,2),mean),type="l",lty=2)
matplot(RES$FTr0[,1:3],type="l",add=T)

matplot(apply(RES$STr,c(1,2),mean),type="l",lty=2)
matplot(RES$STr0[,1:3],type="l",add=T)



matplot(Y[,CL==0],type="l",col=2)
CL1 <- apply(RES$AB[,2,],1,major.vote)




matplot(Y,type="l",col=CL+1) # yeeeeeeee
newobsB <- CL==0 #this shall always be 0


BET    <- RES$AB[newobsB,2,]
psmBET <- comp.psm(t(BET))
image(psmBET)

clB <- minbinder(psmBET)$cl
clV <- minVI(psmBET)$cl
cl_hat <- CL
cl_hat[CL==0] <- clB+J # I add J to separate from the original tr labels
matplot(Y,type="l",col=cl_hat) #yeeeeeeeeeeeee






plot(RES$USL)
############################################################
plot(ts(RES$aDP))
acf(RES$aDP)
hist(RES$aDP,breaks = 100)
###############################################################
RES$AK[[2000]]

matplot(RES$FTr0,type="l",lty=1)
matplot(RES$FTr0+3*RES$STr0,type="l",add=T)
matplot(RES$FTr0-3*RES$STr0,type="l",add=T)
library(tidyverse)
avg.Fte1 <- map_dfc( RES$FTe,~.x[,1])
avg.Fte1 <- rowMeans(avg.Fte1)
avg.Fte2 <- map_dfc(RES$FTe,~.x[,2])
avg.Fte2 <- rowMeans(avg.Fte2)

avg.Ste1 <- map_dfc( RES$STe,~.x[,1])
avg.Ste1 <- rowMeans(avg.Ste1)

avg.Ste2 <- map_dfc(RES$STe,~.x[,2])
avg.Ste2 <- rowMeans(avg.Ste2)




matplot(Y,type="l",col="white")
matplot(RES$FTr0,type="l")
matplot(RES$FTr0+3*RES$STr0,type="l",add=T)
matplot(RES$FTr0-3*RES$STr0,type="l",add=T)
matplot(avg.Fte1,type="l",add=T)
matplot(avg.Fte2,type="l",add=T)
matplot(avg.Fte1-3*avg.Ste1,type="l",add=T)
matplot(avg.Fte2-3*avg.Ste2,type="l",add=T)
matplot(avg.Fte1+3*avg.Ste1,type="l",add=T)
matplot(avg.Fte2+3*avg.Ste2,type="l",add=T)


dim(X)
dim(Y)
table(mX$Var1)

mX <- reshape2::melt(X)
mY <- reshape2::melt(Y)

mX <- cbind(mX, c=rep(1:3,c(1000,2000,3000)) )
mY <- cbind(mY, c=rep(1:5,rep(2000,5)) )

nrow(mY)
nrow(mX)


py <- ggplot(mY)+geom_line(aes(x=Var1,y=value,col=factor(c),group=Var2))+theme_bw()
py


px <- ggplot(mX %>% filter(c==2))+geom_line(aes(x=Var1,y=value,
  col=factor(c),group=Var2))+theme_bw()
px



ggplot()

Tr1 <- reshape2::melt(RES$FTr0)
Tr2 <- reshape2::melt(RES$FTr0+3*RES$STr0)
Tr3 <- reshape2::melt(RES$FTr0-3*RES$STr0)
Tr4 <- cbind(a=Tr2,b=Tr3)


minX <- min(mX)
maxX <- max(mX)
px0 <- ggplot(mX)+geom_line(aes(x=Var1,y=value,
  col=factor(c),group=Var2),alpha=.5,lwd=1)+theme_bw()+
  xlab("x")+ylab("y")+theme(legend.position="none")+
  xlim(minX,maxX)
px0


px1 <- ggplot()+geom_line(data=Tr1,
  aes(x=Var1, y=value, group=Var2,col=factor(Var2)))+theme_bw()+
  geom_ribbon(data=Tr4,aes(x=a.Var1,ymin=a.value,ymax=b.value,group=a.Var2,fill=a.Var2),alpha=.3)+geom_line(data=Tr1,
    aes(x=Var1, y=value, group=Var2,col=factor(Var2)),lwd=1.2)+
xlab("x")+ylab("y")+theme(legend.position="none")+
  xlim(minX,maxX)

library(patchwork)
px0+px1





ggplot()+geom_line(data=mX,aes(x=Var1,y=value,
  col=factor(c),group=Var2),alpha=.1,lwd=1)+

  geom_line(data=Tr1,
  aes(x=Var1, y=value, group=Var2,col=factor(Var2)))+theme_bw()+
  geom_ribbon(data=Tr4,aes(x=a.Var1,ymin=a.value,ymax=b.value,group=a.Var2,fill=a.Var2),alpha=.5)+geom_line(data=Tr1,
    aes(x=Var1, y=value, group=Var2,col=factor(Var2)),lwd=1.2)+
  xlab("x")+ylab("y")+theme(legend.position="none")+
  xlim(minX,maxX)






matplot(avg.Fte1,type="l",add=T)
matplot(avg.Fte2,type="l",add=T)
matplot(avg.Fte1-3*avg.Ste1,type="l",add=T)
matplot(avg.Fte2-3*avg.Ste2,type="l",add=T)
matplot(avg.Fte1+3*avg.Ste1,type="l",add=T)
matplot(avg.Fte2+3*avg.Ste2,type="l",add=T)



