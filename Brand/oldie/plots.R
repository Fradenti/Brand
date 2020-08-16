library(mvnfast)
library(mcclust)
library(mcclust.ext)
library(fda)
library(MCMCpack)
library(patchwork)
library(gganimate)
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


X <- cbind(X,
           cos(x)-2+rnorm(t,0,.25),
           cos(x-2.5)-2+rnorm(t,0,.25))
matplot(X,type="l")

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

# matplot(apply(RES$FTr,c(1,2),mean),type="l",lty=2)
# matplot(RES$FTr0[,1:3],type="l",add=T)

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
matplot(X,type="l")

matplot(-cos(x)+rnorm(t,0,.25), add=T,type="l")


mX <- reshape2::melt(X)
mY <- reshape2::melt(Y)

mX <- cbind(mX, c=c(rep(1:3,c(1000,2000,3000)),rep(1,100),rep(2,100) ))
mY <- cbind(mY, c=rep(1:5,rep(2000,5)) )


nrow(mY)
nrow(mX)
matplot(X,type="l")


py <- ggplot(mY)+geom_line(aes(x=Var1,y=value,col=factor(c),group=Var2))+theme_bw()
py


px <- ggplot(mX %>% filter(c==2))+geom_line(aes(x=Var1,y=value,
  col=factor(c),group=Var2))+theme_bw()
px



# gettin serious ----------------------------------------------------------

Tr1 <- reshape2::melt(RES$FTr0)
Tr2 <- reshape2::melt(RES$FTr0+3*RES$STr0)
Tr3 <- reshape2::melt(RES$FTr0-3*RES$STr0)
Tr4 <- cbind(a=Tr2,b=Tr3)


minX <- min(mX)
maxX <- max(mX)
px0 <- ggplot(mX)+geom_line(aes(x=Var1,y=value,
  col=factor(c),group=Var2),alpha=.25,lwd=.75)+theme_bw()+
  xlab("x")+ylab("y")+theme(legend.position="none")+
  ggtitle("Training Set")+#ylim(-3.5,1.5)+
  scale_color_manual(values=c(1:3))+
  scale_y_continuous(breaks=c(-3:3))
px0


px0Y <- ggplot(mY)+geom_line(aes(x=Var1,y=value,
                                col=factor(c),group=Var2),alpha=.25,lwd=.75)+theme_bw()+
  xlab("x")+ylab("y")+theme(legend.position="none")+ggtitle("Training Set - Data")+
  scale_y_continuous(breaks=c(-3:3))+
  scale_color_manual(values=rep("royalblue",5))+ggtitle("Test Set")#+ylim(-3.2,1.5)
px0Y

px0+px0Y
px0/px0Y
ggsave("TestAndTraining.png",width = 8,height = 6)
ggsave("TestAndTraining.pdf",width = 8,height = 6)


px1a <- ggplot(mX)+geom_line(aes(x=Var1,y=value,
                                col=factor(c),group=Var2),alpha=.25,lwd=.75)+theme_bw()+
  xlab("x")+ylab("y")+theme(legend.position="none")+
  scale_y_continuous(breaks=c(-3:3))+
  ggtitle("Training Set")+#ylim(-1.5,1.5)+
  scale_color_manual(values=c(1:3))
px1a

px1b <- ggplot()+geom_line(data=Tr1,
  aes(x=Var1, y=value, group=Var2,col=factor(Var2)))+theme_bw()+
  scale_y_continuous(breaks=c(-3:3))+ylim(-3,1.8)+
  geom_ribbon(data=Tr4,aes(x=a.Var1,ymin=a.value,ymax=b.value,group=a.Var2,fill=a.Var2),alpha=.3)+geom_line(data=Tr1,
    aes(x=Var1, y=value, group=Var2,col=factor(Var2)),lwd=1.2)+
  xlab("x")+ylab("y")+theme(legend.position="none")+
  scale_color_manual(values=c(1:3))+
  scale_fill_manual(values=c(1:3))+#ylim(-3.2,1.8)+
  ggtitle("Training Set - Robust Mean and Standard Deviation")#+ylim(-4,1.5)

px1c <- ggplot()+theme_void()+  scale_y_continuous(breaks=c(-3:3))


px1a/px1c
ggsave("OnlyTrain.png",width = 8,height = 6)
ggsave("OnlyTrain.pdf",width = 8,height = 6)

(px1a)/px1b
ggsave("output_gif/TrainExtract.png",width = 8,height = 6)
ggsave("output_gif/TrainExtract.pdf",width = 8,height = 6)


v <- apply(RES$AB[,2,]==0,1,function(x) mean(x==0))
mY2  <- cbind(mY,b=rep(v,rep(100,100)))
px1Y <- ggplot(mY2)+geom_line(aes(x=Var1,y=value,
                                 col=factor(c),group=Var2),alpha=.25,lwd=.75)+theme_bw()+
  scale_y_continuous(breaks=c(-3:3))+
  scale_color_manual(values=rep("royalblue",5))+ggtitle("Test Set")+
  xlab("x")+ylab("y")+theme(legend.position="none")+ggtitle("Training Set - Data")+
  geom_line(aes(x=Var1,y=value,
                group=Var2),col=(ifelse(mY2$b==0,"darkblue","red")),alpha=.25,lwd=.75)+
  gganimate::transition_layers(from_blank = F)+gganimate::enter_fade()+gganimate::exit_fade()
  #+ylim(-3.2,1.5)
px1Y
gganimate::anim_save(animation = px1Y,filename = "/Users/frali/Documents/GitHub/Brand-package/Brand/output_gif/TestProb.gif")


v <- apply(RES$AB[,2,]==0,1,function(x) mean(x==0))
mY2  <- cbind(mY,b=rep(v,rep(100,100)))
px2Y <- ggplot(mY2)+geom_line(aes(x=Var1,y=value,
  col=factor(c),group=Var2),alpha=0,lwd=.75)+theme_bw()+
  scale_y_continuous(breaks=c(-3:3))+
  scale_color_manual(values=rep("transparent",5))+
  ggtitle("Test Set")+
  xlab("x")+ylab("y")+theme(legend.position="none")+ggtitle("Training Set - Data")+
  geom_line(aes(x=Var1,y=value,
    group=Var2),col="royalblue",alpha=.25,lwd=.75)+
  geom_line(aes(x=Var1,y=value,
    group=Var2),col=(ifelse(mY2$b==0,"darkblue","red")),alpha=.25,lwd=.75)+
  geom_line(aes(x=Var1,y=value,
    group=Var2),col=(ifelse(mY2$b==0,"transparent","red")),alpha=(ifelse(mY2$b==0,0,.25)),lwd=.75)+
  geom_line(aes(x=Var1,y=value,
    group=Var2),col=(ifelse(mY2$c==4,"darkorange",ifelse(mY2$c==5,"purple","transparent"))),alpha=(ifelse(mY2$b==0,0,.25)),lwd=.75)+
  gganimate::transition_layers(from_blank = F,
                               keep_layers = c(Inf,0,0,0,Inf),
                               layer_length = 20,
                               transition_length = 5)+
  gganimate::enter_fade()+gganimate::exit_fade()

px2Y
gganimate::anim_save(animation = px2Y,filename = "TestProb2.gif")


  plot(3,col="darkorange")
# ggplot()+geom_line(data=mX,aes(x=Var1,y=value,
#   col=factor(c),group=Var2),alpha=.1,lwd=1)+
#
#   geom_line(data=Tr1,
#   aes(x=Var1, y=value, group=Var2,col=factor(Var2)))+theme_bw()+
#   geom_ribbon(data=Tr4,aes(x=a.Var1,ymin=a.value,ymax=b.value,group=a.Var2,fill=a.Var2),alpha=.5)+geom_line(data=Tr1,
#     aes(x=Var1, y=value, group=Var2,col=factor(Var2)),lwd=1.2)+
#   xlab("x")+ylab("y")+theme(legend.position="none")+
#   xlim(minX,maxX)
#




#
# matplot(avg.Fte1,type="l",add=T)
# matplot(avg.Fte2,type="l",add=T)
# matplot(avg.Fte1-3*avg.Ste1,type="l",add=T)
# matplot(avg.Fte2-3*avg.Ste2,type="l",add=T)
# matplot(avg.Fte1+3*avg.Ste1,type="l",add=T)
# matplot(avg.Fte2+3*avg.Ste2,type="l",add=T)
#


















#
#
#
#
#
#
#
#
#
#
#
#
#
# # -------------------------------------------------------------------------
#
# y1 <- mvtnorm::rmvnorm(500,c(4,4),sigma = matrix(c(1,.75,.75,1),2,2))
# y2 <- mvtnorm::rmvnorm(500,c(0,0),sigma = matrix(c(1,-.75,-.75,1),2,2))
# y3 <- mvtnorm::rmvnorm(500,c(-4,4),sigma = matrix(c(1,0,0,1),2,2))
# x1 <- mvtnorm::rmvnorm(500,c(4,4),sigma = matrix(c(1,.75,.75,1),2,2))
# x2 <- mvtnorm::rmvnorm(500,c(0,0),sigma = matrix(c(1,-.75,-.75,1),2,2))
# x4 <- mvtnorm::rmvnorm(500,c(-4,4),sigma = matrix(c(1,0,0,1),2,2))
# x4 <- mvtnorm::rmvnorm(500,c(-4,-4),sigma = matrix(c(1,0,0,1),2,2))
# x5 <- mvtnorm::rmvnorm(500,c(4,-4),sigma = matrix(c(1,.5,.5,1),2,2))
#
# lab <- rep(1:8,rep(500,8))
# l2  <- rep(1:2,c(1500,2500))
#
# X <- as_tibble(cbind(rbind(y1,y2,y3,x1,x2,x3,x4,x5),lab,l2))
# X <- X %>% mutate(l2=ifelse(l2==1,"Training set","Test set"))
# X <- X %>% mutate(l2=factor(l2,levels = c("Training set","Test set")))
#
#
# ggplot(data = X)+theme_bw()+xlab("x")+ylab("y")+
#   geom_point(aes(x=V1,y=V2))+facet_wrap(~factor(l2)+


