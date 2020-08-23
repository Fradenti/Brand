library(ggplot2)
library(gganimate)
set.seed(1234)
X <- mvtnorm::rmvnorm(500,c(0,0),sigma=matrix(c(1,.7,.7,1),2,2))
L <- matrix(runif(100,-5,3),50,2)
Y <- rbind(X,L)


IL <- sort(apply(L,1,function(x) sum((colMeans(Y)-x)^2)),index=T,decreasing = T)

  ggplot()+
    theme_bw()+xlab("x")+ylab("y")+
    geom_point(aes(x=Y[,1],y=Y[,2]),col="royalblue",alpha=.5)+
    geom_point(aes(x=mean(Y[,1]),mean(Y[,2])),size=5,col="darkblue")+
    stat_ellipse(aes(x=Y[,1],y=Y[,2]),col="darkblue")
  ggsave("/home/fra/Documents/GitHub/adaptdaBNP/output/MBC2_presentation/figures/MCD1.pdf",height = 3,width = 6)



  ggplot()+
    theme_bw()+xlab("x")+ylab("y")+
    geom_point(aes(x=Y[,1],y=Y[,2]),col="royalblue",alpha=.5)+
    geom_point(aes(x=mean(Y[,1]),mean(Y[,2])),size=5,col="darkblue")+
    stat_ellipse(aes(x=Y[,1],y=Y[,2]),col="darkblue")+
    geom_point(aes(x=L[IL$ix[1:30],1],L[IL$ix[1:30],2]),pch="X",col="red",size=4)
    ggsave("/home/fra/Documents/GitHub/adaptdaBNP/output/MBC2_presentation/figures/MCD2.pdf",height = 3,width = 6)


    ggplot()+
      theme_bw()+xlab("x")+ylab("y")+
      geom_point(aes(x=Y[,1],y=Y[,2]),col="royalblue",alpha=.5)+
      geom_point(aes(x=mean(Y[,1]),mean(Y[,2])),size=5,col="darkblue")+
      stat_ellipse(aes(x=Y[,1],y=Y[,2]),col="darkblue")+
      geom_point(aes(x=L[IL$ix[1:30],1],L[IL$ix[1:30],2]),pch="X",col="red",size=4)+
    #geom_point(aes(x=Y[,1],y=Y[,2]),col="red",alpha=.5)+
    geom_point(aes(x=mean(X[,1]),mean(X[,2])),size=5,col="red")+
    stat_ellipse(aes(x=X[,1],y=X[,2]),col="darkred")
    ggsave("/home/fra/Documents/GitHub/adaptdaBNP/output/MBC2_presentation/figures/MCD3.pdf",height = 3,width = 6)
