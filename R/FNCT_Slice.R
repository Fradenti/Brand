

a.AND.b_sigmaT_train <- function(sigma.est.train, # matrix TxG
                                 vg,t,G){
  ab      <- array(NA,c(t,G,2))
  ab[,,1] <- 2 + (sigma.est.train)^2/vg
  ab[,,2] <- ((sigma.est.train)^2/vg+1)*sigma.est.train

  return(ab)
}


extract_robust_tr_v2 <-
  Vectorize(function(single_func_set,
                     x = x, basis_evaluated,
                     h_MCD = .75) {
    smooth_representation  <-   smooth.basis(
      argvals = x,
      y = as.matrix(single_func_set),
      fdParobj = basis
    )
    basis_coef_matrix <- t(smooth_representation$fd$coefs)
    robust_procedure <- rrcov::CovMrcd(x = basis_coef_matrix, alpha = h_MCD)
    # regularized version, suitable for high dimensional data

    effe.train    <- basis_evaluated%*%robust_procedure$center

    dif <- apply(single_func_set[,robust_procedure$best],2,function(x)  x - effe.train)

    sigma2 <- apply(dif, 1, function(k)
      sum(k ^ 2) / (length(k) - 1))



    list(center=effe.train,sigma2=sigma2,basis=robust_procedure$center)
  }, vectorize.args = "single_func_set")


#' Functional Brand
#'
#' @param Y xxx
#' @param prior xxx
#' @param L xxx
#' @param nsim xxx
#' @param thinning xxx
#' @param burn_in xxx
#' @param fixed_alphaDP xxx
#' @param kappa xxx
#' @param verbose xxx
#' @param learning_type xxx
#'
#' @return model
#' @export
#'
Brand_fct <- function(Y,
                      prior,
                      L,
                      # L numero possibili gruppi
                      nsim,
                      thinning,
                      burn_in,
                      fixed_alphaDP = F,
                      kappa = .5,
                      verbose = 1,
                      learning_type = c("transductive", "inductive")) {

  #####################################################################
  ############
  basis       <- prior$basis
  beta.train  <- prior$beta.train
  sigma.train0<- prior$sigma.train # now this is a matrix of functions
  effe.train0 <- basis%*%beta.train
  vg          <- prior$vg
  KappaG      <- prior$KappaG
  n           <- ncol(Y)
  t           <- nrow(Y)
  G           <- ncol(beta.train)
  K           <- nrow(beta.train) # numero basi
  ############
  R           <- rowSums(basis^2)
  alphaD      <- rep(1,G)
  OneOverR    <- 1/R
  # HYPERPRIOR PARAM for STOCHASTIC TRAINING
  A_B_IG_Train_Stochastic <- a.AND.b_sigmaT_train(sigma.est.train = sigma.train0,vg = vg,
                                                  t = t,G = G)
  ###############################################
  # Containers
  EFFE.TEST   <- vector("list",length = nsim)#array( NA,dim = c(t,L,nsim))
  SIGMA.TEST  <- vector("list",length = nsim)
  ###
  PROB        <- matrix(NA,nsim, G+1)
  OMEGA       <- vector("list",length = nsim)
  TAU         <- vector("list",length = nsim)
  AK          <- vector("list",length = nsim)
  ALPHA.BETA  <- array( NA,dim = c(n,2,nsim))
  AlphaDP     <- numeric(nsim)
  USLICE      <- matrix(NA,n,nsim)
  ###
  if (learning_type %in% c("transductive", "training_stochastic")) {
    EFFE.TRAIN  <- array(NA, dim = c(t, G, nsim))
    SIGMA.TRAIN <- array(NA, dim = c(t, G, nsim))
  }
  ###
  ##############################################
  ### Hyperparameters
  aDir     <- prior$aDir
  aDP      <- prior$aDP
  a_H      <- prior$a_H  # prior varianza gruppi nuovi
  b_H      <- prior$b_H
  a_tau    <- prior$a_tau  # varianza dei betini
  b_tau    <- prior$b_tau  # varianza dei betini
  a_alpha  <- prior$a_alphaDP
  b_alpha  <- prior$b_alphaDP
  s_tau    <- prior$s_tau
  # inizialization
  pidir         <- c(MCMCpack::rdirichlet(1,aDir))
  alphabeta     <- matrix(NA,n,2)
  alphabeta[,1] <- sample(0:G,n,replace = T)
  alphabeta[,2] <- ifelse(alphabeta[,1]>0,0,sample(L))
  u             <- stats::rbeta(n = L, shape1 = 1, shape2 = aDP)
  omega         <- StickBreaker_cpp(V = u)
  effe.test     <- matrix(NA,t,L)
  sigma.test    <- matrix(NA,t,L)
  effe.train    <- effe.train0 #matrix(NA,t,G)
  sigma.train   <- sigma.train0 #matrix(NA,t,G)
  tau           <- numeric(L)
  ak            <- numeric(L)
  for(k in 1:L){
    sigma.test[,k] <- MCMCpack::rinvgamma(t,a_H,b_H)
    tau[k]         <- MCMCpack::rinvgamma(1,a_tau,b_tau)
    ak[k]          <- stats::rnorm(1,0,sqrt(s_tau))
    effe.test[,k]  <- stats::rnorm(t,0,sd=sqrt(tau[k]*R))
  }
  ZETA <- alphaneta_to_zeta(alphabeta,n = n,G = G)
  #####################################################################################
  g2 <- function(j,G) ifelse(j<=G, (1-kappa)/(G+1),
                             (1-kappa)/(G+1) *
                               (((G+1)*kappa)/(G*kappa+1)) ^ (j-G-1)  )
  #####################################################################################
  if (verbose) {
    cat("MCMC progress:\n")
    utils::flush.console()
    pbar <- utils::txtProgressBar(min = 1,
                           max = nsim*thinning + burn_in,
                           style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }

  for(sim in 1:(nsim*thinning + burn_in)){

    ############################################################################
    Ui    <- stats::runif(n, 0, g2(ZETA,G = G,kappa = kappa))
    L.z   <- threshold.slice(min(Ui),kappa,G)
    #if(length(L.z)==0){ L.z <- 1 }
    G.z  <- max(alphabeta[,2])
    L    <- max(c(L.z, G.z))

    xi   <- g2(1:(L), G = G,kappa = kappa)
    PL   <- c(1:(L))

    # old labels = G
    L_new = L - G
    #############################################################################
    pidir   <- c(Update.pi(alphabeta = alphabeta, aDir = aDir, G = G))
    #######################################################################
    u       <- UPD_Sticks_Beta_cpp(AB = alphabeta,L_new = L_new,alphaDP = aDP)
    if(L==1){
      omega <- 1
    }else{
      omega <- StickBreaker_cpp(u)
    }

    log_pitilde <- log_pitilde.maker(pidir, omega, G)  # pitilde contains L = G + L_new components
    #######################################################################
    tau        <- Update_tau_l(a_tau = a_tau,b_tau = b_tau,
                               effe_test = as.matrix(effe.test),
                               L_new = L_new,
                               t = t,
                               OOR = OneOverR,
                               a_l = ak)
    ak         <- Update_a_l(effe_test = as.matrix(effe.test),L_new = L_new,tau_l = tau,
                             R = R,s = s_tau,t = t,OneOverR = OneOverR)
    #####
    effe.test  <- Update_effe_test_t_cpp(Y = Y,
                                         alphabeta =  alphabeta, L_new =  L_new,
                                         R =  R,
                                         ak = ak,
                                         t = t,
                                         tau_k = tau,
                                         sigma_test = as.matrix(sigma.test))

    sigma.test <- Update_sigma_test_t_cpp(Y = Y, alphabeta = alphabeta, t = t,
                                          asig = a_H,bsig =  b_H,L_new =  L_new,
                                          effe_test = effe.test)
    #######################################################################
    if (learning_type %in% c("transductive", "training_stochastic")) {

      effe.train  <- Update_effe_train_t_cpp(
        Y = Y,
        G = G,
        alphabeta =  alphabeta,
        t = t,
        f_bar_g = effe.train0,
        sigma_train = sigma.train,
        KAPPAG = KappaG
      )

      sigma.train <- Update_sigma_train_t_cpp(
        Y = Y,
        alphabeta = alphabeta,
        t = t,
        effe_train = effe.train,
        G = G,
        a_priorG = A_B_IG_Train_Stochastic[, , 1],
        b_priorG = A_B_IG_Train_Stochastic[, , 2]
      )
    }
    #######################################################################
    effe_all  <- cbind(effe.train, effe.test)
    sigma_all <- cbind(sigma.train,sigma.test)
    ###############################################################################
    ZETA <- Upd_ZETA_t_SLICE(Y = Y,effe_ALL = effe_all,
                             sigma_ALL = sigma_all,
                             Uitilde = Ui,
                             xitilde = xi,
                             log_pitilde = log_pitilde,
                             G = G,n = n,L_new = L_new,poss_lab = PL)
    ########################################################################
    # prime 4 righe forse lasciabili per questioni di output piu leggero
    #ind2         <- ZETA[ZETA>G]
    #u.ind        <- unique(sort(ind2))-G
    #effe.test    <- effe.test[,u.ind]
    #sigma.test   <- sigma.test[,u.ind]
    #ZETA[ZETA>G] <- as.numeric(factor(ZETA[ZETA>G]))+G

    alphabeta <- zeta_to_alphabeta(ZETA,n,G)
    ############################################################################################
    if (fixed_alphaDP == F){
      aDP <-  Update_concentration_parameter(alphabeta = alphabeta,
                                             a_alpha = a_alpha, b_alpha = b_alpha,
                                             aDP = aDP)
    }
    ############################################################################################

    if (sim > burn_in && ((sim - burn_in) %% thinning == 0)) {
      rr                <- floor((sim - burn_in)/thinning);
      PROB[rr,]         <- pidir
      ALPHA.BETA[,,rr]  <- alphabeta
      if (learning_type %in% c("transductive", "training_stochastic")) {
        EFFE.TRAIN[, , rr]  <- effe.train
        SIGMA.TRAIN[, , rr] <- sigma.train
      }
      OMEGA[[rr]]       <- omega
      EFFE.TEST[[rr]]   <- effe.test
      SIGMA.TEST[[rr]]  <- sigma.test
      TAU[[rr]]         <- tau
      AK[[rr]]          <- ak
      AlphaDP[rr]       <- aDP
      USLICE[,rr]       <- Ui
    }
    ################################################################
    ################################################################
    if (verbose) {
      ipbar <- ipbar + 1
      utils::setTxtProgressBar(pbar, ipbar)
    }
  }

  # If training fixed I only report the output from robust information extraction
  if (learning_type %in% c("inductive", "training_fixed")) {
    EFFE.TRAIN  <- effe.train
    SIGMA.TRAIN <- sigma.train
  }
  out <-  list(
    P   =  PROB,
    AB  =  ALPHA.BETA,
    O   =  OMEGA,
    TAU =  TAU,
    AK  =  AK,
    aDP =  AlphaDP,
    FTr =  EFFE.TRAIN,
    STr =  SIGMA.TRAIN,
    FTr0 =  effe.train0,
    STr0 =  sigma.train0,
    FTe =  EFFE.TEST,
    STe =  SIGMA.TEST,
    USL = USLICE,
    prior=prior)
  return(out)
}

