library(LaplacesDemon)
library(MCMCpack)

major.vote <- function(x){
  as.numeric(names(sort(table(x),decreasing = TRUE))[1]) 
}


Update.pi <- function(alphabeta, aDir, G) {
  alpha        <- alphabeta[, 1]
  ns           <- table(factor(alpha, levels = 0:G))
  # tabulate??
  nj           <- c(ns[-1], ns[1])
  aDirtilde    <- nj + aDir
  return(MCMCpack::rdirichlet(1, aDirtilde))
}

################################ robust estimates extractor

log_pitilde.maker <- function(pi, omega, G) {
  return(c(log(pi[-(G + 1)]), log(pi[G + 1]) + log(omega)))
}


StartingEst <- function(X, categ, raw_MCD, h_MCD, G) {
  
  if (h_MCD == 1) {
    #standard non-robust methods are computed
    estimates_from_train <-
      mclust::covw(X = X, Z = mclust::unmap(categ))
    xbar_j               <- estimates_from_train$mean
    S2_j                 <- estimates_from_train$S
  } else {
    robust_estimates_from_train <-
      lapply(unique(sort(categ)), function(group)
        rrcov::CovMcd(x = X[categ == group, ], alpha = h_MCD))
    
    if (raw_MCD) {
      xbar_j <-
        sapply(1:G, function(g)
          robust_estimates_from_train[[g]]$center)
      S2_j   <-
        array(
          sapply(1:G, function(g)
            robust_estimates_from_train[[g]]$cov) / nrow(X) * (nrow(X) - 1),
          dim = c(p, p, G)
        )
    } else{
      xbar_j <-
        sapply(1:G, function(g)
          robust_estimates_from_train[[g]]$raw.center)
      S2_j   <-
        array(
          sapply(1:G, function(g)
            robust_estimates_from_train[[g]]$raw.cov) / nrow(X) * (nrow(X) -
                                                                     1),
          dim = c(p, p, G)
        )
    }
  }
  return(list(xbar_j, S2_j))
}

alphaneta_to_zeta <- function(AB, n, G) {
  Z           <- AB[, 1]
  Z[AB[, 2] > 0] <- AB[AB[, 2] > 0, 2] + G
  return(Z)
}

zeta_to_alphabeta <- function(Z, n, G) {
  AB   <- matrix(0, n, 2)
  ind1 <- which(Z <= G)
  AB[ind1, 1]  <- Z[ind1]
  AB[-ind1, 2] <- Z[-ind1] - G
  return(AB)
}


###################################################################################

Brand_mlvt <- function(Y,
                                             X,
                                             categ,
                                             prior,
                                             L,
                                             h_MCD = .75,
                                             # percentage of untrimmed observations for computing robust prior mean vector and scatter from train set
                                             raw_MCD = F,
                                             NSIM,
                                             thinning,
                                             burn_in,
                                             verbose = T,
                                             fixed_alphaDP = F,
                                             kappa = .75,
                                             light = F,
                                             learning_type = c("transductive", "inductive")) {
  #####################################################################################
  Y <-
    data.matrix(Y) # for avoiding problems in dealing with data.frame objects
  X <- data.matrix(X)
  #####################################################################################
  n      <- nrow(Y) # sample size test data
  p      <- ncol(Y)
  G      <- length(unique(categ)) # # known categories
  #####################################################################################
  SERob  <-
    StartingEst(
      X = X,
      categ = categ,
      raw_MCD = raw_MCD,
      h_MCD = h_MCD,
      G=G
    )
  xbar_j <- SERob[[1]]
  S2_j   <- SERob[[2]]
  #####################################################################################
  ### Hyperparameters
  aDir     <- prior$aDir
  aDP      <- prior$aDP
  m_H      <- prior$m_H
  k_H      <- prior$k_H
  v_H      <- prior$v_H
  S_H      <- prior$S_H
  k_g      <- prior$k_g
  v_g      <- prior$v_g
  a_alpha  <- prior$a_alphaDP
  b_alpha  <- prior$b_alphaDP
  #####################################################################################
  ### Containers
  Alpha_Beta   <- array(NA, dim = c(n, 2, NSIM))
  MU_new       <- vector("list", length = NSIM)
  SIGMA_new    <- vector("list", length = NSIM)
  
  if (learning_type %in% c("transductive", "training_stochastic")) {
    MU_train     <- array(NA, dim = c(G, p, NSIM))
    SIGMA_train  <- array(NA, dim = c(p, p, G, NSIM))
  }
  
  PiDir        <- matrix(NA, NSIM, G + 1) # G+1-th ? pi0
  Omega        <- vector("list", length = NSIM)
  AlphaDP      <- numeric(NSIM)
  U_Slice      <- matrix(NA, n, NSIM)
  #####################################################################################
  # inizialization
  pidir         <- c(rdirichlet(1, aDir))
  alphabeta     <- matrix(0, n, 2)
  alphabeta[, 1] <- sample(0:G, n, replace = T)
  u             <- rbeta(n = L,
                         shape1 = 1,
                         shape2 = aDP)
  omega         <- StickBreaker_cpp(V = u)
  alphabeta[alphabeta[, 1] == 0, 2] <-
    as.numeric(factor(sample(
      1:L,
      size =
        sum(alphabeta[, 1] ==
              0),
      replace = T,
      prob = omega
    )))
  Sigma         <- replicate(L, riwish(v = v_H, S = S_H))
  mu            <-
    sapply(1:L, function(ind)
      mvtnorm::rmvnorm(
        n = 1,
        mean = m_H,
        sigma = Sigma[, , ind] / k_H
      ))
  mu.train      <- xbar_j
  Sigma.train   <- S2_j
    # #####################################################################################
  g2 <- function(j,G) ifelse(j<=G, (1-kappa)/(G+1), 
                             (1-kappa)/(G+1) *
                               (((G+1)*kappa)/(G*kappa+1)) ^ (j-G-1)  )
  
  ZETA <- alphaneta_to_zeta(alphabeta, n = n, G = G)
  # #####################################################################################
  if (verbose) {
    cat("MCMC progress:\n")
    flush.console()
    pbar <- txtProgressBar(min = 1,
                           max = NSIM * thinning + burn_in,
                           style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  # #####################################################################################
  for (sim in 1:(NSIM * thinning + burn_in)) {
    
    
    Ui    <- runif(n, 0, 100*g2(ZETA,G = G))/100
    L.z   <- G + 1 + floor(
      (log(Ui) - (log( (G * kappa + 1) / (G + 1) ) + 
                    log( (1 - kappa) / (G * kappa + 1) )))/
        log(((G * kappa + kappa) / (G * kappa + 1)) ))
    if(length(L.z)==0){ L.z <- 1 }
    G.z  <- max(alphabeta[,2])
    L    <- max(c(L.z, G.z))
    xi   <- g2(1:(L), G = G)
    PL   <- c(1:(L))
    
    # old labels = G
    L_new = L - G
    ##################################################################################
    pidir      <-
      c(Update.pi(
        alphabeta = alphabeta,
        aDir = aDir,
        G = G
      ))
    ###################################################################################
    u          <-
      UPD_Sticks_Beta_cpp(AB = alphabeta,
                          L_new = L_new,
                          alphaDP = aDP)
    if (L == 1) {
      omega <- 1
    } else{
      omega <- StickBreaker_cpp(u)
    }
    
    # omegatilde <- Omegatilde.maker(pidir,omega,G)
    log_pitilde <- log_pitilde.maker(pidir, omega, G)
    ###################################################################################
    ####################################################################################
    TH         <- Update_Theta_cpp(
      Y = Y,
      AB = alphabeta,
      m_H = m_H,
      k_H = k_H,
      v_H = v_H,
      S_H = S_H,
      L_new = L_new,
      p = p
    )
    mu         <- TH[[1]]
    Sigma      <- TH[[2]]
    #####################################################################################
    if (learning_type %in% c("transductive", "training_stochastic")) {
      THETA_g    <- Update_Theta_cpp_TRAIN(
        Y = Y,
        AB = alphabeta,
        MU_g = xbar_j,
        k_g = k_g,
        v_g = v_g,
        SIGMA_g =  S2_j * (v_g - p - 1),
        G = G,
        p = p
      )
      
      mu.train    <- THETA_g[[1]]
      Sigma.train <- THETA_g[[2]]
    }
    mu_all <- cbind(mu.train, mu)
    Sigma_all <- abind::abind(Sigma.train, Sigma, along = 3)
    
    ZETA <-
      Upd_Zeta_cpp(
        Y = Y,
        mu_all = mu_all,
        Sigma_all = Sigma_all,
        Uitilde = Ui,
        xitilde = xi,
        log_pitilde = log_pitilde,
        G = G,
        n = n,
        L_new = L_new,
        poss_lab = 1:(L_new + G)
      )
    
    ind2 <- ZETA[ZETA > G]
    u.ind <- unique(sort(ind2)) - G
    mu    <- mu[, u.ind]
    Sigma <- Sigma[, , u.ind]
    ZETA[ZETA > G] <- as.numeric(factor(ZETA[ZETA > G])) + G
    
    alphabeta <- zeta_to_alphabeta(ZETA, n, G)
    
    if (fixed_alphaDP  & sim == 1) {
      aDP      <- aDP
    }else if (fixed_alphaDP == F){
      beta0  <- alphabeta[,2]
      uz     <- length(unique(beta0[beta0>0]))
      n_nov  <-  sum(alphabeta[,1]==0)
      if(uz==0){ # if no obs in novelty, sample from prior
        aDP <- rgamma(1,a_alpha,b_alpha)
      }else{
      eta    <- rbeta(1, aDP + 1, n_nov)
      Q      <- (a_alpha + uz - 1) / (n_nov * (b_alpha - log(eta)))
      pi_eta <- Q / (1 + Q)
      aDP    <-
        ifelse(
          runif(1) < pi_eta,
          rgamma(1, a_alpha + uz,   b_alpha - log(eta)),
          rgamma(1, a_alpha + uz - 1, b_alpha - log(eta))
        )
    }}
    ############################################################################################
    
    if (sim > burn_in &&
        ((sim - burn_in) %% thinning == 0)) {
      rr                 <- floor((sim - burn_in) / thinning)
      
      PiDir[rr, ]         <- pidir
      Alpha_Beta[, , rr]   <- alphabeta
      Omega[[rr]]        <- omega
      MU_new[[rr]]       <- t(mu)
      SIGMA_new[[rr]]    <- Sigma
      if (learning_type %in% c("transductive", "training_stochastic")) {
        MU_train[, , rr]     <- t(mu.train)
        SIGMA_train[, , , rr] <- Sigma.train
      }
      AlphaDP[rr]        <- aDP
      U_Slice[, rr]       <- Ui
    }
    ################################################################
    ################################################################
    #   if (sim%%(verbose.step*thinning) == 0) {
    #     cat(paste("Sampling iteration: ", sim, " out of ",NSIM*thinning + burn_in,"\n",round(aDP,3),"\n"))}
    
    if (verbose) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
  }
  
  # If training fixed I only report the output from robust information extraction
  if (learning_type %in% c("inductive", "training_fixed")) {
    MU_train     <- t(mu.train)
    SIGMA_train <- Sigma.train
  }
  
  if (light) {
    out <-  list(
      P         = PiDir,
      AB        = Alpha_Beta,
      O         = Omega,
      Mu        = MU_new,
      Sigma     = SIGMA_new,
      Mu.train  = MU_train,
      Sig.train = SIGMA_train,
      x_b       = xbar_j,
      s_b       = S2_j,
      alphaDP   = AlphaDP,
      Uslice    = U_Slice
    )
    
  } else{
    out <-  list(
      P         = PiDir,
      AB        = Alpha_Beta,
      O         = Omega,
      Mu        = MU_new,
      Sigma     = SIGMA_new,
      Mu.train  = MU_train,
      Sig.train = SIGMA_train,
      x_b       = xbar_j,
      s_b       = S2_j,
      alphaDP   = AlphaDP,
      Uslice    = U_Slice,
      prior     = prior
    )
  }
  
  
  return(out)
}
