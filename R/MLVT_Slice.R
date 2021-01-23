



StartingEst <- function(X, categ, raw_MCD, h_MCD, G) {

  p <- ncol(X)

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



###################################################################################

#' Multivariate Brand
#'
#' @param Y data
#' @param X data
#' @param categ xxx
#' @param prior xxx
#' @param L xxx
#' @param h_MCD xxx
#' @param raw_MCD xxx
#' @param NSIM xxx
#' @param thinning xxx
#' @param burn_in xxx
#' @param verbose xxx
#' @param fixed_alphaDP xxx
#' @param kappa xxx
#' @param light xxx
#' @param learning_type xxx
#'
#' @return model
#' @export
#'
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
  Y <- data.matrix(Y) # for avoiding problems in dealing with data.frame objects
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
  LZ           <- numeric(NSIM)
  U_Slice      <- matrix(NA, n, NSIM)
  #####################################################################################
  # inizialization
  pidir         <- c(MCMCpack::rdirichlet(1, aDir))
  alphabeta     <- matrix(0, n, 2)
  alphabeta[, 1] <- sample(0:G, n, replace = T)
  u             <- stats::rbeta(n = L,
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
  Sigma         <- replicate(L, MCMCpack::riwish(v = v_H, S = S_H))
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


  ZETA <- alphaneta_to_zeta(alphabeta, n = n, G = G)
  # #####################################################################################
  if (verbose) {
    cat("MCMC progress:\n")
    utils::flush.console()
    pbar <- utils::txtProgressBar(min = 1,
                           max = NSIM * thinning + burn_in,
                           style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  # #####################################################################################
  for (sim in 1:(NSIM * thinning + burn_in)) {


    Ui    <- stats::runif(n, 0, g2(ZETA,G = G,kappa = kappa))
    L.z   <- threshold.slice(min(Ui),kappa,G)
    #if(length(L.z)==0){ L.z <- 1 }
    G.z  <- max(alphabeta[,2])
    L    <- max(c(L.z, G.z))

    xi   <- g2(1:(L), G = G,kappa = kappa)
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

    #ind2 <- ZETA[ZETA > G]
    #u.ind <- unique(sort(ind2)) - G
    #mu    <- mu[, u.ind]
    #Sigma <- Sigma[, , u.ind]
    #ZETA[ZETA > G] <- as.numeric(factor(ZETA[ZETA > G])) + G

    alphabeta <- zeta_to_alphabeta(ZETA, n, G)

    ############################################################################################
    if (fixed_alphaDP == F){
      aDP <-  Update_concentration_parameter(alphabeta = alphabeta,
                                             a_alpha = a_alpha, b_alpha = b_alpha,
                                             aDP = aDP)
    }
    ############################################################################################

    if (sim > burn_in &&
        ((sim - burn_in) %% thinning == 0)) {
      rr                 <- floor((sim - burn_in) / thinning)
      LZ[rr]             <- L.z
      PiDir[rr, ]        <- pidir
      Alpha_Beta[, , rr] <- alphabeta
      Omega[[rr]]        <- omega
      MU_new[[rr]]       <- t(mu)
      SIGMA_new[[rr]]    <- Sigma
      if (learning_type %in% c("transductive", "training_stochastic")) {
        MU_train[, , rr]     <- t(mu.train)
        SIGMA_train[, , , rr] <- Sigma.train
      }
      AlphaDP[rr]         <- aDP
      U_Slice[, rr]       <- Ui
    }
    ################################################################
    ################################################################
    #   if (sim%%(verbose.step*thinning) == 0) {
    #     cat(paste("Sampling iteration: ", sim, " out of ",NSIM*thinning + burn_in,"\n",round(aDP,3),"\n"))}

    if (verbose) {
      ipbar <- ipbar + 1
      utils::setTxtProgressBar(pbar, ipbar)
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
      Uslice    = U_Slice,
      LZ        = LZ
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
      prior     = prior,
      LZ        = LZ

    )
  }


  return(out)
}
