#' Gibbs sampler for functional BRAND
#'
#' The function runs a Gibbs sampler for the functional brand model.
#'
#' @param Y a matrix containing n functions (columns) of length t (rows)
#' @param prior a list with prior information needed to run the model
#' @param L the starting number of novelty components
#' @param nsim the number of saved iterations
#' @param thinning the thinning period
#' @param burn_in the number of discarded iterations
#' @param fixed_alphaDP logical, should the concentration parameter be held fixed?
#' @param kappa the parameter defining the
#' @param verbose logical, should the iterations be printed?
#' @param learning_type Inductive or Trunsductive learning
#'
#'
#'
#'
#' @return A list containing
#'     \itemize{
#'     \item \code{...}
#'     \item \code{...}
#'     }
#'
#' @export
functional_brand <- function( Y, prior, L = 10, nsim = 1000, thinning = 1, burn_in = 1000, fixed_alphaDP = F, kappa = .5, verbose = 1, learning_type = c("transductive", "inductive") ) {
  # Extraction of Prior Information -----------------------------------------
  basis        <- prior$basis
  beta_train   <- prior$beta_train
  sigma_train0 <- prior$sigma_train

# Basis functions ---------------------------------------------------------
  effe_train0  <- basis %*% beta_train
  vg           <- prior$vg
  KappaG       <- prior$KappaG
  R            <- rowSums(basis ^ 2)
  OneOverR     <- 1 / R

# Dimensions --------------------------------------------------------------
  n            <- ncol(Y)
  t            <- nrow(Y)
  J            <- ncol(beta_train)
  K            <- nrow(beta_train)

# Hyperpriors for Stochastic training -------------------------------------
  alphaD      <- rep(1, J)
  A_B_IG_Train_Stochastic <- a_b_sigmat_train(sigma_est_train = sigma_train0,
                                              vg = vg,
                                              t = t,
                                              G = J)

# Containers --------------------------------------------------------------
  EFFE_TEST   <- vector("list", length = nsim)
  SIGMA_TEST  <- vector("list", length = nsim)
  ###
  PROB        <- matrix(NA, nsim, J + 1)
  OMEGA       <- vector("list", length = nsim)
  TAU         <- vector("list", length = nsim)
  AK          <- vector("list", length = nsim)
  ALPHA_BETA  <- array(NA, dim = c(n, 2, nsim))
  AlphaDP     <- numeric(nsim)
  USLICE      <- matrix(NA, n, nsim)
  ###
  #if (learning_type %in% c("transductive", "training_stochastic")) {
    EFFE_TRAIN  <- array(NA, dim = c(t, J, nsim))
    sigma_train <- array(NA, dim = c(t, J, nsim))
  #}
  ###
  ##############################################
  ### Hyperparameters
  aDir     <- prior$aDir
  aDP      <- prior$aDP
  a_H      <- prior$a_H    # prior varianza gruppi nuovi
  b_H      <- prior$b_H
  a_tau    <- prior$a_tau  # varianza dei betini
  b_tau    <- prior$b_tau  # varianza dei betini
  a_alpha  <- prior$a_alphaDP
  b_alpha  <- prior$b_alphaDP
  s_tau    <- prior$s_tau
  # inizialization
  pidir         <- c(MCMCpack::rdirichlet(1, aDir))
  alphabeta     <- matrix(NA, n, 2)
  alphabeta[, 1] <- sample(0:J, n, replace = T)
  alphabeta[, 2] <- ifelse(alphabeta[, 1] > 0, 0, sample(L))
  u             <- stats::rbeta(n = L, shape1 = 1, shape2 = aDP)
  omega         <- StickBreaker_cpp(V = u)
  effe_test     <- matrix(NA, t, L)
  sigma_test    <- matrix(NA, t, L)
  effe_train    <- effe_train0 #matrix(NA,t,J)
  sigma_train   <- sigma_train0 #matrix(NA,t,J)
  tau           <- numeric(L)
  ak            <- numeric(L)

  for (k in 1:L) {
    sigma_test[, k] <- MCMCpack::rinvgamma(t, a_H, b_H)
    tau[k]          <- MCMCpack::rinvgamma(1, a_tau, b_tau)
    ak[k]           <- stats::rnorm(1, 0, sqrt(s_tau))
    effe_test[, k]  <- stats::rnorm(t, 0, sd = sqrt(tau[k] * R))
  }
  ZETA <- alphabeta_to_zeta(alphabeta, n = n, J = J)

  if (verbose) {
    cat("MCMC progress:\n")
    flush.console()
    pbar <- txtProgressBar(min = 1,
      max = nsim * thinning + burn_in,
      style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }

  for (sim in 1:(nsim * thinning + burn_in)) {
    ############################################################################
    Ui    <- runif(n, 0, g2(ZETA, J = J))
    L_z   <- J + 1 + floor((log(Ui) - (log( (J * kappa + 1) / (J + 1) ) +
             log( (1 - kappa) / (J * kappa + 1) ))) /
             log(((J * kappa + kappa) / (J * kappa + 1))))

    if (length(L_z) == 0) {
      L_z <- 1
    }

    J_z  <- max(alphabeta[, 2])
    L    <- max(c(L_z, J_z))
    xi   <- g2(1:L, J = J)
    PL   <- c(1:(J + L))
    #############################################################################
    pidir   <- c(update_pi(alphabeta = alphabeta, aDir = aDir, J = J ))
    #######################################################################
    u  <-  UPD_Sticks_Beta_cpp(AB = alphabeta, L = L, alphaDP = aDP)
    if (L == 1) {
      omega <- 1
    } else{
      omega <- StickBreaker_cpp(u)
    }
    log_omegatilde <- log_omegatilde_maker(pidir, omega, J)
    #######################################################################
    tau        <- Update_tau_l(
      a_tau = a_tau,
      b_tau = b_tau,
      effe_test = as.matrix(effe_test),
      L = L,
      t = t,
      OOR = OneOverR,
      a_l = ak    )
    ak         <- Update_a_l(effe_test = as.matrix(effe_test), L = L, tau_l = tau,
                             R = R, s = s_tau, t = t ,OneOverR = OneOverR)
    #######################################################################
    effe_test  <- Update_effe_test_t_cpp(
      Y = Y,
      alphabeta =  alphabeta,
      L =  L,
      R =  R,
      ak = ak,
      t = t,
      tau_k = tau,
      sigma_test = as.matrix(sigma_test)
    )

    sigma_test <-
      Update_sigma_test_t_cpp(
        Y = Y,
        alphabeta = alphabeta,
        t = t,
        asig = a_H,
        bsig =  b_H,
        L = L,
        effe_test = effe_test
      )
    #######################################################################
    if (learning_type %in% c("transductive", "training_stochastic")) {
      effe_train  <- Update_effe_train_t_cpp(
        Y = Y,
        G = J,
        alphabeta =  alphabeta,
        t = t,
        f_bar_g = effe_train0,
        sigma_train = sigma_train,
        KAPPAG = KappaG
      )

      sigma_train <- Update_sigma_train_t_cpp(
        Y = Y,
        alphabeta = alphabeta,
        t = t,
        effe_train = effe_train,
        G = J,
        a_priorG = A_B_IG_Train_Stochastic[, , 1],
        b_priorG = A_B_IG_Train_Stochastic[, , 2]
      )
    }
    #######################################################################
    effe_all  <- cbind(effe_train, effe_test)
    sigma_all <- cbind(sigma_train, sigma_test)
    ###############################################################################
    ZETA <- Upd_ZETA_t_SLICE(
      Y = Y,
      effe_ALL = effe_all,
      sigma_ALL = sigma_all,
      Uitilde = Ui,
      xitilde = xi,
      log_omegatilde = log_omegatilde,
      J = J,
      n = n,
      L = L,
      poss_lab = PL
    )
    ########################################################################
    ind2         <- ZETA[ZETA > J]
    u.ind        <- unique(sort(ind2)) - J
    effe_test    <- effe_test[, u.ind]
    sigma_test   <- sigma_test[, u.ind]
    ZETA[ZETA > J] <- as.numeric(factor(ZETA[ZETA > J])) + J
    alphabeta <- zeta_to_alphabeta(ZETA, n, J)
    ###############################################################################
    if (fixed_alphaDP  & sim == 1) {
      aDP      <- aDP
    } else if (fixed_alphaDP == F) {
      beta0  <- alphabeta[, 2]
      uz     <- length(unique(beta0[beta0 > 0]))
      eta    <- stats::beta(1, aDP + 1, n)
      Q      <- (a_alpha + uz - 1) / (n * (b_alpha - log(eta)))
      pi_eta <- Q / (1 + Q)
      aDP    <-
        ifelse(
          runif(1) < pi_eta,
          rgamma(1, a_alpha + uz,   b_alpha - log(eta)),
          rgamma(1, a_alpha + uz - 1, b_alpha - log(eta))
        )
    }
    ###############################################################################
    if (sim > burn_in && ((sim - burn_in) %% thinning == 0)) {
      rr                <- floor((sim - burn_in) / thinning)

      PROB[rr, ]         <- pidir
      ALPHA_BETA[, , rr]  <- alphabeta
      if (learning_type %in% c("transductive", "training_stochastic")) {
        EFFE_TRAIN[, , rr]  <- effe_train
        sigma_train[, , rr] <- sigma_train
      }
      OMEGA[[rr]]       <- omega
      EFFE_TEST[[rr]]   <- effe_test
      SIGMA_TEST[[rr]]  <- sigma_test
      TAU[[rr]]         <- tau
      AK[[rr]]          <- ak
      AlphaDP[rr]       <- aDP
      USLICE[, rr]       <- Ui
    }
    ################################################################
    if (verbose) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
  }

  # If training fixed I only report the output from robust information extraction
  if (learning_type %in% c("inductive", "training_fixed")) {
    EFFE_TRAIN  <- effe_train
    sigma_train <- sigma_train
  }
  out <-  list(
    P   =   PROB,
    AB  =   ALPHA_BETA,
    O   =   OMEGA,
    TAU =   TAU,
    AK  =   AK,
    aDP =   AlphaDP,
    FTr =   EFFE_TRAIN,
    STr =   sigma_train,
    FTr0 =  effe_train0,
    STr0 =  sigma_train0,
    FTe =   EFFE_TEST,
    STe =   SIGMA_TEST,
    USL =   USLICE
  )
  return(out)
}
