
# Update Dirichlet weights ------------------------------------------------

update_pi <- function(alphabeta,aDir,J){
  alpha        <- alphabeta[,1]
  ns           <- table(factor(alpha,levels = 0:J))
  nj           <- c(ns[-1],ns[1])
  aDirtilde    <- nj+aDir
  return( MCMCpack::rdirichlet(1,aDirtilde))
}


# Sigma(t) parameters -----------------------------------------------------

a_b_sigmat_train <- function(sigma_est_train, # matrix TxG
  vg,t,G){
  ab        <- array(NA, c(t, G, 2))
  ab[, , 1] <- 2 + (sigma_est_train) ^ 2 / vg
  ab[, , 2] <- ((sigma_est_train) ^ 2 / vg  + 1) * sigma_est_train
  return(ab)
}


# BNP weights -------------------------------------------------------------

omegatilde_maker <- function(pi, omega, J) {
  return(c ( pi[-(J + 1)], pi[J + 1] * omega) )
}

log_omegatilde_maker <- function(pi, omega, J) {
  return( c( log( pi[-(J + 1)] ), log( pi[J + 1] ) + log( omega ) ) )
}


# Clean Indexes -----------------------------------------------------------

alphabeta_to_zeta <- function(AB, n, J) {
  Z              <- AB[, 1]
  Z[AB[, 2] > 0] <- AB[AB[, 2] > 0, 2] + J
  return(Z)
}

zeta_to_alphabeta <- function(Z, n, J) {
  AB           <- matrix(0, n, 2)
  ind1         <- which(Z <= J)
  AB[ind1, 1]  <- Z[ind1]
  AB[-ind1, 2] <- Z[-ind1] - J
  return(AB)
}



# Slice g -----------------------------------------------------------------

#' Auxiliary Slice function
#'
#' @param j the membership index
#' @param J the number of known components
#' @param kappa the parameter driving the decay of the slice sticks
#'
#' @return the value of the slice sticks at point j
#'
g2 <- function(j, J, kappa=.75){
  ifelse(j <= J, (1 - kappa) / (J + 1),
                 (1 - kappa) / (J + 1) * (((J + 1) * kappa) /
                 (J * kappa + 1)) ^ (j - J - 1))

}
