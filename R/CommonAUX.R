# Common Auxiliary Functions


#' Major Vote Computer
#'
#' @param x xxx
#'
#' @return major label
#' @export
#'
#' @examples
#' major.vote(  c(1,1,1,2,2,2,3,3,3,4,4,4,4) )
#'
major.vote <- function(x){
  as.numeric(names(sort(table(x),decreasing = TRUE))[1])
}


Update.pi <- function(alphabeta,aDir,G){
  alpha        <- alphabeta[,1]
  ns           <- table(factor(alpha,levels = 0:G))
  nj           <- c(ns[-1],ns[1])
  aDirtilde    <- nj+aDir
  return(MCMCpack::rdirichlet(1,aDirtilde))
}


g2 <- function(j,G,kappa) ifelse(j<=G, (1-kappa)/(G+1),
                           (1-kappa)/(G+1) *
                             (((G+1)*kappa)/(G*kappa+1)) ^ (j-G-1)  )

threshold.slice <- function(Ui, kappa, G){
  G + 1 + floor(
    (log(Ui) - (log( (G * kappa + 1) / (G + 1) ) +
                  log( (1 - kappa) / (G * kappa + 1) )))/
      log(((G * kappa + kappa) / (G * kappa + 1)) ))
}


log_pitilde.maker <- function(pi, omega, G) {
  return(c(log(pi[-(G + 1)]), log(pi[G + 1]) + log(omega)))

}

alphaneta_to_zeta <- function(AB,n,G){
  Z           <- AB[,1]
  Z[AB[,2]>0] <- AB[AB[,2]>0,2] + G
  return(Z)
}

zeta_to_alphabeta <- function(Z,n,G){
  AB   <- matrix(0,n,2)
  ind1 <- which(Z<=G)

  if(length(ind1)==0){
    AB[,2] <- Z - G
  }else{
    AB[ind1,1]  <- Z[ind1]
    AB[-ind1,2] <- Z[-ind1] - G
  }


  return(AB)
}




Update_concentration_parameter <- function(alphabeta, a_alpha, b_alpha, aDP){

  beta0  <- alphabeta[,2]
  uz     <- length(unique(beta0[beta0>0]))
  n_nov  <- sum(alphabeta[,1]==0)
  if(uz==0){ # if no obs in novelty, sample from prior
    aDP <- stats::rgamma(1,a_alpha,b_alpha)
  }else{
    eta    <- stats::rbeta(1, aDP + 1, n_nov)
    Q      <- (a_alpha + uz - 1) / (n_nov * (b_alpha - log(eta)))
    pi_eta <- Q / (1 + Q)
    aDP    <-
      ifelse(
        stats::runif(1) < pi_eta,
        stats::rgamma(1, a_alpha + uz,     b_alpha - log(eta)),
        stats::rgamma(1, a_alpha + uz - 1, b_alpha - log(eta))
      )
  }
  return(aDP)
}
