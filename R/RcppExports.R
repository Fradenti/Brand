# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

PSM <- function(inds) {
    .Call(`_brand_PSM`, inds)
}

LogSumExp <- function(logX) {
    .Call(`_brand_LogSumExp`, logX)
}

dmvnrm_arma <- function(x, mean, sigma, logd = FALSE) {
    .Call(`_brand_dmvnrm_arma`, x, mean, sigma, logd)
}

UPD_Sticks_Beta_cpp <- function(AB, L_new, alphaDP) {
    .Call(`_brand_UPD_Sticks_Beta_cpp`, AB, L_new, alphaDP)
}

StickBreaker_cpp <- function(V) {
    .Call(`_brand_StickBreaker_cpp`, V)
}

Update_Theta_cpp <- function(Y, AB, m_H, k_H, v_H, S_H, L_new, p) {
    .Call(`_brand_Update_Theta_cpp`, Y, AB, m_H, k_H, v_H, S_H, L_new, p)
}

Update_Theta_cpp_TRAIN <- function(Y, AB, MU_g, k_g, SIGMA_g, v_g, G, p) {
    .Call(`_brand_Update_Theta_cpp_TRAIN`, Y, AB, MU_g, k_g, SIGMA_g, v_g, G, p)
}

Upd_Zeta_cpp <- function(Y, mu_all, Sigma_all, Uitilde, xitilde, log_pitilde, G, n, L_new, poss_lab) {
    .Call(`_brand_Upd_Zeta_cpp`, Y, mu_all, Sigma_all, Uitilde, xitilde, log_pitilde, G, n, L_new, poss_lab)
}

Upd_alphabeta_t_cpp <- function(Y, effe_train, effe_test, sigma_train, sigma_test, pidir, omega, G, n, L, t, poss_lab) {
    .Call(`_brand_Upd_alphabeta_t_cpp`, Y, effe_train, effe_test, sigma_train, sigma_test, pidir, omega, G, n, L, t, poss_lab)
}

DN <- function(X, M, S) {
    .Call(`_brand_DN`, X, M, S)
}

Update_effe_test_t_cpp <- function(Y, alphabeta, L_new, t, tau_k, R, sigma_test, ak) {
    .Call(`_brand_Update_effe_test_t_cpp`, Y, alphabeta, L_new, t, tau_k, R, sigma_test, ak)
}

Update_effe_train_t_cpp <- function(Y, alphabeta, G, t, f_bar_g, sigma_train, KAPPAG) {
    .Call(`_brand_Update_effe_train_t_cpp`, Y, alphabeta, G, t, f_bar_g, sigma_train, KAPPAG)
}

Update_sigma_test_t_cpp <- function(Y, alphabeta, effe_test, L_new, t, asig, bsig) {
    .Call(`_brand_Update_sigma_test_t_cpp`, Y, alphabeta, effe_test, L_new, t, asig, bsig)
}

Update_sigma_train_t_cpp <- function(Y, alphabeta, effe_train, G, t, a_priorG, b_priorG) {
    .Call(`_brand_Update_sigma_train_t_cpp`, Y, alphabeta, effe_train, G, t, a_priorG, b_priorG)
}

Update_tau_l <- function(a_tau, b_tau, effe_test, L_new, t, OOR, a_l) {
    .Call(`_brand_Update_tau_l`, a_tau, b_tau, effe_test, L_new, t, OOR, a_l)
}

Update_a_l <- function(effe_test, L_new, tau_l, R, s, t, OneOverR) {
    .Call(`_brand_Update_a_l`, effe_test, L_new, tau_l, R, s, t, OneOverR)
}

Upd_ZETA_t_SLICE <- function(Y, effe_ALL, sigma_ALL, log_pitilde, Uitilde, xitilde, G, n, L_new, poss_lab) {
    .Call(`_brand_Upd_ZETA_t_SLICE`, Y, effe_ALL, sigma_ALL, log_pitilde, Uitilde, xitilde, G, n, L_new, poss_lab)
}

