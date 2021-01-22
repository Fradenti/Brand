#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// General functions


// [[Rcpp::export]]
arma::mat PSM(arma::mat inds){
  int nsim = inds.n_rows;
  int n= inds.n_cols;
  arma::mat PSM(n,n), D(n,n); 
  D.eye(); PSM.zeros();
  
  for(int i=0; i<n; i++){
    for(int j=i+1; j<n; j++){
      arma::colvec Z = (inds.col(i)-inds.col(j));
      arma::uvec success = find(Z==0);
      PSM(i,j) = success.n_elem;
      PSM(j,i) = PSM(i,j);
    }
  }
  
  return(PSM/nsim + D);
  
}
// [[Rcpp::export]]
double LogSumExp(arma::colvec logX){
  double a = max(logX);
  return(  a + log(accu( exp( logX-a ) )));
}


const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
double dmvnrm_arma(arma::mat x,
                   arma::rowvec mean,
                   arma::mat sigma,
                   bool logd = false) {
  //int n = x.n_rows;
  int xdim = x.n_cols;
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  arma::vec z = rooti * arma::trans( x - mean) ;
  out      = constants - 0.5 * arma::sum(z%z) + rootisum;

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}



////////////// BRAND MULTIVARIATE////////

// [[Rcpp::export]]
arma::colvec UPD_Sticks_Beta_cpp(arma::mat AB,
                                 int L_new,
                                 double alphaDP){

  arma::colvec beta_lab           = AB.col(1);
  arma::uvec   inds               = find(beta_lab>0);
  arma::colvec beta_lab_newgroups = beta_lab.elem(inds);

  arma::colvec Sticks(L_new);

  if(beta_lab_newgroups.n_rows==0){

    for(int l=0; l<(L_new); l++){  //for(int l=0; l<(L_new-1); l++){
      Sticks[l] = R::rbeta(1 ,
                           alphaDP  );  }
  }else{
    for(int l=0; l<(L_new); l++){ //   for(int l=0; l<(L_new-1); l++){
      Sticks[l] = R::rbeta(1 + accu( beta_lab_newgroups == (l+1)),
                           alphaDP + accu(beta_lab_newgroups > (l+1)) );
    }
  }

  return Sticks;
}


// [[Rcpp::export]]
arma::colvec StickBreaker_cpp(arma::colvec V) {
  int N = V.size();
  arma::colvec pi(N), mV2(N);
  arma::colvec mV = 1-V;
  mV2   = exp(arma::cumsum(log(mV)));
  mV2   = arma::shift(mV2,+1);
  mV2(0) = 1;
  return(V%mV2);
}



// [[Rcpp::export]]
List Update_Theta_cpp(arma::mat Y,
                      arma::mat AB,
                      arma::colvec m_H, double k_H, double v_H, arma::mat S_H,
                      int L_new, int p){


  arma::mat    mu_new(p,L_new);
  arma::cube   Sig_new(p,p,L_new);
  arma::colvec beta_lab = AB.col(1);

  for(int l=0; l<L_new; l++){
    arma::uvec ind_l = find( beta_lab == (l+1) ) ;
    arma::mat  Y_l   = Y.rows(ind_l);
    int n_l          = (Y_l.n_rows);

    double k_ = n_l + k_H;
    double v_ = n_l + v_H;

    arma::colvec m_(p);
    arma::mat    S_(p,p);
    // Rcout << "ok1 \n";



    if(n_l>1){
      //Rcout << Y_l;
      arma::colvec my = mean(Y_l,0).t();
      //Rcout << my;
      m_ = (k_H * m_H + my * n_l) / k_;
      // Rcout << "ok2 \n";

      S_ =  S_H + cov(Y_l) * (n_l-1) +
        (k_H * n_l)/(k_H+n_l) * (my-m_H)*((my-m_H).t());
      //Rcout << (my-m_H)*((my-m_H).t());
      // Rcout << S_;
      // Rcout << "ok3 \n";
    }else if(n_l == 1 ){

      arma::colvec my = Y_l.t();
      m_ = (k_H * m_H + my * n_l) / k_;
      S_ =  S_H + (k_H* n_l)/(k_H+n_l) * (my-m_H)*((my-m_H).t());

    }else{
      m_ = m_H;
      S_ = S_H;
    }



    Sig_new.slice(l) = arma::iwishrnd( (S_), v_);

    mu_new.col(l)    = arma::mvnrnd(m_, Sig_new.slice(l) /k_);

  }

  return List::create(
    _["mu"]      = mu_new,
    _["Sigma"]   = Sig_new);
}













// [[Rcpp::export]]
List Update_Theta_cpp_TRAIN(arma::mat Y,
                            arma::mat AB,
                            arma::mat MU_g,  double k_g,
                            arma::cube SIGMA_g, double v_g,
                            int G, int p){


  arma::mat    Xg(p,G);
  arma::cube   S2g(p,p,G);
  arma::colvec alpha_lab = AB.col(0);

  for(int g=0; g<G; g++){
    arma::uvec ind_g = find( alpha_lab == (g+1) ) ;
    arma::mat  Y_g   = Y.rows(ind_g);
    int n_g          = (Y_g.n_rows);

    double k_ = n_g + k_g;
    double v_ = n_g + v_g;

    arma::colvec m_(p);
    arma::mat    S_(p,p);

    arma::colvec   m_g = MU_g.col(g);
    arma::mat      S_g = SIGMA_g.slice(g);

    if(n_g>1){
      arma::colvec my = mean(Y_g,0).t();
      m_ = (k_g * m_g + my * n_g) / k_;

      S_ =  S_g + cov(Y_g) * (n_g-1) +
        (k_g * n_g)/(k_g+n_g) * (my-m_g)*((my-m_g).t());

    }else if(n_g == 1 ){

      arma::colvec my = Y_g.t();
      m_ = (k_g * m_g + my * n_g) / k_;
      S_ =  S_g + (k_g* n_g)/(k_g+n_g) * (my-m_g)*((my-m_g).t());

    }else{
      m_ = m_g;
      S_ = S_g;
    }


    // Rcout << "ok4 \n";

    S2g.slice(g) = arma::iwishrnd( (S_), v_);
    // Rcout << "ok5 \n";

    Xg.col(g)    = arma::mvnrnd(m_, S2g.slice(g) /k_);
    // Rcout << "ok6 \n";

  }

  return List::create(
    _["mu"]      = Xg,
    _["Sigma"]   = S2g);
}



// [[Rcpp::export]]
arma::colvec Upd_Zeta_cpp(
    arma::mat Y,
    arma::mat mu_all,
    arma::cube Sigma_all,
    arma::colvec Uitilde,
    arma::colvec xitilde,
    arma::colvec log_pitilde,
    int G, int n,  int L_new,
    arma::colvec poss_lab){
  int toto = L_new+G;
  arma::colvec res_Z(n);  res_Z.fill(0);
  arma::colvec likel(toto);

  for(int ii=0; ii<n; ii++) {

    for(int h=0; h<(toto); h++){
      if(Uitilde[ii] < xitilde[h]){
        likel(h) =  dmvnrm_arma(Y.row(ii), mu_all.col(h).t() ,Sigma_all.slice(h),true) - log(xitilde[h]);
      }else{
        likel(h) = - arma::datum::inf;
      }
    }
    likel = log_pitilde + likel ;
    arma::colvec r = exp( likel - LogSumExp(likel) );
    res_Z[ii] = RcppArmadillo::sample(poss_lab, 1, TRUE, r)[0];

  }
  return(res_Z);
}

////////////// BRAND FUNCTIONAL////////

// [[Rcpp::export]]
arma::mat Upd_alphabeta_t_cpp(arma::mat Y,
                              arma::mat effe_train,
                              arma::mat effe_test,
                              arma::mat sigma_train,
                              arma::mat sigma_test,
                              arma::colvec pidir,
                              arma::colvec omega,
                              int G, int n,  int L, int t,
                              arma::colvec poss_lab){
  arma::mat res_AB(n,2);
  res_AB.fill(0);
  double pi0 = pidir(G);
  pidir.shed_row(G); // removed last element
  //Rcout << "ok1";
  arma::colvec likold(G);
  arma::colvec liknew(L);

  for(int  ii=0; ii<n; ii++) {
    for(int j=0; j<G ; j++){
      likold(j) = arma::accu(log(normpdf(Y.col(ii), effe_train.col(j), sqrt(sigma_train.col(j) ) )));
    }
    for(int l=0; l<L; l++){

      liknew(l) = arma::accu(log(normpdf(Y.col(ii), effe_test.col(l),  sqrt(sigma_test.col(l))  )));

    }
    //Rcout << "ok2";

    arma::colvec oldprob = log(pidir) + likold;
    arma::colvec newprob = log(pi0) + log(omega) + liknew;

    arma::colvec r = join_cols(oldprob, newprob);
    r = exp( r - max(r));
    int ind = RcppArmadillo::sample(poss_lab, 1, TRUE, r)[0];
    //Rcout << "ok3";

    if (ind <= G) {
      res_AB(ii,0) = ind;
    } else{
      res_AB(ii,1) = ind - G;
    }
  }
  return(res_AB);
}


// [[Rcpp::export]]
arma::colvec DN(arma::colvec X, arma::colvec M, arma::colvec S){
  return(arma::normpdf(X,M,S));
}


/*
 // [[Rcpp::export]]
 arma::mat Upd_alphabeta_t_cpp_SLICE(arma::mat Y,
 arma::mat effe_train,
 arma::mat effe_test,
 arma::mat sigma_train,
 arma::mat sigma_test,
 arma::colvec pidir,
 arma::colvec omega,
 arma::colvec Ui,
 arma::colvec xi,
 int G, int n,  int L,
 arma::colvec poss_lab){
 arma::mat res_AB(n,2);
 res_AB.fill(0);
 double pi0 = pidir(G);
 pidir.shed_row(G); // removed last element
 //Rcout << "ok1";
 arma::colvec likold(G);
 arma::colvec liknew(L);

 for(int ii=0; ii<n; ii++) {
 for(int j=0; j<G; j++){

 likold(j) = arma::accu(log(normpdf( Y.col(ii),
 effe_train.col(j),
 sqrt(sigma_train.col(j)))
 ));

 }
 for(int l=0; l<L; l++){
 if(Ui[ii] < xi[l]){
 liknew(l) = arma::accu(log(normpdf(Y.col(ii), effe_test.col(l), sqrt(sigma_test.col(l)))));
 }else{
 liknew(l) = -arma::datum::inf;
 }
 }

 arma::colvec oldprob = log(pidir) + likold;
 arma::colvec newprob = log(pi0)   + log(omega) + liknew;

 //  Rcout<< oldprob;
 //  Rcout<< newprob;
 arma::colvec r = join_cols(oldprob, newprob);
 r = exp( r - LogSumExp(r) );
 int ind = RcppArmadillo::sample(poss_lab, 1, TRUE, r)[0];
 //Rcout << "ok3";

 if (ind <= G) {
 res_AB(ii,0) = ind;
 } else{
 res_AB(ii,1) = ind - G;
 }
 }

 return(res_AB);
 }
 */

// [[Rcpp::export]]
arma::mat Update_effe_test_t_cpp(arma::mat Y, arma::mat alphabeta,
                                 int L_new, int t, arma::colvec tau_k, arma::colvec R,
                                 arma::mat sigma_test,arma::colvec ak){

  arma::mat effe_test_new(t,L_new);
  arma::colvec beta_lab = alphabeta.col(1);
  arma::colvec k_(t), m_(t), Sy(t);



  for(int l=0; l<L_new; l++){

    arma::uvec ind_l = find( beta_lab == (l+1) ) ;
    arma::mat  Y_l   = Y.cols(ind_l);

    arma::colvec tau_k_times_R = (tau_k[l]*R);

    int n_l          = (Y_l.n_cols);

    if(n_l>1){

      Sy  = arma::sum(Y_l,1);
      k_  = n_l * pow(sigma_test.col(l),-1) + 1/tau_k_times_R;
      //m_  = (Sy % pow(sigma_test.col(l),-1) + ak[l]/tau_k[l]) % pow(k_,-1);

      // post cap's review
      m_                    = (Sy % pow(sigma_test.col(l),-1) + ak[l]/tau_k_times_R) % pow(k_,-1);
      ////////////////////
      effe_test_new.col(l)  = m_ + as<arma::vec>(rnorm(t))%arma::sqrt(pow(k_,-1));


    }else if( n_l==1){

      Sy  = (Y_l);
      k_  = n_l * pow(sigma_test.col(l),-1) + 1/tau_k_times_R;
      //m_  = (Sy % pow(sigma_test.col(l),-1)+ ak[l]/tau_k[l]) % pow(k_,-1);
      m_  = (Sy % pow(sigma_test.col(l),-1)+ ak[l]/tau_k_times_R) % pow(k_,-1);
      effe_test_new.col(l)  = m_ + as<arma::vec>(rnorm(t))%arma::sqrt(pow(k_,-1));


    }else{
      effe_test_new.col(l) = ak[l] + as<arma::vec>(rnorm(t)) % arma::sqrt(tau_k_times_R);
    }

  }
  return(effe_test_new);
}



// [[Rcpp::export]]
arma::mat Update_effe_train_t_cpp(arma::mat Y, arma::mat alphabeta,
                                  int G, int t,
                                  arma::mat f_bar_g,
                                  arma::mat sigma_train,
                                  double KAPPAG){

  arma::mat    effe_train_new(t,G);
  arma::colvec alpha_lab = alphabeta.col(0);
  arma::colvec k_(t), m_(t), Sy(t);



  for(int g=0; g<G; g++){

    arma::uvec ind_g = find( alpha_lab == (g+1) ) ;
    arma::mat  Y_g   = Y.cols(ind_g);


    int n_g          = (Y_g.n_cols);

    if(n_g>1){

      Sy  = arma::sum(Y_g,1);
      k_  = n_g * pow( sigma_train.col(g), -1 ) + 1 / KAPPAG;

      m_  = (Sy % pow(sigma_train.col(g),-1) + f_bar_g.col(g)/KAPPAG) % pow(k_,-1);
      ////////////////////
      effe_train_new.col(g)  = m_ + as<arma::vec>(rnorm(t))%arma::sqrt(pow(k_,-1));


    }else if( n_g==1){

      Sy  = (Y_g);
      k_  = n_g * pow(sigma_train.col(g),-1) +  1 / KAPPAG;;
      m_  = (Sy % pow(sigma_train.col(g),-1) + f_bar_g.col(g)/KAPPAG) % pow(k_,-1);
      effe_train_new.col(g)  = m_ + as<arma::vec>(rnorm(t)) % arma::sqrt( pow( k_,- 1));


    }else{
      effe_train_new.col(g) = f_bar_g.col(g) + as<arma::vec>(rnorm(t)) * pow(KAPPAG,.5);
    }

  }
  return(effe_train_new);
}


// [[Rcpp::export]]
arma::mat Update_sigma_test_t_cpp(arma::mat Y, arma::mat alphabeta,
                                  arma::mat effe_test,
                                  int L_new, int t,
                                  double asig, double bsig){

  arma::mat sigma_test_new(t,L_new);
  arma::colvec beta_lab = alphabeta.col(1);

  for(int l=0; l<L_new; l++){

    arma::uvec ind_l = find( beta_lab == (l+1) ) ;
    arma::mat  Y_l   = Y.cols(ind_l);

    int n_l          = (Y_l.n_cols);


    if(n_l>1){

      arma::mat fake_FL = arma::repelem(effe_test.col(l),1,n_l);
      arma::mat diff2 = pow((Y_l - fake_FL),2);

      for( int small_t=0; small_t<t; small_t++){
        double B = accu( diff2.row(small_t) );
        double bpost  = .5  * B + bsig;
        double apost  = n_l * 1/2 + asig;
        sigma_test_new(small_t,l) = 1/ rgamma(1,apost,1/(bpost))[0];
      }

    }else if( n_l==1){

      //  Rcout << "2.0";

      arma::colvec diff2 = pow( (Y_l- effe_test.col(l)) ,2 ) ;
      //  Rcout << "2.1";

      for( int small_t=0; small_t<t; small_t++){

        double bpost  = .5  *  diff2(small_t) + bsig;
        //      Rcout << "2.2";

        double apost  = 1/2 + asig;
        sigma_test_new(small_t,l) = 1/ rgamma(1,apost,1/(bpost))[0];
      }
      //   Rcout << "2";

    }else{
      for( int small_t=0; small_t<t; small_t++){
        sigma_test_new(small_t,l) = 1/ rgamma(1,asig,1/(bsig))[0];
      }
      //   Rcout << "3";
    }
  }
  return(sigma_test_new);
}



// [[Rcpp::export]]
arma::mat Update_sigma_train_t_cpp(arma::mat Y, arma::mat alphabeta,
                                   arma::mat effe_train,
                                   int G, int t,
                                   arma::mat a_priorG, arma::mat b_priorG){

  arma::mat sigma_train_new(t,G);
  arma::colvec alpha_lab = alphabeta.col(0);

  for(int g=0; g<G; g++){

    arma::uvec ind_g = find( alpha_lab == (g+1) ) ;
    arma::mat  Y_g   = Y.cols(ind_g);

    int n_g          = (Y_g.n_cols);

    if(n_g>1){

      arma::mat fake_FG = arma::repelem(effe_train.col(g),1,n_g);
      arma::mat diff2 = pow((Y_g - fake_FG),2);

      for( int small_t=0; small_t<t; small_t++){
        double B = accu( diff2.row(small_t) );
        double bpost  = .5  * B   + b_priorG(small_t,g);
        double apost  = n_g * 1/2 + a_priorG(small_t,g);
        sigma_train_new(small_t,g) = 1/ rgamma(1,apost,1/(bpost))[0];
      }

    }else if( n_g==1){


      arma::colvec diff2 = pow( (Y_g - effe_train.col(g)) ,2 ) ;

      //  Rcout << "2.1";

      for( int small_t=0; small_t<t; small_t++){


        double B = accu( diff2.row(small_t) );
        double bpost  = .5  * B   + b_priorG(small_t,g);
        double apost  =  1/2 + a_priorG(small_t,g);
        sigma_train_new(small_t,g) = 1/ rgamma(1,apost,1/(bpost))[0];

      }
      //   Rcout << "2";

    }else{
      for( int small_t=0; small_t<t; small_t++){
        sigma_train_new(small_t,g) = 1/ rgamma(1, a_priorG(small_t,g),1/( b_priorG(small_t,g)))[0];
      }
      //   Rcout << "3";
    }
  }
  return(sigma_train_new);
}





// [[Rcpp::export]]
arma::colvec Update_tau_l(double a_tau, double b_tau,
                          arma::mat effe_test, int L_new, int t,
                          arma::colvec OOR, arma::colvec a_l){
  arma::colvec new_tau_l(L_new);
  double apost = a_tau + t*.5;
  int oldL = effe_test.n_cols;

  if(oldL>L_new){

    for(int l=0;l<L_new;l++){

      double bpost = 1/(accu( pow( (effe_test.col(l)-a_l[l]),2 ) % OOR )*.5 + b_tau);
      new_tau_l[l] =      1/ rgamma(1,apost,
                                    bpost)[0];
    }
  }else{
    for(int l=0;l<oldL;l++){

      double bpost = 1/( accu( pow( (effe_test.col(l)-a_l[l]),2 ) % OOR )*.5 + b_tau);
      new_tau_l[l] =      1/ rgamma(1,apost,
                                    bpost)[0];
    }
    for(int l=oldL;l<L_new;l++){
      new_tau_l[l] =      1/ rgamma(1,a_tau,
                                    b_tau)[0];
    }
  }
  return(new_tau_l);
}

// [[Rcpp::export]]
arma::colvec Update_a_l(arma::mat effe_test, int L_new,
                        arma::colvec tau_l,
                        arma::colvec R, double s, int t,
                        arma::colvec OneOverR){

  arma::colvec a_l_new(L_new);
  int oldL = effe_test.n_cols;
  if(oldL>L_new){

    for(int l=0;l<L_new;l++){
      double SumY        = accu(effe_test.col(l) % OneOverR);
      //Rcout << SumY;
      double k_          = accu(OneOverR)/(tau_l[l]) + 1/(s);
      //Rcout << k_;
      double m_          = ( SumY / (tau_l[l]) )/k_;
      a_l_new[l]   =   m_ + (rnorm(1))[0] * pow(1/k_ ,.5);
    }
  }else{
    for(int l=0;l<oldL;l++){
      double SumY        = accu(effe_test.col(l)%OneOverR);
      //Rcout << SumY;
      double k_          = accu(OneOverR)/(tau_l[l]) + 1/(s);
      //Rcout << k_;
      double m_          = ( SumY/(tau_l[l]) )/k_;
      a_l_new[l]   =   m_ + (rnorm(1))[0] * pow(1/k_ ,.5);
    }
    for(int l=oldL;l<L_new;l++){
      a_l_new[l]   =   0 + (rnorm(1))[0] * pow(s,.5);
    }

  }
  return(a_l_new);
}











// [[Rcpp::export]]
arma::mat Upd_ZETA_t_SLICE(arma::mat Y,
                           arma::mat effe_ALL,
                           arma::mat sigma_ALL,
                           arma::colvec log_pitilde,
                           arma::colvec Uitilde,
                           arma::colvec xitilde,
                           int G, int n,  int L_new,
                           arma::colvec poss_lab){
  int toto = L_new+G;
  arma::colvec res_Z(n);  res_Z.fill(0);
  arma::colvec likel(toto);

  for(int ii=0; ii<n; ii++) {

    for(int l=0; l<toto; l++){
      if(Uitilde[ii] < xitilde[l]){
        likel(l) =  arma::accu(log(normpdf(Y.col(ii), effe_ALL.col(l), sqrt(sigma_ALL.col(l))))) - log(xitilde[l]);
      }else{
        likel(l) =  -arma::datum::inf;
      }
    }

    likel = log_pitilde + likel;

    arma::colvec r = exp( likel - LogSumExp(likel) );
    res_Z[ii] = RcppArmadillo::sample(poss_lab, 1, TRUE, r)[0];

  }

  return(res_Z);
}

