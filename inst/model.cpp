// Space time 
#include <TMB.hpp>
//#include <atomic_math.hpp>  //for D_lgamma


// Variance
template<class Type>
Type var( array<Type> vec ){
  Type vec_mod = vec - (vec.sum()/vec.size());
  Type res = pow(vec_mod, 2.0).sum() / vec.size();
  return res;
}

// square
template<class Type>
Type square(Type x){
  return pow(x,2.0);
}

// sqrt
template<class Type>
Type sqrt(Type x){
  return pow(x,0.5);
}

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}



// Bernoulli
template<class Type>
Type dbern(Type x, Type prob, int give_log=false){
  Type logres;
  if( x==0 ) logres = log( 1-prob );
  if( x==1 ) logres = log( prob );
  if(give_log) return logres; else return exp(logres);
}


// dlnorm -- carefull log(x) is important
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  //return 1/(sqrt(2*M_PI)*sd) * exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// dzinflognorm
template<class Type>
Type dzinflognorm(Type x, Type meanlog, Type encounter_prob, Type log_notencounter_prob, Type sdlog, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dlnorm( x, meanlog - square(sdlog)/2, sdlog, false );
    if(give_log==true) Return = log(encounter_prob) + dlnorm( x, meanlog - square(sdlog)/2, sdlog, true );
  } 
  return Return;
}


// dzinfgamma, shape = 1/CV^2, scale = mean*CV^2
template<class Type>
Type dzinfgamma(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dgamma( x, pow(cv,-2), posmean*pow(cv,2), false );
    if(give_log==true) Return = log(encounter_prob) + dgamma( x, pow(cv,-2), posmean*pow(cv,2), true );
  } 
  return Return;
}


// dzinfnorm
template<class Type>
Type dzinfnorm(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dnorm( x, posmean, posmean*cv, false );
    if(give_log==true) Return = log(encounter_prob) + dnorm( x, posmean, posmean*cv, true );
  }
  return Return;
}


template<class Type>
Type objective_function<Type>::operator() (){
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;
  
  
  ///////////////////////////////////
  // Model configurations and indices
  ///////////////////////////////////
  
  DATA_FACTOR( Options_vec );
  // Slot 0 - Prior (DEPRECATED): prior on random effects (0=SPDE_GMRF; 1=ICAR_GMRF)
  // Slot 1 - Alpha: Alpha (neighborhood)
  // Slot 2 - IncludeDelta (DEPRECATED): Include hyperdistribution for delta
  // Slot 3 - IncludeEta (DEPRECATED): Include hyperdistribution for eta
  // Slot 4 - SE: Output S_x, total_abundance in ADREPORT (0=no, 1=yes); needed for unbiased map estimates using bias.correct
  // Slot 5 - DataSource: Data sources feeding the model (1: scientific_commercial, 2: scientific_only, 3: commercial_only)
  // Slot 6 - DataObs: Observation model (1: zero-inflated gamma, 2: zero-inflated lognormal, 3: lognormal)
  // Slot 7 - SamplingProcess: Sampling of commercial data contributes to likelihood (0: no, 1: yes)
  // Slot 8 - zero.infl_model (DEPRECATED): Parameterization of zero-inflation
  // Slot 9 - commercial_obs (DEPRECATED): Commercial observations contributes to likelihood (O: no, 1: yes)
  // Slot 10 - b_constraint: Constraint on preferential sampling parameter (b_constraint=1: b > 0 by specifying exp(b), b_constraint=2: b left free)
  // Slot 11 - catchability_random: Relative catchability parameters are considered as a random effect (0: no, 1: yes)
  // Slot 12 - cov_samp_process: Account for covariates in the sampling process (0: no, 1: yes)
  // Slot 13 - const_kappa (DEPRECATED): kappa parameters is constant in time (0: time varying, 1: constant)
  // Slot 14 - const_target (MODIFY NAME): PS parameter is considered as a Gaussian random effect (0:no, fixed effect; 1: yes, Gaussian random effect)
  // Slot 15 - const_q1 (DEPRECATED): zero-inflation parameter q1 is constant over time (0: time varying, 1:constant)
  // Slot 16 - const_k (DEPRECATED): constant catchability over time (0: time varying, 1: constant)
  // Slot 17 const_Sigma (DEPRECATED): constant variance of observation over time (0: time varying, 1: constant)
  // Slot 18 - ref_data (DEPRECATED): old way to define which data source is the reference level
  // Slot 19 - const_tau (DEPRECATED): constant marginal variance of the random effects over time (0: time varying, 1: constant)
  // Slot 20 - const_spphab: constant species-habitat relationship over time (0: time varying, 1: constant)
  // Slot 21 - biomassAR1: temporal AR1 process in the latent field of biomass (0: no temporal correlation, 1: AR1 temporal correlation)
  // Slot 22 - sampAR1: temporal AR1 process in the sampling intensity (0: no temporal correlation, 1: AR1 temporal correlation)
  // Slot 23 - anisotropy: the spatial random effects account for spatial anisotropy (0: no, 1: yes)
  
  // Note: lot of deprecated stuff to keep the model relatively parcimonious
  // But I might come back to those ideas, so I left these options
  
  // Indices
  DATA_INTEGER(n_x);
  DATA_INTEGER(n_t);
  DATA_INTEGER(n_S);
  DATA_INTEGER(n_p);
  DATA_INTEGER(n_eta);
  int t,i,x,p;
  
  
  ///////////////////////////////////////////////////////////////////
  // ---------------------- Biomass field -------------------------//
  ///////////////////////////////////////////////////////////////////
  
  
  ////////
  // Data
  ////////
  
  // SPDE objects
  DATA_STRUCT(spde,spde_t); // standard spde objects
  DATA_STRUCT(spde_aniso,spde_aniso_t); // anisotropic spde objects

  // Prediction stuff
  DATA_MATRIX( Cov_x_pred ); // Covariates
  DATA_IMATRIX(Aix_ij_pred); // Matrix and weight vector to relate mesh nodes to datapoints
  DATA_VECTOR(Aix_w_pred);
  
  
  /////////////
  // Parameters (random and fixed for the latent field only)
  /////////////
  // Intercepts
  PARAMETER_MATRIX(beta_j);
  PARAMETER(beta_j0);
  PARAMETER_VECTOR(beta_j0year);
  PARAMETER_VECTOR(beta_j0season);
  
  // Random effects
  PARAMETER_ARRAY(deltainput_x);
  PARAMETER_ARRAY(epsiloninput_x);
  PARAMETER_VECTOR(logkappa_S);
  PARAMETER_VECTOR(logta_S);
  PARAMETER_VECTOR(logtaAR1_S);
  PARAMETER(rho_epsilon);
  PARAMETER_VECTOR(ln_Hdelta_input);
  
  
  /////////////////
  // derived values
  /////////////////
  // likelihood
  vector<Type> jnll_comp(5);
  jnll_comp.setZero();
  
  // random effect matrices
  matrix<Type> delta_x( n_x , n_t );
  matrix<Type> epsilon_x( n_x , n_t );
  array<Type> delta_p(n_p, n_t);
  array<Type> epsilon_p(n_p, n_t);
  
  delta_x.setZero();
  epsilon_x.setZero();
  delta_p.setZero();
  epsilon_p.setZero();
  
  
  // parameters definition
  vector<Type> logtau_S( n_S );
  vector<Type> logtauAR1_S( n_S );
  vector<Type> MargSD_S( n_S );
  vector<Type> MargSDAR1_S( n_S );
  vector<Type> Range_S( n_S );
  
  for( int s=0; s<n_S; s++){
    
    logtau_S(s) = logta_S(s) - logkappa_S(s);
    logtauAR1_S(s) = logtaAR1_S(s) - logkappa_S(s);
    Range_S(s) = sqrt(8) / exp(logkappa_S(s));
    MargSD_S(s) = 1 / sqrt(4*M_PI) / exp(logtau_S(s)) / exp(logkappa_S(s));
    MargSDAR1_S(s) = 1 / sqrt(4*M_PI) / exp(logtauAR1_S(s)) / exp(logkappa_S(s));
    
  }
  
  // Anisotropy elements
  matrix<Type> H_delta(2,2);
  H_delta(0,0) = exp(ln_Hdelta_input(0));
  H_delta(1,0) = ln_Hdelta_input(1);
  H_delta(0,1) = ln_Hdelta_input(1);
  H_delta(1,1) = (1+ln_Hdelta_input(1)*ln_Hdelta_input(1)) / exp(ln_Hdelta_input(0));
  
  // introduce precision matrix
  SparseMatrix<Type> Q;
  
  for(t=0; t<n_t; t++){
    
    // transform random effect
    for( int x=0; x<n_x; x++){
      
      if(Options_vec(21)==0) delta_x(x,t) = deltainput_x(x,t) / exp(logtau_S(0));
      if(Options_vec(21)==1) epsilon_x(x,t) = epsiloninput_x(x,t) / exp(logtauAR1_S(0));
      
    }
    
    // Precision matrices (perhaps the 2 lines below could be moved outside the loop, should try...)
    if(Options_vec(23)==1) Q = Q_spde(spde_aniso, exp(logkappa_S(0)), H_delta );
    if(Options_vec(23)==0) Q = Q_spde(spde, exp(logkappa_S(0)) );

    if(Options_vec(2)==1){
    
      if(Options_vec(21)==0) jnll_comp(3) += GMRF(Q)(deltainput_x.col(t));
      if(Options_vec(21)==1 & t == 0) jnll_comp(3) +=  GMRF(Q)(epsiloninput_x.col(t));
      if(Options_vec(21)==1 & t > 0) jnll_comp(3) += GMRF(Q)(epsiloninput_x.col(t) - rho_epsilon*epsiloninput_x.col(t-1));
      
    }
    
    // Projection for plotting
    for( int Arow=0; Arow<Aix_ij_pred.rows(); Arow++ ){
      
      int p = Aix_ij_pred(Arow,0);
      int x = Aix_ij_pred(Arow,1);
      if(Options_vec(21)==0) delta_p(p,t) += Aix_w_pred(Arow) * delta_x(x,t);
      if(Options_vec(21)==1) epsilon_p(p,t) += Aix_w_pred(Arow) * epsilon_x(x,t);
    
    }
    
  }
  
  //-------------------------------
  // Prediction of the latent field
  //-------------------------------
  
  matrix<Type> S_p(n_p,n_t);
  matrix<Type> linpredS_p(n_p,n_t);
  linpredS_p.setZero();
  
  vector<Type> beta_j_t(Cov_x_pred.row(0).size());
  Type sum_beta_j;
  
  for(t=0; t<n_t; t++){
    
    if(Options_vec(20)==1){
      beta_j_t = beta_j.row(0);
    }else{
      beta_j_t = beta_j.row(t);
    }
    linpredS_p.col(t) = Cov_x_pred * beta_j_t;
    
  }
  
  for(int t=0; t<n_t; t++){
    for(int p=0; p<n_p; p++){
      if(Options_vec(24)==0) S_p(p,t) = exp(beta_j0 + beta_j0year(t) + beta_j0season(t) +  linpredS_p(p,t) + delta_p(p,t) + epsilon_p(p,t) );
      if(Options_vec(24)==1) S_p(p,t) = plogis(beta_j0 + beta_j0year(t) + beta_j0season(t) +  linpredS_p(p,t) + delta_p(p,t) + epsilon_p(p,t) );
    }
  }
  
  
  ///////////////////////////////////////////////////////////////////
  // ---------------------- Commercial data -----------------------//
  ///////////////////////////////////////////////////////////////////
  
  if(Options_vec(5) == 3 | Options_vec(5) == 1){
    
    ////////
    // Data
    ////////
    
    // indices
    DATA_INTEGER(n_com_i);
    DATA_INTEGER(n_ipp);
    
    // commercial observations stuff
    DATA_MATRIX( Cov_x_com );
    DATA_ARRAY( c_com_x );       	// Response (count) for each observation i (commercial data)
    DATA_VECTOR( y_com_i );       	// Response (0:not surveyed, 1:surveyed) for each site (commercial data)
    DATA_FACTOR( b_com_i );
    DATA_FACTOR( t_com_i );
    DATA_VECTOR( q2_com ); 
    DATA_VECTOR( weights_com );
    DATA_MATRIX( Cov_fb );  //design matrix for sampling of commercial data
    DATA_IMATRIX(Aix_ij_com);
    DATA_VECTOR(Aix_w_com);
    
    // Point process stuff
    DATA_IMATRIX(Aix_ij_ipp);
    DATA_VECTOR( Aix_w_ipp );
    DATA_MATRIX( Cov_x_ipp );
    DATA_MATRIX( Cov_fb_ipp );
    DATA_MATRIX( Cov_fb_pred );
    DATA_VECTOR(W);
    
    //////////////
    // Parameters 
    //////////////
    PARAMETER_VECTOR( q1_com );
    PARAMETER_VECTOR( logSigma_com );
    PARAMETER_VECTOR( k_com );
    PARAMETER_VECTOR( logSigma_catch );
    PARAMETER_VECTOR( logMean_catch );

    /////////////////
    // derived values
    /////////////////
    vector<Type> delta_com_i(n_com_i);
    vector<Type> epsilon_com_i(n_com_i);
    vector<Type> S_com_i(n_com_i);
    delta_com_i.setZero();
    epsilon_com_i.setZero();
    Type k_com_i;

    // Projection from mesh nodes to commercial data points
    for( int Arow=0; Arow<Aix_ij_com.rows(); Arow++ ){

      i = Aix_ij_com(Arow,0);
      x = Aix_ij_com(Arow,1);
      if(Options_vec(21)==0) delta_com_i(i) += Aix_w_com(Arow) * delta_x(x,t_com_i(i));
      if(Options_vec(21)==1) epsilon_com_i(i) += Aix_w_com(Arow) * epsilon_x(x,t_com_i(i));

    }

    ////////////////
    // Latent field
    ///////////////
    matrix<Type> linpredS_com(n_com_i,n_t);
    linpredS_com.setZero();

    for(t=0; t<n_t; t++){

      if(Options_vec(20)==1){
        beta_j_t = beta_j.row(0);
      }else{
        beta_j_t = beta_j.row(t);
      }
      linpredS_com.col(t) = Cov_x_com * beta_j_t;

    }

    for(int i=0; i<n_com_i; i++){

      if(Options_vec(24)==0) S_com_i(i) = exp(beta_j0 + beta_j0year(t_com_i(i)) + beta_j0season(t_com_i(i)) +  linpredS_com(i,t_com_i(i)) + delta_com_i(i) + epsilon_com_i(i)); //
      if(Options_vec(24)==1){
        k_com_i = k_com(b_com_i(i));
        S_com_i(i) = plogis(k_com_i + beta_j0 + beta_j0year(t_com_i(i)) + beta_j0season(t_com_i(i)) +  linpredS_com(i,t_com_i(i)) + delta_com_i(i) + epsilon_com_i(i)); //
      }

    }

    /////////////////////
    // Observation model
    ////////////////////

    if(Options_vec(24)==0){

      // Zero-inflated lognormal (log-link Poisson)
      vector<Type> Sigma_catch = exp(logSigma_catch);
      vector<Type> Sigma_com = exp(logSigma_com);
      vector<Type> q1_com_i(n_com_i);
      Type Sigma_com_i;
      vector<Type>  E_com(n_com_i);
      vector<Type>  encounterprob_com(n_com_i);
      vector<Type>  log_notencounterprob_com(n_com_i);

      for(int i=0; i<n_com_i; i++){

        if( !isNA(y_com_i(i)) ){

          q1_com_i(i) = q1_com(b_com_i(i));
          Sigma_com_i = exp(logSigma_com(b_com_i(i)));
          k_com_i = k_com(b_com_i(i));

          E_com(i) = q2_com(0) * k_com_i * S_com_i(i);

          encounterprob_com(i) = ( 1.0 - exp(-1 * E_com(i) * exp(q1_com_i(i))) );
          log_notencounterprob_com(i) = -1 * E_com(i) * exp(q1_com_i(i));

          if( Options_vec(6)==1 ) jnll_comp(1) -= weights_com(0) * dzinfgamma(y_com_i(i), E_com(i)/encounterprob_com(i), encounterprob_com(i), log_notencounterprob_com(i), Sigma_com_i, true);
          if( Options_vec(6)==2 ) jnll_comp(1) -= weights_com(0) * dzinflognorm(y_com_i(i), log(E_com(i))-log(encounterprob_com(i)), encounterprob_com(i), log_notencounterprob_com(i), Sigma_com_i, true);
          if( Options_vec(6)==3 & y_com_i(i) > 0 ) jnll_comp(1) -= weights_com(0) * dlnorm(y_com_i(i), log(E_com(i)), Sigma_com_i, true);

        }

      }


      if(Options_vec(11) == 1) jnll_comp(1) -= dnorm(log(k_com),logMean_catch(0),Sigma_catch(0),true).sum();

      REPORT( q1_com );
      REPORT( Sigma_com );
      REPORT( encounterprob_com );
      REPORT( log_notencounterprob_com );
      REPORT( E_com );

    }
    
    if(Options_vec(24)==1){


      for(int i=0; i<n_com_i; i++){

        if( !isNA(y_com_i(i)) ){

          // Bernoulli model
          // Type k_com_i;
          Type E_com;
          // k_com_i = k_com(b_com_i(i));
          E_com = q2_com(0) * S_com_i(i); // k_com_i *
          jnll_comp(1) -= weights_com(0) * dbern(y_com_i(i), E_com, true);

        }
      }
    }

    REPORT( delta_com_i );
    REPORT( S_com_i );
    REPORT( epsilon_com_i );

    //////////////////////////////////////
    // Sampling process of commercial data
    /////////////////////////////////////
    
    if( Options_vec(7) == 1){
      
      // Parameters
      PARAMETER_MATRIX(beta_fb);
      PARAMETER_VECTOR(beta_fb0);
      PARAMETER_MATRIX(beta_fb0year);
      PARAMETER_MATRIX(beta_fb0season);
      PARAMETER_MATRIX( par_b );
      PARAMETER_VECTOR(logSigma_targ);
      PARAMETER_VECTOR(Mean_targ);
      PARAMETER_MATRIX( par_byear );
      PARAMETER_VECTOR(logSigma_targyear);
      PARAMETER_MATRIX( par_bseason );
      PARAMETER_VECTOR(logSigma_targseason);
      PARAMETER_MATRIX( par_bseasonyear );
      PARAMETER_VECTOR( logSigma_targseasonyear );
      PARAMETER_ARRAY( etainput_x );
      PARAMETER_ARRAY( psiinput_x );
      PARAMETER(rho_psi);
      PARAMETER(rho_b);
      PARAMETER_VECTOR(ln_Heta_input);
      
      // Derived values
      Type b;
      array<Type> eta_x( n_x , n_t , n_eta );
      vector<Type> etainput2_x(n_x);

      array<Type> psi_x( n_x , n_t , n_eta );
      vector<Type> psiinput2_x(n_x);
      vector<Type> psiinput3_x(n_x);

      vector<Type> psimean_x( n_x );
      vector<Type> eta_i( n_com_i );
      vector<Type> psi_i( n_com_i );
      array<Type> eta_p( n_p, n_t , n_eta );
      array<Type> psi_p( n_p, n_t , n_eta );
      array<Type> eta_ipp( n_ipp, n_t , n_eta );
      array<Type> psi_ipp( n_ipp, n_t , n_eta );

      vector<Type> linpredR_i( n_com_i );
      vector<Type> lambda_i( n_com_i );
      vector<Type> log_lambda_i( n_com_i );
      eta_x.setZero();
      psi_x.setZero();
      eta_i.setZero();
      psi_i.setZero();
      eta_p.setZero();
      psi_p.setZero();
      eta_ipp.setZero();
      psi_ipp.setZero();
      Type fact_S;

      // Anisotropy elements
      matrix<Type> H_eta(2,2);
      H_eta(0,0) = exp(ln_Heta_input(0));
      H_eta(1,0) = ln_Heta_input(1);
      H_eta(0,1) = ln_Heta_input(1);
      H_eta(1,1) = (1+ln_Heta_input(1)*ln_Heta_input(1)) / exp(ln_Heta_input(0));

      // Latent field
      for(t=0; t<n_t; t++){

        for( int f=0; f<n_eta; f++){

          if( Options_vec(0)==0 ){

            for(int s=0; s<n_x; s++){

              eta_x(s,t,f) = etainput_x(s,t,f) / exp(logtau_S(f+1));
              etainput2_x(s) = etainput_x(s,t,f);

              // Auto regression of random spatial effect
              if(Options_vec(22)==1){

                psiinput2_x(s) = psiinput_x(s,t,f);
                psiinput3_x(s) = psiinput_x(s,t-1,f);
                psi_x(s,t,f) = psiinput_x(s,t,f) / exp(logtauAR1_S(f+1));

              }

            }

            if(Options_vec(23)==1) Q = Q_spde(spde_aniso, exp(logkappa_S(f+1)), H_eta );
            if(Options_vec(23)==0) Q = Q_spde(spde, exp(logkappa_S(f+1)) );

            if(Options_vec(3)==1) jnll_comp(4) += GMRF(Q)(etainput2_x); // latent field for sampling
            if(Options_vec(22)==1){
              if(Options_vec(23)==1) Q = Q_spde(spde_aniso, exp(logkappa_S(f+1)), H_eta );
              if(Options_vec(23)==0) Q = Q_spde(spde, exp(logkappa_S(f+1)) );
              if(t == 0) jnll_comp(4) += GMRF(Q)(psiinput2_x); // latent field for sampling
              if(t > 0) jnll_comp(4) += GMRF(Q)(psiinput2_x - rho_psi*psiinput3_x); // latent field for sampling
            }

          }

        }

      }

      // Projection for plotting
      for(t=0; t<n_t; t++){
        for( int Arow=0; Arow<Aix_ij_pred.rows(); Arow++ ){
          for( int f=0; f<n_eta; f++){
            int p = Aix_ij_pred(Arow,0);
            int x = Aix_ij_pred(Arow,1);
            eta_p(p,t,f) += Aix_w_pred(Arow) * eta_x(x,t,f);
            if(Options_vec(22)==1) psi_p(p,t,f) += Aix_w_pred(Arow) * psi_x(x,t,f);
          }
        }
      }

      for( int Arow=0; Arow<Aix_ij_com.rows(); Arow++ ){
        i = Aix_ij_com(Arow,0);
        x = Aix_ij_com(Arow,1);
        eta_i(i) += Aix_w_com(Arow) * eta_x(x,t_com_i(i),b_com_i(i));
        if(Options_vec(22)==1) psi_i(i) += Aix_w_com(Arow) * psi_x(x,t_com_i(i),b_com_i(i));
      }

      for(t=0; t<n_t; t++){
        for( int Arow=0; Arow<Aix_ij_ipp.rows(); Arow++ ){
          for( int f=0; f<n_eta; f++){
            int p = Aix_ij_ipp(Arow,0);
            int x = Aix_ij_ipp(Arow,1);
            eta_ipp(p,t,f) += Aix_w_ipp(Arow) * eta_x(x,t,f);
            if(Options_vec(22)==1) psi_ipp(p,t,f) += Aix_w_ipp(Arow) * psi_x(x,t,f);
          }
        }
      }

      for(int i=0; i<n_com_i; i++){

        if(Options_vec(12)==1){
          linpredR_i = Cov_fb * beta_fb(b_com_i(i));
        }else{
          linpredR_i(i) = 0;
        }

        if(Options_vec(10)==1){

          b=exp(par_b(t_com_i(i),b_com_i(i))) + exp(par_bseason(t_com_i(i),b_com_i(i))) + exp(par_byear(t_com_i(i),b_com_i(i))) + exp(par_bseasonyear(t_com_i(i),b_com_i(i)));

        }

        if(Options_vec(10)==2){

          b=par_b(t_com_i(i),b_com_i(i)) + par_bseason(t_com_i(i),b_com_i(i)) + par_byear(t_com_i(i),b_com_i(i)) + par_bseasonyear(t_com_i(i),b_com_i(i));

        }

        if(Options_vec(24)==0) fact_S = b*(log(S_com_i(i)));
        if(Options_vec(24)==1) fact_S = b*(logit(S_com_i(i)));

        log_lambda_i(i) = beta_fb0(b_com_i(i)) + beta_fb0year(t_com_i(i),b_com_i(i)) + beta_fb0season(t_com_i(i),b_com_i(i)) + fact_S + linpredR_i(i) + eta_i(i) +  psi_i(i); //
        jnll_comp(2) -= weights_com(0) * ( log_lambda_i(i) ); // log density of poisson process (see Diggle (2013))          }
        // For standardisation constant, see below

      }


      /////////////
      // Prediction
      /////////////
      array<Type> lambda_p( n_p , n_t , n_eta );
      matrix<Type> linpredS_p( n_p , n_t );
      matrix<Type> linpredR_p( n_p , n_eta );

      linpredS_p.setZero();

      matrix<Type> par_bAR1( n_t , n_eta );
      vector<Type> vec_par_b(n_t);
      vector<Type> vec_par_byear(n_t);
      vector<Type> vec_par_bseason(n_t);
      vector<Type> vec_par_bAR1(n_t);

      for(t=0; t<n_t; t++){
        if(Options_vec(20)==1){
          beta_j_t = beta_j.row(0);
        }else{
          beta_j_t = beta_j.row(t);
        }
        linpredS_p.col(t) = Cov_x_pred * beta_j_t;
      }

      for( int f=0; f<n_eta; f++){

        linpredR_p.col(f) = Cov_fb_pred * beta_fb.col(f);

        for(int t=0; t<n_t; t++){

          for(int p=0; p<n_p; p++){

            S_p(p,t) = exp(beta_j0 + beta_j0year(t) + beta_j0season(t) + linpredS_p(p,t) + delta_p(p,t) + epsilon_p(p,t));

            if(Options_vec(12)!=1) linpredR_p(p,f) = 0;

            if(Options_vec(10)==1){

              b = exp(par_b(t,f)) + exp(par_bseason(t,f)) + exp(par_byear(t,f)) + exp(par_bseasonyear(t,f));

            }

            if(Options_vec(10)==2){

              b = par_b(t,f) + par_bseason(t,f) + par_byear(t,f) + par_bseasonyear(t,f);

            }

            if(Options_vec(24)==0) fact_S = b*(log(S_com_i(i)));
            if(Options_vec(24)==1) fact_S = b*(logit(S_com_i(i)));

            lambda_p(p,t,f) = exp(  beta_fb0(f) + beta_fb0year(t,f) + beta_fb0season(t,f) + fact_S + linpredR_p(p) + eta_p(p,t,f) + psi_p(p,t,f));
            jnll_comp(2) -= Type(1)-lambda_p(p,t,f);

          }

          if(t==0) par_bAR1(t,f) = par_bseasonyear(t,f);
          if(t>0) par_bAR1(t,f) = par_bseasonyear(t,f) - rho_b * par_bseasonyear(t-1,f);

        }

        if(Options_vec(14) == 0){

          vec_par_b = par_b.col(f);
          vec_par_byear = par_byear.col(f);
          vec_par_bseason = par_bseason.col(f);
          vec_par_bAR1 = par_bAR1.col(f);

          jnll_comp(2) -= dnorm(vec_par_b,Mean_targ(f),exp(logSigma_targ(0)),true).sum();
          jnll_comp(2) -= dnorm(vec_par_byear,0,exp(logSigma_targyear(0)),true).sum();
          jnll_comp(2) -= dnorm(vec_par_bseason,0,exp(logSigma_targseason(0)),true).sum();
          jnll_comp(2) -= dnorm(vec_par_bAR1,0,exp(logSigma_targseasonyear(0)),true).sum();

        }

      }

      // Outputs
      REPORT( lambda_p );
      REPORT(beta_fb0);
      REPORT(beta_fb0year);
      REPORT(beta_fb0season);
      REPORT(linpredR_p);
      REPORT( beta_fb );
      REPORT( par_b );
      REPORT( par_bseason );
      REPORT( par_byear );
      REPORT( par_bseasonyear );
      REPORT( eta_p );
      REPORT( psi_p );
      
    }
    
  }
  
  //////////////////////////////////// Scientific /////////////////////////////////////////
  
  if(Options_vec(5) == 2 | Options_vec(5) == 1){
    
    ////////
    // Data
    ////////
    
    DATA_INTEGER(n_sci_i);
    DATA_MATRIX( Cov_x_sci );
    DATA_VECTOR( y_sci_i );
    DATA_VECTOR( q2_sci );
    DATA_FACTOR( t_sci_i );
    
    //////////////
    // Parameters 
    //////////////
    PARAMETER(logSigma_sci);
    PARAMETER_VECTOR(q1_sci);
    PARAMETER_VECTOR( k_sci );
    DATA_IMATRIX(Aix_ij_sci);
    DATA_VECTOR(Aix_w_sci);
    
    /////////////////
    // derived values
    /////////////////

    // Delta
    vector<Type> delta_sci_i(n_sci_i);
    vector<Type> epsilon_sci_i(n_sci_i);

    delta_sci_i.setZero();
    epsilon_sci_i.setZero();

    Type k_sci_i;

    // Predicted densities
    vector<Type> S_sci_i( n_sci_i );
    for( int Arow=0; Arow<Aix_ij_sci.rows(); Arow++ ){

      i = Aix_ij_sci(Arow,0);
      x = Aix_ij_sci(Arow,1);
      if(Options_vec(21)==0) delta_sci_i(i) += Aix_w_sci(Arow) * delta_x(x,t_sci_i(i));
      if(Options_vec(21)==1) epsilon_sci_i(i) += Aix_w_sci(Arow) * epsilon_x(x,t_sci_i(i));

    }

    ////////////////
    // Latent field
    ///////////////
    matrix<Type> linpredS_sci(n_sci_i,n_t);
    linpredS_sci.setZero();
    for(t=0; t<n_t; t++){

      if(Options_vec(20)==1){
        beta_j_t = beta_j.row(0);
      }else{
        beta_j_t = beta_j.row(t);
      }
      linpredS_sci.col(t) = Cov_x_sci * beta_j_t;

    }

    for(int i=0; i<n_sci_i; i++){

      if(Options_vec(24)==0) S_sci_i(i) = exp(beta_j0 + beta_j0year(t_sci_i(i)) + beta_j0season(t_sci_i(i)) + linpredS_sci(i,t_sci_i(i)) + delta_sci_i(i) + epsilon_sci_i(i)); // beta_j0intra(t_sci_i(i)) +
      if(Options_vec(24)==1){
        k_sci_i = k_sci(0);
        S_sci_i(i) = plogis(k_sci_i + beta_j0 + beta_j0year(t_sci_i(i)) + beta_j0season(t_sci_i(i)) + linpredS_sci(i,t_sci_i(i)) + delta_sci_i(i) + epsilon_sci_i(i)); // beta_j0intra(t_sci_i(i)) +
      }

    }


    /////////////////////
    // Observation model
    ////////////////////

    if(Options_vec(24)==0){

      // Zero-inflated lognormal (log-link Poisson)
      Type Sigma_sci = exp(logSigma_sci);
      vector<Type>  E_sci(n_sci_i);
      vector<Type>  encounterprob_sci(n_sci_i);
      vector<Type>  log_notencounterprob_sci(n_sci_i);

      for(int i=0; i<n_sci_i; i++){
        if( !isNA(y_sci_i(i)) ){

          E_sci(i) = k_sci(0) * q2_sci(0) * S_sci_i(i);

          encounterprob_sci(i) = ( 1.0 - exp(-1 * E_sci(i) * exp(q1_sci(0)) ));
          log_notencounterprob_sci(i) = -1 * E_sci(i) * exp(q1_sci(0));

          if( Options_vec(6)==1 ) jnll_comp(0) -= dzinfgamma(y_sci_i(i), E_sci(i)/encounterprob_sci(i), encounterprob_sci(i), log_notencounterprob_sci(i), Sigma_sci, true);
          if( Options_vec(6)==2 ) jnll_comp(0) -= dzinflognorm(y_sci_i(i), log(E_sci(i))-log(encounterprob_sci(i)), encounterprob_sci(i), log_notencounterprob_sci(i), Sigma_sci, true);
          if( Options_vec(6)==3 & y_sci_i(i) > 0 ) jnll_comp(0) -= dlnorm(y_sci_i(i), log(E_sci(i)), Sigma_sci, true);

        }
      }

      //////////
      // Outputs
      //////////
      REPORT( q1_sci );
      REPORT( Sigma_sci );

      REPORT( E_sci );
      REPORT( encounterprob_sci );
      REPORT( log_notencounterprob_sci );

    }


    if(Options_vec(24)==1){

      for(int i=0; i<n_sci_i; i++){
        if( !isNA(y_sci_i(i)) ){

          // Bernoulli model
          // Type k_sci_i;
          Type E_sci;

          // k_sci_i = k_sci(0);
          E_sci = q2_sci(0) * S_sci_i(i); // k_sci_i *
          jnll_comp(0) -= dbern(y_sci_i(i), E_sci, true);

        }
      }

    }

    REPORT( delta_sci_i);
    REPORT( S_sci_i);
    REPORT( epsilon_sci_i );
    
  }
  
  // Total objective
  Type jnll;
  int n_jnll = jnll_comp.size();
  for(int i=0; i<n_jnll; i++){
    jnll += jnll_comp(i);
  }
  
  vector<Type> total_abundance(n_t);
  for(int t=0; t<n_t; t++){
    total_abundance(t) = S_p.col(t).sum();
  }
  
  // Reporting
  REPORT( S_p );
  REPORT( total_abundance );
  ADREPORT( total_abundance );
  REPORT( Range_S );
  REPORT( MargSD_S );
  REPORT( MargSDAR1_S );
  REPORT( beta_j );
  REPORT( delta_p );
  REPORT( delta_x );
  REPORT( epsilon_p );
  REPORT( epsilon_x );
  
  REPORT( jnll_comp );
  
  if(Options_vec(4)==1){
    ADREPORT(S_p);
  }
  
  
  return jnll;
}
