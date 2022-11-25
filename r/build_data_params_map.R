#####################################################################
## Build input data for fitting the model: 'Data', 'Params' and 'Map'
#####################################################################
# Baptiste Alglave

# Option_vec : model configuration
# Data : list of data (catch data, covariates, Option_vec, spde objects, etc.)
# Params : list of parameters
# Map : parameters that are fixed to a specific value

Options_vec = c( 'Prior'=0, # DELETE
                 'Alpha'=Alpha, # DELETE
                 'IncludeDelta'=1, # DELETE
                 'IncludeEta'=1,  # DELETE
                 'SE'=1,  # DELETE
                 'DataSource' = Estimation_model_i,
                 'DataObs' = DataObs,  # 1 : zinfgamma, 2 : zinflognormal
                 'SamplingProcess' = Samp_process, # 1: sampling process is activated, else : it is ignored
                 'zero.infl_model' = 2, # type of Delta model : default is 2
                 'commercial_obs' = 1, # 1 : commercial observations are considered in the likelihood
                 'b_constraint' = b_constraint, # 1 : b > 0 | 2 : b is free
                 'catchability_random' = T,
                 'cov_samp_process' = 0, # DELETE
                 'const_kappa' = 1, # DELETE
                 'const_target' = 1, # DELETE
                 'const_q1' = 1, # DELETE
                 'const_k' = 1, # DELETE
                 'const_Sigma' = 1, # DELETE
                 'ref_data' = ifelse(ref_data=="com",1,0),
                 'const_tau' = 1, # DELETE
                 'const_spphab' = 1, # DELETE
                 'biomassAR1' = ifelse(biomassAR1==T,1,0),
                 'sampAR1' = ifelse(sampAR1==T,1,0),
                 'anisotropy' = ifelse(anisotropy==T,1,0),
                 'lf_link' = lf_link # latent field link function. 0: log link, 1: logit link
)

n_com.fleets = length(unique(b_com_i))

if( Options_vec['Prior']==0 ){
  
  deltainput_x = matrix(0,nrow = mesh$n, ncol = length(time.step_df$t)) # when SPDE approach
  
  etainput_x = array(0,c(mesh$n,length(time.step_df$t),n_com.fleets))
  
}

## Data & Params
Map = list()
Random = c()
if(Estimation_model_i == 1){   # Integrated model (scientific_commercial)
  
  Data = list( "Options_vec"=Options_vec,
               
               "n_x" = mesh$n, # number of cells
               "n_t" = nrow(time.step_df), # number of time steps
               "n_p" = nrow(Cov_x_pred),
               "n_ipp" = nrow(Cov_x_pred),
               "n_com_i" = length(y_com_i), # number of commercial samples
               "n_sci_i" = length(y_sci_i), # number of scientific samples
               "n_eta" = n_com.fleets, # number of vessels
               "n_S" = 1+n_com.fleets, # number of range and marginal variance parameters
               
               "Cov_x_com"=Cov_x_com, # covariate of the biomass field
               "Cov_x_sci"=Cov_x_sci,
               
               "Cov_fb"=if(is.null(xfb_x)){as.matrix(rep(1,nrow(Cov_x_com)))}else{xfb_x}, # covariate for the commercial sampling equation
               
               "c_com_x"=as.matrix(c_com_x[,,1]), # matrix of fleet counts at each time step in each cells
               "y_com_i"=y_com_i, # commercial CPUE
               "b_com_i" = (as.numeric(b_com_i)-1), # fleet index
               "t_com_i" = t_com_i-1, # time step for ith observation of the commercial data
               "q2_com" = rep(1,length(unique(b_com_i))), # commercial catchability
               "weights_com" = weights_com, # weighting factor for commercial data
               
               "y_sci_i"=y_sci_i, # scientific CPUE
               "t_sci_i" = t_sci_i-1, # time step for ith observation of the scientific data
               "q2_sci" =  1,  # scientific catchability
               
               "spde"=spde, # SPDE objects
               "spde_aniso" =list(),
               "M0"=spde$M0,
               "M1"=spde$M1,
               "M2"=spde$M2,
               
               "Aix_ij_com"=Aix_ij_com,
               "Aix_w_com"=Aix_w_com,
               
               "Aix_ij_sci"=Aix_ij_sci,
               "Aix_w_sci"=Aix_w_sci,
               
               "Cov_fb_mesh"=if(is.null(xfb_x)){as.matrix(rep(1,mesh$n))}else{xfb_x},
               
               "Cov_x_ipp"=Cov_x_pred,
               "Aix_ij_ipp"=Aix_ij_pred,
               "Aix_w_ipp"=Aix_w_pred,
               
               "Cov_fb_ipp"=if(is.null(xfb_x)){as.matrix(rep(1,mesh_pred$n))}else{xfb_x},
               "Cov_x_pred"=Cov_x_pred,
               "Cov_fb_pred"=if(is.null(xfb_x)){as.matrix(rep(1,nrow(Cov_x_pred)))}else{xfb_x},
               "Aix_ij_pred"=Aix_ij_pred,
               "Aix_w_pred"=Aix_w_pred,
               "W"=W
  )
  
  Params = list("beta_j0"=0,
                "beta_j0season"=rep(0,nrow(time.step_df)),
                "beta_j0year"=rep(0,nrow(time.step_df)),
                "beta_j"=array(0,c(nrow(time.step_df),ncol(Data$Cov_x_com))), # linear predictor for abundance 
                "beta_fb0"=rep(0,n_com.fleets),
                "beta_fb0season"=matrix(0,nrow = nrow(time.step_df),ncol = n_com.fleets),
                "beta_fb0year"=matrix(0,nrow = nrow(time.step_df),ncol = n_com.fleets),
                "beta_fb"=matrix(0,nrow = ncol(Data$Cov_fb), ncol = n_com.fleets), # additionnal linear predictor for sampling intensity 
                "par_b"=array(0,c(nrow(time.step_df),n_com.fleets)), # link between abundance and sampling intensity
                "par_bseason"=array(0,c(nrow(time.step_df),n_com.fleets)), # link between abundance and sampling intensity
                "par_byear"=array(0,c(nrow(time.step_df),n_com.fleets)), # link between abundance and sampling intensity
                "par_bseasonyear"=array(0,c(nrow(time.step_df),n_com.fleets)), # link between abundance and sampling intensity
                "logta_S"=rep(0,1+n_com.fleets), # parameters of SPDE object (see below)
                "logtaAR1_S"=rep(0,1+n_com.fleets), # parameters of SPDE object (see below)
                "logkappa_S"=rep(0,1+n_com.fleets), 
                "deltainput_x"=deltainput_x, # input for random noise
                "epsiloninput_x"=deltainput_x,
                "etainput_x"=etainput_x,
                "psiinput_x"=etainput_x,
                "logSigma_com"=rep(0,n_com.fleets),
                "logSigma_sci"=log(1),
                
                "q1_com"=rep(0,n_com.fleets),
                "q1_sci"=rep(0,n_survey),
                
                "k_com" = rep(1,c(n_com.fleets)),
                "k_sci" = 1,
                
                "logSigma_catch" = 0,
                "logMean_catch" = 0,
                "logSigma_targ" = 0,
                "logSigma_targyear" = 0,
                "logSigma_targseason" = 0,
                "logSigma_targseasonyear" = 0,
                "Mean_targ" = rep(0,n_com.fleets),
                "rho_epsilon" = 0,
                "rho_psi" = 0,
                "rho_b" = 0,
                "ln_Hdelta_input" = rep(0, 2),
                "ln_Heta_input" = rep(0, 2)
  )
  
  
  
}else if(Estimation_model_i == 2){ # scientific model (scientific_only)
  
  Data = list( "Options_vec"=Options_vec,
               
               "n_x" = mesh$n, # number of cells
               "n_t" = nrow(time.step_df), # number of time steps
               "n_p" = nrow(Cov_x_pred),
               "n_sci_i" = length(y_sci_i), # number of scientific samples
               "n_eta" = n_com.fleets,
               "n_S" = 1, # number of range and marginal variance parameters
               
               # covariate of the biomass field
               "Cov_x_sci"=Cov_x_sci,
               
               "y_sci_i"=y_sci_i, # scientific CPUE
               "t_sci_i" = (t_sci_i-1), # time step for ith observation of the scientific data
               "q2_sci" =  1,  # scientific catchability
               
               "spde"=spde, # SPDE objects
               "spde_aniso" =list(),
               "M0"=spde$M0,
               "M1"=spde$M1,
               "M2"=spde$M2,
               
               "Aix_ij_sci"=Aix_ij_sci,
               "Aix_w_sci"=Aix_w_sci,
               
               "Cov_x_pred"=Cov_x_pred,
               "Aix_ij_pred"=Aix_ij_pred,
               "Aix_w_pred"=Aix_w_pred,
               "W"=W
  )
  
  Params = list("beta_j0"=0,
                "beta_j0season"=rep(0,nrow(time.step_df)),
                "beta_j0year"=rep(0,nrow(time.step_df)),
                
                "beta_j"=array(0,c(nrow(time.step_df),ncol(Data$Cov_x_com))), # linear predictor for abundance 
                
                "logta_S"=rep(0,1), # parameters of SPDE object (see below)
                "logtaAR1_S"=rep(0,1), # parameters of SPDE object (see below)
                
                "logkappa_S"=rep(0,1), 
                "deltainput_x"=deltainput_x, # input for random noise
                "epsiloninput_x"=deltainput_x,
                
                "logSigma_sci"=log(1),
                
                "q1_sci"=rep(0,n_survey),
                
                "k_sci" = 1,
                
                "rho_epsilon" = 0,
                "ln_Hdelta_input" = rep(0, 2)
  )
  
  Map[["k_sci"]] <- factor(NA)
  
  
}else if(Estimation_model_i == 3){ # commercial model (commercial_only)
  
  Data = list( "Options_vec"=Options_vec,
               
               "n_x" = mesh$n, # number of cells
               "n_t" = nrow(time.step_df), # number of time steps
               "n_p" = nrow(Cov_x_pred),
               "n_ipp" = nrow(Cov_x_pred),
               "n_com_i" = length(y_com_i), # number of commercial samples
               "n_eta" = n_com.fleets, # number of vessels
               "n_S" = 1+n_com.fleets, # number of range and marginal variance parameters
               
               "Cov_x_com"=Cov_x_com, # covariate of the biomass field
               "Cov_fb"=if(is.null(xfb_x)){as.matrix(rep(1,nrow(Cov_x_com)))}else{xfb_x}, # covariate for the commercial sampling equation
               
               "c_com_x"=as.matrix(c_com_x[,,1]), # matrix of fleet counts at each time step in each cells
               "y_com_i"=y_com_i, # commercial CPUE
               "b_com_i" = (as.numeric(b_com_i)-1), # fleet index
               "t_com_i" = (t_com_i-1), # time step for ith observation of the commercial data
               "q2_com" = rep(1,length(unique(b_com_i))), # commercial catchability
               "weights_com" = 1, # weighting factor for commercial data
               
               "spde"=spde, # SPDE objects
               "spde_aniso" =list(),
               "M0"=spde$M0,
               "M1"=spde$M1,
               "M2"=spde$M2,
               
               "Aix_ij_com"=Aix_ij_com,
               "Aix_w_com"=Aix_w_com,
               
               "Cov_fb_mesh"=if(is.null(xfb_x)){as.matrix(rep(1,mesh$n))}else{xfb_x},
               
               "Cov_x_ipp"=Cov_x_pred,
               "Aix_ij_ipp"=Aix_ij_pred,
               "Aix_w_ipp"=Aix_w_pred,
               
               "Cov_fb_ipp"=if(is.null(xfb_x)){as.matrix(rep(1,mesh_pred$n))}else{xfb_x},
               "Cov_x_pred"=Cov_x_pred,
               "Cov_fb_pred"=if(is.null(xfb_x)){as.matrix(rep(1,nrow(Cov_x_pred)))}else{xfb_x},
               "Aix_ij_pred"=Aix_ij_pred,
               "Aix_w_pred"=Aix_w_pred,
               "W"=W
  )
  
  Params = list("beta_j0"=0,
                "beta_j0season"=rep(0,nrow(time.step_df)),
                "beta_j0year"=rep(0,nrow(time.step_df)),
                # "beta_j0intra"=rep(0,nrow(time.step_df)),
                "beta_j"=array(0,c(nrow(time.step_df),ncol(Data$Cov_x_com))), # linear predictor for abundance 
                "beta_fb0"=rep(0,n_com.fleets),
                "beta_fb0season"=matrix(0,nrow = nrow(time.step_df),ncol = n_com.fleets),
                "beta_fb0year"=matrix(0,nrow = nrow(time.step_df),ncol = n_com.fleets),
                # "beta_fb0intra"=rep(0,nrow(time.step_df)),
                "beta_fb"=matrix(0,nrow = ncol(Data$Cov_fb), ncol = n_com.fleets), # additionnal linear predictor for sampling intensity 
                "par_b"=array(0,c(nrow(time.step_df),n_com.fleets)), # link between abundance and sampling intensity
                "par_bseason"=array(0,c(nrow(time.step_df),n_com.fleets)), # link between abundance and sampling intensity
                "par_byear"=array(0,c(nrow(time.step_df),n_com.fleets)), # link between abundance and sampling intensity
                "par_bseasonyear"=array(0,c(nrow(time.step_df),n_com.fleets)), # link between abundance and sampling intensity
                "logta_S"=rep(0,1+n_com.fleets), # parameters of SPDE object (see below)
                "logtaAR1_S"=rep(0,1+n_com.fleets), # parameters of SPDE object (see below)
                "logkappa_S"=rep(0,1+n_com.fleets), 
                "deltainput_x"=deltainput_x, # input for random noise
                "epsiloninput_x"=deltainput_x,
                "etainput_x"=etainput_x,
                "psiinput_x"=etainput_x,
                "logSigma_com"=rep(0,n_com.fleets),
                "q1_com"=rep(0,n_com.fleets),
                "k_com" = rep(1,c(n_com.fleets)),
                "logSigma_catch" = 0,
                "logMean_catch" = 0,
                "logSigma_targ" = 0,
                "logSigma_targyear" = 0,
                "logSigma_targseason" = 0,
                "logSigma_targseasonyear" = 0,
                "Mean_targ" = rep(0,n_com.fleets),
                "rho_epsilon" = 0,
                "rho_psi" = 0,
                "rho_b" = 0,
                "ln_Hdelta_input" = rep(0, 2),
                "ln_Heta_input" = rep(0, 2)
  )
  
}

## Random effect of the latent field
if(biomassAR1 == T){
  
  Map[["deltainput_x"]] <- factor(rep(NA,length(deltainput_x)))
  Random <- c(Random,"epsiloninput_x")
  Map[["logtaAR1_S"]] <- c(1) # carefull : Map of 'biomassAR1' must be before 'sampAR1'
  Map[["logta_S"]] <- c(NA)
  
}else{
  
  Map[["rho_epsilon"]] <- factor(NA)
  Map[["epsiloninput_x"]] <- factor(rep(NA,length(deltainput_x)))
  Random <- c(Random,"deltainput_x")
  Map[["logtaAR1_S"]] <- c(NA)
  Map[["logta_S"]] <- c(1) # carefull : Map of 'biomassAR1' must be before 'sampAR1'
  
}

if( "spde_aniso" %in% names(Data) ) Data[['spde_aniso']] = list("n_s"=MeshList_aniso$anisotropic_spde$n.spde, "n_tri"=nrow(MeshList_aniso$anisotropic_spde$mesh$graph$tv), "Tri_Area"=MeshList_aniso$Tri_Area, "E0"=MeshList_aniso$E0, "E1"=MeshList_aniso$E1, "E2"=MeshList_aniso$E2, "TV"=MeshList_aniso$TV-1, "G0"=MeshList_aniso$anisotropic_spde$param.inla$M0, "G0_inv"=INLA::inla.as.dgTMatrix(solve(MeshList_aniso$anisotropic_spde$param.inla$M0)) )
if(anisotropy == F) Map[["ln_Hdelta_input"]] <- factor(rep(NA, 2))

# constant species habitat relationship
Params[["beta_j"]] = matrix(rep(0,ncol(Cov_x_com)),nrow=1)
Map[["beta_j"]] = factor(matrix(rep(NA,ncol(Cov_x_com)),nrow=nrow(Params[["beta_j"]])))

## Map and Random which are common to model 1 and 3 (integrated and commercial models)
if(Estimation_model_i %in% c(1,3)){
  
  # constant kappa
  Map[["logkappa_S"]] <- c(0,rep(1,n_com.fleets))
  
  ## Targeting
  Map[["par_bseason"]] <- factor(array(NA,c(nrow(time.step_df),n_com.fleets)))
  Map[["par_byear"]] <- factor(array(NA,c(nrow(time.step_df),n_com.fleets)))
  Map[["par_bseasonyear"]] <- factor(array(NA,c(nrow(time.step_df),n_com.fleets)))
  Map[["logSigma_targyear"]] <- factor(NA)
  Map[["logSigma_targseason"]] <- factor(NA)
  Map[["logSigma_targseasonyear"]] <- factor(NA)
  Map[["rho_b"]] <- factor(NA)
  
  
  ## Covariates effects
  if(is.null(xfb_x) & Samp_process==1){
    Map[["beta_fb"]]=factor(matrix(NA,nrow = ncol(Data$Cov_fb), ncol = n_com.fleets)) # additionnal linear predictor for sampling intensity 
  }

  
  ## Eliminate linkeage of density and sampling intensity
  if( EM=="fix_b" ){
    
    Map[["par_b"]] = factor(array(NA,c(nrow(time.step_df),n_com.fleets)))
    Map[["logSigma_targ"]] <- factor(NA)
    Map[["Mean_targ"]] <- factor(rep(NA,n_com.fleets))
    Map[["par_bseason"]] <- factor(array(NA,c(nrow(time.step_df),n_com.fleets)))
    Map[["par_byear"]] <- factor(array(NA,c(nrow(time.step_df),n_com.fleets)))
    Map[["logSigma_targyear"]] <- factor(NA)
    Map[["logSigma_targseason"]] <- factor(NA)
    
  }
  
  ## Auto-regressive process
  if(sampAR1 == T){
    
    Map[["etainput_x"]] <- factor(array(NA,c(nrow(etainput_x),ncol(etainput_x),n_com.fleets)))
    Random <- c(Random,"etainput_x","psiinput_x")
    Map[["logta_S"]] <- c(Map[["logta_S"]],rep(NA,(length(Params$logta_S)-1)))
    Map[["logtaAR1_S"]] <- c(Map[["logtaAR1_S"]],seq(1:(length(Params$logta_S)-1)))
    
  }else{
    
    Map[["rho_psi"]] <- factor(NA)
    Map[["psiinput_x"]] <- factor(array(NA,c(nrow(etainput_x),ncol(etainput_x),n_com.fleets)))
    Random <- c(Random,"etainput_x")
    Map[["logta_S"]] <- c(Map[["logta_S"]],seq(1:(length(Params$logta_S)-1)))
    Map[["logtaAR1_S"]] <- c(Map[["logtaAR1_S"]],rep(NA,(length(Params$logta_S)-1)))
    
  }
  
  ## Anisotropy
  if(anisotropy == F) Map[["ln_Heta_input"]] <- factor(rep(NA, 2))
  
  ## No sampling process
  if(Samp_process == 0){
    
    Params[["etainput_x"]] = NULL
    Params[["psiinput_x"]] = NULL
    Params[["beta_fb"]] = NULL
    Params[["beta_fb0"]] = NULL
    Params[["beta_fb0year"]] = NULL
    Params[["beta_fb0season"]] = NULL
    Params[["par_b"]] = NULL
    Params[["logSigma_targ"]] = NULL
    Params[["Mean_targ"]] = NULL
    Params[["par_byear"]] = NULL
    Params[["logSigma_targyear"]] = NULL
    Params[["par_bseason"]] = NULL
    Params[["logSigma_targseason"]] = NULL
    Params[["par_bseasonyear"]] = NULL
    Params[["logSigma_targseasonyear"]] = NULL
    Params[["par_bseasonyear"]] = NULL
    Params[["rho_psi"]] = NULL
    Params[["rho_b"]] = NULL
    Params[["ln_Heta_input"]] = NULL
    
    if(T %in% str_detect(Random,"etainput_x")) Random = Random[-which(Random == "etainput_x")]
    if(T %in% str_detect(Random,"psiinput_x")) Random = Random[-which(Random == "psiinput_x")]
    if(T %in% str_detect(Random,"par_b")) Random = Random[-which(Random == "par_b")]
    
    Map[["logta_S"]][2:(n_com.fleets+1)]=NA
    Map[["logtaAR1_S"]][2:(n_com.fleets+1)]=rep(NA,n_com.fleets)
    Map[["logkappa_S"]][2:(n_com.fleets+1)]=rep(NA,n_com.fleets)
    Map[["beta_fb0"]]=NULL
    Map[["beta_fb0season"]]=NULL
    Map[["beta_fb0year"]]=NULL
    Map[["beta_fb"]]=NULL
    Map[["par_b"]]=NULL
    Map[["par_bseason"]]=NULL
    Map[["par_byear"]]=NULL
    Map[["par_bseasonyear"]]=NULL
    Map[["etainput_x"]] = NULL
    Map[["psiinput_x"]] = NULL
    Map[["logSigma_targ"]] = NULL
    Map[["logSigma_targseason"]] = NULL
    Map[["logSigma_targseasonyear"]] = NULL
    Map[["logSigma_targyear"]] = NULL
    Map[["Mean_targ"]] = NULL
    Map[["rho_psi"]] = NULL
    Map[["rho_b"]] = NULL
    Map[["rho_psi"]] = NULL
    Map[["ln_Heta_input"]] = NULL
    
  }

  if("logSigma_targ" %in% names(Params)) Map[["logSigma_targ"]] = factor(rep(NA,length(Params$logSigma_targ))) # To delete
  if("Mean_targ" %in% names(Params)) Map[["Mean_targ"]] = factor(rep(NA,length(Params$Mean_targ)))
  
  Map[["logSigma_catch"]] = factor(rep(NA,length(Params$logSigma_catch)))
  Map[["logMean_catch"]] = factor(rep(NA,length(Params$logMean_catch)))

  ## Fix reference level of the first commercial fleet
  Map[["k_com"]] <- seq(1:(length(Params$k_com)))
  Map[["k_com"]][1] <- NA # reference level is the first fleet
  Map[["k_com"]] <- factor(Map[["k_com"]])
  
  if(Samp_process == 1){
    
    Map[["beta_fb0"]] = factor(rep(NA,n_com.fleets))
    Map[["beta_fb0season"]] = matrix(NA,nrow = nrow(time.step_df),ncol = n_com.fleets)
    Map[["beta_fb0season"]] = factor(Map[["beta_fb0season"]])
    
  }
  
}

## Set scientific data as reference data
if(Estimation_model_i == 1 & ref_data == "sci"){
  
  Map[["k_com"]] <- seq(1:(length(Params$k_com)))
  Map[["k_com"]] <- factor(Map[["k_com"]])
  
  Map[["k_sci"]] <- factor(NA)
  
}

## Intercept of latent field is yearly and seasonnal
Map[["beta_j0year"]] = as.numeric(time.step_df$Year) - (min(as.numeric(time.step_df$Year)-1))
Map[["beta_j0year"]][which(Map[["beta_j0year"]] == 1)] = NA
Map[["beta_j0year"]] = factor(Map[["beta_j0year"]])

Map[["beta_j0season"]] = as.numeric(time.step_df[,1]) - (min(as.numeric(time.step_df[,1])-1))
Map[["beta_j0season"]][which(Map[["beta_j0season"]] == month_ref)] = NA
Map[["beta_j0season"]] = factor(Map[["beta_j0season"]])


## If scientific model
if(Estimation_model_i == 2){
  if(str_detect(Version,"com_x_sci_data")) Random = c("epsiloninput_x")

  if(nrow(time.step_df) == 1){
    
    # Map[["epsiloninput_x"]] <- factor(rep(NA,length(deltainput_x)))
    Map[["beta_j0year"]] = factor(rep(NA,nrow(time.step_df)))
    Map[["beta_j0season"]] = factor(rep(NA,nrow(time.step_df)))
    Map[["ln_Hdelta_input"]] <- factor(rep(NA, 2))
    Map[["logkappa_S"]] <- NULL
    # Map[["deltainput_x"]] <- NULL
    # Map[["logtaAR1_S"]] <- factor(NA)
    Map[["logta_S"]] <- factor(NA)
    Map[["k_sci"]] <- factor(NA)
    
  }
  
}

## If only one time step, fix the auto-regressive parameter
if(nrow(time.step_df) == 1) Map[["rho_epsilon"]] <- factor(NA)

## Map of logkappa_S, logta_S, logtaAR1_S
if(Estimation_model_i != 2) Map[["logkappa_S"]]=factor(Map[["logkappa_S"]])
if(!is.null(Map[["logta_S"]])) Map[["logta_S"]]=factor(Map[["logta_S"]])
if(!is.null(Map[["logtaAR1_S"]])) Map[["logtaAR1_S"]]=factor(Map[["logtaAR1_S"]])

## Presence-absence framework
if(lf_link == 1){
  
  if(Estimation_model_i %in% c(1,3)){
    Params$k_com <- rep(0,length(Params$k_com))
    Map$k_com <- NULL
    Map[["q1_com"]] <- factor(rep(NA,length(Params$q1_com)))
    Map[["logSigma_com"]] <- factor(rep(NA,length(Params$logSigma_com)))
  }
  
  if(Estimation_model_i %in% c(1,2)){
    Params$k_sci <- rep(0,length(Params$k_sci))
    Map[["k_sci"]] <- factor(NA)
    Map[["q1_sci"]] <- factor(NA)
    Map[["logSigma_sci"]] <- factor(NA)
  }
  
  if(Estimation_model_i %in% 3) Map$k_com <- factor(c(NA,1:(length(Params$k_com)-1)))
  
}
