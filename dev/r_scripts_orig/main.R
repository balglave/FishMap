#----------------------------------------------
## Toy example to run the spatio-temporal model
#----------------------------------------------

## Packages
library(dplyr)
library(INLA)
library(ggplot2)
library(raster)
library(sf)
library(stringr)
library(tidyr)
library(TMB)

## Load data
#-----------
species <- "Solea_solea" # species of interest
fleet <- c("OTB_DEF_>=70_0","OTB_CEP_>=70_0","OTT_DEF_>=70_0") # Fleet to filter
# The first one will be taken as reference level in the model

data_folder <- paste0("data-raw/",species,"/")

fitted_data <- "biomass" # "biomass" "presabs"

# Scientific data
n_survey <- 1 # number of surveys
load(paste0(data_folder,"survey_data.Rdata")) # scientific observations
scientific_observation <- "CPUE" # 'CPUE' or 'Density'

load(paste0(data_folder,"bathy_sci.Rdata")) # covariates related to the scientific observations
survey_data$bathy_sci <- bathy_sci

survey_data_0 <- survey_data %>% ungroup %>% dplyr::select(-layer)

# 'VMS x logbook' data
load(paste0(data_folder,"vmslogbook_data.Rdata")) # 'vms x logbooks' observations

select_aggreg_level <- paste(fleet,collapse = "|")
vmslogbook_data <- vmslogbook_data %>%
  filter(str_detect(LE_MET_level6,select_aggreg_level))

vmslogbook_data$LE_MET_level6 <- factor(as.character(vmslogbook_data$LE_MET_level6),levels = fleet)

load(paste0(data_folder,"bathy_com.Rdata"))
vmslogbook_data$bathy_com <- bathy_com

vmslogbook_data_0 <- vmslogbook_data %>% ungroup %>% dplyr::select(-layer)

## Time series
#-------------
year_start <- 2018
year_end <- 2018
year_vec <- year_start:year_end
month_start <- 11
month_end <- 11
month_vec <- month_start:month_end
time_step <- "Month" # Month or Quarter

if(time_step == "Month"){

  time.step_df <- expand.grid(month_vec,year_vec)
  colnames(time.step_df) <- c("Month","Year")

  time.step_df <- time.step_df %>%
    arrange(Year,Month) %>%
    mutate(Month = ifelse(Month < 10,paste0("0",Month),Month)) %>%
    mutate(Year_Month = paste0(Year,"_",Month)) %>%
    mutate(t = 1:nrow(time.step_df))
  time.step_df$Year <- as.character(time.step_df$Year)
  time.step_df$Month <- as.character(time.step_df$Month)

}else if(time_step == "Quarter"){

  time.step_df <- expand.grid(1:4,all_years)
  colnames(time.step_df) <- c("Quarter","Year")

  time.step_df <- time.step_df %>%
    arrange(Year,Quarter) %>%
    mutate(Quarter = ifelse(Quarter < 10,paste0("0",Quarter),Quarter)) %>%
    mutate(Year_Quarter = paste0(Year,"_",Quarter)) %>%
    mutate(t = 1:nrow(time.step_df))
  time.step_df$Year <- as.character(time.step_df$Year)
  time.step_df$Quarter <- as.character(time.step_df$Quarter)

}


## Configure spatial domain
#--------------------------
grid_xmin <- -6
grid_xmax <- 0
grid_ymin <- 42
grid_ymax <- 48
grid_limit <- extent(c(grid_xmin,grid_xmax,grid_ymin,grid_ymax))

grid_projection <- "+proj=longlat +datum=WGS84"

resol <- 0.05 # resolution of the discretization grid
create_mesh <- "from_shapefile"
# from_shapefile: the mesh will be more regular on the grid
# from_data: the mesh will be denser in the areas where there are data

# Mesh parameterization
k <- 0.25
Alpha <- 2

load(paste0(data_folder,"study_domain.Rdata"))

## Load domain / mesh / spde object
source("dev/r_scripts_orig/domain_mesh_spde.R")

## Shape scientific data
source("dev/r_scripts_orig/shape_sci_data_st.R")

## Shape commercial data
source("dev/r_scripts_orig/shape_vmslogbook_data_st.R")

if(fitted_data=="presabs"){
  y_com_i[which(y_com_i > 0)] <- 1
  y_sci_i[which(y_sci_i > 0)] <- 1
  lf_link <- 1 # logit link
}

load(paste0(data_folder,"bathy_pred.Rdata"))
cov_x_pred <- matrix(data = bathy_pred[1:nrow(loc_x)], ncol = 1)

########################################################################################
##################### After this point --> need to be cleaned ##########################
########################################################################################

## Fit model
#-----------

## Model configuration
SE <- 1 # run ADREPORT - 0: no, 1: yes
data_source <- 1 # 1: integrated model (scientific + commercial data), 2: scientific model, 3: commercial model
data_obs <- 2 # observation model for biomass - 1: zero_inflated gamma, 2: zero-inflated lognormal, 3: lognormal
samp_process <- 0 # Sampling process - 0: no sampling process, 1: inhomogeneous Poisson point process
b_constraint <- 2 # put constraint on b parameters - 1: b are positive, 2: no constraints
const_spphab <- 1 # Species-habitat relationship - 1: constant in time
cov_samp_process <- 0 # covariate in the sampling process - 0: none, 1: covariate in the samplign process
biomass_temporal <- 1 # Account for temporal correlation in biomass - 0: no, 1: AR1
sampling_temporal <- 0 # Account for temporal correlation in sampling process - 0: no, 1: AR1
anisotropy <- 0 # Account for anisotropy
lf_link <- 0 # link function of the latent field - 0: log (biomass data), 1: logit (presence absence data)
ref_data <- "com" # reference data - "com": commercial (default) or "sci": scientific
EM <- "est_b" # "est_b": b is estimated, "fix_b": b is fixed
month_ref <- 1 # reference month (here January)
cov_samp_process <- 0 # 0: no covariate in the sampling process, 1: put covariate in the sampling process
compute_sd <- F

xfb_x <- NULL # TO DELETE
weights_com <- 1 # TO DELETE

## Build Data, Params and Map objects for model fitting
source("dev/r_scripts_orig/build_data_params_map.R")

# Add link to path
fixwinpath <- function(){
  PATH <- Sys.getenv("PATH")
  PATH <- paste0(R.home(), "/bin/x64;", PATH)
  PATH <- paste0("c:/Rtools/mingw64/bin;", PATH)
  Sys.setenv(PATH=PATH)
}
fixwinpath()
shell("where g++")
shell("where gdb")
shell("where Rterm")

## Model compilation
TMB::compile("inst/model.cpp","-O1 -g",DLLFLAGS="")
dyn.load( dynlib("inst/model") )

## Fit model

# ## Debugging
# source("r/function/MakeADFun_windows_debug.R")
# MakeADFun_windows_debug(cpp_name = "inst/model",  data=Data, parameters=Params,  random=Random)
# TMB::gdbsource("inst/model.R",interactive = T) ## Non-interactive
# dyn.unload( dynlib( "inst/model" ) )
# #-----------

start_time <- Sys.time()

obj = MakeADFun( data=Data, parameters=Params,  random=Random, map = Map, silent = TRUE,hessian = T )
obj$fn( obj$par )

# Parameters boundary for optimization
Lower <- -50
Upper <- 50

if(T %in% str_detect(names(obj$par),"rho_")){ # constraints on the bounds of rho

  Lower <- rep(-50,length(obj$par))
  Upper <- rep(50,length(obj$par))

  Lower[which(str_detect(names(obj$par),"rho_"))] <- -0.99
  Upper[which(str_detect(names(obj$par),"rho_"))] <- 0.99

}

opt = nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, lower=Lower, upper=Upper, control=list(trace=1,maxit=200))
opt[["diagnostics"]] = data.frame( "Param"=names(obj$par), "Lower"=-Inf, "Est"=opt$par, "Upper"=Inf, "gradient"=obj$gr(opt$par) )
report = obj$report() # output values
converge=opt$convergence # convergence test
if(compute_sd) SD = sdreport(obj,bias.correct=F,ignore.parm.uncertainty=T) # compute standard deviation

dyn.unload( dynlib( "inst/model" ) )

end_time <- Sys.time()
end_time - start_time

## Plot model outputs
#--------------------

if(nrow(time.step_df)==1){

  pred_df <- cbind(loc_x[,c("long","lati")],S_x=report$S_p[1:nrow(loc_x),]) %>%
    pivot_longer(cols = starts_with("S_x."),names_to = "t", values_to = "S_x") %>%
    mutate(t = as.numeric(str_replace(t,"S_x.",""))) %>%
    inner_join(time.step_df)

  pred_sf <- st_as_sf(pred_df,coords = c("long","lati"))

  pred_plot <- ggplot(pred_sf)+
    geom_sf(aes(col=S_x))+
    scale_color_distiller(palette = "Spectral")+
    facet_wrap(.~Year_Month)


}else if(nrow(time.step_df)>1){

  pred_df <- cbind(loc_x[,c("long","lati")],S_x=report$S_p[1:nrow(loc_x),])

  pred_sf <- st_as_sf(pred_df,coords = c("long","lati"))

  pred_plot <- ggplot(pred_sf)+
    geom_sf(aes(col=S_x))+
    scale_color_distiller(palette = "Spectral")

}

x11();plot(pred_plot)

for(i in 1:3){
  eta_df <- cbind(loc_x[,c("long","lati")],eta=report$lambda_p[1:nrow(loc_x),,i]) %>%
    pivot_longer(cols = starts_with("eta."),names_to = "t", values_to = "eta") %>%
    mutate(t = as.numeric(str_replace(t,"eta.",""))) %>%
    inner_join(time.step_df)

  eta_sf <- st_as_sf(eta_df,coords = c("long","lati"))

  eta_plot <- ggplot(eta_sf)+
    geom_sf(aes(col=log(eta)))+
    scale_color_distiller(palette = "Spectral")+
    facet_wrap(.~Year_Month)

  x11();plot(eta_plot)
}
