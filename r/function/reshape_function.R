#####################################
## Function to fit model on real data
#####################################

## Scientific data
#-----------------
#' @title reshape_sci_data()
#' 
#' @param list_grid_vms_survey : commercial data (discretized)
#' @param loc_x dataframe with grid points and corresponding covariates
#' @param scale_factor for rescaling area units when using number of indiv. per unit of surface (e.g. indiv/km^2 --> indiv/ha^2)
#' @param scientific_observation unit of scientific observations ("CPUE", "Density")
#' 
#' @return y_sci_i : scientific observations
#' @return index_sci_i : cells sampled
#' @return t_sci_i : time step

reshape_sci_data <- function(survey_data,
                             loc_x,
                             scale_factor = 0.01,
                             scientific_observation = "CPUE"){
  
  # Observations must intersect study domain
  toto <- SpatialPointsDataFrame(coords=survey_data[,c("long","lati")],
                                 data=data.frame(key = 1:nrow(survey_data)),
                                 proj4string=crs(grid_projection))
  toto_2 <- over(toto,study_domain_2)
  
  survey_data <- survey_data %>%
    dplyr::filter(!is.na(long)) %>%
    mutate(Sci.obs_spp = ifelse(scientific_observation == "Density", CatchWgt_spp/SweptArea_km2,CatchWgt_spp/HaulDur)) %>%
    arrange(cell)
  
  n_samp_sci <- nrow(survey_data)
  
  # density data
  y_sci_i <- survey_data$Sci.obs_spp
  
  # corresponding cell
  index_sci_i <- survey_data$cell
  
  # time step
  t_sci_i <- survey_data$t
  
  list_res <- list(y_sci_i = y_sci_i,
                   index_sci_i = index_sci_i,
                   t_sci_i = t_sci_i)
  
  return(list_res)

}



## Commercial data
#-----------------
#' @title reshape_com_data()
#' 
#' @param list_grid_vms_survey commercial data (discretized)
#' @param loc_x dataframe with grid points and corresponding covariates
#' @param select_aggreg_level aggrgation level (eg métier)
#' @param cluster_df vessels cluster
#' @param ref_vessel reference vessel (fixed catchability)
#' @param cov_fb fishing behavior covariates
#' @param com_sample If T, take randomly n_samp_com samples without replacing in the commercial data
#' @param filter_cluster the cluster to filter
#' @param Kfold_valid if T, make Kfold validation (i.e divide the dataset in folds, separate fitting data and validation data)
#' @param nb_K number of folds
#' @param validation_K Kfold iterator (Kfold dataset)
#' @param gridpoint_cov grid points with covariates
#' @param stratidied.KFold.Fact factors that are used to build a stratified Kfold
#' @param multiple_cluster If T, cluster data in several fleets  
#' 
#' @return y_com_i : commercial observations
#' @return c_com_x : number of points in a cell
#' @return index_com_i : cells sampled
#' @return b_com_i : fleet/cluster corresponding to observation in y_com_i
#' @return VE_i : vessel index
#' @return xfb_x : covariate of preferential sampling
#' @return t_com_i : time step
#' @return ptvms_wgs84_fit : fitted dataset
#' @return ptvms_wgs84_valid : full dataset
#' @return ptvms_wgs84_full : validation dataset
#' 

reshape_com_data <- function(ptvms_wgs84,
                             loc_x,
                             select_aggreg_level,
                             cluster_df,
                             ref_vessel,
                             cov_fb,
                             com_sample = F, 
                             n_samp_com,
                             filter_cluster = NA,
                             Kfold_valid = F,
                             nb_K = NULL,
                             validation_K = NULL,
                             gridpoint_cov = NA,
                             stratidied.KFold.Fact = NA,
                             multiple_cluster=F,
                             multiple_vessels=F){
  
  # Observations must intersect study domain
  toto <- SpatialPointsDataFrame(coords=ptvms_wgs84[,c("long","lati")],
                                 data=data.frame(key = 1:nrow(ptvms_wgs84)),
                                 proj4string=crs(grid_projection))
  toto_2 <- over(toto,study_domain_2)
  ptvms_wgs84 <- ptvms_wgs84 %>% dplyr::select(-ID)
  ptvms_wgs84_full <- cbind(ptvms_wgs84,ID = toto_2[,1]) %>%
    filter(!is.na(ID))
  
  # Add cell index (cell) and cluster (clust) to commercial data
  ptvms_wgs84_full <- ptvms_wgs84_full %>%
    left_join(cluster_df,by=c("VE_REF"="aggreg_level")) %>%
    arrange(cell)
  
  if(multiple_cluster == F) ptvms_wgs84_full$clust <- 1
  
  ptvms_wgs84_full$clust <- factor(ptvms_wgs84_full$clust)
  
  # Add vessel index (VE_i) to commercial data
  
  if(multiple_vessels == T){
    df_vessel <- data.frame(VE_REF_gr = levels(factor(ptvms_wgs84_full$VE_REF_gr)))
    df_vessel$VE_REF_gr <- as.character(df_vessel$VE_REF_gr)
    df_vessel$VE_REF_gr <- c(df_vessel$VE_REF_gr[which(df_vessel$VE_REF_gr == ref_vessel)],df_vessel$VE_REF_gr[which(df_vessel$VE_REF_gr != ref_vessel)])
    df_vessel <- df_vessel %>% mutate(vessel_index = 1:nrow(df_vessel))
    ptvms_wgs84_full <- inner_join(ptvms_wgs84_full,df_vessel,by="VE_REF_gr")
  }
  
  if(multiple_vessels == F) ptvms_wgs84_full$vessel_index <- 1
  
  ptvms_wgs84_full$vessel_index <- factor(ptvms_wgs84_full$vessel_index)
  
  if(com_sample == T){
    # modify size of commercial data
    nrow_vms <- nrow(ptvms_wgs84_full)
    ptvms_wgs84_full <- ptvms_wgs84_full[sample(x = 1:nrow_vms, size = n_samp_com,replace = F),]
  }
  
  if(!is.na(filter_cluster)){
    ptvms_wgs84_full <- filter(ptvms_wgs84_full,clust == filter_cluster)
  }
  
  if(Kfold_valid == T){
    
    # titi <- gridpoint_cov[,c("layer","long","lati",colnames(gridpoint_cov)[which(colnames(gridpoint_cov) %in% stratidied.KFold.Fact)])] %>%
    #   mutate(bathy = cut(bathy,breaks = bathy_breaks))
    # 
    # titi <- data.table(titi)
    # tata <- data.table(ptvms_wgs84_full)
    # toto <- merge(titi,tata,by = c("layer","long","lati"),allow.cartesian=TRUE)
    # 
    # toto <- as.data.frame(toto)
    
    ptvms_wgs84_full <- ptvms_wgs84_full %>%
      dplyr::arrange(across(all_of(stratidied.KFold.Fact))) %>%
      mutate(Kfold = rep_len(1:nb_K,nrow(ptvms_wgs84_full))) %>%
      # dplyr::select(-all_of(colnames(gridpoint_cov)[which(colnames(gridpoint_cov) %in% stratidied.KFold.Fact)])) %>%
      dplyr::arrange(cell)
    
    ptvms_wgs84_fit <- ptvms_wgs84_full %>%
      filter(Kfold != validation_K)
    
    ptvms_wgs84_valid <- ptvms_wgs84_full %>%
      filter(Kfold == validation_K)
    
  }else{
    ptvms_wgs84_fit <- ptvms_wgs84_full
    ptvms_wgs84_valid <- NULL
  }
  
  n_samp_com <- nrow(ptvms_wgs84_fit)
  
  if(n_samp_com > 0){
    # nb of fishing operation in a cell
    # df_c_com_x <- ptvms_wgs84_fit %>%
    #   dplyr::select("layer","clust","t") %>%
    #   dplyr::group_by(layer,clust,t) %>%
    #   count() %>%
    #   dplyr::select("layer","clust","t","n") %>%
    #   pivot_wider(id_cols = c("layer","t"), names_from = "clust",values_from = "n",names_prefix = "fleet_") %>%
    #   full_join(loc_x[,c("layer","cell")],by = "layer") %>%
    #   arrange(cell)
    c_com_x <- array(NA,dim = c(nrow(loc_x),nrow(time.step_df),length(levels(ptvms_wgs84_fit$clust))))
    # for(clust_i in 1:as.numeric(levels(ptvms_wgs84_fit$clust))){
    #   df_c_com_x <- ptvms_wgs84_fit %>%
    #     filter(clust == clust_i) %>%
    #     dplyr::select("layer","t") %>%
    #     dplyr::group_by(layer,t) %>%
    #     count() %>%
    #     arrange(t) %>%
    #     pivot_wider(id_cols = c("layer"), names_from = "t",values_from = "n") %>%
    #     full_join(loc_x[,c("layer","cell")],by = "layer") %>%
    #     arrange(cell) %>%
    #     ungroup() %>%
    #     dplyr::select(-cell,-layer)
    #   for(y in 1:nrow(time.step_df)){
    #     if(paste(y) %in% colnames(df_c_com_x)) vec_c_com_x <- df_c_com_x[,paste(y)]
    #     if(!paste(y) %in% colnames(df_c_com_x)) vec_c_com_x <- rep(NA,nrow(df_c_com_x))
    #     if(y==1) mat_c_com_x <- vec_c_com_x
    #     if(y>1) mat_c_com_x <- cbind(mat_c_com_x,vec_c_com_x)
    #   }
    #   colnames(mat_c_com_x) <- paste(1:12)
    #   c_com_x[,,clust_i] <- as.matrix(mat_c_com_x)
    # }
    
    # c_com_x <- NULL
    c_com_x[is.na(c_com_x)] <- 0
    
    # CPUE data
    y_com_i <- ptvms_wgs84_fit$CPUE_spp
    
    # corresponding layer
    index_com_i <- ptvms_wgs84_fit$cell
    
    # clusters
    b_com_i <- ptvms_wgs84_fit$clust
    
    # vessel index
    VE_i <- ptvms_wgs84_fit$vessel_index
    
    # fishing-behaviour covariates
    xfb_x <- ptvms_wgs84_fit[,cov_fb]
    
    # time_step
    if(is.null(ptvms_wgs84_fit$t)){
      ptvms_wgs84_fit$t <- rep(1,nrow(ptvms_wgs84_fit))
    }else{
      t_com_i <- ptvms_wgs84_fit$t
    }
  }
  list_res <- list(y_com_i = y_com_i, c_com_x = c_com_x, index_com_i = index_com_i, b_com_i = b_com_i, VE_i = VE_i, xfb_x = xfb_x, t_com_i = t_com_i, ptvms_wgs84_fit = ptvms_wgs84_fit, ptvms_wgs84_valid = ptvms_wgs84_valid,ptvms_wgs84_full = ptvms_wgs84_full)
  return(list_res)
}

