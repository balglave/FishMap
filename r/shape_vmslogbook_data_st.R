#----------------------------
## Shape 'VMS x logbook' data
#----------------------------

## Join commercial data with time dataframe
vmslogbook_data_1 <- vmslogbook_data_0 %>%
  inner_join(time.step_df)

## Convert commercial data into sf object to intersect with gridpolygon
vmslogbook_data_sf <- st_as_sf(vmslogbook_data_1,
                               coords = c("long","lati"),
                               crs = grid_projection)

## Intersect 'VMS x logbook' data with spatial domain
vmslogbook_data_sf <- vmslogbook_data_sf[st_intersects(vmslogbook_data_sf,gridpolygon_sf) %>% lengths > 0,]
vmslogbook_data_2 <- st_join(vmslogbook_data_sf,gridpolygon_sf) %>%
  mutate(long = sf::st_coordinates(.)[,1],
         lati = sf::st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  dplyr::select(-geometry)

# CPUE data
y_com_i <- vmslogbook_data_2$CPUE_spp

# time step
t_com_i <- vmslogbook_data_2$t

# clusters
b_com_i <- vmslogbook_data_2$f

# Covariates
Cov_x_com <- matrix(data = vmslogbook_data_2$bathy_com, ncol = 1)

## To delete
c_com_x <- array(NA,dim = c(nrow(loc_x),nrow(time.step_df),length(unique(vmslogbook_data$f))))
c_com_x[which(is.na(c_com_x))] <- 0

## Mesh objects for commercial data
A <- inla.spde.make.A(mesh, loc=as.matrix(vmslogbook_data_2[,c("long","lati")] ))
A <- as( A, "dgTMatrix" )
Aix_ij_com <- cbind(A@i,A@j)
Aix_w_com <- A@x

