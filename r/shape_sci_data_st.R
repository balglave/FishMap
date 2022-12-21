#-----------------------
## Shape scientific data
#-----------------------

survey_data_0$Month <- as.character(survey_data_0$Month)
survey_data_0$Year <- as.character(survey_data_0$Year)

## Join scientific data with time step data frame
survey_data_1 <- survey_data_0 %>%
  mutate(Month = ifelse(Month %in% paste0(1:9),
                        paste0("0",Month),paste0(Month))) %>%
  inner_join(time.step_df)

## Convert commercial data into sf object to intersect with gridpolygon
survey_data_sf <- st_as_sf(survey_data_1,
                          coords = c("long","lati"),
                          crs = grid_projection) %>%
  ungroup()

## Discretize scientific data
survey_data_sf <- survey_data_sf[st_intersects(survey_data_sf,gridpolygon_sf) %>% lengths > 0,]
survey_data_2 <- st_join(survey_data_sf,gridpolygon_sf) %>%
  mutate(long = sf::st_coordinates(.)[,1],
         lati = sf::st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  dplyr::select(-geometry)

if(scientific_observation == "Density") survey_data_2$Sci.obs_spp <- Sci.obs_spp$CatchWgt_spp/Sci.obs_spp$SweptArea_km2
if(scientific_observation == "CPUE") survey_data_2$Sci.obs_spp <- survey_data_2$CatchWgt_spp/survey_data_2$HaulDur

y_sci_i <- survey_data_2$Sci.obs_spp # scientific observation
t_sci_i <- survey_data_2$t # time step
cov_x_sci <- matrix(data = survey_data_2$bathy_sci, ncol = 1) # covariates

## Mesh objects for scientific data
A <- inla.spde.make.A(mesh, loc=as.matrix(survey_data_2[,c("long","lati")] ))
A <- as( A, "dgTMatrix")
Aix_ij_sci <- cbind(A@i,A@j)
Aix_w_sci <- A@x
