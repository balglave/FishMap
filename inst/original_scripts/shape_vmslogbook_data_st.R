#----------------------------
## Shape 'VMS x logbook' data
#----------------------------
#' @importFrom dplyr inner_join mutate select filter group_by count arrange full_join ungroup
#' @importFrom INLA inla.spde.make.A
#' @importFrom sf st_as_sf st_intersects st_join st_coordinates
#' @importFrom tidyr pivot_wider

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
  mutate(long = st_coordinates(.)[,1],
         lati = st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  dplyr::select(-geometry)

# CPUE data
y_com_i <- vmslogbook_data_2$CPUE_spp

# time step
t_com_i <- vmslogbook_data_2$t

# clusters
b_com_i <- vmslogbook_data_2$f

# Number of fishing points per cell
c_com_x <- array(NA,
                 dim = c(nrow(loc_x),
                         nrow(time.step_df),
                         length(unique(vmslogbook_data_2$f))))

for(f_i in 1:length(unique(vmslogbook_data_2$f))){

  for(t_i in 1:length(unique(vmslogbook_data_2$t))){

    df_c_com_x <- vmslogbook_data_2 %>%
      filter(f == f_i & t == t_i) %>%
      dplyr::select("layer","t") %>%
      group_by(layer,t) %>%
      count() %>%
      arrange(t) %>%
      pivot_wider(id_cols = c("layer"), names_from = "t",values_from = "n") %>%
      full_join(loc_x[,c("layer","cell")],by = "layer") %>%
      inner_join(loc_x[,c("layer","cell")]) %>%
      arrange(cell) %>%
      ungroup() %>%
      dplyr::select(-cell,-layer)
    c_com_x[,t_i,f_i] <- as.matrix(df_c_com_x)

  }

}

c_com_x[is.na(c_com_x)] <- 0

## Mesh objects for commercial data
A <- inla.spde.make.A(mesh, loc=as.matrix(vmslogbook_data_2[,c("long","lati")] ))
A <- as( A, "dgTMatrix" )
Aix_ij_com <- cbind(A@i,A@j)
Aix_w_com <- A@x
