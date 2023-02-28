#-------------------------------
## Build spatial domain and mesh
#-------------------------------
## Build regular grid objects
# Create raster of the grid
grid <- raster(grid_limit,res = resol,crs = grid_projection)

# Create spatialpolygon based on raster
gridpolygon <- rasterToPolygons(grid)
gridpolygon$layer <- c(1:length(gridpolygon$layer))
gridpolygon_sf <- st_as_sf(gridpolygon)

# Create spatial points based on raster
raster_to_point <- rasterToPoints(grid)
datapoint <- SpatialPointsDataFrame(coords=raster_to_point,
                                    data=data.frame(layer = 1:nrow(raster_to_point)),
                                    proj4string=crs(grid_projection))
datapoint_sf <- st_as_sf(datapoint)
datapoint_sf_2 <- datapoint_sf[st_intersects(datapoint_sf,study_domain) %>% lengths > 0,]

# Cross grid with study domain
datapoint_sf_2 <- st_join(datapoint_sf_2,study_domain)
datapoint_2 <- datapoint_sf_2 %>%
  mutate(long = sf::st_coordinates(.)[,1],
         lati = sf::st_coordinates(.)[,2]) %>%
  as.data.frame %>%
  dplyr::select(-geometry)
loc_x <- datapoint_2 %>% mutate(cell = 1:nrow(datapoint_2))

## Build mesh
# mesh boundary
bound <- inla.nonconvex.hull(unique(as.matrix(datapoint_2[,c("long","lati")])),convex=-0.05)
bound2 <- inla.nonconvex.hull(unique(as.matrix(datapoint_2[,c("long","lati")])),convex=-0.2)

# create mesh
if(create_mesh == "from_shapefile"){
  
  mesh <- inla.mesh.2d(
    loc=as.matrix(datapoint_2[,c("long","lati")]), ## provide locations or domain
    boundary=list(bound,bound2),
    max.edge=c(1/k, 2/k), ## mandatory
    cutoff=0.1/k,crs = inla.CRS(projargs = crs(grid_projection))) ## good to have >0
  
}else if(create_mesh == "from_data"){
  
  ptvms_wgs84_sf <- st_as_sf(vmslogbook_data,coords = c("long","lati"),crs=grid_projection)
  ptvms_wgs84_sf_2 <- ptvms_wgs84_sf[st_intersects(ptvms_wgs84_sf,study_domain_sf) %>% lengths > 0,]
  ptvms_wgs84_sf_2 <- st_join(ptvms_wgs84_sf_2,study_domain_sf)
  
  mesh <- inla.mesh.2d(
    as.matrix(st_coordinates(ptvms_wgs84_sf_2)), ## provide locations or domain
    boundary=list(bound,bound2),
    max.edge=c(1/k, 2/k), ## mandatory
    cutoff=0.1/k)
  
}

## Anisotropic objects
anisotropic_mesh <- mesh
Dset = 1:2
# Triangle info
TV = anisotropic_mesh$graph$tv       # Triangle to vertex indexing
V0 = anisotropic_mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
V1 = anisotropic_mesh$loc[TV[,2],Dset]
V2 = anisotropic_mesh$loc[TV[,3],Dset]
E0 = V2 - V1                      # E = edge for each triangle
E1 = V0 - V2
E2 = V1 - V0
Tri_Area = rep(NA, nrow(E0))
crossprod_fn = function(Vec1,Vec2) abs(det( rbind(Vec1,Vec2) ))
for(i in 1:length(Tri_Area)) Tri_Area[i] = crossprod_fn( E0[i,],E1[i,] )/2   # T = area of each triangle
anisotropic_spde = INLA::inla.spde2.matern(anisotropic_mesh, alpha=2)
MeshList_aniso <- list(anisotropic_spde = anisotropic_spde,
                       Tri_Area = Tri_Area,
                       E0 = E0,
                       E1 = E1,
                       E2 = E2,
                       TV = TV)

## SPDE objects
spde <- (inla.spde2.matern(mesh, alpha=Alpha)$param.inla)[c("M0","M1","M2")]

## Prediction mesh and related spde objects
mesh_pred <- inla.mesh.create( loc_x[,c('long','lati')] )
W <- rep(1,nrow(loc_x))
W <- c(W,rep(0,mesh_pred$n - nrow(loc_x)))

A <- inla.spde.make.A(mesh, loc=as.matrix(loc_x[,c("long","lati")] ))
A <- as( A, "dgTMatrix" )
Aix_ij_pred <- cbind(A@i,A@j)
Aix_ij_pred <- Aix_ij_pred[1:(nrow(loc_x)*3),]
Aix_w_pred <- A@x