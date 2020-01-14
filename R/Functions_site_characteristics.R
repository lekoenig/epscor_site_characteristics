
## Functions called in epscor_stream_characteristics project
## Last updated 11 December 2019 by LE Koenig


## Function to snap points to nearest point along lines:
  # sf native functions don't seem to work quite how we want; see full post here: https://stackoverflow.com/questions/51292952/snap-a-point-to-the-closest-point-on-a-line-segment-using-sf
  st_snap_points = function(x, y, max_dist = 1000) {
    
    if (inherits(x, "sf")) n = nrow(x)
    if (inherits(x, "sfc")) n = length(x)
    
    out = do.call(c,
                  lapply(seq(n), function(i) {
                    nrst = st_nearest_points(st_geometry(x)[i], y)
                    nrst_len = st_length(nrst)
                    nrst_mn = which.min(nrst_len)
                    if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                    return(st_cast(nrst[nrst_mn], "POINT")[2])
                  })
    )
    return(out)
  }

## Function to match COMID with flowline data by VPU hydroregion:

  #  This function subsets the national flowline data by vpu to speed up processing
  
  link_comid2 <-  function (points, CRS, vpu, maxDist){
    
    flowline.sub <- flowline[which(flowline$VPUID==vpu), ]
    flowline.sub <- sf::st_transform(flowline.sub,crs=5070)
    
    # read in sample points
    points <- sf::st_as_sf(points,coords=c(2,3),crs=CRS)
    points <- sf::st_transform(points,crs=5070)
    
    # find flowline reaches within maxDist of sample point and choose closest flowline
    join.dat <- sf::st_join(points,flowline.sub,sf::st_is_within_distance,dist=maxDist)
    dist <- st_distance(points,flowline.sub[which(flowline.sub$COMID %in% join.dat$COMID),])
    choose.dist <- which.min(dist)
    
    if(!length(dist))  {join.dat2 <- join.dat} else {
      join.dat2 <- join.dat[choose.dist,]}
    join.dat2$near_dist_m <- if(!length(dist)) {NA} else {as.numeric(dist[choose.dist])}
    join.dat2 <- join.dat2[,c("site","COMID","GNIS_NAME","SLOPE","REACHCODE","FTYPE","FCODE","VPUID",
                              "StreamCalc","StreamOrde","AreaSqKM","TotDASqKM","QE_MA","near_dist_m")]
    
    # snap point to closest line and save new coordinates:
    if(!is.na(join.dat2$COMID)){
      snap.pt <- st_snap_points(points,flowline.sub,max_dist = 200)
      # Transform new point to NAD83 and output coords:
      snap.pt <- snap.pt %>% st_transform(.,crs=4269) %>%
        st_as_sf(.) %>%
        mutate(Lat_SnapNHD = st_coordinates(.)[2],
               Lon_SnapNHD = st_coordinates(.)[1],
               Datum_SnapNHD = "NAD83")
    }
    points2 <- points %>% st_transform(.,crs=4269) %>% 
      mutate(Lat = st_coordinates(.)[2],
             Lon = st_coordinates(.)[1],
             Datum = "NAD83")
    if(is.na(join.dat2$COMID)){
      join.dat2$Lat = points2$Lat
      join.dat2$Lon = points2$Lon
      join.dat2$Datum = points2$Datum
    } else {
      join.dat2$Lat = snap.pt$Lat_SnapNHD
      join.dat2$Lon = snap.pt$Lon_SnapNHD
      join.dat2$Datum = snap.pt$Datum_SnapNHD
    }
    
    return(join.dat2)
  }
  
  
## Modify function above to match site with NHDHR data:
  
  link_comid3 <-  function (points, CRS, maxDist){
    
    flowline.sub <- sf::st_transform(nhdhr_dataset,crs=5070) %>%
                    # drop M dimension:
                    sf::st_zm(.,drop=TRUE,what="ZM")
    
    # read in sample points
    points <- sf::st_as_sf(points,coords=c(2,3),crs=CRS)
    points <- sf::st_transform(points,crs=5070)
    
    # find flowline reaches within maxDist of sample point and choose closest flowline
    join.dat <- sf::st_join(points,flowline.sub,sf::st_is_within_distance,dist=maxDist)
    dist <- st_distance(points,flowline.sub[which(flowline.sub$NHDPlusID %in% join.dat$NHDPlusID),])
    choose.dist <- which.min(dist)
    
    if(!length(dist))  {join.dat2 <- join.dat} else {
      join.dat2 <- join.dat[choose.dist,]}
    join.dat2$near_dist_m <- if(!length(dist)) {NA} else {as.numeric(dist[choose.dist])}
    join.dat2 <- join.dat2[,c("site","NHDPlusID","GNIS_Name","Slope","ReachCode.x","FType","FCode","VPUID.x",
                              "StreamCalc","StreamOrde","AreaSqKm","TotDASqKm","near_dist_m")]
    names(join.dat2)[names(join.dat2) == "ReachCode.x"] <- "ReachCode"
    names(join.dat2)[names(join.dat2) == "VPUID.x"] <- "VPUID"

    # snap point to closest line and save new coordinates:
    if(!is.na(join.dat2$NHDPlusID)){
      snap.pt <- st_snap_points(points,flowline.sub,max_dist = 200)
      # Transform new point to NAD83 and output coords:
      snap.pt <- snap.pt %>% st_transform(.,crs=4269) %>%
                 st_as_sf(.) %>%
                 mutate(Lat_SnapNHD = st_coordinates(.)[2],
                        Lon_SnapNHD = st_coordinates(.)[1],
                        Datum_SnapNHD = "NAD83")
    }
    points2 <- points %>% st_transform(.,crs=4269) %>% 
                          mutate(Lat = st_coordinates(.)[2],
                                 Lon = st_coordinates(.)[1],
                                 Datum = "NAD83")
    if(is.na(join.dat2$NHDPlusID)){
      join.dat2$Lat = points2$Lat
      join.dat2$Lon = points2$Lon
      join.dat2$Datum = points2$Datum
    } else {
      join.dat2$Lat = snap.pt$Lat_SnapNHD
      join.dat2$Lon = snap.pt$Lon_SnapNHD
      join.dat2$Datum = snap.pt$Datum_SnapNHD
    }

    return(join.dat2)
  }
  

## Function to estimate river slope over a given river reach:
  
  est_slope <- function(site,nhdplusid,HRflowlines,reach.length,dem){
    
    site <- site %>% st_transform(.,crs=5070)
    fromnode.up.fline <- HRflowlines$FromNode[which(HRflowlines$NHDPlusID == nhdplusid)]
    up.fline.row <- which(HRflowlines$ToNode == fromnode.up.fline & HRflowlines$StreamOrde == HRflowlines$StreamOrde[which(HRflowlines$NHDPlusID == nhdplusid)])
    rvr.reach <- HRflowlines[c(up.fline.row,which(HRflowlines$NHDPlusID == nhdplusid)),] %>% 
                 st_transform(.,crs=5070) %>% st_zm(.,drop=TRUE,what="ZM")
    
    rvr.reach.pts <- st_sample(x=rvr.reach,type="regular",size=sum(as.numeric(st_length(rvr.reach)))) %>%
                         st_cast(.,"POINT")
    prox.node <- which.min(st_distance(site,rvr.reach.pts))
    upstr.reach <- rvr.reach.pts[c(1:prox.node)] 
    reach.sub <- upstr.reach[c(which.min(abs(reach.length - as.numeric(st_distance(site,upstr.reach)))):length(upstr.reach))]
    
    top.elev <- raster::extract(dem,reach.sub[1] %>% st_sf(.) %>% st_transform(.,crs=4269))
    upstr.elev <- raster::extract(dem,reach.sub %>% st_sf(.) %>% st_transform(.,crs=4269))
    site.elev <- raster::extract(dem, site %>% st_sf(.) %>% st_transform(.,crs=4269))
    reach.length.meters <- as.numeric(st_distance(site,reach.sub[1]))
    reach.slope <- (top.elev - site.elev)/reach.length.meters
    reach.slope2 <- (max(upstr.elev) - site.elev)/reach.length.meters
    reach.slope3 <- (max(upstr.elev) - min(upstr.elev))/reach.length.meters
    return(reach.slope)
    
  }
  
  
  est_slope_flowlines <- function(site,nhdplusid,flowlines,reach.length,dem){
    
    site <- site %>% st_transform(.,crs=5070)
    flowlines <- flowlines %>% st_transform(.,crs=5070)
    snap.fline <- sf::st_join(site,flowlines,sf::st_is_within_distance,dist=10)
    joined.fline <- flowlines[which(flowlines$arcid %in% snap.fline$arcid),]
    
    near.flines <- flowlines[which(flowlines$arcid == 38151),]
    near.flines.pts <- st_sample(x=near.flines,type="regular",size=sum(as.numeric(st_length(near.flines)))) %>% st_cast(.,"POINT")
    
    joined.fline.pts <- st_sample(x=joined.fline,type="regular",size=sum(as.numeric(st_length(joined.fline)))) %>% st_cast(.,"POINT")
    prox.node <- which.min(st_distance(site,joined.fline.pts))
    joined.fline.pts <- joined.fline.pts[c(prox.node:length(joined.fline.pts))]
    
    up.fline.pts <- c(near.flines.pts,joined.fline.pts)
    
    upstr.reach <- up.fline.pts

    top.elev <- raster::extract(dem,upstr.reach[1] %>% st_sf(.) %>% st_transform(.,crs=4269))
    upstr.elev <- raster::extract(dem,upstr.reach %>% st_sf(.) %>% st_transform(.,crs=4269))
    site.elev <- raster::extract(dem, site %>% st_sf(.) %>% st_transform(.,crs=4269))
    reach.length.meters <- as.numeric(st_distance(site,upstr.reach[1]))
    reach.slope <- (top.elev - site.elev)/reach.length.meters
    reach.slope2 <- (max(upstr.elev) - site.elev)/reach.length.meters
    reach.slope3 <- (max(upstr.elev) - min(upstr.elev))/reach.length.meters
    return(reach.slope)
    
  }
  

  
  
  
