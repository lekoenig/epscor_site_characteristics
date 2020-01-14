## The objectives of this script are to:
##     Snap EPSCoR sensor locations to nearest NHDV2 COMID and/or NHDHR_beta NHDPLUSID
##     Run spatial analyses and output site summary metrics for each site/surrounding watershed
## Lauren E Koenig
## Last updated 5 January 2020


# Load packages:
  library(dplyr)         # general data cleaning/aggregating
  library(sf)            # used for geospatial analyses
  library(mapview)       # plot spatial objects
  #devtools::install_github("USGS-R/nhdplusTools")  # install Blodgett nhdplusTools ()
  library(nhdplusTools)  # USGS OWI package for interfacing with NHDV2 and NHDHR
  #library(raster)        # work with raster data; raster overrides some base functions, so we will just call raster with ::raster instead
  
# Call functions here:
source("./R/Functions_site_characteristics.R")
# Please don't print those HR NHDPlusID's in scientific notation:
options(scipen = 999)

###################################################################    
###        Snap EPSCoR sensor locations to NHD Flowlines        ###
###################################################################
  
  ## Load National Hydrography Data (NHDV2):
  
    # Access seamless national database - https://www.epa.gov/waterdata/nhdplus-national-data 
    # warning: this dataset is very memory-intensive
    
    # Indicate where NHD national data is stored (nhdplusTools):
    #nhdplusTools::nhdplus_path("/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb")

    # Stage the national data so that it is easier to access (nhdplusTools). Note that this takes *a lot of memory* and so once this line has been run, can also call the saved RDS file:
    ##staged_data <- nhdplusTools::stage_national_data(include = "flowline",
    ##                                                 output_path = "./output/spatial")  
    # Read in flowline data and filter out coastlines:
    ##flowline <- readRDS(staged_data$flowline)
    ##names(flowline)[1:10]
    ##flowline <- flowline[-which(flowline$FTYPE=="Coastline"),]
    ##saveRDS(flowline,file="./output/spatial/nhdplus_flowlines_omitCoastlines.rds")
    ##rm(staged_data)
  
    # If the above code chunk has already been run, read in the file:
    flowline <- readRDS("./output/spatial/nhdplus_flowlines_omitCoastlines.rds")

    
  ## Define EPSCoR sensor locations:
    
    # NH EPSCoR:
    epscor.sites <- data.frame(site = c("LMP","GOF","BDC","DCF","SBM","WHB","HBF","MCQ","BEF","TPB"),
                               Lat_WGS84 = c(43.104007,42.948108,43.09331,43.134728,43.1704,43.122246,43.954853,42.964753,44.061732,43.317758),
                               Lon_WGS84 = c(-70.96261,-71.463267,-70.98904,-71.183921,-71.217305,-71.004924,-71.722534,-71.477969,-71.294576,-71.167497))
    # Include CT River (can comment out this chunk if omitted):
    siteNumber  <-'01193050'
    Conn <- dataRetrieval::readNWISsite(siteNumber)
    Conn <- Conn[,c("site_no","station_nm","dec_lat_va","dec_long_va","dec_coord_datum_cd")] %>% 
      sf::st_as_sf(coords = c("dec_long_va", "dec_lat_va"), crs = 4269) %>%
      sf::st_transform(.,crs=4326) %>%
      mutate(site = "Conn",
             Lat_WGS84 = sf::st_coordinates(.)[2],
             Lon_WGS84 = sf::st_coordinates(.)[1]) %>%
      sf::st_drop_geometry(.)
    
    # List of sites (datum = WGS84, crs = 4326):
    site.list <- rbind(epscor.sites,Conn[,c("site","Lat_WGS84","Lon_WGS84")]) %>%
                 sf::st_as_sf(coords = c("Lon_WGS84","Lat_WGS84"),crs = 4326)
    
    
  ## Spatial join - EPSCoR sensor locations to NHDV2 flowlines:
    
    # Define crs:
    CRS.def <- st_crs(flowline)
    
    # Define VPU region (the code chunk below can be modified to account for sites spanning multiple VPU's, but for now is focused on NH EPSCoR extent in VPU01):
    VPUID <- "01"
    
    # Format site list:
    sites <- site.list %>% st_transform(.,crs=CRS.def) %>%
             mutate(X = st_coordinates(.)[,1],
                    Y = st_coordinates(.)[,2])
    
    # Note that this join may take awhile to run (outputs row number as it runs to monitor progress)
    # max distance to nearest COMID is currently set at 200 m
    nhd.dat.ls <- list()
    for(i in 1:length(sites$site)){
      pts <- as.data.frame(sites)[i,c("site","X","Y")]
      vpu <- VPUID
      dat <- try(link_comid2(points = pts,CRS = CRS.def,vpu=vpu,maxDist = 200))
      nhd.dat.ls[[i]] <- dat
      print(i)
    }
    
    nhd.dat <- do.call("rbind",nhd.dat.ls)
    nhd.dat <- nhd.dat %>% st_drop_geometry()
    saveRDS(nhd.dat,"./output/epscor_nhdv2.rds")
    
    
  ## Spatial join - EPSCoR sensor locations to NHDHR (beta) flowlines:
    
    # Load HR flowlines:
    #nhdhrflowline_0106 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0106_HU4_GDB.gdb/",
    #                              layer = "NHDFlowline")
    #nhdhrflowline_0107 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0107_HU4_GDB.gdb/",
    #                             layer = "NHDFlowline")
    #nhdhrflowline_0108 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0108_HU4_GDB.gdb/",
    #                              layer = "NHDFlowline")
    nhdhrflowline_0106 <- read_sf(dsn = "C:/Users/lak17007/Documents/Spatial_data/NHDPlus/NHDPlusHR/NHDPLUS_H_0106_HU4_GDB.gdb",
                                  layer = "NHDFlowline")
    nhdhrflowline_0107 <- read_sf(dsn = "C:/Users/lak17007/Documents/Spatial_data/NHDPlus/NHDPlusHR/NHDPLUS_H_0107_HU4_GDB.gdb",
                                  layer = "NHDFlowline")
    nhdhrflowline_0108 <- read_sf(dsn = "C:/Users/lak17007/Documents/Spatial_data/NHDPlus/NHDPlusHR/NHDPLUS_H_0108_HU4_GDB.gdb",
                                  layer = "NHDFlowline")
    nhdhr_flowline <- rbind(nhdhrflowline_0106,
                            nhdhrflowline_0107,
                            nhdhrflowline_0108)
    # Load HR VAA:
    #nhdhrvaa_0106 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0106_HU4_GDB.gdb/",
    #                         layer = "NHDPlusFlowlineVAA")
    #nhdhrvaa_0107 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0107_HU4_GDB.gdb/",
    #                         layer = "NHDPlusFlowlineVAA")
    #nhdhrvaa_0108 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0108_HU4_GDB.gdb/",
    #                         layer = "NHDPlusFlowlineVAA")
    nhdhrvaa_0106 <- read_sf(dsn = "C:/Users/lak17007/Documents/Spatial_data/NHDPlus/NHDPlusHR/NHDPLUS_H_0106_HU4_GDB.gdb",
                             layer = "NHDPlusFlowlineVAA")
    nhdhrvaa_0107 <- read_sf(dsn = "C:/Users/lak17007/Documents/Spatial_data/NHDPlus/NHDPlusHR/NHDPLUS_H_0107_HU4_GDB.gdb",
                             layer = "NHDPlusFlowlineVAA")
    nhdhrvaa_0108 <- read_sf(dsn = "C:/Users/lak17007/Documents/Spatial_data/NHDPlus/NHDPlusHR/NHDPLUS_H_0108_HU4_GDB.gdb",
                             layer = "NHDPlusFlowlineVAA")
    nhdhr_vaa <- rbind(nhdhrvaa_0106,
                       nhdhrvaa_0107,
                       nhdhrvaa_0108)
    # Combine flowline and vaa tables for HR data:
    nhdhr_dataset <- left_join(nhdhr_flowline,nhdhr_vaa,by="NHDPlusID")
    
   # Join site locations with NHDHR data:
    nhdhr.dat.ls <- list()
    for(i in 1:length(sites$site)){
      pts <- as.data.frame(sites)[i,c("site","X","Y")]
      dat <- try(link_comid3(points = pts,CRS = CRS.def,maxDist = 200))
      nhdhr.dat.ls[[i]] <- dat
      print(i)
    }
    
    nhdhr.dat <- do.call("rbind",nhdhr.dat.ls)
    nhdhr.dat <- nhdhr.dat %>% st_drop_geometry()
    saveRDS(nhdhr.dat,"./output/epscor_nhdhr.rds")
    rm(flowline)
    rm(nhdhrflowline_0106)
    rm(nhdhrflowline_0107)
    rm(nhdhrflowline_0108)
    rm(nhdhrvaa_0106)
    rm(nhdhrvaa_0107)
    rm(nhdhrvaa_0108)

###################################################################    
###            Compare delineated watershed boundaries          ###
###################################################################
    
    # Load boundaries from UNH WQAL (Michelle Shattuck):
    ws.dat <- data.frame(site = nhdhr.dat$site,
                         site.elev.meters.wqal = NA,
                         site.elev.meters.streamstats = NA,
                         mean.ws.slope.deg.wqal = NA,
                         mean.ws.slope.deg.streamstats = NA,
                         stream.slope = NA,
                         TotArea_km2_wqal = NA,
                         TotArea_km2_streamstats = NA,
                         nlcd2016_PctDev_wqal = NA,
                         nlcd2016_PctDev_streamstats = NA,
                         nlcd2016_PctAg_wqal = NA,
                         nlcd2016_PctAg_streamstats = NA,
                         nlcd2016_PctForest_wqal = NA,
                         nlcd2016_PctForest_streamstats = NA,
                         nlcd2016_PctWetland_wqal = NA,
                         nlcd2016_PctWetland_streamstats = NA,
                         popden2010_wqal = NA,
                         popden2010_streamstats = NA)
    
    # Load NLCD2016 land cover data:
      
      # unzip folder that contains nlcd2016 land cover data for hydroregion VPU01
      #unzip("./data/landuse/nlcd2016_vpu01.zip", exdir = "./data/landuse/")
      unzip("./data/landuse/nlcd2016_vpu01.zip", exdir = "./data/landuse")
      # load raster:
      nlcd <- raster::raster("./data/landuse/nlcd2016_vpu01/nlcd2016_01.tif")
      #raster::plot(nlcd)
      
    # Use FedData package to download NLCD land cover data:
      # load vpu01 boundary to use as a template for clipping conus-level land cover data:
      ##unzip("./data/watershed_boundaries/HUC2_boundary/VPU01.zip", exdir = "./data/watershed_boundaries/HUC2_boundary/")
      # load shapefile:
      ##vpu01.shp <- st_read("./data/watershed_boundaries/HUC2_boundary/VPU01_boundary.shp")
      ##st_crs(vpu01.shp)           
      # remove files to save memory:
      ##file.remove(grep(list.files(path="./data/watershed_boundaries/HUC2_boundary/",full.names = T), pattern=".zip", inv=T, value=T),recursive=T)
      
      # load nlcd data from FedData (only accepts 2011, 2006, and 2001 for the year):
      ##nlcd.dat <- FedData::get_nlcd(vpu01.shp,
      ##                 label = "vpu01",
      ##                 year = 2011)
    
    # Load Pop density 2010 data:
      unzip("./data/popden2010/pden2010_block_vpu01.zip", exdir = "./data/popden2010")
      # load pop den raster:
      pop.raster <- raster::raster("./data/popden2010/pden2010_block_vpu01/pden2010_60m.tif")
      
      

    for(i in 1:length(nhdhr.dat$site)){
  
      ## 1. WQAL
        ## Load shapefile:
        # unzip folder that contains the site shapefile:
        if(nhdhr.dat$site[i]!="Conn"){
        unzip(paste("./data/watershed_boundaries/WQAL_Shattuck/",nhdhr.dat$site[i],"_shp.zip",sep=""), exdir = "./data/watershed_boundaries/WQAL_Shattuck")
        # load shapefile:
        wqal.shp <- st_read(paste("./data/watershed_boundaries/WQAL_Shattuck/",nhdhr.dat$site[i],".shp",sep="")) %>%
                    st_transform(.,crs = 4269)
        # remove files to save memory:
        file.remove(grep(list.files(path="./data/watershed_boundaries/WQAL_Shattuck",full.names = T), pattern=".zip", inv=T, value=T),recursive=T)
        
        ## Calculate watershed area:
        ws.dat[i,"TotArea_km2_wqal"] <- st_area(wqal.shp)/(10^6)
        
        ## Calculate watershed elevation/topography data:
        # load 1/3 arc-second (~10 m) elevation data from the National Elevation Dataset (NED):
        ws.dem <- FedData::get_ned(template = wqal.shp,res = "13",
                                   label = as.character(nhdhr.dat$site[i]),
                                   raw.dir = "./data/RAW/NED",
                                   extraction.dir = "./data/EXTRACTIONS/NED")
        # convert dem raster to slope raster:
        ws.slope <- raster::terrain(x = ws.dem, opt="slope", unit="degrees", neighbors=8)
        # get slope raster values and calculate mean
        mean.ws.slope <- raster::getValues(ws.slope) %>% mean(.,na.rm=T)
        # extract elevation at sample location:
        pt <- nhdhr.dat[i,] %>% st_as_sf(coords=c("Lon","Lat"),crs=4269)
        site.elev <- raster::extract(ws.dem,pt)
        
        # estimate streambed slope over a 200 m reach:
        # if(!is.na(nhdhr.dat$NHDPlusID[i])){
        # reach.slope <- est_slope(site = pt, nhdplusid = nhdhr.dat$NHDPlusID[i],
        #                        HRflowlines = nhdhr_dataset,reach.length = 200, dem = ws.dem)
        # } else {
        # reach.slope <- NA
        # }
    
        # remove files to save memory:
        unlink("./data/RAW", recursive = TRUE)  
        unlink("./data/EXTRACTIONS", recursive = TRUE)  
        
        ws.dat[i,"mean.ws.slope.deg.wqal"] <- mean.ws.slope
        ws.dat[i,"site.elev.meters.wqal"] <- site.elev
        
        ## Calculate mean watershed NLCD classes
        # Crop and mask nlcd raster before summarizing land cover data (don't actually need to crop raster if just tabulating areas):
        nlcd.cr <- raster::crop(x = nlcd, y = wqal.shp %>% st_transform(.,crs = st_crs(nlcd)))
        nlcd.fr <- raster::rasterize(x = wqal.shp %>% st_transform(.,crs = st_crs(nlcd)),y = nlcd.cr)
        nlcd.lr <- raster::mask(x = nlcd.cr, mask = nlcd.fr)
        
        # Extract land cover data within watershed boundaries:
        lndcvr.wqal <- raster::extract(x = nlcd.lr,
                                       y = wqal.shp %>% st_transform(.,crs = st_crs(nlcd.lr)),
                                       method = "simple")
        # Summarize proportion of each cover type for each site:
        sum.lndcvr.wqal <- lapply(lndcvr.wqal,
                                  function(x) {prop.table(table(x))})
        sum.lndcvr.wqal <- do.call("rbind",sum.lndcvr.wqal)
        
        # Populate land cover results (note legend here: https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend):
        ws.dat[i,"nlcd2016_PctDev_wqal"] <- sum(sum.lndcvr.wqal[which(colnames(sum.lndcvr.wqal) %in% c("21","22","23","24"))])*100
        ws.dat[i,"nlcd2016_PctAg_wqal"] <- sum(sum.lndcvr.wqal[which(colnames(sum.lndcvr.wqal) %in% c("81","82"))])*100
        ws.dat[i,"nlcd2016_PctForest_wqal"] <- sum(sum.lndcvr.wqal[which(colnames(sum.lndcvr.wqal) %in% c("41","42","43"))])*100
        ws.dat[i,"nlcd2016_PctWetland_wqal"] <- sum(sum.lndcvr.wqal[which(colnames(sum.lndcvr.wqal) %in% c("90","95"))])*100
      
        ## Calculate mean population density for each watershed:
        # Crop and mask pop density raster before summarizing pop density data (don't actually need to crop raster if just tabulating areas):
        pop.cr <- raster::crop(x = pop.raster, y = wqal.shp %>% st_transform(.,crs = st_crs(pop.raster)))
        pop.fr <- raster::rasterize(x = wqal.shp %>% st_transform(.,crs = st_crs(pop.raster)),y=pop.cr)
        pop.lr <- raster::mask(x = pop.cr,mask=pop.fr)
        
        # Extract pop density data within watershed boundaries:
        pop.wqal <- raster::extract(x = pop.lr,
                                    y = wqal.shp %>% st_transform(.,crs = st_crs(pop.lr)),
                                    method = "simple")
 
        # Summarize population density for each site:
        meanpop.wqal <- do.call("rbind",lapply(pop.wqal,function(x) {mean(x)}))
        ws.dat[i,"popden2010_wqal"]<- meanpop.wqal
        
        
        }
      
      
      ## 2. StreamStats
        ## Load shapefile:
        # unzip folder that contains the site shapefile:
        unzip(paste("./data/watershed_boundaries/StreamStats/",nhdhr.dat$site[i],".zip",sep=""), exdir = "./data/watershed_boundaries/StreamStats")
        # load shapefile:
        streamstats.shp <- st_read(paste("./data/watershed_boundaries/StreamStats/","globalwatershed.shp",sep="")) %>%
                           st_transform(.,crs = 4269)
        # remove files to save memory:
        file.remove(grep(list.files(path="./data/watershed_boundaries/StreamStats",full.names = T), pattern=".zip", inv=T, value=T),recursive=T)
        
        ## Calculate watershed area:
        ws.dat[i,"TotArea_km2_streamstats"] <- st_area(streamstats.shp)/(10^6)
        
        ## Calculate watershed elevation/topography data:
        # load 1/3 arc-second elevation data from the National Elevation Dataset (NED):
        ws.dem <- FedData::get_ned(template = streamstats.shp,res = "13",
                                   label = as.character(nhdhr.dat$site[i]),
                                   raw.dir = "./data/RAW/NED",
                                   extraction.dir = "./data/EXTRACTIONS/NED")
        # convert dem raster to slope raster:
        ws.slope <- raster::terrain(x = ws.dem, opt="slope", unit="degrees", neighbors=8)

        # get slope raster values and calculate mean
        mean.ws.slope <- raster::getValues(ws.slope) %>% mean(.,na.rm=T)
        # extract elevation at sample location:
        pt <- nhdhr.dat[i,] %>% st_as_sf(coords=c("Lon","Lat"),crs=4269)
        site.elev <- raster::extract(ws.dem,pt)
        
        # estimate streambed slope over a 200 m reach:
        if(!is.na(nhdhr.dat$NHDPlusID[i])){
          reach.slope <- est_slope(site = pt, nhdplusid = nhdhr.dat$NHDPlusID[i],
                                   HRflowlines = nhdhr_dataset,reach.length = 200, dem = ws.dem)
        } else if(nhdhr.dat$site[i] == "BDC"){
          # Note that slope estimation along BDC flowline is pretty much manually calculated:
          LMP.flowlines <- st_read("./data/flowlines/LMP/LMP_stream_network.shp") %>% st_transform(.,4269)
          reach.slope <- est_slope_flowlines(site=pt,flowlines = LMP.flowlines,dem = ws.dem)
        } else {
          reach.slope <- NA
        }
        
        # remove files to save memory:
        unlink("./data/RAW", recursive = TRUE)  
        unlink("./data/EXTRACTIONS", recursive = TRUE)  
        
        ws.dat[i,"mean.ws.slope.deg.streamstats"] <- mean.ws.slope
        ws.dat[i,"site.elev.meters.streamstats"] <- site.elev
        ws.dat[i,"stream.slope"] <- reach.slope
        
        ## Calculate mean watershed NLCD classes
        # Crop and mask nlcd raster before summarizing land cover data (don't actually need to crop raster if just tabulating areas):
        nlcd.cr <- raster::crop(x = nlcd, y = streamstats.shp %>% st_transform(.,crs = st_crs(nlcd)))
        nlcd.fr <- raster::rasterize(x = streamstats.shp %>% st_transform(.,crs = st_crs(nlcd)),y = nlcd.cr)
        nlcd.lr <- raster::mask(x = nlcd.cr, mask = nlcd.fr)
   
        # Extract land cover data within watershed boundaries:
        lndcvr.streamstats <- raster::extract(x = nlcd.lr,
                                              y = streamstats.shp %>% st_transform(.,crs = st_crs(nlcd.lr)),
                                              method = "simple")
        # Summarize proportion of each cover type for each site:
        sum.lndcvr.streamstats <- lapply(lndcvr.streamstats,
                                  function(x) {prop.table(table(x))})
        sum.lndcvr.streamstats <- do.call("rbind",sum.lndcvr.streamstats)
        
        # Populate land cover results (note legend here: https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend):
        ws.dat[i,"nlcd2016_PctDev_streamstats"] <- sum(sum.lndcvr.streamstats[which(colnames(sum.lndcvr.streamstats) %in% c("21","22","23","24"))])*100
        ws.dat[i,"nlcd2016_PctAg_streamstats"] <- sum(sum.lndcvr.streamstats[which(colnames(sum.lndcvr.streamstats) %in% c("81","82"))])*100
        ws.dat[i,"nlcd2016_PctForest_streamstats"] <- sum(sum.lndcvr.streamstats[which(colnames(sum.lndcvr.streamstats) %in% c("41","42","43"))])*100
        ws.dat[i,"nlcd2016_PctWetland_streamstats"] <- sum(sum.lndcvr.streamstats[which(colnames(sum.lndcvr.streamstats) %in% c("90","95"))])*100
        
        
        ## Calculate mean population density for each watershed:
        # Crop and mask pop density raster before summarizing pop density data (don't actually need to crop raster if just tabulating areas):
        pop.cr <- raster::crop(x = pop.raster, y = streamstats.shp %>% st_transform(.,crs = st_crs(pop.raster)))
        pop.fr <- raster::rasterize(x = streamstats.shp %>% st_transform(.,crs = st_crs(pop.raster)),y=pop.cr)
        pop.lr <- raster::mask(x = pop.cr,mask=pop.fr)
        
        # Extract pop density data within watershed boundaries:
        pop.streamstats <- raster::extract(x = pop.lr,
                                    y = streamstats.shp %>% st_transform(.,crs = st_crs(pop.lr)),
                                    method = "simple")
        
        # Summarize population density for each site:
        meanpop.streamstats <- do.call("rbind",lapply(pop.streamstats,function(x) {mean(x)}))
        ws.dat[i,"popden2010_streamstats"]<- meanpop.streamstats
        
        
    print(i)
        
    }
    
   
    # Evaluate calculated watershed data and save to file:
    print(ws.dat %>% as_tibble(.))
    write.csv(ws.dat,"./output/epscor_landcover_compare.csv",row.names = FALSE)
    
    
    
    # remove extracted land use and pop density raster files to save memory:
    file.remove(grep(list.files(path="./data/landuse/nlcd2016_vpu01",full.names = T), pattern=".zip", inv=T, value=T),recursive=T)
    file.remove(grep(list.files(path="./data/landuse",full.names = T), pattern=".zip", inv=T, value=T),recursive=T)
    
    file.remove(grep(list.files(path="./data/popden2010/pden2010_block_vpu01/",full.names = T), pattern=".zip", inv=T, value=T),recursive=T)

    
    ## 3. Join StreamCat covariates:
    # data downloaded from EPA StreamCat by HydroRegion
    
    StreamCat.gathered <- NA
      
    # Load zip folder that contains StreamCat data
    folder <- paste("./data/StreamCat")
    zips <- list.files(folder,full.names = TRUE)
    for(j in 1:length(zips)){
      StreamCat.dat <- do.call("rbind", lapply(zips[j],function(x) {
        df <- read.csv(unz(x, unzip(x, list=TRUE)$Name), header = TRUE,sep = ",")}))
      if(j == 1){
        StreamCat.dat2 <- StreamCat.dat
      } else {
        StreamCat.dat2 <- left_join(StreamCat.dat2,StreamCat.dat,by="COMID")
      }
    }
    
    StreamCat.all <- StreamCat.dat2[,c("COMID","WsAreaSqKm.x","BFIWs","BFICat","PctImp2011Ws",
                                          "PctImp2011Cat","Precip8110Ws","Precip8110Cat","Tmean8110Ws",
                                          "Tmean8110Cat","RunoffWs","RunoffCat","ClayWs","ClayCat",
                                          "SandWs","SandCat","OmWs","OmCat","PermWs","PermCat",
                                          "RckDepWs","RckDepCat","WtDepWs","WtDepCat",
                                          "HUDen2010Ws","HUDen2010Cat","PopDen2010Ws","PopDen2010Cat")]
        
    StreamCat.gather <- nhdhr.dat
    # Add COMID to data set that we're compiling:
    for(b in 1:length(nhdhr.dat$site)){
    StreamCat.gather$COMID[b] <- nhd.dat$COMID[which(nhd.dat$site==nhdhr.dat$site[b])]
    }
    # Re-arrange columns:
    StreamCat.gather <- StreamCat.gather[,c("site","COMID","NHDPlusID","GNIS_Name","Slope","ReachCode","FType",
                          "FCode","VPUID","StreamCalc","StreamOrde","AreaSqKm","TotDASqKm",
                          "near_dist_m","Lat","Lon","Datum")]
    # Join StreamCat.gather with StreamCat covariates:
    StreamCat.gather2 <- left_join(StreamCat.gather,StreamCat.all,by="COMID")
    
    # Save compiled NHD and StreamCat covariates to file:
    write.csv(StreamCat.gather2,"./output/epscor_site_covariates.csv",row.names = FALSE)
    
  
    
    
    
  