## The objectives of this script are to:
##     Snap EPSCoR sensor locations to nearest NHDV2 COMID and/or NHDHR_beta NHDPLUSID
##     Run spatial analyses and output site summary metrics for each site/surrounding watershed
## Lauren E Koenig
## Last updated 11 December 2019


# Load packages:
  library(dplyr)         # general data cleaning/aggregating
  library(sf)            # used for geospatial analyses
  library(mapview)       # plot spatial objects
  #devtools::install_github("USGS-R/nhdplusTools")  # install Blodgett nhdplusTools ()
  library(nhdplusTools)  # USGS OWI package for interfacing with NHDV2 and NHDHR

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
    nhdplusTools::nhdplus_path("/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb")

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
    nhdhrflowline_0106 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0106_HU4_GDB.gdb/",
                                  layer = "NHDFlowline")
    nhdhrflowline_0107 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0107_HU4_GDB.gdb/",
                                  layer = "NHDFlowline")
    nhdhrflowline_0108 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0108_HU4_GDB.gdb/",
                                  layer = "NHDFlowline")
    nhdhr_flowline <- rbind(nhdhrflowline_0106,
                            nhdhrflowline_0107,
                            nhdhrflowline_0108)
    # Load HR VAA:
    nhdhrvaa_0106 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0106_HU4_GDB.gdb/",
                             layer = "NHDPlusFlowlineVAA")
    nhdhrvaa_0107 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0107_HU4_GDB.gdb/",
                             layer = "NHDPlusFlowlineVAA")
    nhdhrvaa_0108 <- read_sf(dsn = "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusHR/NHDPLUS_H_0108_HU4_GDB.gdb/",
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

        
###################################################################    
###            Compare delineated watershed boundaries          ###
###################################################################
    
    
    
    
    
    
    
  