#### Summarize spatial characteristics for NH EPSCoR sites

The objective of this script is to join the NH EPSCoR site locations to the National Hydrography Dataset and to calculate summary metrics for each site/upstream watershed.    

The data folder contains the watershed boundaries for each EPSCoR site from Michelle Shattuck and UNH WQAL as well as those derived from USGS StreamStats.  

The epscor_site_characteristics.R script will:  
- link each EPSCoR site to the nearest NHDV2 COMID and NHDHR NHDPlusID  
- load watershed boundaries from UNH WQAL and USGS StreamStats  
- calculate upstream watershed area  
- estimate mean watershed slope and streambed slope from 10-m DEM, National Elevation Dataset (NED)  
- summarize land use characteristics from the NLCD 2016 data set (https://www.mrlc.gov/data/nlcd-2016-land-cover-conus)  
- summarize population density characteristics from the U.S. block-level population density dataset 2010 (https://www.sciencebase.gov/catalog/item/57753ebee4b07dd077c70868)  
  
