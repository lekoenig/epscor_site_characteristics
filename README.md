#### Summarize spatial characteristics for NH EPSCoR sites

The objective of this script is to join the NH EPSCoR site locations to the National Hydrography Dataset and to calculate summary metrics for each site/upstream watershed.    

The data folder contains the watershed boundaries for each EPSCoR site from Michelle Shattuck and UNH WQAL and as well as those derived from USGS StreamStats.  

The epscor_site_characteristics.R script will:  
- link each EPSCoR site to the nearest NHDV2 COMID and NHDHR NHDPlusID  
- load watershed boundaries from UNH WQAL and USGS StreamStats  
- calculate upstream watershed area  
  
