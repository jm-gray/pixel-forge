LSR_LANDSAT_8 - Land Surface Reflectance - L8 OLI/TIRS
LSR_LANDSAT_ETM_COMBINED - Landsat 7 Enhanced Thematic Mapper Plus (1999 - present) - Land Surface Reflectance
LSR_LANDSAT_TM - Landsat 4-5 Thematic Mapper (1982 - present) - Land Surface Reflectance


urmia_simple <- shapefile("/Volumes/research/fer/jmgray2/NIP/GIS/UrmiaBasinEarthExplorerSimplified/UrmiaBasinEarthExplorerSimplified.shp")

earthexplorer_search(usgs_eros_username="joshgray",usgs_eros_password="redpoint1", datasetName="LSR_LANDSAT_8", "lowerLeft"=list(latitude=75,longitude=-135), "upperRight"=list(latitude=90,longitude=-120), startDate="2006-01-01",endDate="2007-12-01", includeUnknownCloudCover=T,minCloudCover=0,maxCloudCover=100, place_order = F, verbose=T)

l8_results <- earthexplorer_search(usgs_eros_username="joshgray",usgs_eros_password="redpoint9", datasetName="LSR_LANDSAT_8", sp=urmia_simple, startDate="2015-1-1", endDate="2017-1-20",includeUnknownCloudCover=T,minCloudCover=0,maxCloudCover=50,place_order=F,verbose=T)

earthexplorer_search(usgs_eros_username="joshgray",usgs_eros_password="redpoint9",datasetName="GLS2005",lowerLeft=list(latitude=75,longitude=-135),upperRight=list(latitude=90,longitude=-120),startDate="2006-01-01",endDate="2007-12-01",includeUnknownCloudCover=T,minCloudCover=0,maxCloudCover=100,place_order=F,verbose=T)
