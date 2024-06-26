#########################################
### Process CMIP6 historical run data ###
#########################################

# input is the netcdf file AFTER it has been processed by the cdo scripts
# output is the austral average covariate for 1980-2010

# 2024-06-20

process_hist <- function(path, project = T){
  
  # historical files have data from 1850 to 2014
  # 165 years * 12 months = 1980 layers
  # trim to 1980 to 2010
  # layer 1561 should be 1980
  # layer 1932 is 2010.12.16 
  hist <- stack(path)
  
  cat("object contains", raster::nlayers(hist), "data layers\n") # check it contains 1980 layers at start
  
  hist <- subset(hist, 1561:1932) # trim to 1980 to 2010
  
  # then generate austral summer index - October to March for remaining 31 years
  # 31 = 12 * 31
  index <- c(seq(1, 372, 12),
             seq(2, 372, 12),
             seq(3, 372, 12),
             seq(10, 372, 12),
             seq(11, 372, 12),
             seq(12, 372, 12))
  index <- sort(index)
  
  # process...
  hist <- subset(hist, index) # subset raster stack of 372 months to only austral summer
  hist <- calc(hist, mean) # calculate seasonal mean
  
  # project to stereographic if requested
  prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

    if (project == T){
    # crop and project
    hist <- raster::crop(hist, raster::extent(-180, 180, -90, -30)) # crop to southern hemisphere
    
    hist <- raster::projectRaster(hist,
                                  res = 100000,
                                  method = "bilinear",
                                  crs = "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    return(hist)
  } else {
    return(hist)
  }
}

# ends