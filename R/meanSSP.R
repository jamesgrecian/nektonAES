#############################################
### Function to process SSP netcdf files ###
############################################

# input is simply the path to the netcdf file for a specific covariate and CMIP6 model
# returns austral summer mean projected if required

# Generate austral summer mean of product from CMIP6 SSP runs
meanSSP <- function(path, project = T){
  
  # SSP run files have data from 2015-2071
  # 86 years * 12 months = 1032 layers
  # trim from 2071-2100
  # layer 673 should be 2071
  
  dat <- raster::stack(path)
  
  cat("assuming object contains 1032 layers; object contains", raster::nlayers(dat), "data layers\n")
  
  dat <- subset(dat, 673:1032)
  
  # then generate austral summer index - October to March for remaining 30 years
  index <- c(seq(1, 360, 12),
             seq(2, 360, 12),
             seq(3, 360, 12),
             seq(10, 360, 12),
             seq(11, 360, 12),
             seq(12, 360, 12))
  index <- sort(index)
  
  # process...
  dat <- subset(dat, index) # subset raster stack of 360 months to the 180 austral summer months
  dat <- calc(dat, mean) # calculate seasonal mean
  
  if (project == T){
    # crop and project
    dat <- raster::crop(dat, raster::extent(-180, 180, -90, -30)) # crop to southern hemisphere
    
    dat <- raster::projectRaster(dat,
                                  res = 100000,
                                  method = "bilinear",
                                  crs = "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    return(dat)
  } else {
    return(dat)
  }
}

# ends
