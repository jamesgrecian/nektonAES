#####################################
### Process CMIP6 SSP245 run data ###
#####################################

# input is the netcdf file AFTER it has been processed by the cdo scripts
# output is the austral average covariate for 2070-2100

# 2024-07-01

process_ssp <- function(path, project = T){
  
  # ssp245 files have data from 2015 to 2100
  # 86 years * 12 months = 1032 layers
  # trim to 2071 to 2100
  # layer 673 should be 2071
  # layer 1032 is 2100.12.16 
  r <- stack(path)
  
  cat("object contains", raster::nlayers(r), "data layers\n") # check it contains 1980 layers at start
  
  r <- subset(r, 673:1032) # trim to 1980 to 2010

  # then generate austral summer index - October to March for remaining 30 years
  index <- c(seq(1, 360, 12),
             seq(2, 360, 12),
             seq(3, 360, 12),
             seq(10, 360, 12),
             seq(11, 360, 12),
             seq(12, 360, 12))
  index <- sort(index)
  
  # process...
  r <- subset(r, index) # subset raster stack of 360 months to only austral summer
  r <- calc(r, mean) # calculate seasonal mean
  
  # project to stereographic if requested
  prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  if (project == T){
    # crop and project
    r <- raster::crop(r, raster::extent(-180, 180, -90, -30)) # crop to southern hemisphere
    
    r <- raster::projectRaster(r,
                                  res = 100000,
                                  method = "bilinear",
                                  crs = "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    return(r)
  } else {
    return(r)
  }
}

# include a short version for files that only contain 2065-2100
process_ssp_short <- function(path, project = T){
  
  # some ssp245 files only have data from 2065 to 2100
  # 36 years * 12 months = 432 layers
  # trim to 2071 to 2100
  # layer 73 should be 2071
  # layer 1032 is 2100.12.16 
  r <- stack(path)
  
  cat("object contains", raster::nlayers(r), "data layers\n") # check it contains 1980 layers at start
  
  r <- subset(r, 73:432) # trim to 1980 to 2010
  
  # then generate austral summer index - October to March for remaining 30 years
  index <- c(seq(1, 360, 12),
             seq(2, 360, 12),
             seq(3, 360, 12),
             seq(10, 360, 12),
             seq(11, 360, 12),
             seq(12, 360, 12))
  index <- sort(index)
  
  # process...
  r <- subset(r, index) # subset raster stack of 360 months to only austral summer
  r <- calc(r, mean) # calculate seasonal mean
  
  # project to stereographic if requested
  prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  if (project == T){
    # crop and project
    r <- raster::crop(r, raster::extent(-180, 180, -90, -30)) # crop to southern hemisphere
    
    r <- raster::projectRaster(r,
                               res = 100000,
                               method = "bilinear",
                               crs = "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    return(r)
  } else {
    return(r)
  }
}

# ends