######################################
### Function to return AES polygon ###
######################################

# 2024-05-01

# Adapted from Reisinger marine predator AES code
# alternatively provide fixed quantile `q` and ignore `prob`
# input is a raster and the probability threshold

aes_poly <- function(r, prob = NULL, q = NULL){
  
  if(is.null(q)){
    # Calculate quantiles
    q <- quantile(r, p = prob)
  }
  if(is.null(prob)){
    q <- q
  }
  
  # Threshold raster based on quantiles
  r[r < q] <- NA
  r[!is.na(r)] <- 1
  
  # Convert to polygons
  c <- rasterToPolygons(r, dissolve = T)
  r_poly <- c |> st_as_sf()
  
  area_thresh <- units::set_units(250*250, km^2) # threshold for the smallest individual polygon area allowed
  r_poly_dropped <- smoothr::drop_crumbs(r_poly, area_thresh)
  r_poly_filled <- smoothr::fill_holes(r_poly_dropped, area_thresh)
  r_poly_smooth <- smoothr::smooth(r_poly_filled, method = "ksmooth", smoothness = 10)
  
  return(r_poly_smooth)
}

# ends