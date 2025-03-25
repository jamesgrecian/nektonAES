##############################################################################
### Helper function to generate CMIP6 model by scenario nekton AES outputs ###
##############################################################################

# 2024-07-24

# input is a focal CMIP6 model and scenario of interest
# output is a tibble of predictions averaged across guilds

futureAES <- function(model, scenario){
  fn <- list.files("data/ensemble model outputs", pattern = scenario, full.names = T)
  x <- lapply(fn, readRDS)      # load raster stacks as a list
  x <- lapply(x, subset, model) # subset each stack to only focal GCM
  x <- stack(x)                 # convert list to stack
  
  names(x) <- c("species01", "species02", "species03", "species04", "species05",
                "species06", "species07", "species08", "species09", "species10",
                "species11", "species12", "species13", "species14", "species15",
                "species16", "species17", "species18", "species19", "species20",
                "species21", "species22", "species23", "species24", "species25",
                "species26", "species27", "species28")
  
  # should we average across all species?
  # or should we break into groups?
  fish <- subset(x, c(1:13))
  krill <- subset(x, 14)
  squid <- subset(x, c(15:28))
  
  # average across groups for equal weighting
  all_sp <- stack(mean(fish), krill, mean(squid))
  all_sp <- mean(all_sp)
  
  return(all_sp)
}

# ends
