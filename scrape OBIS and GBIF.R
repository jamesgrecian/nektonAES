################################################
### Source occurrence data from OBIS and GBIF ###
################################################

# 2023-04-23

# Supplement existing data with Pleuragramma and Bathylagus data from OBIS and GBIF data sets

# require packages
require(tidyverse)

##############################
### Bathylagus antarcticus ###
##############################

### download OBIS data set ###
dat_obis <- robis::occurrence("Bathylagus antarcticus")

# format dates
dat_obis <- dat_obis |> mutate(eventDate = ymd_hms(eventDate, tz = "UTC")) 
dat_obis <- dat_obis |> mutate(date_start = as.POSIXct(date_start / 1000, origin = "1970-01-01"),
                               date_end = as.POSIXct(date_end / 1000, origin = "1970-01-01"))
dat_obis <- dat_obis |> mutate(eventDate = coalesce(eventDate, date_start)) # if there is a startDate use that instead...
dat_obis <- dat_obis |> mutate(new_ymd = make_date(year, month, day)) # there are also year, month and day columns
dat_obis <- dat_obis |> mutate(eventDate = coalesce(eventDate, new_ymd)) # combine those and use those if eventDate is a NA
# this leaves 57 NAs but they aren't real dates so need to be dropped
dat_obis <- dat_obis |> filter(year(eventDate) >= 1955) # remove prior to 1955 following the Liu 2023 paper
dat_obis <- dat_obis |> filter(decimalLatitude < -40) # filter to same latitudinal boundary as others
dat_obis <- dat_obis |> filter(month(eventDate) %in% c(10, 11, 12, 1, 2, 3)) # only summer records

# now pull columns of interest
dat_obis <- dat_obis |> dplyr::select("species", "eventDate", "decimalLongitude", "decimalLatitude", "basisOfRecord")

### download GBIF data ###
dat_gbif <- rgbif::occ_data(scientificName = "Bathylagus antarcticus", hasCoordinate = T, limit = 25000)
dat_gbif <- dat_gbif$data

# format dates
source("R/prepGBIF.R")
dat_gbif <- prepGBIF(dat_gbif) # format dates
dat_gbif <- dat_gbif |> mutate(eventDate = ymd_hms(eventDate, tz = "UTC")) # 12112 failed to parse
dat_gbif <- dat_gbif |> mutate(new_ymd = make_date(year, month, day))
dat_gbif <- dat_gbif |> mutate(eventDate = coalesce(eventDate, new_ymd)) # leaves 18 NAs that aren't real dates and need to be dropped

dat_gbif <- dat_gbif |> filter(year(eventDate) > 1955)
dat_gbif <- dat_gbif |> filter(decimalLatitude < -40)
dat_gbif <- dat_gbif |> filter(month(eventDate) %in% c(10, 11, 12, 1, 2, 3))

# now pull columns of interest
dat_gbif <- dat_gbif |> dplyr::select("scientificName", "eventDate", "decimalLongitude", "decimalLatitude", "basisOfRecord")
names(dat_gbif)[1] <- "species"
dat_gbif$species <- "Bathylagus antarcticus"

### combine ###
dat <- dat_obis |> rbind(dat_gbif)
dat$species <- "Bathylagus antarcticus" # force same species format
dat <- dat |> arrange(eventDate) 
dat <- dat |> mutate(basisOfRecord = fct_recode(basisOfRecord, # format specimen type
                                                "PreservedSpecimen" = "PRESERVED_SPECIMEN",
                                                "HumanObservation" = "HUMAN_OBSERVATION",
                                                "MaterialSample" = "MATERIAL_SAMPLE"))
# identify duplicates following Freer et al. 2023
# same species, lon, lat, year, month and day
dat <- dat |> mutate(eventDate = as.Date(floor_date(eventDate, "day"))) # round date to nearest day
dat <- dat |> dplyr::select(!basisOfRecord)
# thre are some lat lons that are almost exactly the same - a rounding error?
# round to 5 decimal places - equivalent to 1 metre!
dat <- dat |> mutate(decimalLongitude = round(decimalLongitude, 5),
                     decimalLatitude = round(decimalLatitude, 5))

# drop duplicates
dat <- dat |> distinct()

# output
write_csv(dat, "data/bathylagus.csv")

rm(dat)
rm(dat_obis)
rm(dat_gbif)

####################
### Pleuragramma ###
####################

### download OBIS data set ###
dat_obis <- robis::occurrence(scientificname = "Pleuragramma antarctica")

# format dates
source("R/prepOBIS.R")
dat_obis <- prepOBIS(dat_obis) # format dates
dat_obis <- dat_obis |> mutate(eventDate = ymd_hms(eventDate, tz = "UTC")) 

dat_obis <- dat_obis |>
  mutate(date_start = as.POSIXct(date_start / 1000, origin = "1970-01-01"),
         date_end = as.POSIXct(date_end / 1000, origin = "1970-01-01"))

# if eventDate is missing, use date_start
dat_obis <- dat_obis |> mutate(eventDate = coalesce(eventDate, date_start))

# if eventDate and date_start are missing, create a new date from year, month and day  
dat_obis <- dat_obis |> mutate(new_ymd = make_date(year, month, day))
dat_obis <- dat_obis |> mutate(eventDate = coalesce(eventDate, new_ymd))

# filter
dat_obis <- dat_obis |> filter(year(eventDate) >= 1955) # remove prior to 1955 following the Liu 2023 paper
dat_obis <- dat_obis |> filter(decimalLatitude < -40) # filter to same latitudinal boundary as others
dat_obis <- dat_obis |> filter(month(eventDate) %in% c(10, 11, 12, 1, 2, 3)) # only summer records

# now pull columns of interest
dat_obis <- dat_obis |> dplyr::select("species", "eventDate", "decimalLongitude", "decimalLatitude", "basisOfRecord")

### download GBIF data set ###
require(rgbif)
dat_gbif <- rgbif::occ_data(scientificName = "Pleuragramma antarctica", hasCoordinate = T, limit = 25000)
dat_gbif <- dat_gbif$data

dat_gbif <- dat_gbif |>
  mutate(eventDate = case_when(eventDate == "2007-12" ~ "2007-12-01 00:00:00",
                               eventDate == "2006-01-01T23:12Z/2006-01-02T00:00Z" ~ "2006-01-01 00:00:00",
                               eventDate == "1996-02-07T23:54Z/1996-02-08T00:01Z" ~ "1996-02-07 00:00:00",
                               eventDate == "1996-01-19/1996-03-31" ~ "1996-01-19 00:00:00",
                               eventDate == "1989-01" ~ "1989-01-01 00:00:00",
                               eventDate == "1989-02" ~ "1989-02-01 00:00:00",
                               eventDate == "1985-01" ~ "1985-01-01 00:00:00",
                               eventDate == "1959-01-03/1959-01-04" ~ "1959-01-03 00:00:00",
                               eventDate == "1959-01-04/1959-01-05" ~ "1959-01-04 00:00:00",
                               eventDate == "1959-01-01/1959-01-03" ~ "1959-01-01 00:00:00",
                               eventDate == "1959-01-24/1959-01-26" ~ "1959-01-24 00:00:00",
                               eventDate == "1959-01-27/1959-01-29" ~ "1959-01-27 00:00:00",
                               eventDate == "1959-01-29/1959-01-31" ~ "1959-01-29 00:00:00",
                               eventDate == "1959-01-31/1959-02-02" ~ "1959-01-31 00:00:00",
                               eventDate == "1959-01-19/1959-02-02" ~ "1959-01-19 00:00:00",
                               eventDate == "1958-11-11/1958-11-12" ~ "1958-11-13 00:00:00",
                               eventDate == "1940-12" ~ "1949-12-01 00:00:00",
                               .default = as.character(eventDate)))

dat_gbif <- dat_gbif |> mutate(eventDate = ymd_hms(eventDate, tz = "UTC")) 
dat_gbif <- dat_gbif |> mutate(new_ymd = make_date(year, month, day))
dat_gbif <- dat_gbif |> mutate(eventDate = coalesce(eventDate, new_ymd)) # leaves 18 NAs that aren't real dates and need to be dropped

dat_gbif <- dat_gbif |> filter(year(eventDate) > 1955)
dat_gbif <- dat_gbif |> filter(decimalLatitude < -40)
dat_gbif <- dat_gbif |> filter(month(eventDate) %in% c(10, 11, 12, 1, 2, 3))

# now pull columns of interest
dat_gbif <- dat_gbif |> dplyr::select("scientificName", "eventDate", "decimalLongitude", "decimalLatitude", "basisOfRecord")
names(dat_gbif)[1] <- "species"
dat_gbif$species <- "Pleuragramma antarctica"

### combine ###
dat <- dat_obis |> rbind(dat_gbif)
dat$species <- "Pleuragramma antarctica" # force same species format
dat <- dat |> arrange(eventDate) 
dat <- dat |> mutate(basisOfRecord = fct_recode(basisOfRecord, # format specimen type
                                                "PreservedSpecimen" = "PRESERVED_SPECIMEN",
                                                "HumanObservation" = "HUMAN_OBSERVATION",
                                                "MaterialSample" = "MATERIAL_SAMPLE"))
# identify duplicates following Freer et al. 2023
# same species, lon, lat, year, month and day
dat <- dat |> mutate(eventDate = as.Date(floor_date(eventDate, "day"))) # round date to nearest day
dat <- dat |> dplyr::select(!basisOfRecord)
# thre are some lat lons that are almost exactly the same - a rounding error?
# round to 5 decimal places - equivalent to 1 metre!
dat <- dat |> mutate(decimalLongitude = round(decimalLongitude, 5),
                     decimalLatitude = round(decimalLatitude, 5))

# drop duplicates
dat <- dat |> distinct()

# add data point from BAS
foo <- tibble(species = "Pleuragramma antarctica",
              eventDate = "2003-01-12 09:02:00",
              decimalLongitude = -36.15525,
              decimalLatitude = -53.818033)
dat <- rbind(dat, foo)

# output
write_csv(dat, "data/pleuragramma.csv")

rm(dat)
rm(dat_obis)
rm(dat_gbif)

# ends
