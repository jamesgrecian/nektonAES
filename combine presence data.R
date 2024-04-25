###################################################
### Process species presence data for modelling ###
###################################################

# 2024-04-24

# load the 6 csv files - including the new notolepis, pleuragramma and bathylagus files
# format the file structure to be universal: group, species, date, lon, and lat
# drop salps and the other cephalopod species - to speed up process down the line
# append background points

# load libraries
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)

### load raw data ###
mycto <- read_csv("data/myctophids.csv")
krill <- read_csv("data/krill.csv")
cepha <- read_csv("data/cephalopods.csv")
bathy <- read_csv("data/bathylagus.csv")
pleur <- read_csv("data/pleuragramma.csv")
notol <- read_csv("data/notolepis.csv")

# pull columns of interest
mycto <- mycto |> dplyr::select(NAME, eventDate, LON, LAT)
cepha <- cepha |> dplyr::select(ScientificName_accepted, Longitude, Latitude)

krill <- krill |> dplyr::select(DATE, LONGITUDE, LATITUDE, STANDARDISED_KRILL_UNDER_1M2, NUMBER_OF_SALPS_UNDER_1M2)
names(krill)[4:5] <- c("Euphausia superba", "Salpidae")
krill <- krill |> pivot_longer(cols = 4:5, names_to = "species", values_to = "count_under_1m3")
krill <- krill |> dplyr::select(species, DATE, LONGITUDE, LATITUDE, count_under_1m3)
krill <- krill |> drop_na(count_under_1m3) # drop cases where no krill/salps counted # 5.5k empty rows dropped
krill <- krill |> filter(species == "Euphausia superba") # drop salps
krill <- krill |> filter(count_under_1m3 > 0) # drop empty trawls

# only consider austral spring/ summer... October to March
mycto <- mycto |> filter(month(eventDate) %in% c(10, 11, 12, 1, 2, 3)) # Jen already filtered the data...
krill <- krill %>% mutate(DATE = dmy(DATE))
krill <- krill %>% filter(month(DATE) %in% c(10, 11, 12, 1, 2, 3)) # loose about ~500 records (~3.5%)
krill <- krill %>% filter(year(DATE) > 1970) # start at 1975 - looses around 3500 records
# date data not available for cephalopods

# unify column names
mycto$group <- "myctophids"
names(mycto) <- c("species", "date", "lon", "lat", "group")
mycto <- mycto |> dplyr::select("group", "species", "date", "lon", "lat")

krill$group <- krill$species
krill <- krill |> dplyr::select("group", "species", "DATE", "LONGITUDE", "LATITUDE")
names(krill) <- c("group", "species", "date", "lon", "lat")

names(cepha) <- c("species", "lon", "lat")
cepha$group <- "cephalopods"
cepha$date <- NA
cepha <- cepha |> dplyr::select("group", "species", "date", "lon", "lat")

bathy$group <- "fish"
bathy <- bathy |> dplyr::select("group", "species", "eventDate", "decimalLongitude", "decimalLatitude")
names(bathy) <- c("group", "species", "date", "lon", "lat")

pleur$group <- "fish"
pleur <- pleur |> dplyr::select("group", "species", "eventDate", "decimalLongitude", "decimalLatitude")
names(pleur) <- c("group", "species", "date", "lon", "lat")

notol$group <- "fish"
notol$date <- NA
notol <- notol |> dplyr::select("group", "NAME", "date", "LON", "LAT")
names(notol) <- c("group", "species", "date", "lon", "lat")

### combine data frames ###
dat <- rbind(mycto, krill, cepha, bathy, pleur, notol)
dat <- dat |> arrange(group, species, date)
dat <- dat |> filter(lat < -40)

### drop species ###
# with too few points, or ones that based on previous testing won't run due to clustering
dat |> 
  group_by(species) |> 
  tally() |> 
  arrange(desc(n)) |> 
  print(n = Inf)

dat <- dat |> 
  filter(!species %in% c("Batoteuthis skolops", "Chiroteuthis veranyi", "Galiteuthis suhmi", "Grimpoteuthis megaptera",
                         "Histioteuthis miranda", "Illex argentinus", "Loligo gahi", "Lycoteuthis lorigera",
                         "Mastigoteuthis psychrophila", "Moroteuthis knipovitchi", "Parateuthis tunicata", "Semirossia patagonica",
                         "Taningia danae", "Teuthowenia pellucida", "Abraliopsis gilchristi"))

### background points ###

# define stereographic projection around Antarctica
prj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# append xy
dat[c("x", "y")] <- dat %>% st_as_sf(coords = c("lon", "lat")) %>% 
  sf::st_set_crs(4326) %>%
  st_transform(prj) %>%
  st_coordinates()

# generate 10000 background points for each species
source("R/pseudoAbs.R")
pabs <- pseudoAbs(x_min = -180,
                  x_max = 180,
                  y_min = -80,
                  y_max = -40,
                  projection = prj,
                  n = 10000 * length(unique(dat$species)))

# create dummy dataframe replicating structure of original data
# replace locations with pseudo absences coordinates
pseudo <- dat |> 
  group_by(group, species) |>
  dplyr::select(group, species, date) |> 
  slice(rep(1, 10000)) |>
  ungroup()

pseudo$date <- NA
pseudo[c("x", "y")] <- pabs
pseudo[c("lon", "lat")] <- pseudo |> 
  st_as_sf(coords = c("x", "y")) |> 
  sf::st_set_crs(prj) |>
  st_transform(4326) |>
  st_coordinates()

pseudo <- pseudo |> dplyr::select("group", "species", "date", "lon", "lat", "x", "y")

# define 1s and 0s
dat$PresAbs <- 1
pseudo$PresAbs <- 0

# combine original data with pseudo absences
dat <- rbind(dat, pseudo)
dat <- dat %>% dplyr::select("group", "species", "date", "PresAbs", "lon", "lat", "x", "y")

# plot to check
ggplot() + 
  geom_point(aes(x = x, y = y, colour = factor(PresAbs)), data = dat) +
  facet_wrap(~PresAbs)

# save outputted data frame
saveRDS(dat, "data/presence_absence_data_10k.rds")

# ends