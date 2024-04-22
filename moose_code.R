setwd("/home/dlongert/Documents/MDS/block6/Data-589/moose_spatial_project")
library(spatstat)
library(sf)

load("BC_Covariates.Rda")
bc_covars = DATA

load("BC_Parks.Rda")
parks_data = DATA

data = read.csv("moose_dataset.csv", sep = "\t")

# Isolating for observations within British Columbia
data$stateProvince = as.factor(data$stateProvince)
data = data[which(data$stateProvince == "British columbia" |
                  data$stateProvince == "British Columbia (Prov)" |
                  data$stateProvince == "Canada - British Columbia (BC)" |
                  data$stateProvince == "Bc" |
                  data$stateProvince == "British Columbia"),]
data$stateProvince = droplevels(data$stateProvince)

# Filtering out irrelevant columns
data = data[,-c(1:8,11:17, 20:21, 24:25, 27:30, 34:50)]

# Removing rows with missing lat/long coordinates
data = data[!is.na(data$decimalLatitude),]
data = data[!duplicated(data$decimalLatitude),]
rownames(data) = NULL

# Setting map projection
data_sf = st_as_sf(data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
wgs84 = st_crs("+proj=longlat +datum=WGS84 +no_defs")
bc_albers = st_crs("+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Transform coordinates to BC Albers projection
data_bc_albers = st_transform(data_sf, crs = bc_albers)

# Extract the transformed coordinates
bc_x = st_coordinates(data_bc_albers)[, 1]
bc_y = st_coordinates(data_bc_albers)[, 2]

data$x_coord = bc_x
data$y_coord = bc_y

bc_win = as.owin(parks_data$Window)

moose_ppp = ppp(x = bc_x, # X coordinates
                y = bc_y, # Y coordinates
                window = bc_win) # Observation window

# Filtering for points inside the window
inside_window = inside.owin(moose_ppp, w=bc_win)
moose_ppp = moose_ppp[inside_window]

# Initial plot of moose within British Columbia
plot(moose_ppp)

