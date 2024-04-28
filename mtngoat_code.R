setwd("/home/dlongert/Documents/MDS/block6/Data-589/mountain_goat")
suppressMessages(library(spatstat))
suppressMessages(library(sf))
suppressMessages(library(splines))

load("BC_Covariates.Rda")
bc_covars = DATA

load("BC_Parks.Rda")
parks_data = DATA

data = read.csv("mtngoat_dataset.csv", sep = "\t")

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

# Removing rows with missing or duplicated lat/long coordinates
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
data$x_coord = st_coordinates(data_bc_albers)[, 1]
data$y_coord = st_coordinates(data_bc_albers)[, 2]

################################################################################
################ CREATING PPP OBJECT AND BASE PLOT #############################

# Mountain Goat ppp object
goat_ppp = ppp(x = data$x_coord,
               y = data$y_coord,
               window = as.owin(parks_data$Window)) 

# Filtering for points inside the window
goat_ppp = goat_ppp[inside.owin(goat_ppp, w=as.owin(parks_data$Window))]

# Initial plot of Mountain Goats within British Columbia
plot(goat_ppp, pch = 16, cex = 0.6, main = "Spatial distribution of Mountain Goats")

################################################################################
##### ASSESSING HOMOGENEITY AND CLUSTERING PATTERNS IN PIKA OCCURRENCE DATA ####

# Quadrat test to test homogeneity assumption
Q <- quadratcount(goat_ppp,
                  nx = 2,
                  ny = 3)
quadrat.test(Q)

# Mountain goat kernel density plot
lambda_goat = density(goat_ppp, sigma = bw.ppl)
par(mfrow = c(2,1))
plot(lambda_goat, main = "Mountain Goat Kernel Density", ribbon = F)
plot(goat_ppp, pch = 16, cex = 0.6, col = "white", add = T)
plot(goat_ppp, pch = 16, cex = 0.5, col = "black", add = T)


# Ripley's K model to access clustering patterns in Mountain Goats
lambda_pos <- density(goat_ppp,
                      sigma=bw.ppl,
                      positive=TRUE)

E_inhom <- envelope(goat_ppp,
                    Kinhom,
                    simulate = expression(rpoispp(lambda_pos)),
                    correction="border",
                    rank = 1,
                    nsim = 19,
                    fix.n = TRUE)

plot(E_inhom, main = "", lwd = 2, xlim = c(0,220000))

################################################################################
############### CORRELATION PATTERNS IN MOUNTAIN GOATS #########################

# rho_hat for forest and elevation
rho_forest = rhohat(goat_ppp, bc_covars$Forest)
rho_elevation = rhohat(goat_ppp, bc_covars$Elevation)


# plot rho models
par(mfrow = c(1,2))
plot(rho_forest)
plot(rho_elevation, xlim = c(0, 2600))

# assessing multicollinearity in model covariates
cor.im(bc_covars[c(2,3)], use = "pairwise.complete.obs")

################################################################################
################# CONSTRUCTING AND VALIDATING PP MODEL #########################

# Fit with elevation and HFI index as covariates and use LR test to compare 
# fit with null model to validate
fit = suppressWarnings(ppm(goat_ppp ~ Elevation + Forest, data = bc_covars))
fit_null = ppm(goat_ppp ~ 1)
suppressWarnings(anova(fit_null, fit, test = "LRT"))

# both covariates look nonlinear, so fit with GAM and use LR test to compare 
# fit with null model to validate
fit_gam = suppressWarnings(ppm(goat_ppp ~ bs(Elevation,7) + bs(Forest, 8), data = bc_covars, use.gam = TRUE))
suppressWarnings(anova(fit, fit_gam, test = "LRT"))

# Plot GAM model
plot(fit_gam,
     pch = 16,
     se = FALSE,
     superimpose = FALSE,
     log = TRUE,
     n = 500,
     main = "Estimated mountain goat intensity")

