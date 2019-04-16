# convert posterior optimized emissions and uncertainty to netcdf format
# author: Lewis Kunik

## prerequisite scripts:
##	inversion.r
##	make_Vshat.r
##
## output files:
## 	posterior.nc - netcdf file formatted to give posterior emissions
##	as well as posterior uncertainty in lon x lat x time format
##

#load package dependencies
library(ncdf4)

# run dependent scripts
source("config.r")
source("src/post_uncert_func.r") #function to extract posterior gridded uncertainty

# specify start/end from config.r variables
flux_times <- seq(from = flux_start_POSIX, to = flux_end_POSIX-1, by = flux_t_res) #flux_t_res defined in config.r

# load spatial covariance matrix E
sp_cov_file <- paste0(out_path, "sp_cov.rds")
E <- readRDS(sp_cov_file)

# load temporal covariance matrix D
tmp_cov_file <- paste0(out_path, "tmp_cov.rds")
D <- readRDS(tmp_cov_file)

# load prior uncertainty and convert to ntime x ncell matrix
sigma_file <- paste0(out_path, "sigma.rds")
sigma <- readRDS(sigma_file)
sigma_mat <- matrix(sigma, nrow = ntimes, byrow = T) # ntimes defined in config.r

# load model data mismatch matrix
R_file <- paste0(out_path, "R.rds")
R <- readRDS(R_file)
nobs <- nrow(R)

# load HQHtranspose file
HQHt_file <- paste0(out_path, "HQHt.rds")
HQHt <- readRDS(HQHt_file)

# load lonlat_domain grid info
lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)

# load posterior emissions
shat_file <- paste0(out_path, "s_hat.rds")
shat_vec <- readRDS(shat_file)
shat_mat <- matrix(shat_vec, nrow = ntimes, byrow = T)

# ~~~~~~~~~~~~~~~~~~~~~~ set up grid formatting ~~~~~~~~~~~~~~~~~~~~~~~~#

# define smallest possible rectangle to hold the posterior grid
lon <- round(seq(min(lonlat_domain[,1]), max(lonlat_domain[,1]), by = lon_res), round_digs)
lat <- round(seq(min(lonlat_domain[,2]), max(lonlat_domain[,2]), by = lat_res), round_digs)
lonlat_post <- as.matrix(expand.grid(lon, lat)) #get all lon-lat pairs in this new grid

nlon <- length(lon)
nlat <- length(lat)

# get index mapping of non-zeo domain cells in the new rectangular grid
iDomain <- apply(lonlat_domain, FUN = function(x) which((lonlat_post[, 1] == x[1]) &
    (lonlat_post[, 2] == x[2])), MARGIN = 1)

# Define aggregation operator to weight each cell by its fractional area
cellareas <- abs(2 * pi * (6371009^2) * (sin((lonlat_domain[, 2] - lat_res/2) *
		pi/180) - sin((lonlat_domain[, 2] + lat_res/2) * pi/180))/(360/lon_res))
total_area <- sum(cellareas)
W <- matrix(cellareas/total_area, nrow = 1)  #W is the 'aggregation operator'

# set up posterior gridded uncertainty array for each time step
post_uncert_mat <- array(0, dim = c(ntimes, ncells))

# ~~~~~~~~~~~~~~~~~~ obtain posterior gridded uncert ~~~~~~~~~~~~~~~~~~~#

# loop through times to get gridded uncertainties of each time step
for (ii in 1:ntimes) {
	print(paste("extracting gridded uncertainty of timestep", ii))
	# use the provided function to obtain posterior uncertainty mat for this time
	Vshat_time <- get_post_uncert_time(ii, ii, E, D, sigma, R, HQHt, lonlat_domain)
	post_uncert_mat[ii, ] <- W %*% Vshat_time #get the gridded uncertainty for this time
}

# set up 20d arrays with NA values
post_grid_2d <- array(NA, dim = c(ntimes, nrow(lonlat_post)))
uncert_grid_2d <- array(NA, dim = c(ntimes, nrow(lonlat_post)))

# assign non-NA values to the proper spots in arrays
post_grid_2d[, iDomain] <- shat_mat
uncert_grid_2d[, iDomain] <- post_uncert_mat

# set up 3-d arrays with NA values
post_grid_3d <- array(NA, dim = c(nlon, nlat, ntimes))
uncert_grid_3d <- array(NA, dim = c(nlon, nlat, ntimes))

# loop through times to assign lon x lat matrixes to time dimension
for (ii in 1:ntimes) {

		post_grid_3d[, , ii] <- t(matrix(post_grid_2d[ii,], nrow = nlat, byrow = T))
		uncert_grid_3d[, , ii] <- t(matrix(uncert_grid_2d[ii,], nrow = nlat, byrow = T))

}


# ~~~~~~~~~~~~~~~~~ save to netcdf file ~~~~~~~~~~~~~~~~~~~ #

# define dimension variables
time_dim <- ncdim_def("time", "seconds_since_1970_01_01", as.numeric(flux_times), longname = "beginning of time-step, seconds since R epoch: 1970-01-01 00:00:00")
lat_dim <- ncdim_def("lat", "degrees_north", lat, longname="latitude (center of cell)")
lon_dim <- ncdim_def("lon", "degrees_east", lon, longname="longitude (center of cell)")

# define emissions and uncertainty variables
post_emiss_var <- ncvar_def("emiss", flux_units, list(lon_dim, lat_dim, time_dim),
              longname="gridded posterior emissions, optimized")
uncert_var <- ncvar_def("uncertainty", flux_units, list(lon_dim, lat_dim, time_dim),
              longname="gridded covariance-aggregated posterior uncertainty")

# define list of variables to save
post_vars <- list(post_emiss_var, uncert_var)

# open file for writing
nc_filename <- paste0(out_path, "posterior.nc")
nc_post <- nc_create(nc_filename, post_vars)
ncvar_put(nc_post, post_emiss_var, post_grid_3d, count = c(nlon, nlat, ntimes)) #add emissions
ncvar_put(nc_post, uncert_var, uncert_grid_3d, count = c(nlon, nlat, ntimes)) #add uncertainty

nc_close(nc_post) #close netcdf file
