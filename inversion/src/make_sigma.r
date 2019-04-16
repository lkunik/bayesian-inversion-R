# make sigma (gridded prior uncertainty) vector for domain
# author: Lewis Kunik

## prerequisite scripts:
##  none (but prior_uncert.nc must be present in the include/ directory)
##
## output files:
##  sigma.rds - file containing a vector of length (#cells * #times)
##  with values of prior uncertainty for every timestep

# load package dependencies
library(ncdf4)

# run dependent scripts
source("config.r")


# ~~~~~~~~~~~~~~~~~~~~~~~ create time bins ~~~~~~~~~~~~~~~~~~~~~~~#

# establish the markers for the time-bin cutoffs (includes an end-timestamp after
# the last flux hour for the purpose of the cut() function)
time_bins <- seq(from = flux_start_POSIX, to = flux_end_POSIX, by = flux_t_res) #flux_t_res defined in config.r


# ~~~~~~~~~~~~~~~~~~~~~~~ load in prior uncertainty file ~~~~~~~~~~~~~~~~~~~~~~~#

print("Loading prior uncertainty file")

# load ncdf file from outer filepath
nc_sigma <- nc_open(sigma_file)

# get lon, lat, time, and flux from the netcdf file
nc_lat <- round(ncvar_get(nc_sigma, "lat"), round_digs)
nc_lon <- round(ncvar_get(nc_sigma, "lon"), round_digs)
nc_time <- ncvar_get(nc_sigma, "time")
nc_uncert <- ncvar_get(nc_sigma, "uncertainty")
nc_close(nc_sigma)

# convert time field from seconds-since-epoch to POSIX time
class(nc_time) <- c("POSIXt", "POSIXct")
attributes(nc_time)$tzone <- "UTC"

# NOTE: uncertainty vals are intended to have dimensions (lon x lat x time)

# get the number of lon/lat in the uncertainty array
nlon <- length(nc_lon)
nlat <- length(nc_lat)
ntime <- length(nc_time)

# in case the time field has length = 1
if(ntime == 1)
    nc_uncert <- array(nc_uncert, dim = c(nlon, nlat, ntime))

# ~~~~~~~~~~~~~~ obtain avg uncertainties for each timestep ~~~~~~~~~~~~~~ #

# break up emiss times into time bins corresponding to inversion flux times
times_cut_all <- as.POSIXct(cut(nc_time, breaks = time_bins), tz = "UTC")

# filter for NA values
if (any(is.na(times_cut_all))) {
    times_cut <- times_cut_all[-(which(is.na(times_cut_all)))]  #remove NA values
} else {
    times_cut <- times_cut_all  #keep as is
}

# this is the number of timesteps covered by the prior uncertainty
nbins <- length(unique(times_cut))

# Average the uncertainties to obtain mean vals for each time bin
uncert_bins <- array(0, dim = c(nlon, nlat, nbins))
for (ii in 1:nbins) {
    ibin <- which(times_cut_all == unique(times_cut)[ii])  #get indices of timesteps corresponding to this bin
    uncert_arr <- array(nc_uncert[, , ibin], dim = c(nlon, nlat, length(ibin)))  #this is needed in case t_res is hourly
    uncert_bins[, , ii] <- apply(uncert_arr, FUN = mean, MARGIN = c(1, 2))  #aggregate over bin's timesteps
}

# ~~~~~~~~~~~~~~ mask grid to include only domain cells ~~~~~~~~~~~~~~ #

# get lonlat pairs for the sigma file
lonlat_sigma <- expand.grid(nc_lon, nc_lat)

# load in lonlat_domain file specified in config.r
lonlat_domain <- readRDS(lonlat_domain_file)

# note which grid cells are in the domain
iDomain <- apply(lonlat_domain, FUN = function(x) which((lonlat_sigma[, 1] == x[1]) &
    (lonlat_sigma[, 2] == x[2])), MARGIN = 1)

# subset for the domain cells
sigma_mat <- apply(uncert_bins, FUN = function(x) x[iDomain], MARGIN = 3)

# convert to vector format
sigma <- as.vector(sigma_mat)


# ~~~~~~~~~~~~~~~~~~~~ Save prior uncertainty file ~~~~~~~~~~~~~~~~~~~~~~~#

print("Saving formatted prior emission uncertainty to sigma file")

filepath <- paste0(out_path, "sigma.rds")
saveRDS(sigma, filepath)
