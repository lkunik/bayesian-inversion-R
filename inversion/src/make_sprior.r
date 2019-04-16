# format prior emissions into vector of flux values
# author: Lewis Kunik

## prerequisite scripts:
##  none (but prior_emiss.nc must be present in the include/ directory)
##
## output files:
##  sprior.rds - file containing a vector of length (#cells * #times)
##  with values of gridded prior emissions for every timestep

# load package dependencies
library(ncdf4)

# run dependent scripts
source("config.r")


# ~~~~~~~~~~~~~~~~~~~~~~~ create time bins ~~~~~~~~~~~~~~~~~~~~~~~#

# establish the markers for the time-bin cutoffs (includes an end-timestamp after
# the last flux hour for the purpose of the cut() function)
time_bins <- seq(from = flux_start_POSIX, to = flux_end_POSIX, by = flux_t_res) #flux_t_res defined in config.r


# ~~~~~~~~~~~~~~~~~~~~~~~ load in Prior emissions file ~~~~~~~~~~~~~~~~~~~~~~~#

print("Loading prior emissions file")

# load ncdf file from outer filepath
nc_prior <- nc_open(prior_file)

# get lon, lat, time, and flux from the netcdf file
nc_lat <- round(ncvar_get(nc_prior, "lat"), round_digs)
nc_lon <- round(ncvar_get(nc_prior, "lon"), round_digs)
nc_time <- ncvar_get(nc_prior, "time")
nc_emiss <- ncvar_get(nc_prior, "emiss")
nc_close(nc_prior)

# convert time field from seconds-since-epoch to POSIX time
class(nc_time) <- c("POSIXt", "POSIXct")
attributes(nc_time)$tzone <- "UTC"

# NOTE: prior emission are intended to have dimensions (lon x lat x time)

# get the number of lon/lat in the emission array
nlon <- length(nc_lon)
nlat <- length(nc_lat)
ntime <- length(nc_time)

#in case the time field has length = 1
if(ntime == 1)
    nc_emiss <- array(nc_emiss, dim = c(nlon, nlat, ntime))

# ~~~~~~~~~~~~~~ obtain avg emissions for each timestep ~~~~~~~~~~~~~~ #

# break up emiss times into time bins corresponding to inversion flux times
times_cut_all <- as.POSIXct(cut(nc_time, breaks = time_bins), tz = "UTC")

# filter for NA values
if (any(is.na(times_cut_all))) {
    times_cut <- times_cut_all[-(which(is.na(times_cut_all)))]  #remove NA values
} else {
    times_cut <- times_cut_all  #keep as is
}

# this is the number of timesteps covered by the prior emissions
nbins <- length(unique(times_cut))

# Average the emissions to obtain mean vals for each time bin
emiss_bins <- array(0, dim = c(nlon, nlat, nbins))
for (ii in 1:nbins) {
    ibin <- which(times_cut_all == unique(times_cut)[ii])  #get indices of timesteps corresponding to this bin
    emiss_arr <- array(nc_emiss[, , ibin], dim = c(nlon, nlat, length(ibin)))  #this is needed in case t_res is hourly
    emiss_bins[, , ii] <- apply(emiss_arr, FUN = mean, MARGIN = c(1, 2))  #sum over bin's timesteps

}

# ~~~~~~~~~~~~~~ mask grid to include only domain cells ~~~~~~~~~~~~~~ #

# get lonlat pairs for the prior file
lonlat_prior <- expand.grid(nc_lon, nc_lat)

# load in lonlat_domain file specified in config.r
lonlat_domain <- readRDS(lonlat_domain_file)

# note which grid cells are in the domain
iDomain <- apply(lonlat_domain, FUN = function(x) which((lonlat_prior[, 1] == x[1]) &
    (lonlat_prior[, 2] == x[2])), MARGIN = 1)

# subset for the domain cells
sprior_mat <- apply(emiss_bins, FUN = function(x) x[iDomain], MARGIN = 3)

# convert to vector format
sprior <- as.vector(sprior_mat)

# ~~~~~~~~~~~~~~~~~~~~ Save sprior emissions file ~~~~~~~~~~~~~~~~~~~~~~~#

print("Saving formatted prior emissions to sprior file")

filepath <- paste0(out_path, "sprior.rds")
saveRDS(sprior, filepath)
