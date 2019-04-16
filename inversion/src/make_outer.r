# format outer-domain emissions into vector of flux values
# author: Lewis Kunik

## prerequisite scripts:
##  none (but outer_emiss.nc must be present in the include/ directory)
##
## output files:
##  outer.rds - file containing a vector of length (#outer domain cells * #times)
##  with values of outer-domain emissions for every timestep

# load package dependencies
library(ncdf4)

# run dependent scripts
source("config.r")


# ~~~~~~~~~~~~~~~~~~~~~~~ create time bins ~~~~~~~~~~~~~~~~~~~~~~~#

# establish the markers for the time-bin cutoffs (includes an end-timestamp after
# the last flux hour for the purpose of the cut() function)
time_bins <- seq(from = flux_start_POSIX, to = flux_end_POSIX, by = flux_t_res) #flux_t_res defined in config.r



# ~~~~~~~~~~~~~~~~~~~~~~~ load in Outer emissions file ~~~~~~~~~~~~~~~~~~~~~~~#

print("Loading outer emissions file")

# load ncdf file from outer filepath
nc_outer <- nc_open(outer_file)

# get lon, lat, time, and flux from the netcdf file
nc_lat <- round(ncvar_get(nc_outer, "lat"), round_digs)
nc_lon <- round(ncvar_get(nc_outer, "lon"), round_digs)
nc_time <- ncvar_get(nc_outer, "time")
nc_emiss <- ncvar_get(nc_outer, "emiss")
nc_close(nc_outer)

#convert time field from seconds-since-epoch to POSIX time
class(nc_time) <- c("POSIXt", "POSIXct")
attributes(nc_time)$tzone <- "UTC"

# NOTE: outer emissions are intended to have dimensions (lon x lat x time)

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

#filter for NA values
if (any(is.na(times_cut_all))) {
    times_cut <- times_cut_all[-(which(is.na(times_cut_all)))]  #remove NA values
} else {
    times_cut <- times_cut_all  #keep as is
}

# this is the number of timesteps covered by the outer emissions
nbins <- length(unique(times_cut))

# Average the emissions to obtain mean vals for each time bin
emiss_bins <- array(0, dim = c(nlon, nlat, nbins))
for (ii in 1:nbins) {
    ibin <- which(times_cut_all == unique(times_cut)[ii])  #get indices of timesteps corresponding to this bin
    emiss_arr <- array(nc_emiss[, , ibin], dim = c(nlon, nlat, length(ibin)))  #this is needed in case t_res is hourly
    emiss_bins[, , ii] <- apply(emiss_arr, FUN = mean, MARGIN = c(1, 2))  #sum over bin's timesteps
}

# ~~~~~~~~~~~~~~ mask grid to include only domain cells ~~~~~~~~~~~~~~ #

# get lonlat pairs for the large outer-domain emissions file
lonlat_outer_emiss <- expand.grid(nc_lon, nc_lat)

# load in lonlat_outer (subset of the outer-emissions grid) file specified in config.r
lonlat_outer <- readRDS(lonlat_outer_file)

# Efficient way to subset for the outer-domain when both grids are large.
# create a rectangular index map of the loaded emissions, find the inversion's
# "outer domain" NESW bounds, and grab the indices within that inner rectangular.
# NOTE: this requires that the outer-domain bounds are rectangular
iouter_big <- array(1:nrow(lonlat_outer_emiss), dim = c(nlon, nlat))
sub_lon <- unique(lonlat_outer[, 1])
sub_lat <- unique(lonlat_outer[, 2])
imin_lat <- which(nc_lat == min(sub_lat)) #Southern bound
imax_lat <- which(nc_lat == max(sub_lat)) #Northern bound
imin_lon <- which(nc_lon == min(sub_lon)) #Western bound
imax_lon <- which(nc_lon == max(sub_lon)) #Eastern bound

#subset for that inner rectangle
iemiss_outer <- iouter_big[imin_lon:imax_lon, imin_lat:imax_lat]  #indexed as [lon, lat]
iouter <- as.vector(iemiss_outer)

# subset for the far-field domain cells
outer_mat <- apply(emiss_bins, FUN = function(x) x[iouter], MARGIN = 3)

# convert to vector format
outer <- as.vector(outer_mat)

# ~~~~~~~~~~~~~~~~~~~~~~~ save to file ~~~~~~~~~~~~~~~~~~~~~~~#

print("saving formatted outer domain emissions file")
filepath <- paste0(out_path, "outer.rds")
saveRDS(outer, filepath)
