# Create temporal covariance matrix
# author: Lewis Kunik

## prerequisite scripts:
##  none
## output files:
##  tmp_cov.rds - file containing a (#times x #times) matrix of temporal covariance values

# load package dependencies
library(lubridate)

# run dependent scripts
source("config.r")


# ~~~~~~~~~~~~~~~~ Load time domain information ~~~~~~~~~~~~~~~~#

# The flux end time is 1 hour past what we need to capture, so discard last time
times <- seq(from = flux_start_POSIX, to = flux_end_POSIX - 3600, by = flux_t_res) #flux_t_res defined in config.r


# ~~~~~~~~~~~~~~~~ Calculate Xt and D ~~~~~~~~~~~~~~~~#

print("Creating temporal covariance matrix D")

# get time differences in seconds
t_diff <- sapply(times, FUN = function(x) abs(times - x))

# convert to days
Xt <- t_diff * (1/(24 * 60 * 60))

# Calculate D, using temporal correlation parameter defined in config.r
if (lt > 0) {
    D <- exp(-Xt/lt)
} else {
    D <- diag(rep(1, ntimes))
}

# only define correlation between times that have matching hour-of-day
flux_times_hr <- hour(times) #all unique hours of day
for (ii in 1:length(unique(flux_times_hr))) {

    this_hr <- unique(flux_times_hr)[ii]
    ithis_hr <- which(flux_times_hr == this_hr) #which other flux times are at this hour-of-day?
    inotthis_hr <- which(flux_times_hr != this_hr) #which other flux times AREN'T at this hour?
    icoords <- as.matrix(expand.grid(ithis_hr, inotthis_hr)) #define x/y index pairs of nonmatching hours
    D[icoords] <- 0
}

print("saving temporal covariance file")

# ~~~~~~~~~~~~~~~~ save temporal cov file ~~~~~~~~~~~~~~~~#

filepath <- paste0(out_path, "tmp_cov.rds")
saveRDS(D, filepath)
