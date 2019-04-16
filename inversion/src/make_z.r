# calculate co2 enhancements (observed minus background)
# author: Lewis Kunik

## prerequisite scripts:
##  make_receptors.r
##  make_bg.r
##
## output files:
##  z.rds - file containing a vector of measured anthropogenic CO2 enhancement
##  from emissions within the domain

# run dependent scripts
source("config.r")


# load in receptor info
recep_filepath <- paste0(out_path, "receptors.rds")
receptors <- readRDS(recep_filepath)
recep_times <- as.numeric(receptors[, 1])
recep_names <- receptors[, 2]
nobs <- nrow(receptors)

# convert receptor times from seconds-since-epoch to POSIX
class(recep_times) <- c("POSIXt", "POSIXct")
attributes(recep_times)$tzone <- "UTC"

# get each site's indices in receptor list
isites <- lapply(sites, FUN = function(x) grep(x, recep_names))

# obs_file defined in config.r
obs <- readRDS(obs_file)
obs_times <- as.numeric(obs[, 2])

#convert obs times from seconds-since-epoch to POSIX
class(obs_times) <- c("POSIXt", "POSIXct")
attributes(obs_times)$tzone <- "UTC"

# load in the bkgd file
bg <- readRDS(paste0(out_path, "bg.rds"))

# establish the vector for obs
raw_obs <- rep(0, nobs)

# nsites defined in config.r
for (ii in 1:nsites) {

    sitename <- sites[ii]
    isite <- isites[[ii]]  #indices of this site's appearances in receptor list

    if (length(isite) == 0)
        next

    iobs_site <- which(obs[, 1] == sitename)  #indices for this site within the obs file
    co2_site <- as.numeric(obs[iobs_site, 3])  #available co2 vals for this site
    times_site <- obs_times[iobs_site]  #available times for co2 vals at this site

    # now filter obs for the footprints we have
    iz <- which(times_site %in% recep_times[isite])
    raw_obs[isite] <- co2_site[iz]
}

# if toggled, aggregate the subsetted hours into 1 single daily value

if (aggregate_obs) {

    # load aggregated receptor file
    recep_aggr_filepath <- paste0(out_path, "receptors_aggr.rds")
    receptors_aggr <- readRDS(recep_aggr_filepath)
    nobs_new <- nrow(receptors_aggr)

    # load in the obs times from config.r
    # establish daily time-bins for aggregating
    time_bins_daily <- seq(from = obs_start_POSIX, to = obs_end_POSIX + 3600, by = 24 * 3600)

    # cut the original receptor times into 1-day bins
    times_cut_day <- as.POSIXct(cut(recep_times, breaks = time_bins_daily), tz = "UTC")

    # this will hold our new, aggregated obs values
    obs_aggr <- rep(0, nobs_new)
    inew <- 0  #counter

    # ndays defined in config.r
    for (ii in 1:ndays) {

        # get this day's indices within the original receptor list
        iday <- which(times_cut_day == time_bins_daily[ii])

        # loop through all sites to get each site's aggregated obs for the day
        for (jj in 1:nsites) {

            isite <- isites[[jj]]  #find which receps correspond to this site
            idaysite <- isite[which(isite %in% iday)]  #get the overlap of these two - THIS site on THIS day

            # if num obs for this site on this day is below minimum defined in config.r, skip
            if (length(idaysite) < min_agg_obs)
              next

            inew <- inew + 1  #this is a counter, the same way that Hsplit does it

            # fill the new obs vector with your brand new mean value
            obs_aggr[inew] <- mean(raw_obs[idaysite])

        } # end 1:nsites for-loop
    } # end 1:ndays for-loop

    z <- obs_aggr - bg

} else {
    z <- raw_obs - bg
}

# ~~~~~~~~~~~~~~~~~~~ save z to file ~~~~~~~~~~~~~~~~~~~#

# save measured signals
print("saving anthro signals to z.rds")
filepath <- paste0(out_path, "z.rds")
saveRDS(z, filepath)
