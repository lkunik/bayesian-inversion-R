# filter available footprints for the timespan defined in config.r, and for
# available observations in the obs/bkgd files
# author: Lewis Kunik

## prerequisite scripts:
##  none (but obs.rds and background.rds must be present in include directory)
##
## output files:
##  receptors.rds - contains an (n x 2) matrix, col 1 has time of each receptor,
##  col 2 has filepath of each receptor

# load package dependencies
library(lubridate)

# run dependent scripts
source("config.r")

# ~~~~~~~~~~~~~~~~~~~ set up time parameters ~~~~~~~~~~~~~~~~~~~#

# list all desired obs times
master_times_all <- seq(from = obs_start_POSIX, to = obs_end_POSIX, by = obs_t_res)  # obs_t_res defined in config.r

# remove times outside of the desired subset times
isubset <- which(hour(master_times_all) %in% subset_hours_utc)
master_times <- master_times_all[isubset]

# load in obs file
obs <- readRDS(obs_file)

# load in background file
bkgd <- readRDS(bg_file)
ibg_nNA <- which(!is.na(as.numeric(bkgd[,2])))
bkgd_times <- as.numeric(bkgd[ibg_nNA, 1])
class(bkgd_times) <- c("POSIXt", "POSIXct")
attributes(bkgd_times)$tzone <- "UTC"

# set up variables
recep_list_unsorted <- array(NA, dim = c(0, 1))
times_unsorted <- array(NA, dim = c(0, 1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### For each site, determine which of these times we have footprints for
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# refer to sites object in config.r
for (ii in 1:length(sites)) {

    # get the site code
    site <- sites[ii]

    # get path to the site's footprint directory
    site_dir <- paste0(foot_dir, site, "/")

    # Obtain list of foot files in this site's directory
    foot_files <- list.files(site_dir)

    # get the POSIX times of all the available footprints in the site's directory
    foot_times <- as.POSIXct(strptime(substr(foot_files, 1, 10), format = "%Y%m%d%H",
        tz = "UTC"))

    ## now filter the footprint files' times for the subsetted times of day

    # format the desired times into YYYYmmddHH format that exists as substring in the
    # footprint filenames
    times_file_fmt <- format(master_times, format = "%Y%m%d%H", tz = "UTC")

    # filter the DESIRED subsetted times for those which exist in the footprint
    # directory (they should all be there, but this is precautionary)
    isubset_foots <- which(sapply(times_file_fmt, FUN = function(x) any(grepl(x,
        foot_files))))
    # get TF-mask of footprint files that match our DESIRED subsetted times
    footsTF <- foot_times %in% master_times[isubset_foots]

    # IMPORTANT! We don't always have data for the desired subsetted times in our
    # inversion, so we need to check which of these subsetted times we DO have data for
    isite <- which(obs[, 1] == site)  #indices of obs that match this site
    obs_times <- as.numeric(obs[isite, 2])  #get the master list of POSIX times for this site's available data

    #convert obs_times from seconds-since-epoch to POSIX
    class(obs_times) <- c("POSIXt", "POSIXct")
    attributes(obs_times)$tzone <- "UTC"

    # get TF-mask of this site's data, and bkgd data, that match the footprint times
    # we have
    obs_existTF <- (foot_times %in% obs_times) & (foot_times %in% bkgd_times)

    # get the indices of the OVERLAP between footprints in the dir that match desired
    # subsetted times, and footprints in the dir that exist in the DATA
    iOverlap <- which(footsTF & obs_existTF)

    # if no footprints overlap with obs and bkgd, skip to next
    if (length(iOverlap) == 0)
        next

    # add foot files and receptor times to the running lists of receptors/times
    recep_list_unsorted <- c(recep_list_unsorted, paste0(site_dir, foot_files[iOverlap]))
    times_unsorted <- c(times_unsorted, foot_times[iOverlap])



}

#sort the receptor list by time (sites fall into order listed in config.r)
time_list <- sort(times_unsorted)
itime_sort <- order(times_unsorted)
recep_list <- recep_list_unsorted[itime_sort]


# combine the times and filepaths to get the var we want to save
recep_list_w_dates <- cbind(time_list, recep_list)
colnames(recep_list_w_dates) <- c("seconds_since_1970_01_01", "file_path")

# save the receptor file as 2-D array of receptor times and files (times are
# saved as seconds since the R epoch, 1970-01-01 00:00:00)
filepath <- paste0(out_path, "receptors.rds")
saveRDS(recep_list_w_dates, filepath)


# ~~~~~~~~~~~~~~~~ Now aggregate receptor obs if desired ~~~~~~~~~~~~~~~~ #

if (aggregate_obs) {
    # get a list of each site's indices within the receptor list we've created
    isites <- lapply(sites, FUN = function(x) grep(x, recep_list))
    recep_times <- as.numeric(time_list)
    class(recep_times) <- c("POSIXt", "POSIXct")
    attributes(recep_times)$tzone <- "UTC"
    # make the new *aggregated* receptor list, to be used later on

    # prepare to aggregate footprint files based on day
    time_bins_daily <- seq(from = obs_start_POSIX, to = obs_end_POSIX + 3600, by = 24 * 3600)
    times_cut_day <- as.POSIXct(cut(recep_times, breaks = time_bins_daily), tz = "UTC")

    # establish the new receptor list
    new_recep_list <- array(NA, dim = c(0, 2))
    colnames(new_recep_list) <- c("site", "day")

    #ndays and nsites loaded from config.r
    for (ii in 1:ndays) {
        for (jj in 1:nsites) {

            isite <- isites[[jj]]  #find which receps correspond to this site
            iday <- which(times_cut_day == time_bins_daily[ii])  #find which recep times correspond to this day
            idaysite <- isite[which(isite %in% iday)]  #get the overlap of these two - THIS site on THIS day

            # if num obs for this site on this day is below minimum defined in config.r, skip
            if (length(idaysite) < min_agg_obs)
              next

            # save the site/day in order so we can later reference the site/day ordering for
            # the R matrix and obs vectors get the sitename for this site
            sitename <- sites[jj]
            new_recep_list <- rbind(new_recep_list, c(sitename, time_bins_daily[ii]))
        }
    }

    # save the new receplist
    filepath <- paste0(out_path, "receptors_aggr.rds")
    saveRDS(new_recep_list, filepath)

}
