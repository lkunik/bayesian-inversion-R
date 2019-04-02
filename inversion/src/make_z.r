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

make_z <- function(receptor_files,
                   receptor_times,
                   obs_mat,
                   bkgd,
                   aggregated_receptor_mat = NA) {



    nobs <- length(receptor_files)

    # convert receptor times from seconds-since-epoch to POSIX
    class(receptor_times) <- c("POSIXt", "POSIXct")
    attributes(receptor_times)$tzone <- "UTC"

    # get each site's indices in receptor list
    isites <- lapply(sites, FUN = function(x) grep(x, receptor_files))

    obs_times <- as.numeric(obs_mat[, 2])

    #convert obs times from seconds-since-epoch to POSIX
    class(obs_times) <- c("POSIXt", "POSIXct")
    attributes(obs_times)$tzone <- "UTC"


    # establish the vector for obs
    raw_obs <- rep(0, nobs)

    # nsites defined in config.r
    for (ii in 1:nsites) {

        sitename <- sites[ii]
        isite <- isites[[ii]]  #indices of this site's appearances in receptor list

        if (length(isite) == 0)
            next

        iobs_site <- which(obs_mat[, 1] == sitename)  #indices for this site within the obs file
        co2_site <- as.numeric(obs_mat[iobs_site, 3])  #available co2 vals for this site
        times_site <- obs_times[iobs_site]  #available times for co2 vals at this site

        # now filter obs for the footprints we have
        iz <- which(times_site %in% receptor_times[isite])
        raw_obs[isite] <- co2_site[iz]
    }

    # if toggled, aggregate the subsetted hours into 1 single daily value

    if (aggregate_obs) {

        nobs_new <- nrow(aggregated_receptor_mat)

        # load in the obs times from config.r
        y1 <- obs_year_start
        y2 <- obs_year_end

        m1 <- obs_month_start
        m2 <- obs_month_end

        d1 <- obs_day_start
        d2 <- obs_day_end

        h1 <- obs_hour_start
        h2 <- obs_hour_end

        mn1 <- obs_min_start
        mn2 <- obs_min_end

        # establish daily time-bins for aggregating
        time_bins_daily <- seq(from = ISOdatetime(y1, m1, d1, h1, mn1, 0, tz = "UTC"),
            to = ISOdatetime(y2, m2, d2, h2, mn2, 0, tz = "UTC") + 3600, by = 24 * 3600)

        # cut the original receptor times into 1-day bins
        times_cut_day <- as.POSIXct(cut(receptor_times, breaks = time_bins_daily), tz = "UTC")

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

        z <- obs_aggr - bkgd

    } else {
        z <- raw_obs - bkgd
    }

    # ~~~~~~~~~~~~~~~~~~~ save z to file ~~~~~~~~~~~~~~~~~~~#

    # return obs
    z


}
