# assemble background values in inversion file format
# author: Lewis Kunik

## prerequisite scripts:
##  make_receptors.r
##  make_Hs_bkgd.r
##
## output files:
##  bg.rds - file containing a vector of background values corresponding to the
##  sites/times specified in 'receptors.rds'

# run dependent scripts
source("config.r")

make_bg <- function(receptor_files,
                    receptor_times,
                    bkgd_vals,
                    bkgd_times,
                    Hs_outer = NA,
                    Hs_bio = NA,
                    aggregated_receptor_mat = NA){

    class(receptor_times) <- c("POSIXt", "POSIXct")
    attributes(receptor_times)$tzone <- "UTC"
    nobs <- length(receptor_files) #total number of observations

    # grab indices of all sites' footprints/obs in the receptor list
    isites <- lapply(sites, FUN = function(x) grep(x, recep_names))

    class(bkgd_times) <- c("POSIXt", "POSIXct")
    attributes(bkgd_times)$tzone <- "UTC"

    # establish the vector for obs. nsites loaded from config.r
    bg <- rep(0, nobs)
    for (ii in 1:nsites) {

        # get sitename and site indices
        sitename <- sites[ii]
        isite <- isites[[ii]]  #isite corresponds to the idxs of this site in the RECEPTOR LIST (not the obs file)

        # skip if no obs whatsoever
        if (length(isite) == 0)
            next

        # now filter base bkgd for times corresponding to this site's receptors
        iz <- which(bkgd_times %in% recep_times[isite])
        bg_z <- bkgd_vals[iz]

        # get bg for this site and add to the vector
        bg[isite] <- bg_z

    }

    # now, if toggled, aggregate the subsetted hours into 1 single daily value. Do
    # this in the same way that H and R are aggregated based on their obs indexes

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
        times_cut_day <- as.POSIXct(cut(recep_times, breaks = time_bins_daily), tz = "UTC")

        # this will hold our new, aggregated bg values
        bg_aggr <- rep(0, nobs_new)
        inew <- 0  #counter
        for (ii in 1:ndays) {

            # get the idx's of which day this corresponds to within the original receptor list
            iday <- which(times_cut_day == time_bins_daily[ii])

            # loop through all sites to get each site's aggregated bgval for the day
            # nsites loaded from config.r
            for (jj in 1:nsites) {

                isite <- isites[[jj]]  #find which receps correspond to this site
                idaysite <- isite[which(isite %in% iday)]  #get the overlap of these two - THIS site on THIS day

                # if num obs for this site on this day is below minimum defined in config.r, skip
                if (length(idaysite) < min_agg_obs)
                  next

                inew <- inew + 1  #this is a counter, the same way that Hsplit does it

                # fill the new bg vector with your brand new mean background value
                bg_aggr[inew] <- mean(bg[idaysite])
            }
        }

        bg_all <- bg_aggr

    } else {

        bg_all <- bg

    }

    # if including far-field contributions from outside the optimization domain,
    # include these here
    if (include_outer) {

      # load Hs_outer.rds
      Hs_outer_file <- paste0(out_path, "Hs_outer.rds")
      Hs_outer <- readRDS(Hs_outer_file)

      bg_all <- bg_all + Hs_outer #add to background
    }

    # if including contributions from biogenic fluxes, include these here
    if (include_bio) {

      # load Hs_bio.rds
      Hsbio_file <- paste0(out_path, "Hs_bio.rds")
      Hsbio <- readRDS(Hsbio_file)

      bg_all <- bg_all + Hsbio #add to background
    }

    # save measured signals
    #print("saving bkgd vals to bg.rds")
    #filepath <- paste0(out_path, "bg.rds")
    #saveRDS(bg_all, filepath)

    bg_all

}
