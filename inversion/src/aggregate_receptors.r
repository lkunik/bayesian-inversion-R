


aggregate_receptors <- function(receptor_files,
                                receptor_times){


    # get a list of each site's indices within the receptor list we've created
    isites <- lapply(sites, FUN = function(x) grep(x, receptor_files))

    class(receptor_times) <- c("POSIXt", "POSIXct")
    attributes(receptor_times)$tzone <- "UTC"
    # make the new *aggregated* receptor list, to be used later on

    # prepare to aggregate footprint files based on day
    time_bins_daily <- seq(from = ISOdatetime(y1, m1, d1, h1, mn1, 0, tz = "UTC"),
        to = ISOdatetime(y2, m2, d2, h2, mn2, 0, tz = "UTC") + 3600, by = 24 * 3600)

    times_cut_day <- as.POSIXct(cut(receptor_times, breaks = time_bins_daily), tz = "UTC")

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
    #filepath <- paste0(out_path, "receptors_aggr.rds")
    #saveRDS(new_recep_list, filepath)

    #return aggregated receptor list
    new_recep_list

}
