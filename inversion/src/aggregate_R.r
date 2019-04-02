
## aggregate_R Function declaration:
## aggregate the R matrix to daily observations (errors are generally reduced because
## 1 daily representation of mixing ratio is more accurate than individual obs)
## Inputs:
# R - the existing model-data mismatch matrix you wish to aggregate
# site_vec - vector of site strings corresponding to the receptor list
# datetimes - vector of POSIX times corresponding to the receptor list

aggregate_R <- function(R, site_vec, datetimes, site_vec_aggr, dates_aggr){

    nobs <- length(site_vec) #get number of total obs, pre-aggregation
    nobs_aggr <- length(site_vec_aggr) #number of total obs after aggregation
    R_aggr_diag <- rep(0, nobs_aggr) #this will hold the new aggregated R

    # vector of unique dates that appear in the receptor list
    posix_dates <- unique(ISOdate(year(datetimes), month(datetimes),
                          day(datetimes), tz = "UTC"))

    ndays <- length(posix_dates)

    # loop through the sites
    for(site in unique(site_vec_aggr)){

        # get the indices on the diagonal that correspond to this site
        siteTF <- grepl(site, site_vec)
        siteTF_aggr <- grepl(site, site_vec_aggr)

        # list of each date-time's matching days for this site
        isitedays <- lapply(posix_dates, FUN = function(x) which(year(datetimes) %in% year(x) &
                                            month(datetimes) %in% month(x) &
                                            day(datetimes) %in% day(x) & siteTF))

        # same as above, but for the aggregated matrix
        isitedays_aggr <- lapply(posix_dates, FUN = function(x) which(year(dates_aggr) %in% year(x) &
                                            month(dates_aggr) %in% month(x) &
                                            day(dates_aggr) %in% day(x) & siteTF_aggr))


        # loop through all days and aggregate
        for (ii in 1:ndays) {

                idaysite <- isitedays[[ii]] #indices of this site's obs on this day
                iaggr <- isitedays_aggr[[ii]]

                # if num obs for this site on this day is below minimum defined in config.r, skip
                if (length(idaysite) < min_agg_obs)
                  next

                # establish an aggregation operator for this site's obs on this day
                # in general (if no missing data/footprints) it's 1/length(subset_hours_utc)
                W <- rep(0, nobs)
                W[idaysite] <- 1/length(idaysite)
                new_diag <- W %*% R #flatten R into a vector, where only vals from this site/day are non-zero

                R_chunk <- new_diag[idaysite]

                # NOTE: take mean
                R_aggr_diag[iaggr] <- mean(R_chunk)
            }

    } #end sites for-loop

    R_aggregated <- diag(R_aggr_diag)

    #return R_aggr
    R_aggregated
}
