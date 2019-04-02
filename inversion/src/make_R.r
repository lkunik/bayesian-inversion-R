# Assemble the model data mismatch matrix
# author: Lewis Kunik

## prerequisite scripts:
##  make_receptor.s
##  make_z.rds
##
## output files:
##  R.rds - file containing a (#obs x #obs) matrix with values of model data mismatch
## (error covariance of observational constraints)

#load required libraries
library(lubridate)

# run dependent scripts
source("src/add_off_diags_R.r")
source("src/aggregate_R.r")

make_R <- function(receptor_files,
                   receptor_times,
                   lonlat_domain,
                   receptors_aggr = NA) {

    source("config_R_uncert.r")

    nobs <- length(receptor_files)

    # convert receptor times from seconds-since-epoch to POSIX
    class(receptor_times) <- c("POSIXt", "POSIXct")
    attributes(receptor_times)$tzone <- "UTC"


    # get each site's indices in receptor list
    isites <- lapply(sites, FUN = function(x) grep(x, receptor_files))

    ncells <- nrow(lonlat_domain)

    # Model Data Mismatch is made up of multiple parts R = R_part + R_aggr + R_eddy +
    # R_bkgd + R_transPBL + R_transWIND + R_instr + R_bio + R_other


    # ~~~~~~~~~~~~~~~ R_part ~~~~~~~~~~~~~~~#
    R_part <- rep(rmse_part^2, nobs)


    # ~~~~~~~~~~~~~~~ R_aggr ~~~~~~~~~~~~~~~#
    R_aggr <- rep(rmse_aggr^2, nobs)


    # ~~~~~~~~~~~~~~~ R_eddy ~~~~~~~~~~~~~~~#
    R_eddy <- rep(rmse_eddy^2, nobs)


    # ~~~~~~~~~~~~~~~ R_bkgd ~~~~~~~~~~~~~~~#
    R_bkgd <- rep(rmse_bkgd^2, nobs)


    # ~~~~~~~~~~~~~~~ R_transPBL ~~~~~~~~~~~~~~~#
    R_transPBL <- rep(rmse_transPBL^2, nobs)


    # ~~~~~~~~~~~~~~~ R_transWIND ~~~~~~~~~~~~~~~#
    R_transWIND <- rep(rmse_transWIND^2, nobs)


    # ~~~~~~~~~~~~~~~ R_instr ~~~~~~~~~~~~~~~#
    R_instr <- rep(rmse_instr^2, nobs)

    # if extra instrument error for specific sites
    for(ii in 1:nsites){
      extra_site_rmse <- extra_sites_err_rmse[[ii]]
      isite <- isites[[ii]]
      R_instr[isite] = extra_site_rmse^2#R_instr[isite] + extra_site_rmse^2
    }


    # ~~~~~~~~~~~~~~~ R_bio ~~~~~~~~~~~~~~~#
    R_bio <- rep(rmse_bio^2, nobs)


    # ~~~~~~~~~~~~~~~ R_other ~~~~~~~~~~~~~~~#
    R_other <- rep(rmse_other^2, nobs)


    # ~~~~~~~~~~~~~~~ Add up all components of R ~~~~~~~~~~~~~~~#

    R <- diag(R_part + R_eddy + R_aggr + R_bkgd + R_transPBL +
                R_transWIND + R_instr + R_bio + R_other)


    # ~~~~~~ Now add covariance (off-diagonals) ~~~~~~~#
    ncomponents <- nrow(R_arr)

    for(ii in 1:ncomponents){

        rmse_in <- as.numeric(R_arr[ii, "RMSE"])
        cor_scale_in <- as.numeric(R_arr[ii, "correlation_scale"])
        non_decaying_in <- as.logical(R_arr[ii, "non_decaying?"])
        cor_btwn_days_in <- as.logical(R_arr[ii, "correlate_btwn_days?"])

        if(is.na(cor_scale_in) & !non_decaying_in)
          next

        for(site in sites){

            R <- add_off_diags(R, receptor_files, receptor_times, site, rmse_in, cor_scale_in,
                                      non_decaying_in, cor_btwn_days_in)
        }
    }


    # now aggregate into daily observational errors if aggregate_obs is toggled

    if (aggregate_obs) {

        recep_sites_aggr <- receptors_aggr[,1]
        recep_times_aggr <- receptors_aggr[,2]

        class(recep_times_aggr) <- c("POSIXt", "POSIXct")
        attributes(recep_times_aggr)$tzone <- "UTC"

        # turn into diagonal n x n matrix
        R_aggregated <- aggregate_R(R, receptor_files, receptor_times, recep_sites_aggr, recep_times_aggr)

        # save R to file
        #print("saving model data mismatch matrix to R.rds")
        #filepath <- paste0(out_path, "R.rds")
        #saveRDS(R_aggregated, filepath)

        R_return <- R_aggregated

    } else {

        # save R to file
        #print("saving model data mismatch matrix to R.rds")
        #filepath <- paste0(out_path, "R.rds")
        #saveRDS(R, filepath)
        R_return <- R

    }

    #return the R file
    R_return

}
