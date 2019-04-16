# post_uncert_func.r
# author: Lewis Kunik

# This script defines a function that can be used by other scripts to calculate
# time-aggregated grid-scale posterior uncertainty, given start/end timesteps

## input variables:
##  tstart - beginning timestep over which to aggregate Q
##  tstop - ending timestep over which to aggregate Q
##  sp_cov - spatial covariance, matrix of dims (#cells x #cells)
##  tmp_cov - temporal covariance, matrix of dims (#times x #times)
##  sigma - prior uncertainty, vector of length (#times * #cells)
##  R - model data mismatch, matrix of dims (#obs x #obs)
##  HQHt - H * Q * t(H) from inversion equation, matrix of dims (#obs x #obs)
##  lonlat - matrix defining the grid cell-to-lon/lat conversion, number of rows = #cells
##
## output variables:
##  Vshat_bar - grid-scale aggregated posterior error covariance (aggregated over
##  timesteps specified by tstart & tstop), matrix of dims (#cells x #cells)


get_post_uncert_time <- function(tstart, tstop, sp_cov, tmp_cov, sigma, R, HQHt,
    lonlat) {

    # spatial and temporal covariance matrices
    E <- sp_cov
    D <- tmp_cov

    # ensure that D exists in matrix-conformable format, time x time, even if time = 1
    if (ntimes == 1)
        D <- as.matrix(D)

    # sigma must be formatted as num. times x num cells
    sigma_mat <- matrix(sigma, nrow = ntimes, byrow = T)
    nobs <- nrow(R)
    ncells <- nrow(lonlat)

    time_span <- tstop - tstart + 1

    # sigma sub-matrix denoted as S
    S <- sigma_mat[tstart:tstop, ]

    # ensure that S exists in matrix-conformable format, time x ncells, even if time = 1
    if (time_span == 1)
        S <- matrix(S, nrow = time_span, byrow = T)

    # ~~~~~~~~~~~~~~~~~~~~~~~ calculate Qsum ~~~~~~~~~~~~~~~~~~~~~~~#

    # make Qsum - formula shown at this link:
    # https://www.esrl.noaa.gov/gmd/ccgg/carbontracker-lagrange/doc/inversion.html#posterior-uncertainty
    # (equation 5 shows Qsum = (t(S) * D * S) * E, substituting x = t(S) * D.
    x <- t(S) %*% as.matrix(D[tstart:tstop, tstart:tstop])
    Qsum <- (x %*% S) * E  #Multiplication with E is an element-by-element multiplication, not matrix multiplication.


    # ~~~~~~~~~~~~~~~~~~~ calculate and return Vshat_bar ~~~~~~~~~~~~~~~~~ #

    # HQ_path = HQdir
    HQsum <- array(0, dim = c(nobs, ncells))
    # make HQsum
    for (ii in tstart:tstop) {
        HQ_file <- paste("HQ/HQ", formatC(ii, width = filename_width, flag = "0"), ".rds", sep = "")
        HQii <- readRDS(HQ_file)
        HQsum <- HQsum + HQii
    }

    # compute (HQHt + R)^-1
    HQHtRi <- solve(HQHt + R)

    # get k, number of timesteps - comes from config.r constant variables
    k <- (tstop - tstart) + 1
    k2 <- k^2

    # perform calculations from documentation
    HQT_sum_HQHtRi <- t(HQsum) %*% HQHtRi
    a <- HQT_sum_HQHtRi %*% HQsum
    Vshat_bar <- (Qsum - a)/k2

    # return Vshat_bar
    Vshat_bar

}
