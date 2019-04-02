## add_off_diags() Function declaration:
## add off-diagonals to the model data mismatch matrix for a particular site and
## a particular R-component
## Inputs:
# R - the existing model-data mismatch matrix you wish to add off-diagonals to
# site_vec - vector of site strings (or filepaths) corresponding to the receptor list
# datetimes - vector of POSIX times corresponding to the receptor list
# site - string corresponding to the site desired to add off-diagonals for
# rmse - the RMSE value of the particular R component
# correlation_scale - temporal correlation scale (in days) for this component
# non_decaying - boolean. Do you want to treat the correlation as non-decaying?
#                (i.e. rmse is constant down the off-diagonals)
# correlate_btwn_days - boolean. Do you want to extend the correlation of errors
#                       across all days (T) or limit to obs within same day (F)?


add_off_diags <- function(R, site_vec, datetimes, site, rmse, correlation_scale,
                          non_decaying = F, correlate_btwn_days = F){

    # get the indices on the diagonal that correspond to this site
    siteTF <- grepl(site, site_vec)

    # unique list of the dates that appear in the receptor list
    posix_dates <- unique(ISOdate(year(datetimes), month(datetimes),
                          day(datetimes), tz = "UTC"))

    # list of each date-time's matching days for this site
    isitedays <- lapply(posix_dates, FUN = function(x) which(year(datetimes) %in% year(x) &
                                        month(datetimes) %in% month(x) &
                                        day(datetimes) %in% day(x) & siteTF))

    # form the index-pairs matrix (will append to this next)
    match_pairs <- array(NA, dim = c(0,2))
    if(correlate_btwn_days){
      # if correlating between all days, get the site-matching indices
      match_pairs <- expand.grid(which(siteTF), which(siteTF))
    } else{
        # otherwise get those pairs which are only within the given day
        for(ii in 1:length(isitedays)){
            isiteday <- isitedays[[ii]]
            match_pairs <- rbind(match_pairs, expand.grid(isiteday, isiteday))
        }
    }

    # define the diagonal indices, which should NOT be altered
    idiag <- which(apply(match_pairs, FUN = function(x) x[1] == x[2], MARGIN = 1))
    irow <- match_pairs[-idiag, 1] #filter out diagonals
    icol <- match_pairs[-idiag, 2]
    ioff_diag <- cbind(irow, icol)

    if(!non_decaying){
      # get a matrix of nobs x nobs with time differences of receptors, in secs
      t_diff <- sapply(datetimes, FUN = function(x) abs(datetimes - x))
      R_Xt <- t_diff * (1/(24 * 3600)) #convert from seconds to days
      R_Dt <- exp(-R_Xt/correlation_scale) # establish the decay weights
      R_decay <- R_Dt * rmse^2 #mutliply decay weights by the rmse error
      R[ioff_diag] <- R[ioff_diag] + R_decay[ioff_diag] #add to off-diags of R
    } else{
      R[ioff_diag] <- R[ioff_diag] + rmse^2 #apply the unaltered RMSE value to the off-diagonals desired
    }
    #return R
    R
}
