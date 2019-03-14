# convolve outer domain footprints with far-field anthro and bio emissions
# author: Lewis Kunik

## prerequisite scripts:
##  make_receptors.r
##  make_outer.r
##  make_sbio.r
##  Hsplit.r
##
## output files:
##  Hs_outer.rds - far-field emission contributions from outer_emiss.nc file,
##  corresponding to observations from the receptor file
##  Hs_bio.rds - biogenic flux contributions from bio_flux.nc file,
##  corresponding to observations from the receptor file

# run dependent scripts
source("config.r")


# ~~~~~~~~~~~~~~~~ define a function to read sparse H file ~~~~~~~~~~~~~~~~~ #
read_sparse_h_outer <- function(timestep, nobs, ncells) {

    Hi_tmp <- readRDS(paste0("H_outer/H", formatC(timestep, width = 3, flag = "0"), ".rds"))
    # Populate the H-slice matrix (nobs x ncells) with zeros
    Hi <- array(0, dim = c(nobs, ncells))

    # checking if Hi is empty
    if (length(Hi_tmp[, 1]) > 0) {
        # for every value in H file, locate where that exists and put a value there
        iobs <- Hi_tmp[, 1]
        icell <- Hi_tmp[, 2]
        Hi[cbind(iobs, icell)] <- Hi_tmp[, 3]
    }

    # return H timestep
    return(Hi)
}


# ~~~~~~~~~~~~~~~ Load receptor and grid info ~~~~~~~~~~~~~~~~#

#load domain lonlat pairs
lonlat_domain <- readRDS(lonlat_domain_file)

# lonlat_outer is required if this script is running. load outer lonlat pairs
lonlat_outer <- readRDS(lonlat_outer_file)
ncells_outer <- nrow(lonlat_outer)

iDomain <- apply(lonlat_domain, FUN = function(x) which((lonlat_outer[, 1] %in% x[1]) &
    (lonlat_outer[, 2] %in% x[2])), MARGIN = 1)

if (aggregate_obs) {
    receps_aggr_file <- paste0(out_path, "receptors_aggr.rds")
    receps_aggr <- readRDS(receps_aggr_file)
    nobs <- nrow(receps_aggr)
} else {
    receps_file <- paste0(out_path, "receptors.rds")
    receps <- readRDS(receps_file)
    nobs <- nrow(receps)
}

# ~~~~~~~~~~~~~~~ Load emissions files ~~~~~~~~~~~~~~~#

if (include_outer) {
    outer_file <- paste0(out_path, "outer.rds")
    outer_vec <- readRDS(outer_file)
    outer_mat <- matrix(outer_vec, nrow = ntimes, byrow = T)

    # mask out the inner domain cells, so that the outer emiss contribution only
    # considers outer-domain cells
    outer_mat <- outer_mat
    outer_mat[, iDomain] <- 0

    #set up the convolved vector to add to
    Hs_outer <- rep(0, nobs)
}

if (include_bio) {
    sbio_file <- paste0(out_path, "sbio.rds")
    sbio_vec <- readRDS(sbio_file)
    sbio_mat <- matrix(sbio_vec, nrow = ntimes, byrow = T)

    # set up the convolved vector to add to
    Hsbio <- rep(0, nobs)
}

# convolve each H file and add to running total of Hsbio
for (ii in 1:ntimes) {

    # get H far-field file for time = ii
    Hi <- read_sparse_h_outer(ii, nobs, ncells_outer)

    # add outer outer domain contribution, if toggled
    if (include_outer)
        Hs_outer <- Hs_outer + Hi %*% outer_mat[ii, ]

    # add biospheric contributions to observations, if toggled
    if (include_bio)
        Hsbio <- Hsbio + Hi %*% sbio_mat[ii, ]

}  #end ntimes for-loop


# save the files as vectors of enhancements for each observation in obs space
if (include_outer) {

    filepath <- paste0(out_path, "Hs_outer.rds")
    saveRDS(Hs_outer, filepath)

}

if (include_bio) {

    filepath <- paste0(out_path, "Hs_bio.rds")
    saveRDS(Hsbio, filepath)
    
}
