# make HQ blocks and HQH-transpose calculations
# author: Lewis Kunik

## prerequisite scripts:
##  make_receptors.r
##  make_sigma.r
##  make_tmp_cov.r
##  make_sp_cov.r
##  Hsplit.r
##
## output files:
##  HQ/HQ*.rds - multiples files, 1 per inversion timestep, containing HQ matrix
##  slices each with dimension (#obs x #cells)
##  HQHt.rds - file containing a (#obs x #obs) matrix with the result of H * Q * t(H)

# run dependent scripts
source("config.r")


# ~~~~~~~~~~~~~ define a function to read sparse H file ~~~~~~~~~~~~~~~ #
read_sparse_h <- function(timestep, nobs, ncells) {

    Hi_tmp <- readRDS(paste0("H/H", formatC(timestep, width = filename_width, flag = "0"), ".rds"))
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



# ~~~~~~~~~~~~~~~~~~~~~~~ load in required files ~~~~~~~~~~~~~~~~~~~~~~~#

# load in receptor file
if (aggregate_obs) {
    recep_aggr_filepath <- paste0(out_path, "receptors_aggr.rds")
    receptors_aggr <- readRDS(recep_aggr_filepath)
    nobs <- nrow(receptors_aggr)
} else {
    recep_filepath <- paste0(out_path, "receptors.rds")
    receptors <- readRDS(recep_filepath)
    nobs <- nrow(receptors)
}

# load in prior uncertainty
sigma_file <- paste0(out_path, "sigma.rds")
sigma <- readRDS(sigma_file)

# load in temporal covariance file
tmp_cov_file <- paste0(out_path, "tmp_cov.rds")
D_vals <- readRDS(tmp_cov_file)
D <- matrix(D_vals, nrow = ntimes)

# load in spatial covariance file
sp_cov_file <- paste0(out_path, "sp_cov.rds")
E <- readRDS(sp_cov_file)

# load in lonlat_domain file specified in config.r
lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)


# define a 'blank' HQHt which will be added to momentarily
HQHt <- array(0, dim = c(nobs, nobs))

# ~~~~~~~~~~~~ perform calculations ~~~~~~~~~~~~#
# equations are based off of CT-Lagrange documentation:
# https://www.esrl.noaa.gov/gmd/ccgg/carbontracker-lagrange/doc/intro.html#equation2

# format sigma so it's easily accessible by timestep
sigma_mat <- matrix(sigma, nrow = ntimes, byrow = T)

# Loop through timesteps - ntimes defined in config.r
for (ii in 1:ntimes) {

    HDI <- array(0, dim = c(nobs, ncells))

    # Add up all H slices multiplied by elements of D and diagonal sigma matrix
    for (jj in 1:ntimes) {
        # if temporal covariance for this timestep-pair is zero, continue
        if (D[jj, ii] == 0)
            next
        # get I_sigma corresponding to this timestep
        Isig_j <- diag(sigma_mat[jj, ])

        # load in the H file for time = jj
        Hjj <- read_sparse_h(jj, nobs, ncells)
        HDI <- HDI + (D[jj, ii] * Hjj %*% Isig_j) #add to total

    } #end inner timestep for-loop

    # Get the sigma section for this chunk, in diagonal matrix
    Isig_i <- diag(sigma_mat[ii, ])
    HQ <- HDI %*% E %*% Isig_i

    # read in the H file for time = ii
    Hii <- read_sparse_h(ii, nobs, ncells)

    # HQHT is cumulative over all timeslices, add to total
    HQHt <- HQHt + HQ %*% t(Hii)

    # ~~~~~~~~~~~~~~~ save HQ-slice file ~~~~~~~~~~~~~~~#
    filepath <- paste0("HQ/HQ", formatC(ii, width = filename_width, flag = "0"), ".rds")
    print(paste0("writing file ", filepath))
    saveRDS(HQ, filepath)

} #end outer timestep for-loop


# ~~~~~~~~~~~~~~~~~~~~~~~ save HQH^t file ~~~~~~~~~~~~~~~~~~~~~~~#
filepath <- paste0(out_path, "HQHt.rds")
saveRDS(HQHt, filepath)
