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
source("src/read_sparse_h.r")


make_HQ <- function(sigma,
                    tmp_cov,
                    sp_cov,
                    lonlat_domain,
                    n_obs){

    #get info from inputs
    ncells <- nrow(lonlat_domain)
    D <- matrix(tmp_cov, nrow = ntimes)
    E <- matrix(sp_cov, nrow = ncells)

    # define a 'blank' HQHt which will be added to momentarily
    HQHt <- array(0, dim = c(n_obs, n_obs))

    # ~~~~~~~~~~~~ perform calculations ~~~~~~~~~~~~#
    # equations are based off of CT-Lagrange documentation:
    # https://www.esrl.noaa.gov/gmd/ccgg/carbontracker-lagrange/doc/intro.html#equation2

    # format sigma so it's easily accessible by timestep
    sigma_mat <- matrix(sigma, nrow = ntimes, byrow = T)

    path <- "H/"

    # Loop through timesteps - ntimes defined in config.r
    for (ii in 1:ntimes) {

        HDI <- array(0, dim = c(n_obs, ncells))

        # Add up all H slices multiplied by elements of D and diagonal sigma matrix
        for (jj in 1:ntimes) {
            # if temporal covariance for this timestep-pair is zero, continue
            if (D[jj, ii] == 0)
                next
            # get I_sigma corresponding to this timestep
            Isig_j <- diag(sigma_mat[jj, ])

            # load in the H file for time = jj
            H_file <- paste0(path, "H", formatC(jj, width = 3, flag = "0"), ".rds")
            Hjj <- read_sparse_h(H_file, n_obs, ncells)
            HDI <- HDI + (D[jj, ii] * Hjj %*% Isig_j) #add to total

        } #end inner timestep for-loop

        # Get the sigma section for this chunk, in diagonal matrix
        Isig_i <- diag(sigma_mat[ii, ])
        HQ <- HDI %*% E %*% Isig_i

        # read in the H file for time = ii
        H_file <- paste0(path, "H", formatC(ii, width = 3, flag = "0"), ".rds")
        Hii <- read_sparse_h(H_file, n_obs, ncells)

        # HQHT is cumulative over all timeslices, add to total
        HQHt <- HQHt + HQ %*% t(Hii)

        # ~~~~~~~~~~~~~~~ save HQ-slice file ~~~~~~~~~~~~~~~#
        filepath <- paste0("HQ/HQ", formatC(ii, width = 3, flag = "0"), ".rds")
        print(paste0("writing file ", filepath))
        saveRDS(HQ, filepath)

    } #end outer timestep for-loop

    # ~~~~~~~~~~~~~~~~~~~~~~~ save HQH^t file ~~~~~~~~~~~~~~~~~~~~~~~#
    #filepath <- paste0(out_path, "HQHt.rds")
    #saveRDS(HQHt, filepath)

    #return HQHt
    HQHt

}
