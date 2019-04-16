# calculate the reduced chi-squared metric for the inversion
# author: Lewis Kunik

## Chi-squared equation used is taken from Tarantola, 1987. 'Inverse Problem
## Theory: Methods for Data Fitting & Model Parameter Estimation' This version of
## chi-squared calculation attempts to utilize sparse-matrix methods as in
## 'make_HQ.r'

## prerequisite scripts:
## make_sprior.r
## make_sigma.r
## make_sp_cov.r
## make_tmp_cov.r
## make_z.r
## make_R.r
## Hsplit.r
## make_HQ.r
## inversion.r
##
## output files:
## chi_sq.rds - file containing a single value of the reduced chi-squared fit,
##              combining observational and emissions uncertainty constraints

# run dependent scripts
source("config.r")


# ~~~~~~~~~~~~~~~~~~~ load required files ~~~~~~~~~~~~~~~~~~~#

# load z
z_file <- paste0(out_path, "z.rds")
z <- readRDS(z_file)
nobs <- length(z)

# load R
R_file <- paste0(out_path, "R.rds")
R <- readRDS(R_file)

# load sp_cov
sp_cov_file <- paste0(out_path, "sp_cov.rds")
E <- readRDS(sp_cov_file)

# load tmp_cov
tmp_cov_file <- paste0(out_path, "tmp_cov.rds")
D <- as.matrix(readRDS(tmp_cov_file))

# load sigma
sigma_file <- paste0(out_path, "sigma.rds")
sigma <- readRDS(sigma_file)
sigma_mat <- matrix(sigma, nrow = ntimes, byrow = T)

# load s_hat
s_hat_file <- paste0(out_path, "s_hat.rds")
s_hat <- readRDS(s_hat_file)
s_hat_mat <- matrix(s_hat, nrow = ntimes, byrow = T)
shat_avg <- apply(s_hat_mat, FUN = mean, MARGIN = 2)

# load sprior
sprior_file <- paste0(out_path, "sprior.rds")
sprior <- readRDS(sprior_file)
sprior_mat <- matrix(sprior, nrow = ntimes, byrow = T)
nflux <- length(s_hat)

# get grid info
lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)

# define number of degrees of freedom, v
v <- nobs

# ~~~~~~~~~~ calculate Hshat and Hsprior ~~~~~~~~~~#

# making Hshat
Hshat <- rep(0, nobs)
Hsprior <- rep(0, nobs)

#loop through timesteps and convolve fluxes with footprints
for (ii in 1:ntimes) {
    H_tmp <- readRDS(paste0("H/H", formatC(ii, width = filename_width, flag = "0"), ".rds"))
    Hi <- array(0, dim = c(nobs, ncells))

    # checking if Hi_outer is empty
    if (length(H_tmp[, 1]) > 0) {
        # for every value in H file, locate where that exists and put that value there
        iobs <- H_tmp[, 1]
        icell <- H_tmp[, 2]
        Hi[cbind(iobs, icell)] <- H_tmp[, 3]
    }

    #continuously add to the observation-sized vector
    Hshat <- Hshat + Hi %*% s_hat_mat[ii, ]
    Hsprior <- Hsprior + Hi %*% sprior_mat[ii, ]
}



# ~~~~~~~~~~~~~~~~~~~ calculate reduced Chi^2 obs ~~~~~~~~~~~~~~~~~~~#

# calculate z - Hshat
zHshat <- z - Hshat

# take R^-1
R_inv <- solve(R)

# calculate chi-squared
Chi_sq_z <- (t(zHshat) %*% R_inv %*% zHshat)


# ~~~~~~~~~~~~~~~~~~~ calculate Chi^2 flux ~~~~~~~~~~~~~~~~~~~#
shatsp <- s_hat - sprior
shatsp_mat <- matrix(shatsp, nrow = ntimes, byrow = T)
# the following method is similar to the sparse matrix method in make_HQ.r but
# using state vector residuals instead of H, and improvising a Q^-1 from its
# components E, D, and sigma

D_inv <- solve(D)
E_inv <- solve(E)
Chi_sq_s <- 0  #this is a variable to hold the value (shat - sp)^T * Q^-1 * (shat - sp)

# this takes some time, so utilize print statements
print("calculating reduced Chi-squared using sparse-matrix method")

# Loop through timesteps
for (ii in 1:ntimes) {
    print(paste("adding timestep", ii))

    sDI <- rep(0, ncells)

    # Add up all H slices multiplied by elements of D and diagonal sigma matrix
    for (jj in 1:ntimes) {
        # if temporal covariance for this timestep-pair is zero, continue
        if (D_inv[jj, ii] == 0)
            next
        # get the INVERSE of I_sigma corresponding to this timestep NOTE: because I_sigma
        # is diagonal, we can use pieces of its inverse by taking the inverse of the
        # pieces themselves
        Isig_j <- solve(diag(sigma_mat[jj, ]))
        sjj <- shatsp_mat[jj, ]
        sDI <- sDI + (D_inv[jj, ii] * sjj %*% Isig_j)
    }

    # Get the INVERSE of sigma section for this chunk
    Isig_i <- solve(diag(sigma_mat[ii, ]))
    sQi <- sDI %*% E_inv %*% Isig_i

    # Chi_sq_s is cumulative over all timeslices
    sii <- shatsp_mat[ii, ]
    Chi_sq_s <- Chi_sq_s + sQi %*% sii

}


# calculate Chi-square statistic for combined obs and flux parameters
Chi_sq <- (Chi_sq_z + Chi_sq_s)

# get reduced Chi-square, or 'chi-square per degrees of freedom'
Chi_sq_r <- Chi_sq/v


# ~~~~~~~~~~~~~~~~~~~ print & save results ~~~~~~~~~~~~~~~~~~~#

# print & save obs fit
print(paste("### Reduced Chi-squared  =", round(Chi_sq_r, 3)))

filepath <- paste0(out_path, "chi_sq.rds")
saveRDS(Chi_sq_r, filepath)
