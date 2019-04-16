# Compute time-aggregated, grid-scale uncertainty of posterior
# author: Lewis Kunik

# Vshat = (Q_{sum} - (HQ)t_{sum} (HQHt + R)^{-1} (HQ)_{sum}) / k^2
# (where k = num. timesteps)

## prerequisite scripts:
##  make_sp_cov.r
##  make_tmp_cov.r
##  make_sigma.r
##  make_R.r
##  make_HQ.r
##
## output files:
##  Vshat_bar.rds - file containing a (#cells x #cells) matrix of posterior error
##    covariance values
##  perc_unc_red.rds - file containing a single value, percent uncertainty reduction

# load package dependencies
library(lubridate)

# run dependent scripts
source("config.r")


# specify start/end from config.r variables - establish the centers of time steps
flux_times <- seq(from = flux_start_POSIX + flux_t_res/2,
    to = flux_end_POSIX, by = flux_t_res) #flux_t_res defined in config.r

# load spatial covariance matrix E
sp_cov_file <- paste0(out_path, "sp_cov.rds")
E <- readRDS(sp_cov_file)

# load temporal covariance matrix D
tmp_cov_file <- paste0(out_path, "tmp_cov.rds")
D <- readRDS(tmp_cov_file)
# ntimes defined in config.r
if (ntimes == 1)
    D <- as.matrix(D)  #make sure D is in ntime x ntime format even if ntime = 1

# load prior uncertainty and convert to ntime x ncell matrix
sigma_file <- paste0(out_path, "sigma.rds")
sigma <- readRDS(sigma_file)
sigma_mat <- matrix(sigma, nrow = ntimes, byrow = T) # ntimes defined in config.r

# load model data mismatch matrix
R_file <- paste0(out_path, "R.rds")
R <- readRDS(R_file)
nobs <- nrow(R)

# load HQHtranspose file
HQHt_file <- paste0(out_path, "HQHt.rds")
HQHt <- readRDS(HQHt_file)

# load lonlat_domain grid info
lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)

# tstart and tstop are defined in config.r
time_span <- tstop - tstart + 1
S <- sigma_mat[tstart:tstop, ]  #in CT-Lagrange docs, this sub-matrix of sigma is denoted as S

if (time_span == 1)
    S <- matrix(S, nrow = time_span, byrow = T)  #ensure that S is in matrix-conformable format, 1 x ncells

# ~~~~~~~~~~~~~~~~~~~~~~~ calculate and save Qsum ~~~~~~~~~~~~~~~~~~~~~~~#

# make Qsum - formula documented at this link from CT-Lagrange:
# https://www.esrl.noaa.gov/gmd/ccgg/carbontracker-lagrange/doc/inversion.html#posterior-uncertainty
# (equation 5 shows Qsum = (t(S) * D * S) * E, substituting x = t(S) * D.
x <- t(S) %*% D[tstart:tstop, tstart:tstop]
Qsum <- (x %*% S) * E  #Multiplication with E is element-by-element, not matrix multiplication.

# save Qsum (in case desired)
print("saving Qsum to Qsum.rds")
filepath <- paste0(out_path, "Qsum.rds")
saveRDS(Qsum, filepath)


# ~~~~~~~~~~~~~~~~~~~~~~~ calculate and save Vshat_bar ~~~~~~~~~~~~~~~~~~~~~~~#

# sum up HQ blocks
HQsum <- array(0, dim = c(nobs, ncells)) #create empty array to eventually hold sum

# populate HQsum
for (ii in tstart:tstop) {
    HQ_file <- paste0("HQ/HQ", formatC(ii, width = filename_width, flag = "0"), ".rds")
    HQii <- readRDS(HQ_file)
    HQsum <- HQsum + HQii
}

# compute (HQHt + R)^-1
HQHtRi <- solve(HQHt + R)

# get k, number of timesteps - comes from config.r constant variables
k <- (tstop - tstart) + 1
k2 <- k^2

# perform computation from CT-L equation linked above
HQT_sum_HQHtRi <- t(HQsum) %*% HQHtRi
a <- HQT_sum_HQHtRi %*% HQsum
Vshat_bar <- (Qsum - a)/k2

# save files
print("saving V_shat_bar to Vshat_bar.rds")
filepath <- paste0(out_path, "Vshat_bar.rds")
saveRDS(Vshat_bar, filepath)


# ~~~~~~~~~~~~~~~~~~~~ calculate %uncertainty reduction ~~~~~~~~~~~~~~~~~~~~~#

# Define aggregation operator to weight each cell by its fractional area
finecellareas <- abs(2 * pi * (6371009^2) * (sin((lonlat_domain[, 2] - lat_res/2) *
    pi/180) - sin((lonlat_domain[, 2] + lat_res/2) * pi/180))/(360/lon_res))
total_area <- sum(finecellareas)
W <- matrix(finecellareas/total_area, nrow = 1)  #W is the 'aggregation operator'

# multiply uncertainty covariance matrixes by the aggregation operators for
# covariance-included total uncertainties
prior_tot_unc <- sqrt(W %*% (Qsum/k2) %*% t(W))
post_tot_unc <- sqrt(W %*% (Vshat_bar) %*% t(W))
unc_red <- prior_tot_unc - post_tot_unc

if(post_tot_unc < .01){
	post_tot_unc_char <- formatC(post_tot_unc, format = "e", digits = 3)
	prior_tot_unc_char <- formatC(prior_tot_unc, format = "e", digits = 3)
} else{
	post_tot_unc_char <- round(post_tot_unc, 3)
	prior_tot_unc_char <- round(prior_tot_unc, 3)
}

print(paste("prior TOT uncert =", prior_tot_unc_char, flux_units)) #flux_units defined in config.r
print(paste("post TOT uncert =", post_tot_unc_char, flux_units))

perc_red <- (unc_red/prior_tot_unc) * 100
print(paste("percent uncertainty reduction: (prior - posterior)/prior =", round(perc_red,
    2), "%"))

print("saving percent uncertainty reduction to perc_unc_red.rds")
filepath <- paste0(out_path, "perc_unc_red.rds")
saveRDS(perc_red, filepath)


# ~~~~ now determine prior/posterior uncertainty for subsetted-only times ~~~~ #

times_hr <- hour(flux_times)  #this gives you the hour of timestamp at the CENTER of each timestep-bin
isubset <- which(times_hr %in% subset_hours_utc)

if(length(isubset) > 1){

    print("calculating time-subsetted uncertainty reduction")

    # make Qsum_subset for subsetted times - see equation above for how to calculate
    S_subset <- sigma_mat[isubset, ]  #in CT-Lagrange, this sub-matrix of sigma is denoted as S
    x_subset <- t(S_subset) %*% D[isubset, isubset]
    Qsum_subset <- (x_subset %*% S_subset) * E  #Multiplication with E is an element-by-element multiplication, not matrix multiplication.

    # get HQsum for the subset times
    HQsum_subset <- array(0, dim = c(nobs, ncells)) #create empty array to eventually hold sum
    # populate HQsum_subset
    for (ii in isubset) {
        HQ_file <- paste0("HQ/HQ", formatC(ii, width = filename_width, flag = "0"), ".rds")
        HQii <- readRDS(HQ_file)
        HQsum_subset <- HQsum_subset + HQii
    }

    k2_subset <- (length(isubset)^2) # num. times squared

    # make Vshat_bar representing subsetted times
    HQT_sum_HQHtRi_subset <- t(HQsum_subset) %*% HQHtRi
    a_subset <- HQT_sum_HQHtRi_subset %*% HQsum_subset
    Vshat_bar_subset <- (Qsum_subset - a_subset)/k2_subset

    # multiply uncertainty covariance matrixes by the aggregation operators for
    # covariance-included total uncertainties
    prior_subset_unc <- sqrt(W %*% (Qsum_subset/k2_subset) %*% t(W))
    post_subset_unc <- sqrt(W %*% (Vshat_bar_subset) %*% t(W))

    #calculate the percent uncertainty reduction
    subset_perc_red <- ((prior_subset_unc - post_subset_unc)/prior_subset_unc) * 100

    #define subsetted times
    subset_hour_begin <- subset_hours_utc[1]
    subset_hour_end <- tail(subset_hours_utc, 1)

    if ((flux_t_res < 24*3600) & (subset_hour_begin != 0 | subset_hour_end != 23)) {
        print(paste(subset_hour_begin, "-", subset_hour_end, "UTC prior uncert = ", round(prior_subset_unc,
            2), flux_units))
        print(paste(subset_hour_begin, "-", subset_hour_end, "UTC post uncert =", round(post_subset_unc,
            2), flux_units))
        print(paste(subset_hour_begin, "-", subset_hour_end, "UTC percent uncertainty reduction =",
            round(subset_perc_red, 2), "%"))
    }
}
