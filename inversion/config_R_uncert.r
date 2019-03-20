# configuration file for model-data mismatch input parameters

# Description: Model-data mismatch decribes the error/uncertainty introduced
# into the model from both the observations used and the model-linkages from
# observation space to flux space (i.e. transport). Model-data mismatch is
# named here as the 'R' matrix and is broken down into distinct error components.
# Not all components are relevant to all inversion problems (and more may even
# be introduced), but these hit most large sources of observational errors

# run dependent scripts
source("config.r")

# load z, the observational enhancements
z_file <- paste0(out_path, "z.rds")
if(!file.exists(z_file))
    stop("obs file (z.rds) must be created before running")

Hsbio_file <- paste0(out_path, "Hs_bio.rds")
if(include_bio & !file.exists(Hsbio_file))
    stop("biological enhancements file (Hs_bio.rds) must be created before running")

z <- readRDS(z_file) #read observed enhancement file
if (include_bio)
  Hsbio <- readRDS(Hsbio_file) #read biological enhancement file


# define the overall vatrix containing all the details of component-level errors
R_arr <- array(NA, dim = c(0, 5))
colnames(R_arr) <- c("component", "RMSE", "correlation_scale", "non_decaying?", "correlate_btwn_days?")


# ~~~~~~~~~~~~~~~ R_part ~~~~~~~~~~~~~~~ #
# uncertainty in STILT based on the finite number of particles released

rmse_part <- 0.1 #ppm
Lt_part <- NA #correlation time scale, in days
correlate_full_part <- F #if TRUE, apply RMSE to off_diagonals without decay
R_part_info <- c("R_part", rmse_part, Lt_part, correlate_full_part, F)
R_arr <- rbind(R_arr, R_part_info)


# ~~~~~~~~~~~~~~~ R_aggr ~~~~~~~~~~~~~~~ #
# Uncertainty introduced by aggregating heterogeneous
# fluxes into one homogeneous grid cell

rmse_aggr <- 0.167 * mean(z) + 0.25 * mean(z) #percent of mean signal (ppm)
Lt_aggr <- NA #correlation time scale, in days
correlate_full_aggr <- T #if TRUE, apply RMSE to off_diagonals without decay
R_aggr_info <- c("R_aggr", rmse_aggr, Lt_aggr, correlate_full_aggr, F)
R_arr <- rbind(R_arr, R_aggr_info)


# ~~~~~~~~~~~~~~~ R_eddy ~~~~~~~~~~~~~~~ #
# variance in CO2 within the mixed-layer profile reflecting
# turbulence which the model is unable to capture

rmse_eddy <- 0
Lt_eddy <- NA #correlation time scale, in days
correlate_full_eddy <- F #if TRUE, apply RMSE to off_diagonals without decay
R_eddy_info <- c("R_eddy", rmse_eddy, Lt_eddy, correlate_full_eddy, F)
R_arr <- rbind(R_arr, R_eddy_info)


# ~~~~~~~~~~~~~~~ R_bkgd ~~~~~~~~~~~~~~~ #
# uncertainty in background values

rmse_bkgd <- 1.9
Lt_bkgd <- NA #correlation time scale, in days
correlate_full_bkgd <- T #if TRUE, apply RMSE to off_diagonals without decay
R_bkgd_info <- c("R_bkgd", rmse_bkgd, Lt_bkgd, correlate_full_bkgd, F)
R_arr <- rbind(R_arr, R_bkgd_info)


# ~~~~~~~~~~~~~~~ R_transPBL ~~~~~~~~~~~~~~~ #
# Uncertainty in the transport model's estimation of the PBL height

rmse_transPBL <- 0.07 * mean(z)
Lt_transPBL <- NA #correlation time scale, in days
correlate_full_transPBL <- T #if TRUE, apply RMSE to off_diagonals without decay
R_transPBL_info <- c("R_transPBL", rmse_transPBL, Lt_transPBL, correlate_full_transPBL, F)
R_arr <- rbind(R_arr, R_transPBL_info)


# ~~~~~~~~~~~~~~~ R_transWIND ~~~~~~~~~~~~~~~ #
# Uncertainty in the transport model's estimation of trajectories based on wind error

rmse_transWIND <- 0.35 * mean(z)
Lt_transWIND <- 8/24 #correlation time scale, in days
correlate_full_transWIND <- F #if TRUE, apply RMSE to off_diagonals without decay
R_transWIND_info <- c("R_transWIND", rmse_transWIND, Lt_transWIND, correlate_full_transWIND, F)
R_arr <- rbind(R_arr, R_transWIND_info)


# ~~~~~~~~~~~~~~~ R_instr ~~~~~~~~~~~~~~~ #
# Uncertainty resulting from biases and standard errors in the instrumentation

rmse_instr <- 0.25
Lt_instr <- NA #correlation time scale, in days
correlate_full_instr <- F #if TRUE, apply RMSE to off_diagonals without decay
R_instr_info <- c("R_instr", rmse_instr, Lt_instr, correlate_full_instr, F)
R_arr <- rbind(R_arr, R_instr_info)


# ~~~~~~~~~~~~~~~ R_bio ~~~~~~~~~~~~~~~#
# Uncertainty introduced from including biospheric enhancement in measurements
if (include_bio){
    rmse_bio <- 0.25 * mean(Hsbio)
    Lt_bio <- NA #correlation time scale, in days
    correlate_full_bio <- F #if TRUE, apply RMSE to off_diagonals without decay
} else {
    rmse_bio <- 0
    Lt_bio <- NA #correlation time scale, in days
    correlate_full_bio <- F #if TRUE, apply RMSE to off_diagonals without decay
}
R_bio_info <- c("R_bio", rmse_bio, Lt_bio, correlate_full_bio, F)
R_arr <- rbind(R_arr, R_bio_info)

# ~~~~~~~~~~~~~~~ R_other? ~~~~~~~~~~~~~~~#
# Any additional representation/miscellaneous errors

rmse_other <- 0
Lt_other <- NA #correlation time scale, in days
correlate_full_other <- F #if TRUE, apply RMSE to off_diagonals without decay
R_other_info <- c("R_other", rmse_other, Lt_other, correlate_full_other, F)
R_arr <- rbind(R_arr, R_other_info)


# ~~~~~~~~~~~~~~~ additional errors for sites ~~~~~~~~~~~~~~~#
extra_sites_err_rmse <- array(0, dim = c(1, nsites))
colnames(extra_sites_err_rmse) <- sites

# example:
# extra_sites_err_rmse$site1 <- 5
