# script to calculate optimized fluxes

# author: Lewis Kunik


## prerequisite scripts:
##  make_sprior.r
##  make_R.r
##  make_zHsp.r
##  make_HQ.r
##
## output files:
##  s_hat.rds - file containing vector of length (#cells * #times) with
##    values of optimized emissions values across domain

# load package dependencies
library(lubridate)

# run dependent scripts
source("config.r")


# specify start/end times
y1 <- flux_year_start
y2 <- flux_year_end

m1 <- flux_month_start
m2 <- flux_month_end

d1 <- flux_day_start
d2 <- flux_day_end

h1 <- flux_hour_start
h2 <- flux_hour_end

mn1 <- flux_min_start
mn2 <- flux_min_end

# list all desired obs times
flux_times <- seq(from = ISOdatetime(y1, m1, d1, h1, mn1, 0, tz = "UTC") + (flux_t_res/2),
    to = ISOdatetime(y2, m2, d2, h2, mn2, 0, tz = "UTC"), by = flux_t_res) #flux_t_res defined in config.r


# ~~~~~~~~~~~~~ load in necessary files ~~~~~~~~~~~~~~~~~#

# load in obs - Hsp
zHsp_file <- paste0(out_path, "zHsp.rds")
zHsp <- readRDS(zHsp_file)

# load in sprior
sprior_file <- paste0(out_path, "sprior.rds")
sprior <- readRDS(sprior_file)

# load in R
R_file <- paste0(out_path, "R.rds")
R <- readRDS(R_file)

# load in HQH^t file
HQHt_file <- paste0(out_path, "HQHt.rds")
HQHt <- readRDS(HQHt_file)

# load in HQ files
nobs <- nrow(R)




# ~~~~~~~~~~~~~~~~~ combine all the pieces of the inversion ~~~~~~~~~~~~~#
# s_hat = s_prior + (HQ)^T (HQH^T + R)^-1 (z - Hs_prior)

# piece up the equation above
HQHtR <- HQHt + R
HQHtR_inv <- solve(HQHtR)
HQHtR_zHsp <- HQHtR_inv %*% zHsp
additive_cor <- array(0, dim = c(1, 0))

# load in HQ files - ntimes defined in config.r
for (ii in 1:ntimes) {
    HQ_i <- readRDS(paste0("HQ/HQ", formatC(ii, width = 3, flag = "0"), ".rds"))
    HQTi_s <- t(HQ_i) %*% HQHtR_zHsp
    additive_cor <- c(additive_cor, HQTi_s)

}

# add the corrections to the prior here
s_hat <- sprior + additive_cor


# save s_hat
print("saving optimized emissions to s_hat.rds")
filepath <- paste0(out_path, "s_hat.rds")
saveRDS(s_hat, filepath)


# ~~~~~~~~~~~~~~~~~~ grab the average flux for the domain ~~~~~~~~~~~~~~~~#

# Define aggregation operator to weight each cell by its fractional area
finecellareas <- abs(2 * pi * (6371009^2) * (sin((lonlat_domain[, 2] - lat_res/2) *
    pi/180) - sin((lonlat_domain[, 2] + lat_res/2) * pi/180))/(360/lon_res))
total_area <- sum(finecellareas)
W <- matrix(finecellareas/total_area, nrow = 1)  #W is the 'aggregation operator'
isubset <- which(hour(flux_times) %in% subset_hours_utc)

#format prior vector to average over domain cells
sprior_mat <- matrix(sprior, nrow = ntimes, byrow = T)
sprior_avgs <- apply(sprior_mat, FUN = function(x) W %*% x, MARGIN = 1) #weight by cell areas
sprior_avg <- mean(sprior_avgs)

#format posterior vector to average over domain cells
shat_mat <- matrix(s_hat, nrow = ntimes, byrow = T)
shat_avgs <- apply(shat_mat, FUN = function(x) W %*% x, MARGIN = 1) #weight by cell areas
shat_avg <- mean(shat_avgs)

if(shat_avg < .01){
	shat_avg_char <- formatC(shat_avg, format = "e", digits = 3)
	sprior_avg_char <- formatC(sprior_avg, format = "e", digits = 3)
} else{
	shat_avg_char <- round(shat_avg, 3)
	sprior_avg_char <- round(sprior_avg, 3)
}


print(paste("Overall domain-averaged prior emissions:", sprior_avg_char, flux_units))
print(paste("Overall domain-averaged posterior emissions:", shat_avg_char, flux_units))


#if time is subsetted, grab the average for the constrained time period as well
sprior_subset_avg <- mean(sprior_avgs[isubset])
shat_subset_avg <- mean(shat_avgs[isubset])

if(!is.na(shat_subset_avg) &  shat_subset_avg < .01){
	shat_avg_char <- formatC(shat_subset_avg, format = "e", digits = 3)
	sprior_avg_char <- formatC(sprior_subset_avg, format = "e", digits = 3)
} else{
	shat_avg_char <- round(shat_subset_avg, 3)
	sprior_avg_char <- round(sprior_subset_avg, 3)
}

subset_hour_begin <- subset_hours_utc[1]
subset_hour_end <- tail(subset_hours_utc, 1)
if ((flux_t_res < 24*3600) & (subset_hour_begin != 0 | subset_hour_end != 23)) {
    print("-----------------------------")

    print(paste(subset_hour_begin, "-", subset_hour_end, "UTC average of prior emissions:",
        sprior_avg_char, flux_units))

    print(paste(subset_hour_begin, "-", subset_hour_end, "UTC average of posterior emissions:",
        shat_avg_char, flux_units))
}
