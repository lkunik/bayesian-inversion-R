# Create spatial covariance matrix
# author: Lewis Kunik

## prerequisite scripts:
##  none
## output files:
##  sp_cov.rds - file containing a (#cells x #cells) matrix of spatial covariance values

# load package dependencies
library(geosphere)

# run dependent scripts
source("config.r")

# ~~~~~~~~~~~~~~~~ Load grid info ~~~~~~~~~~~~~~~~#

lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)


# ~~~~~~~~~~~~~~~~ Calculate Xs and E ~~~~~~~~~~~~~~~~#

Xs <- array(0, dim = c(ncells, ncells))  #Xs will be a matrix of separation distances between all grid-cell pairs

print("Creating spatial covariance matrix (E)")

for (i in 1:ncells) {

    # get distance between this cell and its pairing of all other cells after it in
    # the array of cells (in meters)
    dist_arr_m <- distCosine(p1 = c(lonlat_domain[i, 1], lonlat_domain[i, 2]), p2 = cbind(lonlat_domain[i:ncells,
        1], lonlat_domain[i:ncells, 2]), r = 6372795)
    dist_arr_km <- dist_arr_m * 0.001  #convert to km (dependent on lengthscale being in km)
    Xs[i:ncells, i] <- Xs[i, i:ncells] <- dist_arr_km  #assign to respective places in matrix
}

# Calculate E (spatial covariance matrix), using spatial correlation parameter ls defined in config.r
if (ls > 0) {
    E <- exp(-Xs/ls)
} else {
    E <- diag(rep(1, ncells))
}

print("saving spatial covariance file")

# ~~~~~~~~~~~~~~~~ save spatial cov file ~~~~~~~~~~~~~~~~#

filepath <- paste0(out_path, "sp_cov.rds")
saveRDS(E, filepath)
