# calculate observed minus prior enhancements
# author: Lewis Kunik

## prerequisite scripts:
##  make_sprior.r
##  Hsplit.r
##  make_z.r
##
## output files:
##  zHsp.rds - file containing a vector of values equal to (z - Hsp)

# run dependent scripts
source("config.r")

# here is the read sparse h function
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




# load in z (obs - bg)
z_file <- paste0(out_path, "z.rds")
z <- readRDS(z_file)
nobs <- length(z)

# load in prior emissions file
sprior_file <- paste0(out_path, "sprior.rds")
sprior_vec <- readRDS(sprior_file)

# format into matrix so different time-steps are easily accessible
sprior_mat <- matrix(sprior_vec, nrow = ntimes, byrow = T)
ncells <- ncol(sprior_mat)

# get Hsprior by convolving footprints
Hsprior <- rep(0, nobs)

# convolve each H file and add to running total of Hsbio
for (ii in 1:ntimes) {

    # grab H from file
    Hi <- read_sparse_h(ii, nobs, ncells)

    # convolve footprints with prior emissions
    Hsprior <- Hsprior + Hi %*% sprior_mat[ii, ]
}

# subtract Hsprior from z
zHsp <- z - Hsprior

# save Hsprior to file
print("saving prior modeled enhancements to Hsprior.rds")
filepath <- paste0(out_path, "Hsprior.rds")
saveRDS(Hsprior, filepath)

# save zHsp to file
print("saving z - Hsp to zHsp.rds")
filepath <- paste0(out_path, "zHsp.rds")
saveRDS(zHsp, filepath)
